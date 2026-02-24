use std::{
    borrow::Cow,
    marker::PhantomData,
    ops::Range,
    path::PathBuf,
    sync::{Arc, OnceLock},
};

use context_error::*;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use uom::ConstZero;

use crate::{
    BoxedIdentifiedPeptideIter, FastaIdentifier, KnownFileFormat, PSM, PSMData,
    PSMFileFormatVersion, PSMMetaData, PSMSource, PeptidoformPresent, SpectrumId, SpectrumIds,
    common_parser::{Location, OptionalColumn, OptionalLocation},
    helper_functions::explain_number_error,
};
use mzcore::{
    chemistry::Chemical,
    csv::{CsvLine, parse_csv},
    molecular_formula,
    ontology::Ontologies,
    quantities::{Tolerance, WithinTolerance},
    sequence::{
        CompoundPeptidoformIon, CrossLinkName, FlankingSequence, Linked, Modification, Peptidoform,
        PeptidoformIon, SequencePosition, SimpleModification, SimpleModificationInner,
        SloppyParsingParameters,
    },
    system::{Mass, MassOverCharge, Time, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid pLink line",
    "This column is not a number but it is required to be a number in this pLink format",
);
static TYPE_ERROR: (&str, &str) = (
    "Invalid pLink peptide type",
    "This column is not a valid paptide type but it is required to be one of 0/1/2/3 in this pLink format",
);

format_family!(
    PLink,
    Linked, PeptidoformPresent, [&V2_3], b',', None;
    required {
        order: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        title: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        z: Charge, |location: Location, _| location.parse::<isize>(NUMBER_ERROR).map(Charge::new::<mzcore::system::e>);
        /// MH+ mass
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<mzcore::system::dalton>);
        /// MH+ mass
        theoretical_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<mzcore::system::dalton>);
        peptide_type: PLinkPeptideType, |location: Location, _| location.parse::<PLinkPeptideType>(TYPE_ERROR);
        peptidoform: PeptidoformIon, |location: Location, ontologies: &Ontologies| {
            match plink_separate(&location, "Invalid pLink peptide", "peptide").map_err(BoxedError::to_owned)? {
                (pep1, Some(pos1), Some(pep2), Some(pos2)) => {
                    let pep1 = Peptidoform::sloppy_pro_forma_inner(&location.base_context(), location.full_line(), pep1, ontologies, &SloppyParsingParameters::default()).map_err(BoxedError::to_owned)?;
                    let pep2 = Peptidoform::sloppy_pro_forma_inner(&location.base_context(), location.full_line(), pep2, ontologies, &SloppyParsingParameters::default()).map_err(BoxedError::to_owned)?;

                    let mut peptidoform = PeptidoformIon::from_vec(vec![pep1.into(), pep2.into()]).unwrap();
                    peptidoform.add_cross_link(
                        (0, SequencePosition::Index(pos1.0.saturating_sub(1) as usize)),
                        (1, SequencePosition::Index(pos2.0.saturating_sub(1) as usize)),
                        Arc::new(SimpleModificationInner::Mass(mzcore::sequence::MassTag::None, Mass::default().into(), None)),
                        CrossLinkName::Name("1".to_string().into_boxed_str()),
                    );
                    Ok(peptidoform)
                }
                (pep1, Some(pos1), None, Some(pos2)) => {
                    let pep = Peptidoform::sloppy_pro_forma_inner(&location.base_context(), location.full_line(), pep1, ontologies, &SloppyParsingParameters::default()).map_err(BoxedError::to_owned)?;

                    let mut peptidoform = PeptidoformIon::from_vec(vec![pep.into()]).unwrap();
                    peptidoform.add_cross_link(
                        (0, SequencePosition::Index(pos1.0.saturating_sub(1) as usize)),
                        (0, SequencePosition::Index(pos2.0.saturating_sub(1) as usize)),
                        Arc::new(SimpleModificationInner::Mass(mzcore::sequence::MassTag::None, Mass::default().into(), None)),
                        CrossLinkName::Name("1".to_string().into_boxed_str()),
                    );
                    Ok(peptidoform)
                }
                (pep1, Some(pos1), None, None) => {
                    let mut pep = Peptidoform::sloppy_pro_forma_inner(&location.base_context(), location.full_line(), pep1, ontologies, &SloppyParsingParameters::default()).map_err(BoxedError::to_owned)?;
                    pep[SequencePosition::Index(pos1.0.saturating_sub(1) as usize)].modifications.push( SimpleModificationInner::Mass(mzcore::sequence::MassTag::None, Mass::default().into(), None).into());

                    Ok(PeptidoformIon::from_vec(vec![pep.into()]).unwrap())
                }
                (pep1, None, None, None) => {
                    let pep = Peptidoform::sloppy_pro_forma_inner(&location.base_context(), location.full_line(), pep1, ontologies, &SloppyParsingParameters::default()).map_err(BoxedError::to_owned)?;

                    Ok(PeptidoformIon::from_vec(vec![pep.into()]).unwrap())
                }
                _ => unreachable!()
            }
        };
        refined_score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        svm_score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        e_value: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        /// Whether this peptide is a target (false) or decoy (true) peptide
        is_decoy: bool, |location: Location, _| Ok(location.as_str() == "1");
        q_value: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        /// The proteins, per protein the first protein, location, second protein, and location
        proteins: Vec<(FastaIdentifier<Box<str>>, Option<u16>, Option<FastaIdentifier<Box<str>>>, Option<u16>)>, |location: Location, _|  {
            location.array('/').filter(|l| !l.as_str().trim().is_empty()).map(|l| {
                let separated = plink_separate(&l, "Invalid pLink protein", "protein").map_err(BoxedError::to_owned)?;

                Ok(((
                    l.full_line()[separated.0.clone()].trim()).parse().map_err(|err|
                        BoxedError::new(BasicKind::Error,
                        "Invalid pLink modification",
                        format!("A pLink protein should be a valid fasta identifier but the number {}", explain_number_error(&err)),
                        l.context().add_highlight((0, separated.0)).to_owned()))?,
                    separated.1.map(|(a, _)| a),
                    separated.2.map(|p| l.full_line()[p.clone()].trim().parse().map_err(|err|
                        BoxedError::new(BasicKind::Error,
                        "Invalid pLink modification",
                        format!("A pLink protein should be a valid fasta identifier but the number {}", explain_number_error(&err)),
                        l.context().add_highlight((0, p)).to_owned()))).transpose()?,
                    separated.3.map(|(a, _)| a)))
            })
            .collect::<Result<Vec<_>, _>>().map(|mut v| {v.shrink_to_fit(); v})
        };
        /// If true this indicates that this cross-link binds two different proteins
        is_different_protein: bool, |location: Location, _| Ok(location.as_str() == "1");
        raw_file_id: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        is_complex_satisfied: bool, |location: Location, _| Ok(location.as_str() == "1");
        /// Whether this fits within the normal filters applied within pLink
        is_filter_in: bool, |location: Location, _| Ok(location.as_str() == "1");
    }
    optional {
        scan_number: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        raw_file: PathBuf, |location: Location, _| {
            let mut p = PathBuf::from(location.get_string());
            p.shrink_to_fit();
            Ok(Some(p))
        };
    }

    #[expect(clippy::similar_names)]
    fn post_process(source: &CsvLine, mut parsed: Self, ontologies: &Ontologies) -> Result<Self, BoxedError<'static, BasicKind>> {
        // Add all modifications
        let pep1 = parsed.peptidoform.peptidoforms()[0].len();
        let pep2 = parsed.peptidoform.peptidoforms().get(1).map_or(0, Peptidoform::len);

        // All modifications with their attachment, and their index (into the full peptidoform, so anything bigger than the first peptide matches in the second)
        for v in source.column("modifications")?.ignore("null").array(';') {
            let v = v.trim();
            let position_start = v.as_str().rfind('(').ok_or_else(||
                BoxedError::new(BasicKind::Error,
                    "Invalid pLink modification",
                    "A pLink modification should follow the format 'Modification[AA](pos)' but the opening bracket '(' was not found",
                    v.context().to_owned()))?;
            let location_start = v.as_str().rfind('[').ok_or_else(||
                BoxedError::new(BasicKind::Error,
                    "Invalid pLink modification",
                    "A pLink modification should follow the format 'Modification[AA](pos)' but the opening square bracket '[' was not found",
                    v.context().to_owned()))?;
            let index = v.full_line()[v.range.start+position_start+1..v.range.end-1].parse::<usize>().map_err(|err|
                BoxedError::new(BasicKind::Error,
                    "Invalid pLink modification",
                    format!("A pLink modification should follow the format 'Modification[AA](pos)' but the position number {}", explain_number_error(&err)),
                    v.context().to_owned()))?;

            let m = Modification::sloppy_modification(v.full_line(), v.range.start..v.range.start+location_start, None, ontologies).map_err(BoxedError::to_owned)?;

            let pos = if index == 0 {
                (0, SequencePosition::NTerm)
            } else if (1..=pep1).contains(&index) {
                (0, SequencePosition::Index(index-1))
            } else if index == pep1+1 {
                (0, SequencePosition::CTerm)
            } else if index == pep1+2 {
                return Err(BoxedError::new(BasicKind::Error,
                    "Invalid modification location",
                    format!("A modification cannot be located here, at location {}, between the two peptides, presumable on the cross-linker", index+1),
                    source.full_context().to_owned(),
                ));
            } else if index == pep1+3 {
                (1, SequencePosition::NTerm)
            } else if (pep1+4..=pep1+4+pep2).contains(&index) {
                (1, SequencePosition::Index(index-(pep1+4)))
            } else if index == pep1+4+pep2+1 {
                (1, SequencePosition::CTerm)
            } else {
                return Err(BoxedError::new(BasicKind::Error,
                    "Invalid modification location",
                    format!("A modification cannot be located here, at location {}, after peptide 2 has ended (maximal valid location {})", index+1,pep1+4+pep2+2),
                    source.full_context().to_owned(),
                ));
            };

            match pos {
                (peptide, SequencePosition::NTerm) => {
                        parsed.peptidoform.peptidoforms_mut()[peptide].add_simple_n_term(m.clone());
                }
                (peptide, SequencePosition::CTerm) => {
                    parsed.peptidoform.peptidoforms_mut()[peptide].add_simple_c_term(m.clone());
                }
                (peptide, index) => {
                    parsed.peptidoform.peptidoforms_mut()[peptide][index]
                            .modifications
                            .push(m.clone().into());
                }
            }
        }

        if parsed.peptide_type != PLinkPeptideType::Common {
            // Find linker based on left over mass
            let left_over = parsed.theoretical_mass
                - molecular_formula!(H 1 :z+1).monoisotopic_mass()
                - parsed
                    .peptidoform
                    .formulas() // TODO: this is the hot function this takes ~51% of runtime when parsing a pLink file
                    .first()
                    .unwrap()
                    .monoisotopic_mass()
                    - if parsed.peptide_type == PLinkPeptideType::Hydrolysed { molecular_formula!(H 2 O 1).monoisotopic_mass() } else { Mass::ZERO };

            let custom_linkers: Vec<_> = ontologies.custom().data().iter().filter(|m|
                    matches!(***m, SimpleModificationInner::Linker{..})).map(|m| (m.formula().monoisotopic_mass(), m.clone())
                ).collect();

            let fitting = &KNOWN_CROSS_LINKERS.get_or_init(|| {
    [
        ontologies.psimod().get_by_index(&34).unwrap(), // Disulfide: M:L-cystine (cross-link)
        ontologies.unimod().get_by_index(&1905).unwrap(), // BS2G: U:Xlink:BS2G[96]
        ontologies.xlmod().get_by_index(&2008).unwrap(), // BS2G heavy: X:BS2G-d4
        ontologies.unimod().get_by_index(&1898).unwrap(), // DSS: U:Xlink:DSS[138]
        ontologies.xlmod().get_by_index(&2227).unwrap(), // Leiker_clv: X:PL
        ontologies.unimod().get_by_index(&1896).unwrap(), // DSSO: U:Xlink:DSSO[158]
        ontologies.unimod().get_by_index(&1899).unwrap(), // DSBU: U:Xlink:BuUrBu[196]
        ontologies.xlmod().get_by_index(&2115).unwrap(), // BAMG: X:BAMG
        ontologies.xlmod().get_by_index(&2002).unwrap(), // DSS heavy: X:DSS-d4
        ontologies.unimod().get_by_index(&2058).unwrap(), // PhoX: U:Xlink:DSPP[210]
        ontologies.xlmod().get_by_index(&2010).unwrap(), // DMTMM: X:1-ethyl-3-(3-Dimethylaminopropyl)carbodiimide hydrochloride
    ]
    .into_iter()
    .map(|m| (m.formula().monoisotopic_mass(), m))
    .collect()}).iter().chain(custom_linkers.iter()).filter(|(mass, _)| Tolerance::<Mass>::Absolute(Mass::new::<mzcore::system::dalton>(0.01)).within(mass, &left_over)).map(|(_, m)| m).collect_vec();

            match fitting.len() {
                0 => return Err(BoxedError::new(BasicKind::Error,"Invalid pLink peptide", format!("The correct cross-linker could not be identified with mass {:.3} Da, if a non default cross-linker was used add this as a custom linker modification.", left_over.value), source.full_context().to_owned())),
                1 => {
                    // Replace 0 mass mod + determine Nterm or side chain
                    for p in parsed.peptidoform.peptidoforms_mut() {
                        let mut n_term = p.get_n_term().to_vec();
                        let mut c_term = p.get_c_term().to_vec();
                        let len = p.len();

                        for (seq_index, seq) in p.sequence_mut().iter_mut().enumerate() {
                            let is_n_term = seq_index == 0;
                            let is_c_term = seq_index == len;
                            let seq_clone = seq.clone();
                            let mut remove = None;
                            for (index, mut m) in seq.modifications.iter_mut().enumerate() {
                                if let Modification::CrossLink{name, linker, ..} = &mut m{
                                    if name == &CrossLinkName::Name("1".to_string().into_boxed_str()) {
                                        *linker = fitting[0].clone();

                                        if is_n_term && m.is_possible(&seq_clone, SequencePosition::NTerm).any_possible() {
                                            remove = Some(index);
                                            n_term.push(m.clone());
                                        } else if is_c_term && m.is_possible(&seq_clone, SequencePosition::CTerm).any_possible() {
                                            remove = Some(index);
                                            c_term.push(m.clone());
                                        }
                                    }
                                } else if Modification::Simple(Arc::new(SimpleModificationInner::Mass(mzcore::sequence::MassTag::None, Mass::default().into(), None))) == *m {
                                    *m = Modification::Simple(fitting[0].clone());
                                }
                            }
                            if let Some(i) = remove {
                                seq.modifications.remove(i);
                            }
                        }
                        p.set_n_term(n_term);
                        p.set_c_term(c_term);
                    }
                },
                _ => return Err(BoxedError::new(BasicKind::Error,"Invalid pLink peptide", "The correct cross-linker could not be identified, there are multiple cross-linkers within the tolerance bounds.", source.full_context().to_owned())),
            }
        }

        let split = parsed.title.rsplitn(6, '.').collect::<Vec<&str>>();
        if split.len() == 6
        {
            parsed.raw_file = Some(split[0].into());
            parsed.scan_number = Some(split[1].parse::<usize>().map_err(|err| BoxedError::new(BasicKind::Error, "Invalid pLink title", format!("The scan number {}", explain_number_error(&err)), Context::none().lines(0, &*parsed.title).add_highlight((0, split[0].len() + 1, split[1].len())).to_owned()))?);

        }
        parsed.peptidoform.shrink_to_fit();
        Ok(parsed)
    }
);

/// The static known cross-linkers
static KNOWN_CROSS_LINKERS: OnceLock<Vec<(Mass, SimpleModification)>> = OnceLock::new();

/// Separate the pLink format of 'pep(pos)-pep(pos)' in all possible combinations
/// # Errors
/// If the format is invalid
fn plink_separate<'a>(
    location: &'a Location<'a>,
    title: &'static str,
    field: &'static str,
) -> Result<
    (
        Range<usize>,
        Option<(u16, Range<usize>)>,
        Option<Range<usize>>,
        Option<(u16, Range<usize>)>,
    ),
    BoxedError<'a, BasicKind>,
> {
    if let Some((peptide1, peptide2)) = location.as_str().split_once(")-") {
        let first_end = peptide1.rfind('(').ok_or_else(||
            BoxedError::new(BasicKind::Error,
                title,
                format!("A pLink {field} should follow the format 'PEP1(pos1)-PEP2(pos2)' but the opening bracket '(' was not found for PEP1"),
                Context::line(Some(location.line.line_index() as u32), location.full_line(), location.range.start, peptide1.len()).to_owned()))?;
        let second_end = peptide2.rfind('(').ok_or_else(||
            BoxedError::new(BasicKind::Error,
                title,
                format!("A pLink {field} should follow the format 'PEP1(pos1)-PEP2(pos2)' but the opening bracket '(' was not found for PEP2"),
                Context::line(Some(location.line.line_index() as u32), location.full_line(), location.range.start+peptide1.len()+2, peptide2.len())))?;

        let pos1 = location.range.start + first_end + 1..location.range.start + peptide1.len();
        let first_index = location.full_line()[pos1.clone()].parse::<u16>().map_err(|err|
            BoxedError::new(BasicKind::Error,
                title,
                format!("A pLink {field} should follow the format 'PEP1(pos1)-PEP2(pos2)' but the position for PEP1 {}", explain_number_error(&err)),
                Context::line_range(Some(location.line.line_index() as u32), location.full_line(), pos1.clone()).to_owned()))?;
        let pos2 = location.range.start + peptide1.len() + 2 + second_end + 1
            ..location.range.start + peptide1.len() + 2 + peptide2.len() - 1;
        let second_index = location.full_line()[pos2.clone()].parse::<u16>().map_err(|err|
            BoxedError::new(BasicKind::Error,
                title,
                format!("A pLink {field} should follow the format 'PEP1(pos1)-PEP2(pos2)' but the position for PEP1 {}", explain_number_error(&err)),
                Context::line_range(Some(location.line.line_index() as u32), location.full_line(), pos2.clone()).to_owned()))?;

        Ok((
            location.range.start..location.range.start + first_end,
            Some((first_index, pos1)),
            Some(
                location.range.start + peptide1.len() + 2
                    ..location.range.start + peptide1.len() + 2 + second_end,
            ),
            Some((second_index, pos2)),
        ))
    } else {
        // rsplit to prevent picking a bracket in the text field, and then reverse for it to make sense to human brains
        let mut split = location.as_str().rsplitn(3, '(').collect_vec();
        split.reverse();

        match split.len() {
            3 => {
                let start = location.range.start;
                let start_pos1 = start + split[0].len() + 1;
                let start_pos2 = start_pos1 + split[1].len() + 1;
                let end = location.range.end;

                let pos1 = start_pos1..start_pos2 - 2;
                let first_index = location.full_line()[pos1.clone()].parse::<u16>().map_err(|err|
                    BoxedError::new(BasicKind::Error,
                        title,
                        format!("A pLink {field} should follow the format 'PEP(pos1)(pos2)' but the first position {}", explain_number_error(&err)),
                        Context::line_range(Some(location.line.line_index() as u32), location.full_line(), pos1.clone()).to_owned()))?;
                let pos2 = start_pos2..end - 1;
                let second_index = location.full_line()[start_pos2..end-1].parse::<u16>().map_err(|err|
                    BoxedError::new(BasicKind::Error,
                        title,
                        format!("A pLink {field} should follow the format 'PEP(pos1)(pos2)' but the second position {}", explain_number_error(&err)),
                        Context::line_range(Some(location.line.line_index() as u32), location.full_line(), start_pos2..end-1).to_owned()))?;

                Ok((
                    start..start_pos1 - 1,
                    Some((first_index, pos1)),
                    None,
                    Some((second_index, pos2)),
                ))
            }
            2 => {
                let start = location.range.start;
                let start_pos1 = start + split[0].len() + 1;
                let end = location.range.end;

                let pos1 = start_pos1..end - 1;
                let first_index = location.full_line()[pos1.clone()].parse::<u16>().map_err(|err|
                    BoxedError::new(BasicKind::Error,
                        title,
                        format!("A pLink {field} should follow the format 'PEP(pos1)(pos2)' but the first position {}", explain_number_error(&err)),
                        Context::line_range(Some(location.line.line_index() as u32), location.full_line(), pos1.clone()).to_owned()))?;

                Ok((start..start_pos1 - 1, Some((first_index, pos1)), None, None))
            }
            1 => Ok((location.range.clone(), None, None, None)),
            _ => unreachable!(),
        }
    }
}

/// The different types of peptides a cross-link experiment can result in
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum PLinkPeptideType {
    #[default]
    /// No cross-linkers
    Common,
    /// A cross-linker, but hydrolysed/monolinker
    Hydrolysed,
    /// A cross-linker binding to the same peptide in a loop
    LoopLink,
    /// A cross-linker binding to a different peptide (although the peptide can be identical)
    IntraLink,
}

impl std::str::FromStr for PLinkPeptideType {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "0" => Ok(Self::Common),
            "1" => Ok(Self::Hydrolysed),
            "2" => Ok(Self::LoopLink),
            "3" => Ok(Self::IntraLink),
            _ => Err(()),
        }
    }
}

impl std::fmt::Display for PLinkPeptideType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Common => "Common",
                Self::Hydrolysed => "Hydrolysed",
                Self::LoopLink => "Loop link",
                Self::IntraLink => "Intra link",
            }
        )
    }
}

impl mzcore::space::Space for PLinkPeptideType {
    fn space(&self) -> mzcore::space::UsedSpace {
        mzcore::space::UsedSpace::stack(size_of::<Self>())
    }
}

/// The only built in version of pLink export
pub const V2_3: PLinkFormat = PLinkFormat {
    version: PLinkVersion::V2_3,
    order: "order",
    title: "title",
    z: "charge",
    mass: "precursor_mh",
    peptide_type: "peptide_type",
    peptidoform: "peptide",
    theoretical_mass: "peptide_mh",
    refined_score: "refined_score",
    svm_score: "svm_score",
    score: "score",
    e_value: "e-value",
    is_decoy: "target_decoy",
    q_value: "q-value",
    proteins: "proteins",
    is_different_protein: "protein_type",
    raw_file_id: "fileid",
    is_complex_satisfied: "iscomplexsatisfied",
    is_filter_in: "isfilterin",
    scan_number: OptionalColumn::NotAvailable,
    raw_file: OptionalColumn::NotAvailable,
};

/// All possible pLink versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum PLinkVersion {
    /// Built for pLink version 2.3.11, likely works more broadly
    #[default]
    V2_3,
}

impl std::fmt::Display for PLinkVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl PSMFileFormatVersion<PLinkFormat> for PLinkVersion {
    fn format(self) -> PLinkFormat {
        match self {
            Self::V2_3 => V2_3,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V2_3 => "v2.3",
        }
    }
}

impl PSMMetaData for PLinkPSM {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptidoform.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::PLink(self.version)
    }

    fn numerical_id(&self) -> Option<usize> {
        Some(self.order)
    }

    fn id(&self) -> String {
        self.order.to_string()
    }

    fn search_engine(&self) -> Option<mzcv::Term> {
        Some(mzcv::term!(MS:1003432|pLink2))
    }

    fn confidence(&self) -> Option<f64> {
        Some(1.0 - self.score)
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<(f64, mzcv::Term)> {
        Some((
            self.score,
            mzcv::term!(MS:1001153|search engine specific score),
        ))
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        None
    }

    fn charge(&self) -> Option<Charge> {
        Some(self.z)
    }

    fn mode(&self) -> Option<Cow<'_, str>> {
        None
    }

    fn retention_time(&self) -> Option<Time> {
        None
    }

    fn scans(&self) -> SpectrumIds {
        self.scan_number.map_or_else(
            || SpectrumIds::FileNotKnown(vec![SpectrumId::Native(self.title.to_string())]),
            |scan_number| {
                self.raw_file.clone().map_or_else(
                    || SpectrumIds::FileNotKnown(vec![SpectrumId::Number(scan_number)]),
                    |raw_file| {
                        SpectrumIds::FileKnown(vec![(
                            raw_file,
                            vec![SpectrumId::Number(scan_number)],
                        )])
                    },
                )
            },
        )
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        Some(MassOverCharge::new::<mzcore::system::thomson>(
            self.mass.value / self.z.to_float().value,
        ))
    }

    fn experimental_mass(&self) -> Option<Mass> {
        Some(self.mass)
    }

    type Protein = FastaIdentifier<Box<str>>;
    fn proteins(&self) -> Cow<'_, [Self::Protein]> {
        Cow::Owned(
            self.proteins
                .iter()
                .flat_map(|p| p.2.as_ref().map_or_else(|| vec![&p.0], |p2| vec![&p.0, p2]))
                .unique()
                .cloned()
                .collect_vec(),
        )
    }

    fn protein_location(&self) -> Option<Range<u16>> {
        if let Ok((_, loc, None, _)) = self.proteins.iter().exactly_one() {
            // find loc in peptide by searching for the first location in the first peptide with a linker and calculating the correct offset from there
            loc.and_then(|loc| {
                self.peptidoform.peptidoforms().first().and_then(|p| {
                    p.sequence()
                        .iter()
                        .position(|s| {
                            s.modifications.iter().any(|m| {
                                m.is_cross_link()
                                    || m.clone().into_simple().is_some_and(|s| {
                                        matches!(*s, SimpleModificationInner::Linker { .. })
                                    })
                            })
                        })
                        .and_then(|pos| {
                            let base = loc.checked_sub(u16::try_from(pos).ok()?)?;
                            Some(base..base + u16::try_from(p.len()).ok()?)
                        })
                })
            })
        } else {
            None
        }
    }

    fn flanking_sequences(&self) -> (&FlankingSequence, &FlankingSequence) {
        (&FlankingSequence::Unknown, &FlankingSequence::Unknown)
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        None
    }

    fn unique(&self) -> Option<bool> {
        None
    }

    fn reliability(&self) -> Option<crate::Reliability> {
        None
    }

    fn uri(&self) -> Option<String> {
        None
    }
}
