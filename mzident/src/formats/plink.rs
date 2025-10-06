use std::{
    borrow::Cow,
    marker::PhantomData,
    ops::Range,
    path::PathBuf,
    sync::{Arc, LazyLock},
};

use context_error::*;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use uom::ConstZero;

use crate::{
    BoxedIdentifiedPeptideIter, FastaIdentifier, FlankingSequence, IdentifiedPeptidoform,
    IdentifiedPeptidoformData, IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion,
    KnownFileFormat, MetaData, PeptidoformPresent, SpectrumId, SpectrumIds,
    common_parser::{Location, OptionalColumn, OptionalLocation},
    csv::{CsvLine, parse_csv},
    helper_functions::explain_number_error,
};
use mzcore::{
    chemistry::Chemical,
    molecular_formula,
    ontology::{CustomDatabase, Ontology},
    quantities::{Tolerance, WithinTolerance},
    sequence::{
        CompoundPeptidoformIon, CrossLinkName, Linked, Modification, Peptidoform, PeptidoformIon,
        SequencePosition, SimpleModification, SimpleModificationInner, SloppyParsingParameters,
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
        title: String, |location: Location, _| Ok(location.get_string());
        z: Charge, |location: Location, _| location.parse::<isize>(NUMBER_ERROR).map(Charge::new::<mzcore::system::e>);
        /// MH+ mass
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<mzcore::system::dalton>);
        /// MH+ mass
        theoretical_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<mzcore::system::dalton>);
        peptide_type: PLinkPeptideType, |location: Location, _| location.parse::<PLinkPeptideType>(TYPE_ERROR);
        peptidoform: PeptidoformIon, |location: Location, _| {
            match plink_separate(&location, "Invalid pLink peptide", "peptide").map_err(BoxedError::to_owned)? {
                (pep1, Some(pos1), Some(pep2), Some(pos2)) => {
                    let pep1 = Peptidoform::sloppy_pro_forma(location.full_line(), pep1, None, &SloppyParsingParameters::default()).map_err(BoxedError::to_owned)?;
                    let pep2 = Peptidoform::sloppy_pro_forma(location.full_line(), pep2, None, &SloppyParsingParameters::default()).map_err(BoxedError::to_owned)?;

                    let mut peptidoform = PeptidoformIon::from_vec(vec![pep1.into(), pep2.into()]).unwrap();
                    peptidoform.add_cross_link(
                        (0, SequencePosition::Index(pos1.0.saturating_sub(1) as usize)),
                        (1, SequencePosition::Index(pos2.0.saturating_sub(1) as usize)),
                        Arc::new(SimpleModificationInner::Mass(Mass::default().into())),
                        CrossLinkName::Name("1".to_string()),
                    );
                    Ok(peptidoform)
                }
                (pep1, Some(pos1), None, Some(pos2)) => {
                    let pep = Peptidoform::sloppy_pro_forma(location.full_line(), pep1, None, &SloppyParsingParameters::default()).map_err(BoxedError::to_owned)?;

                    let mut peptidoform = PeptidoformIon::from_vec(vec![pep.into()]).unwrap();
                    peptidoform.add_cross_link(
                        (0, SequencePosition::Index(pos1.0.saturating_sub(1) as usize)),
                        (0, SequencePosition::Index(pos2.0.saturating_sub(1) as usize)),
                        Arc::new(SimpleModificationInner::Mass(Mass::default().into())),
                        CrossLinkName::Name("1".to_string()),
                    );
                    Ok(peptidoform)
                }
                (pep1, Some(pos1), None, None) => {
                    let mut pep = Peptidoform::sloppy_pro_forma(location.full_line(), pep1, None, &SloppyParsingParameters::default()).map_err(BoxedError::to_owned)?;
                    pep[SequencePosition::Index(pos1.0.saturating_sub(1) as usize)].modifications.push( SimpleModificationInner::Mass(Mass::default().into()).into());

                    Ok(PeptidoformIon::from_vec(vec![pep.into()]).unwrap())
                }
                (pep1, None, None, None) => {
                    let pep = Peptidoform::sloppy_pro_forma(location.full_line(), pep1, None, &SloppyParsingParameters::default()).map_err(BoxedError::to_owned)?;

                    Ok(PeptidoformIon::from_vec(vec![pep.into()]).unwrap())
                }
                _ => unreachable!()
            }
        };
        /// All modifications with their attachment, and their index (into the full peptidoform, so anything bigger than the first peptide matches in the second)
        ptm: Vec<(SimpleModification, u16)>, |location: Location, custom_database: Option<&CustomDatabase>|
            location.ignore("null").array(';').map(|v| {
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
                let position = v.full_line()[v.location.start+position_start+1..v.location.end-1].parse::<u16>().map_err(|err|
                    BoxedError::new(BasicKind::Error,
                        "Invalid pLink modification",
                        format!("A pLink modification should follow the format 'Modification[AA](pos)' but the position number {}", explain_number_error(&err)),
                        v.context().to_owned()))?;

                Ok((Modification::sloppy_modification(v.full_line(), v.location.start..v.location.start+location_start, None, custom_database).map_err(BoxedError::to_owned)?, position))
            }
        ).collect::<Result<Vec<_>,_>>();
        refined_score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        svm_score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        e_value: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        /// Whether this peptide is a target (false) or decoy (true) peptide
        is_decoy: bool, |location: Location, _| Ok(location.as_str() == "1");
        q_value: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        proteins: Vec<(String, Option<u16>, Option<String>, Option<u16>)>, |location: Location, _|  {
            location.array('/').filter(|l| !l.as_str().trim().is_empty()).map(|l| {
                let separated = plink_separate(&l, "Invalid pLink protein", "protein").map_err(BoxedError::to_owned)?;

                Ok((l.full_line()[separated.0].trim().to_string(), separated.1.map(|(a, _)| a), separated.2.map(|p| l.full_line()[p].trim().to_string()), separated.3.map(|(a, _)| a)))
            })
            .collect::<Result<Vec<_>, _>>()
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
        raw_file: PathBuf, |location: Location, _| Ok(Some(location.get_string().into()));
    }

    #[expect(clippy::similar_names)]
    fn post_process(source: &CsvLine, mut parsed: Self, custom_database: Option<&CustomDatabase>) -> Result<Self, BoxedError<'static, BasicKind>> {
        // Add all modifications
        let pep1 = parsed.peptidoform.peptidoforms()[0].len();
        let pep2 = parsed.peptidoform.peptidoforms().get(1).map_or(0, Peptidoform::len);
        for (m, index) in &parsed.ptm {
            let index = *index as usize;
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
                - molecular_formula!(H 1 Electron -1).monoisotopic_mass()
                - parsed
                    .peptidoform
                    .formulas()
                    .first()
                    .unwrap()
                    .monoisotopic_mass()
                    - if parsed.peptide_type == PLinkPeptideType::Hydrolysed { molecular_formula!(H 2 O 1).monoisotopic_mass() } else { Mass::ZERO };

            let custom_linkers = custom_database.map_or(
                Vec::new(),
                |c| c.iter().filter(|(_,_,m)|
                    matches!(**m, SimpleModificationInner::Linker{..})).map(|(_,_,m)| (m.formula().monoisotopic_mass(), m.clone())
                ).collect());

            let fitting = &KNOWN_CROSS_LINKERS.iter().chain(custom_linkers.iter()).filter(|(mass, _)| Tolerance::<Mass>::Absolute(Mass::new::<mzcore::system::dalton>(0.001)).within(mass, &left_over)).map(|(_, m)| m).collect_vec();

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
                                    if name == &CrossLinkName::Name("1".to_string()) {
                                        *linker = fitting[0].clone();

                                        if is_n_term && m.is_possible(&seq_clone, SequencePosition::NTerm).any_possible() {
                                            remove = Some(index);
                                            n_term.push(m.clone());
                                        } else if is_c_term && m.is_possible(&seq_clone, SequencePosition::CTerm).any_possible() {
                                            remove = Some(index);
                                            c_term.push(m.clone());
                                        }
                                    }
                                } else if Modification::Simple(Arc::new(SimpleModificationInner::Mass(Mass::default().into()))) == *m {
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

        if let Some(m) = IDENTIFER_REGEX
            .captures(&parsed.title)
        {
            parsed.raw_file = Some(m.get(1).unwrap().as_str().into());
            parsed.scan_number = Some(m.get(2).unwrap().as_str().parse::<usize>().unwrap());
        }
        Ok(parsed)
    }
);

/// The Regex to match against pLink title fields
static IDENTIFER_REGEX: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"([^/]+)\.(\d+)\.\d+.\d+.\d+.\w+").unwrap());
/// The static known cross-linkers
static KNOWN_CROSS_LINKERS: LazyLock<Vec<(Mass, SimpleModification)>> = LazyLock::new(|| {
    [
        Ontology::Unimod.find_id(1898, None).unwrap(), // DSS: U:Xlink:DSS[138]
        Ontology::Xlmod.find_id(2002, None).unwrap(),  // DSS heavy: X:DSS-d4
        Ontology::Psimod.find_id(34, None).unwrap(),   // Disulfide: M:L-cystine (cross-link)
        Ontology::Unimod.find_id(1905, None).unwrap(), // BS2G: U:Xlink:BS2G[96]
        Ontology::Xlmod.find_id(2008, None).unwrap(),  // BS2G heavy: X:BS2G-d4
        Ontology::Xlmod.find_id(2010, None).unwrap(), // DMTMM: X:1-ethyl-3-(3-Dimethylaminopropyl)carbodiimide hydrochloride
        Ontology::Unimod.find_id(1896, None).unwrap(), // DSSO: U:Xlink:DSSO[158]
    ]
    .into_iter()
    .map(|m| (m.formula().monoisotopic_mass(), m))
    .collect()
});

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
                Context::line(Some(location.line.line_index() as u32), location.full_line(), location.location.start, peptide1.len()).to_owned()))?;
        let second_end = peptide2.rfind('(').ok_or_else(||
            BoxedError::new(BasicKind::Error,
                title,
                format!("A pLink {field} should follow the format 'PEP1(pos1)-PEP2(pos2)' but the opening bracket '(' was not found for PEP2"),
                Context::line(Some(location.line.line_index() as u32), location.full_line(), location.location.start+peptide1.len()+2, peptide2.len())))?;

        let pos1 =
            location.location.start + first_end + 1..location.location.start + peptide1.len();
        let first_index = location.full_line()[pos1.clone()].parse::<u16>().map_err(|err|
            BoxedError::new(BasicKind::Error,
                title,
                format!("A pLink {field} should follow the format 'PEP1(pos1)-PEP2(pos2)' but the position for PEP1 {}", explain_number_error(&err)),
                Context::line_range(Some(location.line.line_index() as u32), location.full_line(), pos1.clone()).to_owned()))?;
        let pos2 = location.location.start + peptide1.len() + 2 + second_end + 1
            ..location.location.start + peptide1.len() + 2 + peptide2.len() - 1;
        let second_index = location.full_line()[pos2.clone()].parse::<u16>().map_err(|err|
            BoxedError::new(BasicKind::Error,
                title,
                format!("A pLink {field} should follow the format 'PEP1(pos1)-PEP2(pos2)' but the position for PEP1 {}", explain_number_error(&err)),
                Context::line_range(Some(location.line.line_index() as u32), location.full_line(), pos2.clone()).to_owned()))?;

        Ok((
            location.location.start..location.location.start + first_end,
            Some((first_index, pos1)),
            Some(
                location.location.start + peptide1.len() + 2
                    ..location.location.start + peptide1.len() + 2 + second_end,
            ),
            Some((second_index, pos2)),
        ))
    } else {
        // rsplit to prevent picking a bracket in the text field, and then reverse for it to make sense to human brains
        let mut split = location.as_str().rsplitn(3, '(').collect_vec();
        split.reverse();

        match split.len() {
            3 => {
                let start = location.location.start;
                let start_pos1 = start + split[0].len() + 1;
                let start_pos2 = start_pos1 + split[1].len() + 1;
                let end = location.location.end;

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
                let start = location.location.start;
                let start_pos1 = start + split[0].len() + 1;
                let end = location.location.end;

                let pos1 = start_pos1..end - 1;
                let first_index = location.full_line()[pos1.clone()].parse::<u16>().map_err(|err|
                    BoxedError::new(BasicKind::Error,
                        title,
                        format!("A pLink {field} should follow the format 'PEP(pos1)(pos2)' but the first position {}", explain_number_error(&err)),
                        Context::line_range(Some(location.line.line_index() as u32), location.full_line(), pos1.clone()).to_owned()))?;

                Ok((start..start_pos1 - 1, Some((first_index, pos1)), None, None))
            }
            1 => Ok((location.location.clone(), None, None, None)),
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
    ptm: "modifications",
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

impl IdentifiedPeptidoformVersion<PLinkFormat> for PLinkVersion {
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

impl MetaData for PLinkData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptidoform.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::PLink(self.version)
    }

    fn id(&self) -> String {
        self.order.to_string()
    }

    fn confidence(&self) -> Option<f64> {
        Some(1.0 - self.score)
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<f64> {
        Some(self.score)
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        None
    }

    fn charge(&self) -> Option<Charge> {
        Some(self.z)
    }

    fn mode(&self) -> Option<&str> {
        None
    }

    fn retention_time(&self) -> Option<Time> {
        None
    }

    fn scans(&self) -> SpectrumIds {
        self.scan_number.map_or_else(
            || SpectrumIds::FileNotKnown(vec![SpectrumId::Native(self.title.clone())]),
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

    fn protein_names(&self) -> Option<Cow<'_, [FastaIdentifier<String>]>> {
        None
    }

    fn protein_id(&self) -> Option<usize> {
        None
    }

    fn protein_location(&self) -> Option<Range<u16>> {
        None
    }

    fn flanking_sequences(&self) -> (&FlankingSequence, &FlankingSequence) {
        (&FlankingSequence::Unknown, &FlankingSequence::Unknown)
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        None
    }
}
