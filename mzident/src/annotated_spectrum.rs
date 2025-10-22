use std::{
    borrow::Cow,
    marker::PhantomData,
    ops::{Not, Range},
    path::Path,
};

use crate::{
    FastaIdentifier, GeneralIdentifiedPeptidoforms, IdentifiedPeptidoform,
    IdentifiedPeptidoformData, KnownFileFormat, MaybePeptidoform, MetaData, SpectrumId,
    SpectrumIds,
};
use mzannotate::prelude::AnnotatedSpectrum;
use mzcore::{
    ontology::CustomDatabase,
    prelude::*,
    sequence::{FlankingSequence, Linked},
    system::isize::Charge,
    system::{Mass, MassOverCharge, Time},
};

use context_error::{BasicKind, BoxedError, Context, CreateError, FullErrorContent};

use itertools::Itertools;

impl MetaData for AnnotatedSpectrum {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        let cpi: CompoundPeptidoformIon = self
            .analytes
            .iter()
            .filter_map(|a| match &a.target {
                mzannotate::mzspeclib::AnalyteTarget::PeptidoformIon(pep) => Some(pep.clone()),
                _ => None,
            })
            .collect();
        cpi.peptidoform_ions()
            .is_empty()
            .not()
            .then_some(Cow::Owned(cpi))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::AnnotatedSpectrum
    }

    fn id(&self) -> String {
        self.description.index.to_string()
    }

    fn confidence(&self) -> Option<f64> {
        self.interpretations
            .iter()
            .filter_map(|i| i.probability)
            .exactly_one()
            .ok()
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<f64> {
        self.interpretations
            .iter()
            .filter_map(|i| i.probability)
            .exactly_one()
            .ok()
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        None
    }

    fn charge(&self) -> Option<Charge> {
        self.analytes
            .iter()
            .filter_map(|a| match &a.target {
                mzannotate::mzspeclib::AnalyteTarget::PeptidoformIon(pep) => Some(pep),
                _ => None,
            })
            .exactly_one()
            .ok()
            .and_then(|p| p.get_charge_carriers().map(MolecularCharge::charge))
    }

    fn mode(&self) -> Option<Cow<'_, str>> {
        self.description
            .precursor
            .first()
            .map(|p| Cow::Owned(p.activation.methods().iter().map(|d| d.name()).join("+")))
    }

    fn retention_time(&self) -> Option<Time> {
        self.description
            .acquisition
            .scans
            .first()
            .map(|s| Time::new::<mzcore::system::time::s>(s.start_time))
    }

    fn scans(&self) -> SpectrumIds {
        self.description
            .params
            .iter()
            .find(|p| {
                p.controlled_vocabulary
                    .is_some_and(|cv| cv.prefix() == "MS")
                    && p.accession == Some(1003203)
            })
            .map_or_else(
                || SpectrumIds::FileNotKnown(vec![SpectrumId::Index(self.description.index)]),
                |rawfile| {
                    SpectrumIds::FileKnown(vec![(
                        rawfile.value.to_string().into(),
                        vec![SpectrumId::Index(self.description.index)],
                    )])
                },
            )
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        self.description
            .precursor
            .first()
            .and_then(|p| p.ions.first())
            .map(|i| MassOverCharge::new::<mzcore::system::mass_over_charge::thomson>(i.mz))
    }

    fn experimental_mass(&self) -> Option<Mass> {
        self.description
            .precursor
            .first()
            .and_then(|p| p.ions.first())
            .and_then(|i| {
                i.charge
                    .map(|c| Mass::new::<mzcore::system::mass::dalton>(i.mz * f64::from(c)))
            })
    }

    fn protein_names(&self) -> Option<Cow<'_, [FastaIdentifier<String>]>> {
        Some(Cow::Owned(
            self.analytes
                .iter()
                .flat_map(|a| &a.proteins)
                .filter_map(|p| p.accession.clone())
                .map(|a| {
                    a.parse()
                        .unwrap_or_else(|_| FastaIdentifier::Undefined(false, a.to_string()))
                })
                .collect(),
        ))
    }

    fn protein_id(&self) -> Option<usize> {
        None
    }

    fn protein_location(&self) -> Option<Range<u16>> {
        None
    }

    fn flanking_sequences(&self) -> (&FlankingSequence, &FlankingSequence) {
        (
            self.analytes
                .iter()
                .flat_map(|a| &a.proteins)
                .map(|p| &p.flanking_sequences.0)
                .exactly_one()
                .unwrap_or(FlankingSequence::UNKNOWN),
            self.analytes
                .iter()
                .flat_map(|a| &a.proteins)
                .map(|p| &p.flanking_sequences.1)
                .exactly_one()
                .unwrap_or(FlankingSequence::UNKNOWN),
        )
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        self.analytes
            .iter()
            .flat_map(|a| &a.proteins)
            .map(|p| p.database_name.as_deref())
            .exactly_one()
            .ok()
            .flatten()
            .map(|name| {
                (
                    name,
                    self.analytes
                        .iter()
                        .flat_map(|a| &a.proteins)
                        .map(|p| p.database_version.as_deref())
                        .exactly_one()
                        .ok()
                        .flatten(),
                )
            })
    }

    fn annotated_spectrum(&self) -> Option<Cow<'_, AnnotatedSpectrum>> {
        Some(Cow::Borrowed(self))
    }

    fn has_annotated_spectrum(&self) -> bool {
        true
    }
}

impl From<AnnotatedSpectrum> for IdentifiedPeptidoform<Linked, MaybePeptidoform> {
    fn from(value: AnnotatedSpectrum) -> Self {
        Self {
            score: value.confidence(),
            local_confidence: value.local_confidence().map(|v| v.to_vec()),
            data: IdentifiedPeptidoformData::AnnotatedSpectrum(value),
            complexity_marker: PhantomData,
            peptidoform_availability_marker: PhantomData,
        }
    }
}

// TODO: this should be done in the mzannotate crate
pub(crate) fn parse_mzspeclib<'a>(
    path: &Path,
    custom_database: Option<&'a CustomDatabase>,
) -> Result<GeneralIdentifiedPeptidoforms<'a>, BoxedError<'static, BasicKind>> {
    mzannotate::mzspeclib::MzSpecLibTextParser::open_file(path, custom_database)
        .map(move |parser| {
            let b: Box<
                dyn Iterator<
                    Item = Result<
                        IdentifiedPeptidoform<Linked, MaybePeptidoform>,
                        BoxedError<'static, BasicKind>,
                    >,
                >,
            > = Box::new(parser.map(move |s| {
                s.map(Into::into)
                    .map_err(|e| e.convert(|_| BasicKind::Error))
            }));
            b
        })
        .map_err(|e| e.convert(|_| BasicKind::Error))
}
