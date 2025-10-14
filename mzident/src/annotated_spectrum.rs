use std::{
    borrow::Cow,
    marker::PhantomData,
    ops::{Not, Range},
    path::Path,
};

use crate::{
    FastaIdentifier, FlankingSequence, GeneralIdentifiedPeptidoforms, IdentifiedPeptidoform,
    IdentifiedPeptidoformData, KnownFileFormat, MaybePeptidoform, MetaData, SpectrumId,
    SpectrumIds,
};
use mzannotate::{mzspeclib::MzSpecLibTextParseError, prelude::AnnotatedSpectrum};
use mzcore::{
    ontology::CustomDatabase,
    prelude::CompoundPeptidoformIon,
    sequence::Linked,
    system::isize::Charge,
    system::{Mass, MassOverCharge, Time},
};

use context_error::{BasicKind, BoxedError, Context, CreateError};

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
        None
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<f64> {
        None
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
            .and_then(|p| p.get_charge_carriers().map(|c| c.charge()))
    }

    fn mode(&self) -> Option<Cow<'_, str>> {
        self.description
            .precursor
            .first()
            .map(|p| Cow::Owned(p.activation.methods().iter().map(|d| d.name()).join("+")))
    }

    fn retention_time(&self) -> Option<Time> {
        Some(Time::new::<mzcore::system::time::s>(
            self.description.acquisition.start_time(),
        ))
    }

    fn scans(&self) -> SpectrumIds {
        SpectrumIds::FileNotKnown(vec![SpectrumId::Index(self.description.index)])
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        self.description.precursor.first().and_then(|p| {
            p.ions
                .first()
                .map(|i| MassOverCharge::new::<mzcore::system::mass_over_charge::thomson>(i.mz))
        })
    }

    fn experimental_mass(&self) -> Option<Mass> {
        self.description.precursor.first().and_then(|p| {
            p.ions.first().and_then(|i| {
                i.charge
                    .map(|c| Mass::new::<mzcore::system::mass::dalton>(i.mz * c as f64))
            })
        })
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
    _custom_database: Option<&'a CustomDatabase>,
) -> Result<GeneralIdentifiedPeptidoforms<'a>, BoxedError<'static, BasicKind>> {
    let path_string = path.to_string_lossy().to_string();
    let path_string_2 = path_string.clone();
    mzannotate::mzspeclib::MzSpecLibParser::new(std::io::BufReader::new(
        std::fs::File::open(path).map_err(|e| {
            BoxedError::new(
                BasicKind::Error,
                "Could not open mzSpecLib file",
                e.to_string(),
                Context::none().source(path_string.clone()),
            )
        })?,
    ))
    .map(move |parser| {
        let b: Box<
            dyn Iterator<
                Item = Result<
                    IdentifiedPeptidoform<Linked, MaybePeptidoform>,
                    BoxedError<'static, BasicKind>,
                >,
            >,
        > = Box::new(parser.map(move |s| {
            s.map(Into::into).map_err(|e| match e {
                MzSpecLibTextParseError::InvalidMzPAF(e)
                | MzSpecLibTextParseError::InvalidProforma(e)
                | MzSpecLibTextParseError::RichError(e) => e,
                _ => BoxedError::new(
                    BasicKind::Error,
                    "Could not parse mzSpecLib spectrum",
                    format!("{e:?}"),
                    Context::none().source(path_string.clone()),
                ),
            })
        }));
        b
    })
    .map_err(|e| {
        BoxedError::new(
            BasicKind::Error,
            "Could not parse mzSpecLib file",
            format!("{e:?}"),
            Context::none().source(path_string_2),
        )
    })
}
