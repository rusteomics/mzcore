use std::borrow::Cow;

use crate::{CVTerm, EMPTY_FASTA_IDENTIFIER, FastaIdentifier, Reliability};
use mzcore::{
    prelude::SequencePosition,
    sequence::{Linear, Peptidoform, SimpleModification},
};
use mzcv::Curie;

/// Generalised access to meta data of identified peptidoforms
pub trait ProteinMetaData {
    /// Get the full sequence, if present
    fn sequence(&self) -> Option<Cow<'_, Peptidoform<Linear>>>;

    /// Get the numerical protein identifier
    fn numerical_id(&self) -> Option<usize>;

    /// Get the identifier
    fn id(&self) -> FastaIdentifier<&str>;

    /// Get the description
    fn description(&self) -> Option<&str>;

    /// Get the species as a curie term
    fn species(&self) -> Option<Curie>;

    /// Get the species as a human readable name
    fn species_name(&self) -> Option<&str>;

    /// Get the search engine that identified this PSM its score and the term to describe the score
    fn search_engine(&self) -> &[(CVTerm, Option<(f64, CVTerm)>)];

    /// A list of all proteins (identifiers) that cannot be separated based on peptide evidence from this main protein
    fn ambiguity_members(&self) -> &[String];

    /// The database that was used for matching optionally with the version of the database
    fn database(&self) -> Option<(&str, Option<&str>)>;

    /// All known modifications that can occur at PSM level on this sequence, with potentially a probability
    fn modifications(&self) -> &[(Vec<(SequencePosition, Option<f64>)>, SimpleModification)];

    /// The coverage of the protein based on the mapped PSMs
    fn coverage(&self) -> Option<f64>;

    /// The gene ontology terms for this protein
    fn gene_ontology(&self) -> &[Curie];

    /// Get the reliability of this PSM
    fn reliability(&self) -> Option<Reliability>;

    /// Get the URI for this PSM
    fn uri(&self) -> Option<&str>;
}

macro_rules! impl_ref {
    ($t:ty) => {
        impl<T: ProteinMetaData> ProteinMetaData for $t {
            fn sequence(&self) -> Option<Cow<'_, Peptidoform<Linear>>> {
                (**self).sequence()
            }

            fn numerical_id(&self) -> Option<usize> {
                (**self).numerical_id()
            }

            fn id(&self) -> FastaIdentifier<&str> {
                (**self).id()
            }

            fn description(&self) -> Option<&str> {
                (**self).description()
            }

            fn species(&self) -> Option<Curie> {
                (**self).species()
            }

            fn species_name(&self) -> Option<&str> {
                (**self).species_name()
            }

            fn search_engine(&self) -> &[(CVTerm, Option<(f64, CVTerm)>)] {
                (**self).search_engine()
            }

            fn ambiguity_members(&self) -> &[String] {
                (**self).ambiguity_members()
            }

            fn database(&self) -> Option<(&str, Option<&str>)> {
                (**self).database()
            }

            fn modifications(
                &self,
            ) -> &[(Vec<(SequencePosition, Option<f64>)>, SimpleModification)] {
                (**self).modifications()
            }

            fn coverage(&self) -> Option<f64> {
                (**self).coverage()
            }

            fn gene_ontology(&self) -> &[Curie] {
                (**self).gene_ontology()
            }

            fn reliability(&self) -> Option<Reliability> {
                (**self).reliability()
            }

            fn uri(&self) -> Option<&str> {
                (**self).uri()
            }
        }
    };
}

impl_ref!(&T);
impl_ref!(std::rc::Rc<T>);
impl_ref!(std::sync::Arc<T>);

impl ProteinMetaData for Box<dyn ProteinMetaData + '_> {
    fn sequence(&self) -> Option<Cow<'_, Peptidoform<Linear>>> {
        (**self).sequence()
    }

    fn numerical_id(&self) -> Option<usize> {
        (**self).numerical_id()
    }

    fn id(&self) -> FastaIdentifier<&str> {
        (**self).id()
    }

    fn description(&self) -> Option<&str> {
        (**self).description()
    }

    fn species(&self) -> Option<Curie> {
        (**self).species()
    }

    fn species_name(&self) -> Option<&str> {
        (**self).species_name()
    }

    fn search_engine(&self) -> &[(CVTerm, Option<(f64, CVTerm)>)] {
        (**self).search_engine()
    }

    fn ambiguity_members(&self) -> &[String] {
        (**self).ambiguity_members()
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        (**self).database()
    }

    fn modifications(&self) -> &[(Vec<(SequencePosition, Option<f64>)>, SimpleModification)] {
        (**self).modifications()
    }

    fn coverage(&self) -> Option<f64> {
        (**self).coverage()
    }

    fn gene_ontology(&self) -> &[Curie] {
        (**self).gene_ontology()
    }

    fn reliability(&self) -> Option<Reliability> {
        (**self).reliability()
    }

    fn uri(&self) -> Option<&str> {
        (**self).uri()
    }
}

/// A basic type to use a placeholder when a PSM format does not contain any protein level information.
#[derive(Clone, Copy, Debug, Default, serde::Deserialize, serde::Serialize)]
pub struct NoProtein {}

impl mzcore::space::Space for NoProtein {
    fn space(&self) -> mzcore::space::UsedSpace {
        mzcore::space::UsedSpace::default()
    }
}

impl ProteinMetaData for NoProtein {
    fn sequence(&self) -> Option<Cow<'_, Peptidoform<Linear>>> {
        None
    }

    fn numerical_id(&self) -> Option<usize> {
        None
    }

    fn id(&self) -> FastaIdentifier<&str> {
        EMPTY_FASTA_IDENTIFIER
    }

    fn description(&self) -> Option<&str> {
        None
    }

    fn species(&self) -> Option<Curie> {
        None
    }

    fn species_name(&self) -> Option<&str> {
        None
    }

    fn search_engine(&self) -> &[(CVTerm, Option<(f64, CVTerm)>)] {
        &[]
    }

    fn ambiguity_members(&self) -> &[String] {
        &[]
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        None
    }

    fn modifications(&self) -> &[(Vec<(SequencePosition, Option<f64>)>, SimpleModification)] {
        &[]
    }

    fn coverage(&self) -> Option<f64> {
        None
    }

    fn gene_ontology(&self) -> &[Curie] {
        &[]
    }

    fn reliability(&self) -> Option<Reliability> {
        None
    }

    fn uri(&self) -> Option<&str> {
        None
    }
}

impl ProteinMetaData for FastaIdentifier<String> {
    fn sequence(&self) -> Option<Cow<'_, Peptidoform<Linear>>> {
        None
    }

    fn numerical_id(&self) -> Option<usize> {
        None
    }

    fn id(&self) -> FastaIdentifier<&str> {
        self.as_str()
    }

    fn description(&self) -> Option<&str> {
        None
    }

    fn species(&self) -> Option<Curie> {
        None
    }

    fn species_name(&self) -> Option<&str> {
        None
    }

    fn search_engine(&self) -> &[(CVTerm, Option<(f64, CVTerm)>)] {
        &[]
    }

    fn ambiguity_members(&self) -> &[String] {
        &[]
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        None
    }

    fn modifications(&self) -> &[(Vec<(SequencePosition, Option<f64>)>, SimpleModification)] {
        &[]
    }

    fn coverage(&self) -> Option<f64> {
        None
    }

    fn gene_ontology(&self) -> &[Curie] {
        &[]
    }

    fn reliability(&self) -> Option<Reliability> {
        None
    }

    fn uri(&self) -> Option<&str> {
        None
    }
}

impl ProteinMetaData for FastaIdentifier<Box<str>> {
    fn sequence(&self) -> Option<Cow<'_, Peptidoform<Linear>>> {
        None
    }

    fn numerical_id(&self) -> Option<usize> {
        None
    }

    fn id(&self) -> FastaIdentifier<&str> {
        self.as_str()
    }

    fn description(&self) -> Option<&str> {
        None
    }

    fn species(&self) -> Option<Curie> {
        None
    }

    fn species_name(&self) -> Option<&str> {
        None
    }

    fn search_engine(&self) -> &[(CVTerm, Option<(f64, CVTerm)>)] {
        &[]
    }

    fn ambiguity_members(&self) -> &[String] {
        &[]
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        None
    }

    fn modifications(&self) -> &[(Vec<(SequencePosition, Option<f64>)>, SimpleModification)] {
        &[]
    }

    fn coverage(&self) -> Option<f64> {
        None
    }

    fn gene_ontology(&self) -> &[Curie] {
        &[]
    }

    fn reliability(&self) -> Option<Reliability> {
        None
    }

    fn uri(&self) -> Option<&str> {
        None
    }
}

impl ProteinMetaData for FastaIdentifier<&str> {
    fn sequence(&self) -> Option<Cow<'_, Peptidoform<Linear>>> {
        None
    }

    fn numerical_id(&self) -> Option<usize> {
        None
    }

    fn id(&self) -> FastaIdentifier<&str> {
        self.clone()
    }

    fn description(&self) -> Option<&str> {
        None
    }

    fn species(&self) -> Option<Curie> {
        None
    }

    fn species_name(&self) -> Option<&str> {
        None
    }

    fn search_engine(&self) -> &[(CVTerm, Option<(f64, CVTerm)>)] {
        &[]
    }

    fn ambiguity_members(&self) -> &[String] {
        &[]
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        None
    }

    fn modifications(&self) -> &[(Vec<(SequencePosition, Option<f64>)>, SimpleModification)] {
        &[]
    }

    fn coverage(&self) -> Option<f64> {
        None
    }

    fn gene_ontology(&self) -> &[Curie] {
        &[]
    }

    fn reliability(&self) -> Option<Reliability> {
        None
    }

    fn uri(&self) -> Option<&str> {
        None
    }
}
