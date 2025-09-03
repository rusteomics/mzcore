//! Python bindings to the rustyms library.

use std::fmt::Debug;
use std::num::NonZeroU16;

use ordered_float::OrderedFloat;
use pyo3::{exceptions::PyValueError, prelude::*, types::PyType};

use rustyms::{
    annotation::AnnotatableSpectrum,
    chemistry::{Chemical, MultiChemical},
    sequence::{IsAminoAcid, Linked},
    system::dalton,
};

/// Function to find all sequences that match a specific mass. The returned sequences are returned in a deterministic order.
///
/// Parameters
/// ----------
/// mass : number
///     The monoisotopic target mass
/// tolerance : number
///     The tolerance in ppm (parts per million)
/// amino_acids : list[AminoAcid]
///     The list of amino acids to use look into TODO to get a default list
/// fixed : list[tuple[SimpleModification, list[AminoAcid]]]
///     The list of fixed modifications, with a list of all allowed positions, if this list is empty the allowed places from the modification are used
/// variable : list[tuple[SimpleModification, list[AminoAcid]]]
///     The list of variable modifications, with a list of all allowed positions, if this list is empty the allowed places from the modification are used
/// base : Peptidoform | None
///     If present, the base sequence that is always assumed to be present, if this has multiple molecular formulas the lowest mass one is used
///
/// Returns
/// -------
/// list[Peptidoform]
///
#[pyfunction]
fn find_isobaric_sets(
    mass: f64,
    tolerance: f64,
    amino_acids: Vec<AminoAcid>,
    fixed: Vec<(SimpleModification, Vec<AminoAcid>)>,
    variable: Vec<(SimpleModification, Vec<AminoAcid>)>,
    base: Option<Peptidoform>,
) -> Vec<Peptidoform> {
    let amino_acids = amino_acids.into_iter().map(|a| a.0).collect::<Vec<_>>();
    let base = base.and_then(|b| b.0.into_simple_linear());
    let fixed = fixed
        .into_iter()
        .map(|(m, r)| {
            (
                m.0,
                (!r.is_empty()).then(|| {
                    rustyms::sequence::PlacementRule::AminoAcid(
                        r.into_iter().map(|a| a.0).collect(),
                        rustyms::sequence::Position::Anywhere,
                    )
                }),
            )
        })
        .collect::<Vec<_>>();
    let variable = variable
        .into_iter()
        .map(|(m, r)| {
            (
                m.0,
                (!r.is_empty()).then(|| {
                    rustyms::sequence::PlacementRule::AminoAcid(
                        r.into_iter().map(|a| a.0).collect(),
                        rustyms::sequence::Position::Anywhere,
                    )
                }),
            )
        })
        .collect::<Vec<_>>();

    rustyms::prelude::find_isobaric_sets(
        rustyms::system::Mass::new::<dalton>(mass),
        rustyms::quantities::Tolerance::new_ppm(tolerance),
        &amino_acids,
        &fixed,
        &variable,
        base.as_ref(),
    )
    .map(|p| Peptidoform(p.into()))
    .collect()
}

/// Mass mode enum.
#[pyclass(eq, eq_int)]
#[derive(Eq, PartialEq)]
enum MassMode {
    Monoisotopic,
    Average,
    MostAbundant,
}

/// Element.
///
/// A chemical element, with its isotopes and their properties.
///
/// Parameters
/// ----------
/// symbol : str
///
#[pyclass]
#[derive(Clone, Copy, Debug)]
pub struct Element(rustyms::chemistry::Element);

#[pymethods]
impl Element {
    #[new]
    fn new(symbol: &str) -> PyResult<Self> {
        rustyms::chemistry::Element::try_from(symbol)
            .map(Element)
            .map_err(|_| PyValueError::new_err("Invalid element symbol."))
    }

    fn __repr__(&self) -> String {
        format!("Element('{}')", self.0)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    /// Get all available isotopes (N, mass, abundance).
    ///
    /// Returns
    /// -------
    /// list[tuple[int, float, float]]
    ///
    fn isotopes(&self) -> Vec<(u16, f64, f64)> {
        self.0
            .isotopes()
            .iter()
            .map(|i| (i.0, i.1.value, i.2))
            .collect()
    }

    /// The mass of the specified isotope of this element (if that isotope exists).
    ///
    /// Parameters
    /// ----------
    /// isotope : int | None
    ///    The isotope number (default: None).
    ///
    /// Returns
    /// -------
    /// float | None
    ///
    #[pyo3(signature = (isotope=None))]
    fn mass(&self, isotope: Option<u16>) -> Option<f64> {
        self.0
            .mass(isotope.and_then(NonZeroU16::new))
            .map(|mass| mass.value)
    }

    /// The average weight of the specified isotope of this element (if that isotope exists).
    ///
    /// Parameters
    /// ----------
    /// isotope : int | None
    ///     The isotope number (default: None).
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[pyo3(signature = (isotope=None))]
    fn average_weight(&self, isotope: Option<u16>) -> Option<f64> {
        self.0
            .average_weight(isotope.and_then(NonZeroU16::new))
            .map(|mass| mass.value)
    }
}

impl std::fmt::Display for Element {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

/// Molecular formula.
///
/// A molecular formula: a selection of elements of specified isotopes together forming a structure.
///
#[pyclass]
#[derive(Clone, Debug)]
pub struct MolecularFormula(rustyms::chemistry::MolecularFormula);

#[pymethods]
impl MolecularFormula {
    /// Create a new molecular formula.
    ///
    /// Parameters
    /// ----------
    /// elements : list[tuple[Element, int | None, int]]
    ///     List of tuples of elements, isotope numbers and counts.
    ///
    /// Returns
    /// -------
    /// MolecularFormula
    ///
    /// Raises
    /// ------
    /// ValueError
    ///    If an isotope is invalid.
    ///
    #[new]
    fn new(elements: Vec<(Element, Option<u16>, i32)>) -> PyResult<Self> {
        let elements = elements
            .iter()
            .map(|(e, i, n)| (e.0, (*i).and_then(NonZeroU16::new), *n))
            .collect::<Vec<_>>();
        let formula = rustyms::chemistry::MolecularFormula::new(elements.as_slice(), &[])
            .ok_or_else(|| PyValueError::new_err("Invalid isotopes"))?;
        Ok(Self(formula))
    }

    fn __repr__(&self) -> String {
        format!("MolecularFormula({})", self.0)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    /// Create a new molecular formula from a ProForma formula notation string.
    ///
    /// Parameters
    /// ----------
    /// proforma : str
    ///
    /// Returns
    /// -------
    /// MolecularFormula
    ///
    #[classmethod]
    fn from_pro_forma(_cls: &Bound<'_, PyType>, proforma: &str) -> PyResult<Self> {
        rustyms::chemistry::MolecularFormula::from_pro_forma(proforma, .., false, false, true, true)
            .map(MolecularFormula)
            .map_err(|e| PyValueError::new_err(format!("Invalid ProForma string: {e}")))
    }

    /// Create a new molecular formula from a PSI-MOD formula notation string.
    ///
    /// Parameters
    /// ----------
    /// psi_mod : str
    ///
    /// Returns
    /// -------
    /// MolecularFormula
    ///
    #[classmethod]
    fn from_psi_mod(_cls: &Bound<'_, PyType>, psi_mod: &str) -> PyResult<Self> {
        rustyms::chemistry::MolecularFormula::from_psi_mod(psi_mod, ..)
            .map(MolecularFormula)
            .map_err(|e| PyValueError::new_err(format!("Invalid PSI-MOD string: {e}")))
    }

    /// Add the given element to this formula (while keeping it ordered and simplified)
    ///
    /// Parameters
    /// ----------
    /// element : Element
    ///     The element to add.
    /// isotope : int | None
    ///     The isotope number of the element to add (default: None).
    /// n : int
    ///     The number of atoms of this element to add (default: 1).
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the element or isotope is invalid.
    ///
    #[pyo3(signature = (element, isotope=None, n=1))]
    fn add(&mut self, element: &Element, isotope: Option<u16>, n: i32) -> PyResult<Option<()>> {
        if self
            .0
            .add((element.0, isotope.and_then(NonZeroU16::new), n))
        {
            Ok(None)
        } else {
            Err(PyValueError::new_err("Invalid element or isotope"))
        }
    }

    /// Get the elements making this formula.
    ///
    /// Returns
    /// -------
    /// list[tuple[Element, int | None, int]]
    ///
    fn elements(&self) -> Vec<(Element, Option<u16>, i32)> {
        self.0
            .elements()
            .iter()
            .map(|(e, i, n)| (Element(*e), (*i).map(NonZeroU16::get), *n))
            .collect()
    }

    /// Create a new molecular formula with the given global isotope modifications.
    fn with_global_isotope_modifications(
        &self,
        substitutions: Vec<(Element, Option<u16>)>,
    ) -> PyResult<Self> {
        self.0
            .with_global_isotope_modifications(
                substitutions
                    .iter()
                    .map(|(e, i)| (e.0, (*i).and_then(NonZeroU16::new)))
                    .collect::<Vec<_>>()
                    .as_slice(),
            )
            .map_or_else(
                || {
                    Err(PyValueError::new_err(
                        "Invalid global isotope modifications",
                    ))
                },
                |formula| Ok(Self(formula)),
            )
    }

    /// Get the number of electrons (the only charged species, any ionic species is saved as that element +/- the correct number of electrons). The inverse of that number is given as the charge.
    ///
    /// Returns
    /// -------
    /// int
    ///
    fn charge(&self) -> isize {
        self.0.charge().value
    }

    /// The mass of the molecular formula of this element, if all element species (isotopes) exists
    ///
    /// Returns
    /// -------
    /// float
    ///
    fn monoisotopic_mass(&self) -> f64 {
        self.0.monoisotopic_mass().value
    }

    /// The average weight of the molecular formula of this element, if all element species (isotopes) exists.
    ///
    /// Returns
    /// -------
    /// float
    ///
    fn average_weight(&self) -> f64 {
        self.0.average_weight().value
    }

    /// The most abundant mass. This is the isotopic species with the highest abundance when the whole isotope
    /// distribution is generated. Because this uses an averagine model it is not very precise in its mass.
    /// Because it has to generate the full isotope distribution this takes more time then other mass modes.
    ///
    /// Returns
    /// -------
    /// float
    ///
    fn most_abundant_mass(&self) -> f64 {
        self.0.most_abundant_mass().value
    }

    /// Get the mass in the given mode.
    ///
    /// Parameters
    /// ----------
    /// mode : MassMode
    ///    The mode to get the mass in.
    ///
    /// Returns
    /// -------
    /// float
    ///
    /// Raises
    /// ------
    /// ValueError
    ///   If the mode is not one of the valid modes.
    ///
    #[pyo3(signature = (mode=&MassMode::Monoisotopic))]
    fn mass(&self, mode: &MassMode) -> PyResult<f64> {
        match mode {
            MassMode::Monoisotopic => Ok(self.monoisotopic_mass()),
            MassMode::Average => Ok(self.average_weight()),
            MassMode::MostAbundant => Ok(self.most_abundant_mass()),
        }
    }

    /// Create a Hill notation from this collections of elements merged with the ProForma notation for specific isotopes.
    ///
    /// Returns
    /// -------
    /// str
    fn hill_notation(&self) -> String {
        self.0.hill_notation()
    }

    /// Create a Hill notation from this collections of elements merged with the ProForma notation for specific isotopes. Using fancy unicode characters for subscript and superscript numbers.
    ///
    /// Returns
    /// -------
    /// str
    fn hill_notation_fancy(&self) -> String {
        self.0.hill_notation_fancy()
    }

    /// Create a Hill notation from this collections of elements encoded in HTML.
    ///
    /// Returns
    /// -------
    /// str
    fn hill_notation_html(&self) -> String {
        self.0.hill_notation_html()
    }

    /// Get information on all sources of ambiguity and multiplicity in formulas.
    ///
    /// Returns
    /// -------
    /// list[str]
    fn ambiguous_labels(&self) -> Vec<String> {
        self.0.labels().iter().map(ToString::to_string).collect()
    }
}

/// A selection of ions that together define the charge of a peptidoform.
#[pyclass]
#[derive(Debug)]
pub struct MolecularCharge(rustyms::chemistry::MolecularCharge);

#[pymethods]
impl MolecularCharge {
    /// Create a charge state with the given ions.
    ///
    /// Parameters
    /// ----------
    /// charge_carriers : list[tuple[int, MolecularFormula]]
    ///   The charge carriers.
    ///
    /// Returns
    /// -------
    /// MolecularCharge
    ///
    #[new]
    fn new(charge_carriers: Vec<(i32, MolecularFormula)>) -> Self {
        Self(rustyms::chemistry::MolecularCharge {
            charge_carriers: charge_carriers
                .iter()
                .map(|(n, mol)| (*n as isize, mol.0.clone()))
                .collect(),
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "MolecularCharge(charge_carriers={})",
            self.0
                .charge_carriers
                .iter()
                .map(|(n, mol)| format!("({n}, {mol})"))
                .collect::<Vec<_>>()
                .join(", ")
        )
    }

    // Create a default charge state with only protons.
    ///
    /// Parameters
    /// ----------
    /// charge : int
    ///    The charge.
    ///
    /// Returns
    /// -------
    /// MolecularCharge
    ///
    #[classmethod]
    fn proton(_cls: &Bound<'_, PyType>, charge: i32) -> Self {
        Self(rustyms::chemistry::MolecularCharge::proton(charge as isize))
    }

    /// List of counts and molecular formulas for the charge carriers.
    ///
    /// Returns
    /// -------
    /// list[tuple[int, MolecularFormula]]
    ///     The charge carriers.
    #[getter]
    fn charge_carriers(&self) -> Vec<(i32, MolecularFormula)> {
        self.0
            .charge_carriers
            .iter()
            .map(|(n, mol)| (*n as i32, MolecularFormula(mol.clone())))
            .collect()
    }
}

/// Amino acid.
///
/// Parameters
/// ----------
/// name : str
///    The name of the amino acid.
///
#[pyclass]
#[derive(Clone, Copy, Debug)]
pub struct AminoAcid(rustyms::sequence::AminoAcid);

#[pymethods]
impl AminoAcid {
    #[new]
    fn new(name: &str) -> PyResult<Self> {
        rustyms::sequence::AminoAcid::try_from(name).map_or_else(
            |_| Err(PyValueError::new_err("Invalid amino acid")),
            |aa| Ok(Self(aa)),
        )
    }

    fn __str__(&self) -> String {
        self.0.pro_forma_definition().to_string()
    }

    fn __repr__(&self) -> String {
        self.to_string()
    }

    /// Molecular formula(s) of the amino acid.
    ///
    /// Returns a list of molecular formulas that are possible for the amino acid symbol.
    ///
    /// Returns
    /// -------
    /// List[MolecularFormula]
    ///
    fn formulas(&self) -> Vec<MolecularFormula> {
        self.0
            .formulas()
            .iter()
            .map(|f| MolecularFormula(f.clone()))
            .collect()
    }

    /// Molecular formula of the amino acid.
    ///
    /// Returns the molecular formula of the amino acid (the first if multiple are possible).
    ///
    /// Returns
    /// -------
    /// List[MolecularFormula]
    ///
    fn formula(&self) -> MolecularFormula {
        MolecularFormula(self.0.formulas().first().unwrap().clone())
    }

    /// Monoisotopic mass(es) of the amino acid.
    ///
    /// Returns a list of monoisotopic masses that are possible for the amino acid symbol.
    ///
    /// Returns
    /// -------
    /// List[float]
    ///
    fn monoisotopic_masses(&self) -> Vec<f64> {
        self.0
            .formulas()
            .iter()
            .map(|f| f.monoisotopic_mass().value)
            .collect()
    }

    /// Monoisotopic mass of the amino acid.
    ///
    /// Returns the monoisotopic mass of the amino acid (the first if multiple are possible).
    ///
    /// Returns
    /// -------
    /// float
    ///
    fn monoisotopic_mass(&self) -> f64 {
        self.0.formulas().first().unwrap().monoisotopic_mass().value
    }
}

impl std::fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{self:?}",)
    }
}

/// Simple amino acid modification.
///
/// Parameters
/// ----------
/// name : str
///   The name of the modification. Any simple modification as allowed in ProForma (no ambiguous or cross-linked modifications).
///
#[pyclass]
#[derive(Clone, Debug)]
pub struct SimpleModification(rustyms::sequence::SimpleModification);

#[pymethods]
impl SimpleModification {
    #[new]
    fn new(name: &str) -> PyResult<Self> {
        match rustyms::sequence::SimpleModificationInner::parse_pro_forma(
            name,
            0..name.len(),
            &mut vec![],
            &mut vec![],
            None,
        ) {
            Ok((modification, _)) => Ok(Self(modification.defined().unwrap())),
            Err(_) => Err(PyValueError::new_err("Invalid modification")),
        }
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("Modification('{}')", self.0)
    }

    /// Molecular formula of the modification.
    ///
    /// Returns
    /// -------
    /// MolecularFormula
    ///
    fn formula(&self) -> MolecularFormula {
        MolecularFormula(self.0.formula())
    }

    /// Monoisotopic mass of the modification.
    ///
    /// Returns
    /// -------
    /// float
    ///
    fn monoisotopic_mass(&self) -> f64 {
        self.0.formula().monoisotopic_mass().value
    }
}

/// Amino acid modification.
///
/// Parameters
/// ----------
/// name : str
///   The name of the modification.
///
#[pyclass]
#[derive(Clone, Debug)]
pub struct Modification(rustyms::sequence::Modification);

#[pymethods]
impl Modification {
    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("Modification('{}')", self.0)
    }

    /// Molecular formula of the modification.
    ///
    /// Returns
    /// -------
    /// MolecularFormula
    ///
    fn formula(&self) -> MolecularFormula {
        MolecularFormula(self.0.formula())
    }

    /// Monoisotopic mass of the modification.
    ///
    /// Returns
    /// -------
    /// float
    ///
    fn monoisotopic_mass(&self) -> f64 {
        self.0.formula().monoisotopic_mass().value
    }
}

/// A theoretical fragment of a peptidoform.
#[pyclass]
#[derive(Debug)]
pub struct Fragment(rustyms::fragment::Fragment);

#[pymethods]
impl Fragment {
    fn __repr__(&self) -> String {
        format!(
            "Fragment(formula='{:?}', charge={}, ion='{}', peptidoform_ion_index={}, peptidoform_index={}, neutral_loss='{:?}')",
            Self::formula(self),
            self.charge(),
            self.ion().0,
            self.peptidoform_ion_index()
                .map_or("-".to_string(), |p| p.to_string()),
            self.peptidoform_index()
                .map_or("-".to_string(), |p| p.to_string()),
            self.neutral_loss(),
        )
    }

    /// The theoretical composition.
    ///
    /// Returns
    /// -------
    /// MolecularFormula | None
    ///
    #[getter]
    fn formula(&self) -> Option<MolecularFormula> {
        self.0.formula.clone().map(MolecularFormula)
    }

    /// The charge.
    ///
    /// Returns
    /// -------
    /// int
    ///
    #[getter]
    fn charge(&self) -> i16 {
        self.0.charge.value as i16
    }

    /// All possible annotations for this fragment.
    ///
    /// Returns
    /// -------
    /// str
    ///
    #[getter]
    fn ion(&self) -> FragmentType {
        FragmentType(self.0.ion.clone())
    }

    /// The peptidoform this fragment comes from, saved as the index into the list of peptidoforms in the overarching crate::PeptidoformIon struct.
    ///
    /// Returns
    /// -------
    /// int | None
    ///
    #[getter]
    const fn peptidoform_index(&self) -> Option<usize> {
        self.0.peptidoform_index
    }

    /// The peptidoform ion this fragment comes from, saved as the index into the list of peptidoform ions in the overarching crate::CompoundPeptidoformIon struct.
    ///
    /// Returns
    /// -------
    /// int | None
    ///
    #[getter]
    const fn peptidoform_ion_index(&self) -> Option<usize> {
        self.0.peptidoform_ion_index
    }

    /// Any neutral losses applied.
    ///
    /// Returns
    /// -------
    /// list[str]
    ///
    #[getter]
    fn neutral_loss(&self) -> Vec<String> {
        self.0
            .neutral_loss
            .iter()
            .map(ToString::to_string)
            .collect()
    }
}

/// All types of fragments.
#[pyclass]
#[derive(Debug)]
pub struct FragmentType(rustyms::fragment::FragmentType);

/// One block in a sequence meaning an amino acid and its accompanying modifications.
#[pyclass]
#[derive(Debug)]
pub struct SequenceElement(rustyms::sequence::SequenceElement<Linked>);

#[pymethods]
impl SequenceElement {
    fn __repr__(&self) -> String {
        format!(
            "SequenceElement(amino_acid='{}', modifications='{:?}', ambiguous='{:?}')",
            self.aminoacid(),
            self.modifications(),
            self.ambiguous()
        )
    }

    /// The amino acid.
    ///
    /// Returns
    /// -------
    /// AminoAcid
    ///
    #[getter]
    const fn aminoacid(&self) -> AminoAcid {
        AminoAcid(self.0.aminoacid.aminoacid())
    }

    /// All present modifications.
    ///
    /// Returns
    /// -------
    /// list[Modification]
    ///
    #[getter]
    fn modifications(&self) -> Vec<Modification> {
        self.0
            .modifications
            .iter()
            .map(|m| Modification(m.clone()))
            .collect()
    }

    /// If this amino acid is part of an ambiguous sequence group `(QA)?` in ProForma
    ///
    /// Returns
    /// -------
    /// int | None
    ///
    #[getter]
    const fn ambiguous(&self) -> Option<std::num::NonZeroU32> {
        self.0.ambiguous
    }
}

/// Fragmentation model enum.
#[pyclass(eq, eq_int)]
#[derive(Eq, PartialEq)]
enum FragmentationModel {
    All,
    CidHcd,
    Etd,
    Ethcd,
    Ead,
    Eacid,
    Uvpd,
}

/// Helper function to match a [`FragmentationModel`] to a rustyms Model.
fn match_model(
    model: &FragmentationModel,
) -> PyResult<rustyms::annotation::model::FragmentationModel> {
    match model {
        FragmentationModel::All => {
            Ok(rustyms::annotation::model::FragmentationModel::all().clone())
        }
        FragmentationModel::CidHcd => {
            Ok(rustyms::annotation::model::FragmentationModel::cid_hcd().clone())
        }
        FragmentationModel::Etd => {
            Ok(rustyms::annotation::model::FragmentationModel::etd().clone())
        }
        FragmentationModel::Ethcd => {
            Ok(rustyms::annotation::model::FragmentationModel::ethcd().clone())
        }
        FragmentationModel::Ead => {
            Ok(rustyms::annotation::model::FragmentationModel::ead().clone())
        }
        FragmentationModel::Eacid => {
            Ok(rustyms::annotation::model::FragmentationModel::eacid().clone())
        }
        FragmentationModel::Uvpd => {
            Ok(rustyms::annotation::model::FragmentationModel::uvpd().clone())
        }
    }
}

/// Parameters for matching theoretical fragments to measured data
///
/// Parameters
/// ----------
/// parameters : MatchingParameters
///     The parameters
///
#[pyclass]
#[derive(Clone, Debug)]
pub struct MatchingParameters(rustyms::annotation::model::MatchingParameters);

#[pymethods]
impl MatchingParameters {
    /// Create default parameters
    #[staticmethod]
    fn new() -> Self {
        Self(rustyms::annotation::model::MatchingParameters::default())
    }

    /// Set the tolerance to a certain ppm value
    #[setter]
    fn tolerance_ppm(&mut self, tolerance: f64) {
        self.0.tolerance = rustyms::quantities::Tolerance::new_ppm(tolerance);
    }

    /// Set the tolerance to a certain absolute thomson value
    #[setter]
    fn tolerance_thomson(&mut self, tolerance: f64) {
        self.0.tolerance =
            rustyms::quantities::Tolerance::new_absolute(rustyms::system::MassOverCharge::new::<
                rustyms::system::mz,
            >(tolerance));
    }
}

/// A position in a sequence
///
/// Parameters
/// ----------
/// position : SequencePosition
///     The position
///
#[pyclass]
#[derive(Clone, Copy, Debug)]
pub struct SequencePosition(rustyms::sequence::SequencePosition);

#[pymethods]
impl SequencePosition {
    /// Create a N-terminal position
    #[staticmethod]
    fn n_term() -> Self {
        Self(rustyms::sequence::SequencePosition::NTerm)
    }

    /// Create a position based on index (0-based indexing)
    #[staticmethod]
    fn index(index: usize) -> Self {
        Self(rustyms::sequence::SequencePosition::Index(index))
    }

    /// Create a C-terminal position
    #[staticmethod]
    fn c_term() -> Self {
        Self(rustyms::sequence::SequencePosition::CTerm)
    }

    /// Check if this is a N-terminal position
    #[getter]
    fn is_n_term(&self) -> bool {
        matches!(self, Self(rustyms::sequence::SequencePosition::NTerm))
    }

    /// Get the index of this position, if it is a terminal position this returns None.
    #[getter]
    fn get_index(&self) -> Option<usize> {
        match self.0 {
            rustyms::sequence::SequencePosition::Index(i) => Some(i),
            _ => None,
        }
    }

    /// Check if this is a C-terminal position
    #[getter]
    fn is_c_term(&self) -> bool {
        matches!(
            self,
            SequencePosition(rustyms::sequence::SequencePosition::CTerm)
        )
    }
}
/// A compound peptidoform ion with all data as provided by ProForma 2.0.
///
/// Parameters
/// ----------
/// proforma : str
///     The ProForma string.
///
#[pyclass]
#[derive(Clone, Debug)]
pub struct CompoundPeptidoformIon(rustyms::sequence::CompoundPeptidoformIon);

#[pymethods]
impl CompoundPeptidoformIon {
    /// Create a new compound peptidoform ion from a ProForma string.
    #[new]
    fn new(proforma: &str) -> Result<Self, BoxedError> {
        rustyms::sequence::CompoundPeptidoformIon::pro_forma(proforma, None)
            .map(CompoundPeptidoformIon)
            .map_err(|e| BoxedError(e.to_owned()))
    }

    /// Create a new compound peptidoform ion from a peptidoform ion.
    #[staticmethod]
    fn from_peptidoform_ion(peptidoform: PeptidoformIon) -> Self {
        Self(peptidoform.0.into())
    }

    /// Create a new compound peptidoform ion from a peptidoform.
    #[staticmethod]
    fn from_peptidoform(peptidoform: Peptidoform) -> Self {
        Self(peptidoform.0.into())
    }

    /// Get all peptidoform ions making up this compound peptidoform.
    ///
    /// Returns
    /// -------
    /// List[PeptidoformIon]
    ///
    #[getter]
    fn peptidoform_ions(&self) -> Vec<PeptidoformIon> {
        self.0
            .peptidoform_ions()
            .iter()
            .map(|p| PeptidoformIon(p.clone()))
            .collect()
    }

    /// Generate the theoretical fragments for this compound peptidoform ion, with the given maximal charge of the fragments,
    /// and the given model. With the global isotope modifications applied.
    ///
    /// Parameters
    /// ----------
    /// max_charge : int
    ///     The maximal charge of the fragments.
    /// model : FragmentationModel
    ///     The model to use for the fragmentation.
    ///
    /// Returns
    /// -------
    /// list[Fragment]
    ///   The theoretical fragments.
    ///
    fn generate_theoretical_fragments(
        &self,
        max_charge: isize,
        model: &FragmentationModel,
    ) -> PyResult<Vec<Fragment>> {
        Ok(self
            .0
            .generate_theoretical_fragments(
                rustyms::system::isize::Charge::new::<rustyms::system::e>(max_charge),
                &match_model(model)?,
            )
            .iter()
            .map(|f| Fragment(f.clone()))
            .collect())
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("CompoundPeptidoformIon({})", self.0)
    }

    fn __len__(&self) -> usize {
        self.0.peptidoform_ions().len()
    }
}

/// A peptidoform ion with all data as provided by ProForma 2.0.
///
/// Parameters
/// ----------
/// proforma : str
///     The ProForma string.
///
#[pyclass]
#[derive(Clone, Debug)]
pub struct PeptidoformIon(rustyms::sequence::PeptidoformIon);

#[pymethods]
impl PeptidoformIon {
    /// Create a new peptidoform ion from a ProForma string. Panics
    #[new]
    fn new(proforma: &str) -> Result<Self, BoxedError> {
        rustyms::sequence::PeptidoformIon::pro_forma(proforma, None)
            .map(PeptidoformIon)
            .map_err(|e| BoxedError(e.to_owned()))
    }

    /// Create a new peptidoform ion from a peptidoform.
    #[staticmethod]
    fn from_peptidoform(peptidoform: Peptidoform) -> Self {
        Self(peptidoform.0.into())
    }

    /// Get all peptidoforms making up this peptidoform ion.
    ///
    /// Returns
    /// -------
    /// List[Peptidoform]
    ///
    #[getter]
    fn peptidoforms(&self) -> Vec<Peptidoform> {
        self.0
            .peptidoforms()
            .iter()
            .map(|p| Peptidoform(p.clone()))
            .collect()
    }

    /// Generate the theoretical fragments for this peptidoform ion, with the given maximal charge of the fragments,
    /// and the given model. With the global isotope modifications applied.
    ///
    /// Parameters
    /// ----------
    /// max_charge : int
    ///     The maximal charge of the fragments.
    /// model : FragmentationModel
    ///     The model to use for the fragmentation.
    ///
    /// Returns
    /// -------
    /// list[Fragment]
    ///   The theoretical fragments.
    ///
    fn generate_theoretical_fragments(
        &self,
        max_charge: isize,
        model: &FragmentationModel,
    ) -> PyResult<Vec<Fragment>> {
        Ok(self
            .0
            .generate_theoretical_fragments(
                rustyms::system::isize::Charge::new::<rustyms::system::e>(max_charge),
                &match_model(model)?,
            )
            .iter()
            .map(|f| Fragment(f.clone()))
            .collect())
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("PeptidoformIon({})", self.0)
    }

    fn __len__(&self) -> usize {
        self.0.peptidoforms().len()
    }
}

/// A peptidoform with all data as provided by ProForma 2.0.
///
/// Parameters
/// ----------
/// proforma : str
///     The ProForma string.
///
#[pyclass]
#[derive(Clone, Debug)]
pub struct Peptidoform(rustyms::sequence::Peptidoform<Linked>);

#[pymethods]
impl Peptidoform {
    /// Create a new peptidoform from a ProForma string.
    #[new]
    fn new(proforma: &str) -> Result<Self, BoxedError> {
        rustyms::sequence::Peptidoform::pro_forma(proforma, None)
            .map(Peptidoform)
            .map_err(|e| BoxedError(e.to_owned()))
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("Peptidoform({})", self.0)
    }

    fn __len__(&self) -> usize {
        self.0.len()
    }

    /// Labile modifications, which will not be found in the actual spectrum.
    ///
    /// Returns
    /// -------
    /// list[Modification]
    ///
    #[getter]
    fn labile(&self) -> Vec<SimpleModification> {
        self.0
            .get_labile()
            .iter()
            .map(|x| SimpleModification(x.clone()))
            .collect()
    }

    /// N-terminal modification.
    ///
    /// Returns
    /// -------
    /// list[Modification]
    ///
    #[getter]
    fn n_term(&self) -> Vec<Modification> {
        self.0
            .get_n_term()
            .iter()
            .map(|m| Modification(m.clone()))
            .collect()
    }

    /// C-terminal modification.
    ///
    /// Returns
    /// -------
    /// list[Modification]
    ///
    #[getter]
    fn c_term(&self) -> Vec<Modification> {
        self.0
            .get_c_term()
            .iter()
            .map(|m| Modification(m.clone()))
            .collect()
    }

    /// Sequence of the peptidoform including modifications.
    ///
    /// Returns
    /// -------
    /// list[SequenceElement]
    ///
    #[getter]
    fn sequence(&self) -> Vec<SequenceElement> {
        self.0
            .sequence()
            .iter()
            .map(|x| SequenceElement(x.clone()))
            .collect()
    }

    /// For each ambiguous modification list all possible positions it can be placed on. Indexed by the ambiguous modification id.
    ///
    /// Returns
    /// -------
    /// list[list[int]]
    ///
    #[getter]
    fn ambiguous_modifications(&self) -> Vec<Vec<SequencePosition>> {
        self.0
            .get_ambiguous_modifications()
            .iter()
            .map(|p| p.iter().map(|p| SequencePosition(*p)).collect())
            .collect()
    }

    /// Stripped sequence, meaning the sequence without any modifications.
    ///
    /// Returns
    /// -------
    /// str
    ///
    #[getter]
    fn stripped_sequence(&self) -> String {
        self.0
            .sequence()
            .iter()
            .map(|x| x.aminoacid.pro_forma_definition())
            .collect()
    }

    /// The precursor charge of the peptidoform.
    ///
    /// Returns
    /// -------
    /// int | None
    ///
    #[getter]
    fn charge(&self) -> Option<isize> {
        self.0.get_charge_carriers().map(|c| c.charge().value)
    }

    /// The adduct ions, if specified.
    ///
    /// Returns
    /// -------
    /// MolecularCharge | None
    ///
    #[getter]
    fn charge_carriers(&self) -> Option<MolecularCharge> {
        self.0
            .get_charge_carriers()
            .map(|c| MolecularCharge(c.clone()))
    }

    /// Get a copy of the peptidoform with its sequence reversed.
    ///
    /// Returns
    /// -------
    /// Peptidoform
    ///
    fn reverse(&self) -> Self {
        Self(self.0.reverse())
    }

    /// Gives the formulas for the whole peptidoform. With the global isotope modifications applied. (Any B/Z will result in multiple possible formulas.)
    ///
    /// Returns
    /// -------
    /// List[MolecularFormula] | None
    ///
    fn formula(&self) -> Option<Vec<MolecularFormula>> {
        self.0.clone().into_linear().map(|p| {
            p.formulas()
                .iter()
                .map(|f| MolecularFormula(f.clone()))
                .collect()
        })
    }

    /// Generate the theoretical fragments for this peptidoform, with the given maximal charge of the fragments, and the given model. With the global isotope modifications applied.
    ///
    /// Parameters
    /// ----------
    /// max_charge : int
    ///     The maximal charge of the fragments.
    /// model : FragmentationModel
    ///     The model to use for the fragmentation.
    ///
    /// Returns
    /// -------
    /// list[Fragment] | None
    ///   The theoretical fragments.
    ///
    fn generate_theoretical_fragments(
        &self,
        max_charge: isize,
        model: &FragmentationModel,
    ) -> Option<Vec<Fragment>> {
        self.0.clone().into_linear().map(|p| {
            p.generate_theoretical_fragments(
                rustyms::system::isize::Charge::new::<rustyms::system::e>(max_charge),
                &match_model(model).unwrap(),
            )
            .iter()
            .map(|f| Fragment(f.clone()))
            .collect()
        })
    }
}

#[pyclass]
struct RawPeak(rustyms::spectrum::RawPeak);

#[pymethods]
impl RawPeak {
    fn __repr__(&self) -> String {
        format!("RawPeak(mz={}, intensity={})", self.mz(), self.intensity())
    }

    /// The m/z value of the peak.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn mz(&self) -> f64 {
        self.0.mz.value
    }

    /// The intensity of the peak.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn intensity(&self) -> f64 {
        self.0.intensity.into_inner()
    }
}

impl std::fmt::Display for RawPeak {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({}, {})", self.mz(), self.intensity())
    }
}

#[pyclass]
/// Represents an annotated peak in a mass spectrometry spectrum.
struct AnnotatedPeak(rustyms::annotation::AnnotatedPeak);

#[pymethods]
impl AnnotatedPeak {
    fn __repr__(&self) -> String {
        format!(
            "AnnotatedPeak(experimental_mz={}, intensity={}, annotation=[{:?}])",
            self.experimental_mz(),
            self.intensity(),
            self.annotation(),
        )
    }

    /// The experimental m/z value of the peak.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn experimental_mz(&self) -> f64 {
        self.0.experimental_mz.value
    }

    /// The intensity of the peak.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn intensity(&self) -> f64 {
        self.0.intensity.into_inner()
    }

    /// All annotations of the peak. Can be empty.
    ///
    /// Returns
    /// -------
    /// list[Fragment]
    ///
    #[getter]
    fn annotation(&self) -> Vec<Fragment> {
        self.0
            .annotation
            .iter()
            .map(|x| Fragment(x.clone()))
            .collect()
    }
}

impl std::fmt::Display for AnnotatedPeak {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "({}, {}, {})",
            self.experimental_mz(),
            self.intensity(),
            self.annotation()
                .iter()
                .map(Fragment::__repr__)
                .collect::<Vec<_>>()
                .join(", ")
        )
    }
}

/// A raw spectrum (meaning not annotated yet)
///
/// Parameters
/// ----------
/// title : str
///     The title of the spectrum.
/// num_scans : int
///     The number of scans.
/// rt : float
///     The retention time.
/// precursor_charge : float
///     The found precursor charge.
/// precursor_mass : float
///     The found precursor mass.
/// mz_array : list[float]
///     The m/z values of the peaks.
/// intensity_array : list[float]
///     The intensities of the peaks.
///
/// Returns
/// -------
/// RawSpectrum
///
#[pyclass]
#[derive(Debug)]
pub struct RawSpectrum(rustyms::spectrum::RawSpectrum);

#[pymethods]
impl RawSpectrum {
    /// Create a new raw spectrum.
    ///
    /// Parameters
    /// ----------
    /// title : str
    ///     The title of the spectrum.
    /// num_scans : int
    ///     The number of scans.
    /// rt : float
    ///     The retention time.
    /// precursor_charge : int
    ///     The found precursor charge.
    /// precursor_mass : float
    ///     The found precursor mass.
    /// mz_array : list[float]
    ///     The m/z values of the peaks.
    /// intensity_array : list[float]
    ///     The intensities of the peaks.
    ///
    /// Returns
    /// -------
    /// RawSpectrum
    ///
    #[new]
    #[pyo3(signature = (title, num_scans, mz_array, intensity_array, rt=None, precursor_charge=None, precursor_mass=None))]
    fn new(
        title: &str,
        num_scans: u64,
        mz_array: Vec<f64>,
        intensity_array: Vec<f64>,
        rt: Option<f64>,
        precursor_charge: Option<isize>,
        precursor_mass: Option<f64>,
    ) -> Self {
        let mut spectrum = rustyms::spectrum::RawSpectrum::default();
        spectrum.title = title.to_string();
        spectrum.num_scans = num_scans;
        spectrum.rt = rt.map(rustyms::system::Time::new::<rustyms::system::s>);
        spectrum.charge =
            precursor_charge.map(rustyms::system::isize::Charge::new::<rustyms::system::e>);
        spectrum.mass = precursor_mass.map(rustyms::system::Mass::new::<dalton>);

        let peaks = mz_array
            .into_iter()
            .zip(intensity_array)
            .map(|(mz, i)| rustyms::spectrum::RawPeak {
                mz: rustyms::system::MassOverCharge::new::<rustyms::system::mz>(mz),
                intensity: OrderedFloat(i),
            })
            .collect::<Vec<_>>();

        spectrum.extend(peaks);
        Self(spectrum)
    }

    fn __repr__(&self) -> String {
        format!(
            "RawSpectrum(title='{}', num_scans={}, rt={}, charge={}, mass={}, spectrum=[{}])",
            self.title(),
            self.num_scans(),
            self.rt().map_or("None".to_string(), |v| v.to_string()),
            self.charge().map_or("None".to_string(), |v| v.to_string()),
            self.mass().map_or("None".to_string(), |v| v.to_string()),
            self.spectrum()
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        )
    }

    /// The title.
    ///
    /// Returns
    /// -------
    /// str
    ///
    #[getter]
    fn title(&self) -> String {
        self.0.title.clone()
    }

    /// The number of scans.
    ///
    /// Returns
    /// -------
    /// int
    ///
    #[getter]
    fn num_scans(&self) -> u64 {
        self.0.num_scans
    }

    /// The retention time.
    ///
    /// Returns
    /// -------
    /// float | None
    ///
    #[getter]
    fn rt(&self) -> Option<f64> {
        self.0.rt.map(|v| v.value)
    }

    /// The found precursor charge.
    ///
    /// Returns
    /// -------
    /// float
    #[getter]
    fn charge(&self) -> Option<isize> {
        self.0.charge.map(|v| v.value)
    }

    /// The found precursor mass in Daltons.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn mass(&self) -> Option<f64> {
        self.0.mass.map(|v| v.get::<dalton>())
    }

    /// The peaks of which this spectrum consists.
    ///
    /// Returns
    /// -------
    /// list[RawPeak]
    ///
    #[getter]
    fn spectrum(&self) -> Vec<RawPeak> {
        self.0.clone().into_iter().map(RawPeak).collect()
    }

    /// Annotate this spectrum with the given peptidoform
    ///
    /// Parameters
    /// ----------
    /// peptidoform : CompoundPeptidoformIon
    ///     The peptidoform to annotate the spectrum with.
    /// model : FragmentationModel
    ///     The model to use for the fragmentation.
    /// parameters : MatchingParameters
    ///     The parameters to use for the matching.
    /// mode : MassMode
    ///    The mode to use for the mass.
    ///
    /// Returns
    /// -------
    /// AnnotatedSpectrum
    ///     The annotated spectrum.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the model is not one of the valid models.
    ///
    #[pyo3(signature = (peptidoform, model, parameters, mode=&MassMode::Monoisotopic))]
    fn annotate(
        &self,
        peptidoform: CompoundPeptidoformIon,
        model: &FragmentationModel,
        parameters: &MatchingParameters,
        mode: &MassMode,
    ) -> PyResult<AnnotatedSpectrum> {
        let rusty_model = match_model(model)?;
        let fragments = peptidoform.0.generate_theoretical_fragments(
            self.0
                .charge
                .unwrap_or_else(|| rustyms::system::isize::Charge::new::<rustyms::system::e>(1)),
            &rusty_model,
        );
        Ok(AnnotatedSpectrum(self.0.annotate(
            peptidoform.0,
            &fragments,
            &parameters.0,
            match mode {
                MassMode::Monoisotopic => rustyms::chemistry::MassMode::Monoisotopic,
                MassMode::Average => rustyms::chemistry::MassMode::Average,
                MassMode::MostAbundant => rustyms::chemistry::MassMode::MostAbundant,
            },
        )))
    }
}

/// An annotated spectrum.
#[pyclass]
#[derive(Debug)]
pub struct AnnotatedSpectrum(rustyms::annotation::AnnotatedSpectrum);

#[pymethods]
impl AnnotatedSpectrum {
    fn __repr__(&self) -> String {
        format!(
            "AnnotatedSpectrum(title='{}', num_scans={}, rt={}, charge={}, mass={}, spectrum=[{}])",
            self.title(),
            self.num_scans(),
            self.rt().map_or("None".to_string(), |v| v.to_string()),
            self.charge().map_or("None".to_string(), |v| v.to_string()),
            self.mass().map_or("None".to_string(), |v| v.to_string()),
            self.spectrum()
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        )
    }

    /// The title.
    ///
    /// Returns
    /// -------
    /// str
    ///
    #[getter]
    fn title(&self) -> String {
        self.0.title.clone()
    }

    /// The number of scans.
    ///
    /// Returns
    /// -------
    /// int
    ///
    #[getter]
    fn num_scans(&self) -> u64 {
        self.0.num_scans
    }

    /// The retention time.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn rt(&self) -> Option<f64> {
        self.0.rt.map(|v| v.value)
    }

    /// The found precursor charge.
    ///
    /// Returns
    /// -------
    /// float
    #[getter]
    fn charge(&self) -> Option<isize> {
        self.0.charge.map(|v| v.value)
    }

    /// The found precursor mass in Daltons.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn mass(&self) -> Option<f64> {
        self.0.mass.map(|v| v.get::<dalton>())
    }

    /// The peaks of which this spectrum consists.
    ///
    /// Returns
    /// -------
    /// list[AnnotatedPeak]
    ///
    #[getter]
    fn spectrum(&self) -> Vec<AnnotatedPeak> {
        self.0.clone().into_iter().map(AnnotatedPeak).collect()
    }
}

/// Python bindings to the rustyms library.
#[pymodule]
#[pyo3(name = "rustyms")]
fn rustyms_py03(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<AminoAcid>()?;
    m.add_class::<AnnotatedPeak>()?;
    m.add_class::<AnnotatedSpectrum>()?;
    m.add_class::<BoxedError>()?;
    m.add_class::<CompoundPeptidoformIon>()?;
    m.add_class::<Element>()?;
    m.add_class::<Fragment>()?;
    m.add_class::<FragmentationModel>()?;
    m.add_class::<FragmentType>()?;
    m.add_class::<MassMode>()?;
    m.add_class::<MatchingParameters>()?;
    m.add_class::<Modification>()?;
    m.add_class::<MolecularCharge>()?;
    m.add_class::<MolecularFormula>()?;
    m.add_class::<Peptidoform>()?;
    m.add_class::<PeptidoformIon>()?;
    m.add_class::<RawPeak>()?;
    m.add_class::<RawSpectrum>()?;
    m.add_class::<SequenceElement>()?;
    m.add_class::<SimpleModification>()?;
    m.add_function(wrap_pyfunction!(find_isobaric_sets, m)?)?;
    Ok(())
}

/// An error with context where it originated
#[pyclass]
#[derive(Debug)]
pub struct BoxedError(custom_error::BoxedError<'static, custom_error::BasicKind>);

impl std::error::Error for BoxedError {}

impl std::fmt::Display for BoxedError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl From<BoxedError> for PyErr {
    fn from(value: BoxedError) -> Self {
        PyValueError::new_err(value)
    }
}
