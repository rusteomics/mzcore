//! Python bindings to the mzcore library.
#![allow(clippy::doc_markdown, clippy::trivially_copy_pass_by_ref)]

use std::fmt::Debug;
use std::num::NonZeroU16;

use pyo3::{exceptions::PyValueError, prelude::*, types::PyType};

use mzalign::AlignScoring;
use mzannotate::prelude::{GlycanFragmention, PeptidoformFragmentation, ToMzPAF};
use mzcore::{
    chemistry::{Chemical, MultiChemical},
    sequence::{IsAminoAcid, Linked, SimpleLinear},
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
                    mzcore::sequence::PlacementRule::AminoAcid(
                        r.into_iter().map(|a| a.0).collect(),
                        mzcore::sequence::Position::Anywhere,
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
                    mzcore::sequence::PlacementRule::AminoAcid(
                        r.into_iter().map(|a| a.0).collect(),
                        mzcore::sequence::Position::Anywhere,
                    )
                }),
            )
        })
        .collect::<Vec<_>>();

    mzcore::prelude::find_isobaric_sets(
        mzcore::system::Mass::new::<dalton>(mass),
        mzcore::quantities::Tolerance::new_ppm(tolerance),
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
pub struct Element(mzcore::chemistry::Element);

#[pymethods]
impl Element {
    #[new]
    fn new(symbol: &str) -> PyResult<Self> {
        mzcore::chemistry::Element::try_from(symbol)
            .map(Element)
            .map_err(|()| PyValueError::new_err("Invalid element symbol."))
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
pub struct MolecularFormula(mzcore::chemistry::MolecularFormula);

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
        let formula = mzcore::chemistry::MolecularFormula::new(elements.as_slice(), &[])
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
        mzcore::chemistry::MolecularFormula::from_pro_forma::<false, false>(proforma, ..)
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
        mzcore::chemistry::MolecularFormula::from_psi_mod(psi_mod, ..)
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
        if let Err(e) = self
            .0
            .add((element.0, isotope.and_then(NonZeroU16::new), n))
        {
            Err(PyValueError::new_err(format!(
                "Invalid element or isotope: {e}"
            )))
        } else {
            Ok(None)
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
pub struct MolecularCharge(mzcore::chemistry::MolecularCharge);

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
        Self(mzcore::chemistry::MolecularCharge {
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
        Self(mzcore::chemistry::MolecularCharge::proton(
            mzcore::system::isize::Charge::new::<mzcore::system::e>(charge as isize),
        ))
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
pub struct AminoAcid(mzcore::sequence::AminoAcid);

#[pymethods]
impl AminoAcid {
    #[new]
    fn new(name: &str) -> PyResult<Self> {
        mzcore::sequence::AminoAcid::try_from(name).map_or_else(
            |()| Err(PyValueError::new_err("Invalid amino acid")),
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
pub struct SimpleModification(mzcore::sequence::SimpleModification);

#[pymethods]
impl SimpleModification {
    #[new]
    fn new(name: &str) -> PyResult<Self> {
        match mzcore::sequence::SimpleModificationInner::parse_pro_forma(
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
pub struct Modification(mzcore::sequence::Modification);

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
pub struct Fragment(mzannotate::fragment::Fragment);

#[pymethods]
impl Fragment {
    fn __repr__(&self) -> String {
        format!(
            "Fragment(formula='{:?}', charge={}, ion='{}', peptidoform_ion_index={}, peptidoform_index={}, neutral_loss='{:?}')",
            Self::formula(self),
            self.charge(),
            self.ion().0, // TODO: this could crash
            self.peptidoform_ion_index()
                .map_or_else(|| "-".to_string(), |p| p.to_string()),
            self.peptidoform_index()
                .map_or_else(|| "-".to_string(), |p| p.to_string()),
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

    /// Get the mzPAF annotation of this fragment.
    ///
    /// Returns
    /// -------
    /// str
    ///
    #[getter]
    fn to_mz_paf(&self) -> String {
        self.0.to_mz_paf_string()
    }

    /// The charge.
    ///
    /// Returns
    /// -------
    /// int
    ///
    #[getter]
    const fn charge(&self) -> i16 {
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
pub struct FragmentType(mzannotate::fragment::FragmentType);

#[pymethods]
impl FragmentType {
    /// The kind of fragment
    ///
    #[getter]
    fn kind(&self) -> FragmentKind {
        self.0.kind().into()
    }

    /// The labels for this fragment.
    ///
    /// /// Returns
    /// -------
    /// list[str | None]
    ///     (optional superscript prefix, main label, optional position label)
    ///
    #[getter]
    fn label(&self) -> (Option<String>, String, Option<String>) {
        let (a, b) = self.0.label();
        let c = self.0.position_label();
        (a, b.to_string(), c)
    }
}

/// The kind of fragment
#[pyclass(eq, eq_int)]
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
#[expect(non_camel_case_types)]
pub enum FragmentKind {
    /// a
    a,
    /// b
    b,
    /// c
    c,
    /// d
    d,
    /// v
    v,
    /// w
    w,
    /// x
    x,
    /// y
    y,
    /// z and z·
    z,
    /// glycan Y fragment, generated by one or more branches broken
    Y,
    /// B or glycan diagnostic ion or Internal glycan fragment, meaning both a B and Y breakages (and potentially multiple of both), resulting in a set of monosaccharides
    B,
    /// Immonium ion
    immonium,
    /// Precursor with amino acid side chain loss
    precursor_side_chain_loss,
    /// Diagnostic ion for a given position
    diagnostic,
    /// Internal ion
    internal,
    /// precursor
    precursor,
    /// unknown fragment
    unknown,
}

impl From<mzannotate::fragment::FragmentKind> for FragmentKind {
    fn from(value: mzannotate::fragment::FragmentKind) -> Self {
        match value {
            mzannotate::fragment::FragmentKind::a => Self::a,
            mzannotate::fragment::FragmentKind::b => Self::b,
            mzannotate::fragment::FragmentKind::c => Self::c,
            mzannotate::fragment::FragmentKind::d => Self::d,
            mzannotate::fragment::FragmentKind::v => Self::v,
            mzannotate::fragment::FragmentKind::w => Self::w,
            mzannotate::fragment::FragmentKind::x => Self::x,
            mzannotate::fragment::FragmentKind::y => Self::y,
            mzannotate::fragment::FragmentKind::z => Self::z,
            mzannotate::fragment::FragmentKind::Y => Self::Y,
            mzannotate::fragment::FragmentKind::B => Self::B,
            mzannotate::fragment::FragmentKind::immonium => Self::immonium,
            mzannotate::fragment::FragmentKind::precursor_side_chain_loss => {
                Self::precursor_side_chain_loss
            }
            mzannotate::fragment::FragmentKind::diagnostic => Self::diagnostic,
            mzannotate::fragment::FragmentKind::internal => Self::internal,
            mzannotate::fragment::FragmentKind::precursor => Self::precursor,
            mzannotate::fragment::FragmentKind::unknown => Self::unknown,
        }
    }
}

/// One block in a sequence meaning an amino acid and its accompanying modifications.
#[pyclass]
#[derive(Debug)]
pub struct SequenceElement(mzcore::sequence::SequenceElement<Linked>);

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
    CID,
    ETD,
    ETciD,
    EAD,
    EAciD,
    UVPD,
}

/// Helper function to match a [`FragmentationModel`] to a mzcore Model.
fn match_model(model: &FragmentationModel) -> mzannotate::annotation::model::FragmentationModel {
    match model {
        FragmentationModel::All => mzannotate::annotation::model::FragmentationModel::all().clone(),
        FragmentationModel::CID => mzannotate::annotation::model::FragmentationModel::cid().clone(),
        FragmentationModel::ETD => mzannotate::annotation::model::FragmentationModel::etd().clone(),
        FragmentationModel::ETciD => {
            mzannotate::annotation::model::FragmentationModel::etcid().clone()
        }
        FragmentationModel::EAD => mzannotate::annotation::model::FragmentationModel::ead().clone(),
        FragmentationModel::EAciD => {
            mzannotate::annotation::model::FragmentationModel::eacid().clone()
        }
        FragmentationModel::UVPD => {
            mzannotate::annotation::model::FragmentationModel::uvpd().clone()
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
pub struct MatchingParameters(mzannotate::annotation::model::MatchingParameters);

#[pymethods]
impl MatchingParameters {
    /// Create default parameters
    #[staticmethod]
    fn new() -> Self {
        Self(mzannotate::annotation::model::MatchingParameters::default())
    }

    /// Set the tolerance to a certain ppm value
    #[setter]
    fn tolerance_ppm(&mut self, tolerance: f64) {
        self.0.tolerance = mzcore::quantities::Tolerance::new_ppm(tolerance);
    }

    /// Set the tolerance to a certain absolute thomson value
    #[setter]
    fn tolerance_thomson(&mut self, tolerance: f64) {
        self.0.tolerance =
            mzcore::quantities::Tolerance::new_absolute(mzcore::system::MassOverCharge::new::<
                mzcore::system::thomson,
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
pub struct SequencePosition(mzcore::sequence::SequencePosition);

#[pymethods]
impl SequencePosition {
    /// Create a N-terminal position
    #[staticmethod]
    const fn n_term() -> Self {
        Self(mzcore::sequence::SequencePosition::NTerm)
    }

    /// Create a position based on index (0-based indexing)
    #[staticmethod]
    const fn index(index: usize) -> Self {
        Self(mzcore::sequence::SequencePosition::Index(index))
    }

    /// Create a C-terminal position
    #[staticmethod]
    const fn c_term() -> Self {
        Self(mzcore::sequence::SequencePosition::CTerm)
    }

    /// Check if this is a N-terminal position
    #[getter]
    const fn is_n_term(&self) -> bool {
        matches!(self.0, mzcore::sequence::SequencePosition::NTerm)
    }

    /// Get the index of this position, if it is a terminal position this returns None.
    #[getter]
    const fn get_index(&self) -> Option<usize> {
        match self.0 {
            mzcore::sequence::SequencePosition::Index(i) => Some(i),
            _ => None,
        }
    }

    /// Check if this is a C-terminal position
    #[getter]
    const fn is_c_term(&self) -> bool {
        matches!(self.0, mzcore::sequence::SequencePosition::CTerm)
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
pub struct CompoundPeptidoformIon(mzcore::sequence::CompoundPeptidoformIon);

#[pymethods]
impl CompoundPeptidoformIon {
    /// Create a new compound peptidoform ion from a ProForma string.
    #[new]
    fn new(proforma: &str) -> Result<Self, BoxedError> {
        mzcore::sequence::CompoundPeptidoformIon::pro_forma(proforma, None)
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
    ) -> Vec<Fragment> {
        self.0
            .generate_theoretical_fragments(
                mzcore::system::isize::Charge::new::<mzcore::system::e>(max_charge),
                &match_model(model),
            )
            .iter()
            .map(|f| Fragment(f.clone()))
            .collect()
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

/// A glycan structure
///
/// Parameters
/// ----------
/// iupac : str
///     The iupac condensed definition.
///
#[pyclass]
#[derive(Clone, Debug)]
pub struct GlycanStructure(mzcore::glycan::GlycanStructure);

#[pymethods]
impl GlycanStructure {
    /// Parse a glycan structure from the iupac condensed defiition.
    #[new]
    fn new(iupac: &str) -> Result<Self, BoxedError> {
        mzcore::glycan::GlycanStructure::from_short_iupac(iupac, 0..iupac.len(), 0)
            .map(GlycanStructure)
            .map_err(|e| BoxedError(e.to_owned()))
    }

    /// Get the full formula for this glycan structure.
    fn formula(&self) -> MolecularFormula {
        MolecularFormula(self.0.formula())
    }

    /// Generate the theoretical fragments for this glycan structure, with the given maximal charge of the fragments,
    /// and the given model.
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
    ) -> Vec<Fragment> {
        let full = self.0.formula();
        self.0
            .clone()
            .determine_positions()
            .generate_theoretical_fragments(
                &match_model(model),
                0,
                0,
                &mut mzcore::chemistry::MolecularCharge::proton(
                    mzcore::system::isize::Charge::new::<mzcore::system::e>(max_charge),
                )
                .into(),
                &full.into(),
                None,
            )
            .iter()
            .map(|f| Fragment(f.clone()))
            .collect()
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("GlycanStructure({})", self.0)
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
pub struct PeptidoformIon(mzcore::sequence::PeptidoformIon);

#[pymethods]
impl PeptidoformIon {
    /// Create a new peptidoform ion from a ProForma string. Panics
    #[new]
    fn new(proforma: &str) -> Result<Self, BoxedError> {
        mzcore::sequence::PeptidoformIon::pro_forma(proforma, None)
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
    ) -> Vec<Fragment> {
        self.0
            .generate_theoretical_fragments(
                mzcore::system::isize::Charge::new::<mzcore::system::e>(max_charge),
                &match_model(model),
            )
            .into_iter()
            .map(Fragment)
            .collect()
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
pub struct Peptidoform(mzcore::sequence::Peptidoform<Linked>);

#[pymethods]
impl Peptidoform {
    /// Create a new peptidoform from a ProForma string.
    #[new]
    fn new(proforma: &str) -> Result<Self, BoxedError> {
        mzcore::sequence::Peptidoform::pro_forma(proforma, None)
            .map(Peptidoform)
            .map_err(|e| BoxedError(e.to_owned()))
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("Peptidoform({})", self.0)
    }

    const fn __len__(&self) -> usize {
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
                mzcore::system::isize::Charge::new::<mzcore::system::e>(max_charge),
                &match_model(model),
            )
            .into_iter()
            .map(Fragment)
            .collect()
        })
    }

    /// Align this peptidoform to another peptidoform. This uses a maximal isobaric step of 4.
    ///
    /// Parameters
    /// ----------
    /// other : Peptidoform
    ///     The other peptidoform to align against
    ///
    /// Returns
    /// -------
    /// Alignment | None
    ///   The resulting alignment or None if at least one of the peptiforms was not simple linear.
    ///
    fn align_4(&self, other: &Self) -> Option<Alignment> {
        self.0
            .clone()
            .into_simple_linear()
            .and_then(|s| {
                other.0.clone().into_simple_linear().map(|o| {
                    mzalign::align::<
                        4,
                        mzcore::sequence::Peptidoform<SimpleLinear>,
                        mzcore::sequence::Peptidoform<SimpleLinear>,
                    >(s, o, AlignScoring::default(), mzalign::AlignType::GLOBAL)
                })
            })
            .map(Alignment)
    }
}

#[pyclass]
struct Alignment(
    mzalign::Alignment<
        mzcore::sequence::Peptidoform<SimpleLinear>,
        mzcore::sequence::Peptidoform<SimpleLinear>,
    >,
);

#[pymethods]
impl Alignment {
    fn __repr__(&self) -> String {
        format!(
            "Alignment(seq_a={}, seq_b={}, start_a={}, start_b={}, path='{}', type={}, maximal_step={})",
            self.0.seq_a(),
            self.0.seq_b(),
            self.0.start_a(),
            self.0.start_b(),
            self.0.short(),
            self.0.align_type(),
            self.0.max_step(),
        )
    }
}

#[pyclass]
/// Represents an annotated peak in a mass spectrometry spectrum.
struct AnnotatedPeak(mzannotate::spectrum::AnnotatedPeak<mzannotate::fragment::Fragment>);

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
    const fn experimental_mz(&self) -> f64 {
        self.0.mz.value
    }

    /// The intensity of the peak.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    const fn intensity(&self) -> f32 {
        self.0.intensity
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
            .annotations
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

/// An annotated spectrum.
#[pyclass]
#[derive(Debug)]
pub struct AnnotatedSpectrum(mzannotate::spectrum::AnnotatedSpectrum);

#[pymethods]
impl AnnotatedSpectrum {
    fn __repr__(&self) -> String {
        format!(
            "AnnotatedSpectrum(spectrum=[{}])",
            self.spectrum()
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(", ")
        )
    }

    /// The peaks of which this spectrum consists.
    ///
    /// Returns
    /// -------
    /// list[AnnotatedPeak]
    ///
    #[getter]
    fn spectrum(&self) -> Vec<AnnotatedPeak> {
        self.0.peaks.iter().cloned().map(AnnotatedPeak).collect()
    }
}

/// Python bindings to the mzcore library.
#[pymodule]
#[pyo3(name = "mzcore")]
fn mzcore_py03(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Alignment>()?;
    m.add_class::<AminoAcid>()?;
    m.add_class::<AnnotatedPeak>()?;
    m.add_class::<AnnotatedSpectrum>()?;
    m.add_class::<BoxedError>()?;
    m.add_class::<CompoundPeptidoformIon>()?;
    m.add_class::<Element>()?;
    m.add_class::<Fragment>()?;
    m.add_class::<FragmentationModel>()?;
    m.add_class::<FragmentKind>()?;
    m.add_class::<FragmentKind>()?;
    m.add_class::<FragmentType>()?;
    m.add_class::<GlycanStructure>()?;
    m.add_class::<GlycanStructure>()?;
    m.add_class::<MassMode>()?;
    m.add_class::<MatchingParameters>()?;
    m.add_class::<Modification>()?;
    m.add_class::<MolecularCharge>()?;
    m.add_class::<MolecularFormula>()?;
    m.add_class::<Peptidoform>()?;
    m.add_class::<PeptidoformIon>()?;
    m.add_class::<SequenceElement>()?;
    m.add_class::<SimpleModification>()?;
    m.add_function(wrap_pyfunction!(find_isobaric_sets, m)?)?;
    Ok(())
}

/// An error with context where it originated
#[pyclass]
#[derive(Debug)]
pub struct BoxedError(context_error::BoxedError<'static, context_error::BasicKind>);

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
