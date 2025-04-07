/// The definition of the position of an ion inside a glycan
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct GlycanPosition {
    /// The depth starting at the amino acid
    pub inner_depth: usize,
    /// The series number (from the ion series terminal)
    pub series_number: usize,
    /// The branch naming
    pub branch: Vec<(GlycanBranchIndex, GlycanBranchMassIndex)>,
    /// The aminoacid index where this glycan is attached
    pub attachment: Option<(AminoAcid, SequencePosition)>,
}
