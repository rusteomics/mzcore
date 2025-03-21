/// All possible neutral losses
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum NeutralLoss {
    /// Gain of a specific formula
    Gain(MolecularFormula),
    /// Loss of a specific formula
    Loss(MolecularFormula),
    /// Loss of a side chain of an amino acid
    SideChainLoss(MolecularFormula, crate::AminoAcid),
}

/// A diagnostic ion, defined in M (not MH+) chemical formula
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub struct DiagnosticIon(pub MolecularFormula);
