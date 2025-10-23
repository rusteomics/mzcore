use std::sync::LazyLock;

use itertools::Itertools;

use mzcore::{
    chemistry::{MultiChemical, NeutralLoss},
    quantities::Tolerance,
    sequence::{BACKBONE, IsAminoAcid},
};

use crate::fragment::{BackboneCFragment, BackboneNFragment, Fragment, FragmentType};

/// Any data structure that can be written to mzPAF
pub trait ToMzPAF {
    /// Create a new mzPAF string from this element
    fn to_mz_paf_string(&self) -> String {
        let mut output = String::new();
        self.to_mz_paf(&mut output).unwrap(); // String writing cannot fail
        output
    }
    /// Write the mzPAF encoding of this element to the writer.
    /// # Errors
    /// When the writer errors.
    fn to_mz_paf(&self, w: impl std::fmt::Write) -> std::fmt::Result;
}

impl ToMzPAF for Fragment {
    /// Write the fragment as a [mzPAF](https://www.psidev.info/mzPAF) string.
    // TODO: figure out a way to handle the fallibility (when used on glycans/cross-linked stuff etc.)
    #[expect(clippy::cognitive_complexity)] // It is a very big function but breaking it up might not benefit readers
    fn to_mz_paf(&self, mut w: impl std::fmt::Write) -> std::fmt::Result {
        if self.auxiliary {
            write!(w, "&")?;
        }
        if let Some(number) = self.peptidoform_ion_index {
            write!(w, "{}@", number + 1)?;
        } else {
            write!(w, "0@")?;
        }
        // Push the ion type info (plus maybe some neutral losses if needed)
        match &self.ion {
            FragmentType::a(pos, variant)
            | FragmentType::b(pos, variant)
            | FragmentType::c(pos, variant)
            | FragmentType::x(pos, variant)
            | FragmentType::y(pos, variant) => write!(
                w,
                "{}{}{}",
                self.ion.kind(),
                pos.series_number,
                if *variant == 0 {
                    String::new()
                } else {
                    format!("{variant:+}H")
                }
            )?,
            FragmentType::z(pos, variant) => write!(
                w,
                "{}{}{}",
                self.ion.kind(),
                pos.series_number,
                if *variant == 1 {
                    String::new()
                } else {
                    format!("{:+}H", variant - 1)
                }
            )?,
            FragmentType::d(pos, aa, distance, variant, label) => {
                if *distance == 0 {
                    write!(
                        w,
                        "d{label}{}{}",
                        pos.series_number,
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )?;
                } else if let Some(loss) = aa
                    .satellite_ion_fragments(
                        pos.sequence_index,
                        self.peptidoform_index.unwrap_or_default(),
                        self.peptidoform_ion_index.unwrap_or_default(),
                    )
                    .and_then(|fragments| {
                        fragments
                            .iter()
                            .find(|f| f.0 == *label)
                            .map(|(_, loss)| loss.clone())
                    })
                {
                    write!(
                        w,
                        "a{}-{loss}{}",
                        pos.series_number,
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )?;
                } else {
                    write!(w, "?",)?;
                }
            }
            FragmentType::v(pos, aa, distance, variant) => {
                if *distance == 0 {
                    write!(
                        w,
                        "v{}{}",
                        pos.series_number,
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )?;
                } else {
                    write!(
                        w,
                        "y{}-{}{}",
                        pos.series_number,
                        aa.formulas()
                            .first()
                            .map(|f| f - LazyLock::force(&BACKBONE))
                            .unwrap_or_default(),
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )?;
                }
            }
            FragmentType::w(pos, aa, distance, variant, label) => {
                if *distance == 0 {
                    write!(
                        w,
                        "w{label}{}{}",
                        pos.series_number,
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )?;
                } else if let Some(loss) = aa
                    .satellite_ion_fragments(
                        pos.sequence_index,
                        self.peptidoform_index.unwrap_or_default(),
                        self.peptidoform_ion_index.unwrap_or_default(),
                    )
                    .and_then(|fragments| {
                        fragments
                            .iter()
                            .find(|f| f.0 == *label)
                            .map(|(_, loss)| loss.clone())
                    })
                {
                    write!(
                        w,
                        "z{}-{loss}{}",
                        pos.series_number,
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )?;
                } else {
                    write!(w, "?",)?;
                }
            }
            FragmentType::Precursor => write!(w, "p")?,
            FragmentType::PrecursorSideChainLoss(_, aa) => {
                write!(w, "p-r[sidechain_{aa}]")?;
            }
            FragmentType::Immonium(_, seq) => write!(
                w,
                "I{}{}",
                seq.aminoacid,
                seq.modifications
                    .iter()
                    .filter_map(|m| m.unimod_name().map(|name| format!("[{name}]")))
                    .join("") // TODO: how to handle ambiguous mods? Maybe store somewhere which where applied for this fragment
            )?,
            FragmentType::Unknown(num) => {
                if let Some(num) = num {
                    write!(w, "?{num}",)?;
                } else if let Some(formula) = &self.formula {
                    write!(w, "f{{{formula}}}",)?;
                } else {
                    write!(w, "?",)?;
                }
            }
            FragmentType::Diagnostic(_)
            | FragmentType::B { .. }
            | FragmentType::BComposition(_, _)
            | FragmentType::Y(_)
            | FragmentType::YComposition(_, _) => {
                if let Some(formula) = &self.formula
                    && formula.charge().value != 0
                {
                    // TODO: better way of storing?
                    write!(w, "f{{{}}}", formula.hill_notation_core())?;
                    if formula.additional_mass() != 0.0 {
                        write!(w, "{:+}", formula.additional_mass())?;
                    }
                    if formula.charge().value != 1 {
                        write!(w, "^{:}", formula.charge().value)?;
                    }
                } else {
                    write!(w, "?",)?;
                }
            }
            FragmentType::Internal(Some(name), a, b) => write!(
                w,
                "m{}:{}{}",
                a.sequence_index + 1,
                b.sequence_index + 1,
                match name {
                    (BackboneNFragment::a, BackboneCFragment::x)
                    | (BackboneNFragment::b, BackboneCFragment::y)
                    | (BackboneNFragment::c, BackboneCFragment::z) => "",
                    (BackboneNFragment::a, BackboneCFragment::y) => "-CO",
                    (BackboneNFragment::a, BackboneCFragment::z) => "-CHNO",
                    (BackboneNFragment::b, BackboneCFragment::x) => "+CO",
                    (BackboneNFragment::b, BackboneCFragment::z) => "-NH",
                    (BackboneNFragment::c, BackboneCFragment::x) => "+CHNO",
                    (BackboneNFragment::c, BackboneCFragment::y) => "+NH",
                }
            )?,
            FragmentType::Internal(None, a, b) => {
                write!(w, "m{}:{}", a.sequence_index + 1, b.sequence_index + 1)?;
            }
        }
        // More losses
        for loss in &self.neutral_loss {
            match loss {
                NeutralLoss::SideChainLoss(_, aa) => {
                    write!(w, "-r[sidechain_{aa}]")?;
                }
                NeutralLoss::Gain(1, mol) => {
                    write!(w, "+{mol}")?;
                }
                NeutralLoss::Loss(1, mol) => {
                    write!(w, "-{mol}")?;
                }
                l => write!(w, "{l}")?,
            }
        }
        // Isotopes: TODO: not handled
        // Charge state
        if self.charge.value != 1 {
            write!(w, "^{}", self.charge.value)?;
        }
        // Deviation
        match self.deviation {
            Some(Tolerance::Absolute(abs)) => write!(w, "/{:.3}", abs.value)?,
            Some(Tolerance::Relative(ppm)) => {
                write!(w, "/{:.3}ppm", ppm.get::<mzcore::system::ratio::ppm>())?;
            }
            None => (),
        }
        // Confidence
        if let Some(confidence) = self.confidence {
            write!(w, "*{confidence}")?;
        }
        Ok(())
    }
}
