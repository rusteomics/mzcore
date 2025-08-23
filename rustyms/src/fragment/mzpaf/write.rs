use std::{fmt::Write, sync::LazyLock};

use itertools::Itertools;

use crate::{
    chemistry::MultiChemical,
    fragment::{BackboneCFragment, BackboneNFragment, Fragment, FragmentType, NeutralLoss},
    quantities::Tolerance,
    sequence::{BACKBONE, IsAminoAcid},
};

impl Fragment {
    /// Write the fragment as a [mzPAF](https://www.psidev.info/mzPAF) string.
    // TODO: figure out a way to handle the fallibility (when used on glycans/cross-linked stuff etc)
    pub fn to_mz_paf(&self) -> String {
        let mut output = String::new();
        if self.auxiliary {
            output.push('&');
        }
        if let Some(number) = self.peptidoform_ion_index {
            write!(&mut output, "{}@", number + 1).unwrap();
        } else {
            write!(&mut output, "0@").unwrap();
        }
        // Push the ion type info (plus maybe some neutral losses if needed)
        match &self.ion {
            FragmentType::a(pos, variant)
            | FragmentType::b(pos, variant)
            | FragmentType::c(pos, variant)
            | FragmentType::x(pos, variant)
            | FragmentType::y(pos, variant) => write!(
                &mut output,
                "{}{}{}",
                self.ion.kind(),
                pos.series_number,
                if *variant == 0 {
                    String::new()
                } else {
                    format!("{variant:+}H")
                }
            )
            .unwrap(),
            FragmentType::z(pos, variant) => write!(
                &mut output,
                "{}{}{}",
                self.ion.kind(),
                pos.series_number,
                if *variant == 1 {
                    String::new()
                } else {
                    format!("{:+}H", variant - 1)
                }
            )
            .unwrap(),
            FragmentType::d(pos, aa, distance, variant, label) => {
                if *distance == 0 {
                    write!(
                        &mut output,
                        "d{label}{}{}",
                        pos.series_number,
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )
                    .unwrap();
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
                        &mut output,
                        "a{}-{loss}{}",
                        pos.series_number,
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )
                    .unwrap();
                } else {
                    write!(&mut output, "?",).unwrap();
                }
            }
            FragmentType::v(pos, aa, distance, variant) => {
                if *distance == 0 {
                    write!(
                        &mut output,
                        "v{}{}",
                        pos.series_number,
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )
                    .unwrap();
                } else {
                    write!(
                        &mut output,
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
                    )
                    .unwrap();
                }
            }
            FragmentType::w(pos, aa, distance, variant, label) => {
                if *distance == 0 {
                    write!(
                        &mut output,
                        "w{label}{}{}",
                        pos.series_number,
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )
                    .unwrap();
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
                        &mut output,
                        "z{}-{loss}{}",
                        pos.series_number,
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )
                    .unwrap();
                } else {
                    write!(&mut output, "?",).unwrap();
                }
            }
            FragmentType::Precursor => write!(&mut output, "p").unwrap(),
            FragmentType::PrecursorSideChainLoss(_, aa) => {
                write!(&mut output, "p-r[sidechain_{aa}]").unwrap();
            }
            FragmentType::Immonium(_, seq) => write!(
                &mut output,
                "I{}{}",
                seq.aminoacid,
                seq.modifications
                    .iter()
                    .filter_map(|m| m.unimod_name().map(|name| format!("[{name}]")))
                    .join("") // TODO: how to handle ambiguous mods? maybe store somewhere which where applied for this fragment
            )
            .unwrap(),
            FragmentType::Unknown(num) => {
                if let Some(num) = num {
                    write!(&mut output, "?{num}",).unwrap();
                } else if let Some(formula) = &self.formula {
                    write!(&mut output, "f{{{formula}}}",).unwrap();
                } else {
                    write!(&mut output, "?",).unwrap();
                }
            }
            FragmentType::Diagnostic(_)
            | FragmentType::B { .. }
            | FragmentType::BComposition(_, _)
            | FragmentType::Y(_)
            | FragmentType::YComposition(_, _) => {
                if let Some(formula) = &self.formula {
                    // TODO: better way of storing?
                    write!(&mut output, "f{{{formula}}}",).unwrap();
                } else {
                    write!(&mut output, "?",).unwrap();
                }
            }
            FragmentType::Internal(Some(name), a, b) => write!(
                &mut output,
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
            )
            .unwrap(),
            FragmentType::Internal(None, a, b) => write!(
                &mut output,
                "m{}:{}",
                a.sequence_index + 1,
                b.sequence_index + 1
            )
            .unwrap(),
        }
        // More losses
        for loss in &self.neutral_loss {
            match loss {
                NeutralLoss::SideChainLoss(_, aa) => {
                    write!(&mut output, "-r[sidechain_{aa}]").unwrap();
                }
                NeutralLoss::Gain(1, mol) => {
                    write!(&mut output, "+{mol}").unwrap();
                }
                NeutralLoss::Loss(1, mol) => {
                    write!(&mut output, "-{mol}").unwrap();
                }
                l => write!(&mut output, "{l}").unwrap(),
            }
        }
        // Isotopes: TODO: not handled
        // Charge state
        if self.charge.value != 1 {
            write!(&mut output, "^{}", self.charge.value).unwrap();
        }
        // Deviation
        match self.deviation {
            Some(Tolerance::Absolute(abs)) => write!(&mut output, "/{:.3}", abs.value).unwrap(),
            Some(Tolerance::Relative(ppm)) => write!(
                &mut output,
                "/{:.3}ppm",
                ppm.get::<crate::system::ratio::ppm>()
            )
            .unwrap(),
            None => (),
        }
        // Confidence
        if let Some(confidence) = self.confidence {
            write!(&mut output, "*{confidence}").unwrap();
        }
        output
    }
}
