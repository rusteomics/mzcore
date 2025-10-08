pub(super) use super::*;
use std::fmt::{Display, Write};

/// Display things and allow the use of fancy non ascii characters
pub trait FancyDisplay: Display {
    /// Equivalent of `.to_string()` but then fancier!
    fn to_fancy_string(&self) -> String;
}

impl FancyDisplay for Gene {
    fn to_fancy_string(&self) -> String {
        const fn to_roman(n: usize) -> &'static str {
            ["0", "Ⅰ", "Ⅱ", "Ⅲ", "Ⅳ", "Ⅴ", "Ⅵ", "Ⅶ", "Ⅷ", "Ⅸ", "Ⅹ"][n]
        }
        let mut f = String::new();

        write!(
            f,
            "Ig{}{}{}{}",
            self.chain.to_fancy_string(),
            self.kind.to_fancy_string(),
            self.number
                .as_ref()
                .map_or_else(String::new, |n| format!("({})", to_roman(*n))),
            if self.number.is_some() && !self.family.is_empty() {
                "-"
            } else {
                ""
            }
        )
        .unwrap();

        let mut first = true;
        let mut last_str = false;
        for element in &self.family {
            if !first && !last_str {
                write!(f, "-").unwrap();
            }
            write!(
                f,
                "{}{}",
                element.0.map(|i| i.to_string()).unwrap_or_default(),
                element.1
            )
            .unwrap();
            last_str = !element.1.is_empty();
            first = false;
        }
        f
    }
}

impl FancyDisplay for ChainType {
    fn to_fancy_string(&self) -> String {
        match self {
            Self::Heavy => "",
            Self::LightKappa => "κ",
            Self::LightLambda => "λ",
            Self::Iota => "ι",
        }
        .to_string()
    }
}

impl FancyDisplay for GeneType {
    fn to_fancy_string(&self) -> String {
        match self {
            Self::V => "V",
            Self::J => "J",
            Self::C(None) => "C",
            Self::C(Some(Constant::A)) => "α",
            Self::C(Some(Constant::D)) => "δ",
            Self::C(Some(Constant::E)) => "ε",
            Self::C(Some(Constant::G)) => "ɣ",
            Self::C(Some(Constant::M)) => "μ",
            Self::C(Some(Constant::O)) => "ο",
            Self::C(Some(Constant::T)) => "τ",
        }
        .to_string()
    }
}
