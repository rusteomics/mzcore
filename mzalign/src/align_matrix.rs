use crate::{AlignScoring, AlignType, MatchType, Piece};

pub(super) struct Matrix {
    value: Vec<Vec<Piece>>,
    a: usize,
    b: usize,
}

impl std::fmt::Debug for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use std::fmt::Write;
        for column in &self.value {
            let mut line_0 = String::new();
            let mut line_1 = String::new();
            for cell in column {
                let top = format!("{}/{} {:2}", cell.step_a, cell.step_b, cell.local_score);
                let bottom = format!(
                    "{} {:3}",
                    match cell.match_type {
                        MatchType::FullIdentity => "FI",
                        MatchType::Gap => "G ",
                        MatchType::IdentityMassMismatch => "IM",
                        MatchType::Isobaric => "I ",
                        MatchType::Rotation => "R ",
                        MatchType::Mismatch => "M ",
                    },
                    cell.score
                );
                write!(&mut line_0, "⎡{top:0$}⎤", top.len().max(bottom.len()))?;
                write!(&mut line_1, "⎣{bottom:0$}⎦", top.len().max(bottom.len()))?;
            }
            writeln!(f, "{line_0}")?;
            writeln!(f, "{line_1}")?;
        }
        Ok(())
    }
}

impl Matrix {
    pub(super) fn new(a: usize, b: usize) -> Self {
        Self {
            value: vec![vec![Piece::default(); b + 1]; a + 1],
            a,
            b,
        }
    }

    #[expect(clippy::cast_possible_wrap)]
    pub(super) fn global_start(&mut self, is_a: bool, scoring: AlignScoring<'_>) {
        let max = if is_a { self.a } else { self.b };
        for index in 0..=max {
            self.value[if is_a { index } else { 0 }][if is_a { 0 } else { index }] = Piece::new(
                match index {
                    0 => 0,
                    _ => {
                        scoring.gap_start as isize + (index as isize) * scoring.gap_extend as isize
                    }
                },
                match index {
                    0 => 0,
                    1 => scoring.gap_start as isize + scoring.gap_extend as isize,
                    _ => scoring.gap_extend as isize,
                },
                MatchType::Gap,
                if is_a { u16::from(index != 0) } else { 0 },
                if is_a { 0 } else { u16::from(index != 0) },
            );
        }
    }

    pub(super) fn trace_path(
        &self,
        ty: AlignType,
        high: (isize, usize, usize),
    ) -> (usize, usize, Vec<Piece>) {
        let mut path = Vec::new();
        let mut high = self.find_end(ty, high);

        // Loop back to left side
        while ty.left.global() || !(high.1 == 0 && high.2 == 0) {
            let value = self.value[high.1][high.2];
            if value.step_a == 0 && value.step_b == 0 || !ty.left.global() && value.score < 0 {
                break;
            }
            high = (
                0,
                high.1 - value.step_a as usize,
                high.2 - value.step_b as usize,
            );
            path.push(value);
        }
        (high.1, high.2, path.into_iter().rev().collect())
    }

    fn find_end(&self, ty: AlignType, high: (isize, usize, usize)) -> (isize, usize, usize) {
        if ty.right.global_a() && ty.right.global_a() {
            (self.value[self.a][self.b].score, self.a, self.b)
        } else if ty.right.global_b() {
            let value = (0..=self.a)
                .map(|v| (v, self.value[v][self.b].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            (value.1, value.0, self.b)
        } else if ty.right.global_a() {
            let value = (0..=self.b)
                .map(|v| (v, self.value[self.a][v].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            (value.1, self.a, value.0)
        } else if ty.right.global() {
            let value_a = (0..=self.a)
                .map(|v| (v, self.value[v][self.b].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            let value_b = (0..=self.b)
                .map(|v| (v, self.value[self.a][v].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            if value_a.1 >= value_b.1 {
                (value_a.1, value_a.0, self.b)
            } else {
                (value_b.1, self.a, value_b.0)
            }
        } else {
            high
        }
    }

    /// # Safety
    /// This function assumes the index to be valid. Not upholding this does an out of bounds unsafe [`Vec::get_unchecked`].
    /// A debug assertion hold up this promise on debug builds.
    pub(super) unsafe fn get_unchecked(&self, index: [usize; 2]) -> &Piece {
        debug_assert!(self.value.len() > index[0]);
        debug_assert!(self.value[index[0]].len() > index[1]);
        unsafe { self.value.get_unchecked(index[0]).get_unchecked(index[1]) }
    }

    /// # Safety
    /// This function assumes the index to be valid. Not upholding this does an out of bounds unsafe [`Vec::get_unchecked_mut`].
    /// A debug assertion hold up this promise on debug builds.
    pub(super) unsafe fn get_unchecked_mut(&mut self, index: [usize; 2]) -> &mut Piece {
        debug_assert!(self.value.len() > index[0]);
        debug_assert!(self.value[index[0]].len() > index[1]);
        unsafe {
            self.value
                .get_unchecked_mut(index[0])
                .get_unchecked_mut(index[1])
        }
    }
}

impl std::ops::Index<[usize; 2]> for Matrix {
    type Output = Piece;
    fn index(&self, index: [usize; 2]) -> &Self::Output {
        assert!(index[0] <= self.a + 1);
        assert!(index[1] <= self.b + 1);
        &self.value[index[0]][index[1]]
    }
}
impl std::ops::IndexMut<[usize; 2]> for Matrix {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        assert!(index[0] <= self.a + 1);
        assert!(index[1] <= self.b + 1);
        &mut self.value[index[0]][index[1]]
    }
}
