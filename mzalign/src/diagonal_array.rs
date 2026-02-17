use std::fmt::{Debug, Write};

/// A diagonal array of limited depth that is implemented as a single continuous slice of memory.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(super) struct DiagonalArray<T, const DEPTH: u16> {
    len: usize,
    data: Box<[T]>,
}

impl<T, const DEPTH: u16> DiagonalArray<T, DEPTH> {
    /// Calculate the index of a given point (along the first axis; n) into the array
    const fn length(n: usize) -> usize {
        Self::length_with_depth(n, DEPTH as usize)
    }

    /// Calculate the index of a given point (along the first axis; n) into the array, with a given depth
    const fn length_with_depth(n: usize, depth: usize) -> usize {
        let mi = if n >= depth { depth } else { n }; // Ord::min is not const
        (mi + 1) * mi / 2 + n.saturating_sub(depth) * depth
    }

    /// Calculate the index of a given point (along the first axis; n) into the array, with a given depth
    fn iter(&self) -> impl Iterator<Item = (&T, [usize; 2])> + '_ {
        let mut index = [0, 0];
        self.data.iter().take(Self::length(self.len)).map(move |v| {
            let i = index;
            if index[0] == index[1] {
                index = [index[0] + 1, 0];
            } else {
                index = [index[0], index[1] + 1];
            }
            (v, i)
        })
    }

    /// # Panics
    /// When the indices are not valid
    fn validate_indices(&self, index: [usize; 2]) -> bool {
        assert!(
            index[0] < self.len,
            "First index {} is outside of diagonal array with length {}",
            index[0],
            self.len
        );
        assert!(
            index[1] <= index[0] || index[1] <= DEPTH as usize,
            "Second index {} is outside of diagonal array with length {} at first index {}",
            index[1],
            self.len,
            index[0],
        );
        true
    }

    /// # Safety
    /// This function assumes the index to be valid. Not upholding this does an out of bounds unsafe [`[T]::get_unchecked`].
    /// A debug assertion hold up this promise on debug builds.
    pub(super) unsafe fn get_unchecked(&self, index: [usize; 2]) -> &T {
        debug_assert!(self.validate_indices(index));
        let index = Self::length(index[0]) + index[1];
        unsafe { self.data.get_unchecked(index) }
    }
}

impl<const DEPTH: u16> DiagonalArray<f64, DEPTH> {
    /// Get the minimal value, ignoring the diagonal
    /// # Panics
    /// If the diagonal array is empty
    pub(crate) fn min<const IGNORE_DIAGONAL: bool>(&self) -> (f64, [usize; 2]) {
        self.iter()
            .filter(|(_, [o, i])| o != i || !IGNORE_DIAGONAL)
            .min_by(|(a, _), (b, _)| a.total_cmp(b))
            .map(|(v, i)| (*v, i))
            .expect("Empty diagonal arrays are not valid")
    }

    // Specifically for neighbour joining.
    // This merges the indicated row+column indices and repacks the data to remove the far row+column this does not trim the underlying boxed slice.
    pub(crate) fn merge_columns(
        mut self,
        close_index: usize,
        far_index: usize,
        distance: f64,
    ) -> Self {
        debug_assert!(self.len > 1);
        debug_assert!(close_index < far_index);
        for row in 0..close_index {
            self[[close_index, row]] =
                (self[[close_index, row]] + self[[far_index, row]] - distance) / 2.0;
        }
        for column in close_index + 1..self.len {
            self[[column, close_index]] = (self[[column, close_index]]
                + self[[column.max(far_index), column.min(far_index)]]
                - distance)
                / 2.0;
        }
        // Remove unneeded cells
        let mut index = [0, 0];
        let mut linear_index = 0;

        for i in 0..Self::length(self.len) {
            if index[0] != far_index && index[1] != far_index {
                self.data[linear_index] = self.data[i];
                linear_index += 1;
            }
            if index[0] == index[1] {
                index = [index[0] + 1, 0];
            } else {
                index = [index[0], index[1] + 1];
            }
        }

        self.len -= 1;
        self
    }
}

impl<T: std::fmt::Display, const DEPTH: u16> DiagonalArray<T, DEPTH> {
    #[allow(dead_code)] // Used for debugging purposes
    pub(crate) fn to_csv(&self) -> String {
        let mut output = String::new();
        let depth = (DEPTH as usize).min(self.len);
        for r in 0..depth {
            for column in 0..self.len {
                if r <= column {
                    write!(&mut output, "\t{:.3}", self[[column, r]]).unwrap();
                } else {
                    write!(&mut output, "\t.").unwrap();
                }
            }
            writeln!(&mut output).unwrap();
        }
        output
    }
}

impl<T: Default + Clone, const DEPTH: u16> DiagonalArray<T, DEPTH> {
    /// Create a new diagonal array of the correct size, with all values initialised to the default value of the type
    pub(super) fn new(len: usize) -> Self {
        Self {
            len,
            data: vec![
                T::default();
                Self::length_with_depth(len, DEPTH.saturating_add(1) as usize)
            ]
            .into(),
        }
    }
}

impl<T, const DEPTH: u16> std::ops::Index<[usize; 2]> for DiagonalArray<T, DEPTH> {
    type Output = T;
    /// Index into the diagonal array
    fn index(&self, index: [usize; 2]) -> &Self::Output {
        assert!(self.validate_indices(index));
        let index = Self::length(index[0]) + index[1];
        &self.data[index]
    }
}

impl<T, const DEPTH: u16> std::ops::IndexMut<[usize; 2]> for DiagonalArray<T, DEPTH> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        assert!(self.validate_indices(index));
        let index = Self::length(index[0]) + index[1];
        &mut self.data[index]
    }
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod tests {
    use super::DiagonalArray;

    #[test]
    fn create() {
        let mut array = DiagonalArray::<i8, 2>::new(2);
        array[[0, 0]] = 1;
        array[[1, 0]] = 2;
        array[[1, 1]] = 3;
        assert_eq!(array[[0, 0]], 1);
        assert_eq!(array[[1, 0]], 2);
        assert_eq!(array[[1, 1]], 3);
    }
}
