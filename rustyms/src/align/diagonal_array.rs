/// A diagonal array of limited depth that is implemented as a single continuous slice of memory.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(super) struct DiagonalArray<T, const DEPTH: u16> {
    len: usize,
    data: Box<[T]>,
}

impl<T, const DEPTH: u16> DiagonalArray<T, DEPTH> {
    /// Calculate the index of a given point (along the first axis; n) into the array
    const fn length(n: usize) -> usize {
        let mi = if n >= DEPTH as usize {
            DEPTH as usize
        } else {
            n
        };
        (mi + 1) * mi / 2 + n.saturating_sub(DEPTH as usize) * DEPTH as usize
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

    /// # Safety
    /// This function assumes the index to be valid. Not upholding this does an out of bounds unsafe [`[T]::get_unchecked_mut`].
    /// A debug assertion hold up this promise on debug builds.
    #[expect(dead_code)]
    pub(super) unsafe fn get_unchecked_mut(&mut self, index: [usize; 2]) -> &mut T {
        debug_assert!(self.validate_indices(index));
        let index = Self::length(index[0]) + index[1];
        unsafe { self.data.get_unchecked_mut(index) }
    }
}

impl<T: Default + Clone, const DEPTH: u16> DiagonalArray<T, DEPTH> {
    /// Create a new diagonal array of the correct size, with all values initialised to the default value of the type
    pub(super) fn new(len: usize) -> Self {
        Self {
            len,
            data: vec![T::default(); Self::length(len)].into(),
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
