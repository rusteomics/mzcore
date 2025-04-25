/// A possibly limited diagonal array that is implemented as a single continuous slice of memory.
/// It consists of a
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) struct DiagonalArray<T> {
    len: usize,
    max_depth: usize,
    data: Box<[T]>,
}

impl<T> DiagonalArray<T> {
    /// Calculate the index of a given point (along the first axis; n) into the array with the given `max_depth` (m)
    const fn length(n: usize, m: usize) -> usize {
        let mi = if n >= m { m } else { n };
        (mi + 1) * mi / 2 + n.saturating_sub(m) * m
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
            index[1] <= index[0] || index[1] <= self.max_depth,
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
        let index = Self::length(index[0], self.max_depth) + index[1];
        unsafe { self.data.get_unchecked(index) }
    }

    /// # Safety
    /// This function assumes the index to be valid. Not upholding this does an out of bounds unsafe [`[T]::get_unchecked_mut`].
    /// A debug assertion hold up this promise on debug builds.
    #[expect(dead_code)]
    pub(super) unsafe fn get_unchecked_mut(&mut self, index: [usize; 2]) -> &mut T {
        debug_assert!(self.validate_indices(index));
        let index = Self::length(index[0], self.max_depth) + index[1];
        unsafe { self.data.get_unchecked_mut(index) }
    }
}

impl<T: Default + Clone> DiagonalArray<T> {
    /// Create a new diagonal array of the correct size, with all values initialised to the default value of the type, with up to and including the depth given in `max_depth`
    pub(super) fn new(len: usize, max_depth: u16) -> Self {
        Self {
            len,
            max_depth: max_depth as usize,
            data: vec![T::default(); Self::length(len, (max_depth as usize).saturating_add(1))]
                .into(),
        }
    }
}

impl<T> std::ops::Index<[usize; 2]> for DiagonalArray<T> {
    type Output = T;
    /// Index into the diagonal array
    fn index(&self, index: [usize; 2]) -> &Self::Output {
        assert!(self.validate_indices(index));
        let index = Self::length(index[0], self.max_depth) + index[1];
        &self.data[index]
    }
}

impl<T> std::ops::IndexMut<[usize; 2]> for DiagonalArray<T> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        assert!(self.validate_indices(index));
        let index = Self::length(index[0], self.max_depth) + index[1];
        &mut self.data[index]
    }
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod tests {
    use super::DiagonalArray;

    #[test]
    fn create() {
        let mut array = DiagonalArray::<i8>::new(2, 2);
        array[[0, 0]] = 1;
        array[[1, 0]] = 2;
        array[[1, 1]] = 3;
        assert_eq!(array[[0, 0]], 1);
        assert_eq!(array[[1, 0]], 2);
        assert_eq!(array[[1, 1]], 3);
    }
}
