//! A reader adapter to be both a lines reader and return a hash of the full file afterwards.
//!
//! The implementation is heavily inspired by the Rust standard library.

use std::io::{BufRead, Read};

/// A [`std::io::BufReader`] inspired design that also calculates the Hash of the read file.
#[derive(Debug)]
pub struct HashBufReader<R, H> {
    /// The underlying reader
    r: R,
    /// The underlying hasher
    hash: H,
    /// The buffer for the buf read part
    buf: Vec<u8>,
    /// The current position in the buffer
    pos: usize,
    /// The current amount of data in the buffer
    filled: usize,
}

impl<R: Read, H: sha2::Digest> Read for HashBufReader<R, H> {
    /// Fill the given buffer with data
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        // This read a certain number of bytes from the underlying reader and updates the hasher
        // with those bytes before returning to the caller.
        let res = self.r.read(buf)?;
        self.hash.update(&buf[..res]);
        Ok(res)
    }
}

impl<R: Read, H: sha2::Digest> BufRead for HashBufReader<R, H> {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        // This either gives the left portion of the buffer or if the buffer is exhausted refills
        // the buffer.
        if self.pos >= self.filled {
            self.pos = 0;
            self.filled = 0;
            match self.r.read(&mut self.buf) {
                Ok(r) => {
                    self.filled = r;
                }
                Err(e) => return Err(e),
            }
        }

        Ok(&self.buf[self.pos..self.filled])
    }

    fn consume(&mut self, amount: usize) {
        // Consume actually reads the data, this is also the place to give this data to the hasher.
        let new_pos = std::cmp::min(self.pos + amount, self.filled);
        self.hash.update(&self.buf[self.pos..new_pos]);
        self.pos = new_pos;
    }
}

impl<R: Read, H: sha2::Digest> HashBufReader<R, H> {
    /// Create a new buffered reader that also calculates the hash of the file or stream.
    pub fn new(r: R) -> Self {
        Self {
            r,
            hash: H::new(),
            buf: vec![0; 8192], // Same default buf size as `std::io::BufReader`
            pos: 0,
            filled: 0,
        }
    }

    /// Create a new buffered reader that also calculates the hash of the file or stream.
    pub fn boxed(r: R) -> HashBufReader<Box<dyn Read>, H>
    where
        R: 'static,
    {
        HashBufReader {
            r: Box::new(r),
            hash: H::new(),
            buf: vec![0; 8192], // Same default buf size as `std::io::BufReader`
            pos: 0,
            filled: 0,
        }
    }

    /// Get the hash of the file or stream, this finalises the hasher and so also takes ownership of the [`HashBufReader`]
    pub fn hash(self) -> Vec<u8> {
        self.hash.finalize().to_vec()
    }

    /// Get the next line. It return `None` when the reader is exhausted, and `Some` otherwise.
    pub fn next_line(&mut self) -> Option<std::io::Result<String>> {
        let mut additional_buf = Vec::new();
        match read_until(self, b'\n', &mut additional_buf) {
            Ok(0) => None,
            Ok(_) => {
                if additional_buf.last().copied() == Some(b'\n') {
                    additional_buf.pop();
                    if additional_buf.last().copied() == Some(b'\r') {
                        additional_buf.pop();
                    }
                }
                Some(
                    String::from_utf8(additional_buf)
                        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e)),
                )
            }
            Err(e) => Some(Err(e)),
        }
    }

    /// Get an iterator over the lines in this buf reader. Note that this method takes a mutable
    /// reference to self instead of ownership as done in the standard library
    /// [`std::io::BufReader::lines`]. This allows the caller to use the iterator and then
    /// afterwards still get the hash of the file.
    pub const fn lines(&mut self) -> Lines<'_, R, H> {
        Lines { reader: self }
    }
}

/// An iterator over the lines of an instance of [`HashBufReader`].
///
/// Please see the documentation of [`HashBufReader::lines`] for more details.
#[derive(Debug)]
pub struct Lines<'a, R, H> {
    reader: &'a mut HashBufReader<R, H>,
}

impl<R: Read, H: sha2::Digest> Iterator for Lines<'_, R, H> {
    type Item = std::io::Result<String>;

    fn next(&mut self) -> Option<Self::Item> {
        self.reader.next_line()
    }
}

/// Taken from the standard library [`std::io::Lines`] implementation but adapted for the stable std library interface
/// # Errors
/// When the underlying reader errors
fn read_until<R: BufRead + ?Sized>(
    r: &mut R,
    delim: u8,
    buf: &mut Vec<u8>,
) -> std::io::Result<usize> {
    let mut read = 0;
    loop {
        let (done, used) = {
            let available = match r.fill_buf() {
                Ok(n) => n,
                Err(ref e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(e) => return Err(e),
            };
            if let Some(i) = available.iter().position(|e| *e == delim) {
                buf.extend_from_slice(&available[..=i]);
                (true, i + 1)
            } else {
                buf.extend_from_slice(available);
                (false, available.len())
            }
        };
        r.consume(used);
        read += used;
        if done || used == 0 {
            return Ok(read);
        }
    }
}
