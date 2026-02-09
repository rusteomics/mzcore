//! Decompression for .Z files.
//! Code was reused from unarc-rs by mkrueger (Mike Krüger) December 2025 licensed under MIT OR Apache-2.
//! https://github.com/mkrueger/unarc-rs
//! It was not added as dependency because only the .Z is needed and the unarc-rs has a lot of supported formats with a lot of dependencies.

use std::io::Read;

pub(crate) enum ArchiveError {
    DecompressionFailed { entry: String, reason: String },
    InvalidHeader,
    IO(std::io::Error),
}

impl From<std::io::Error> for ArchiveError {
    fn from(value: std::io::Error) -> Self {
        Self::IO(value)
    }
}

mod lzw;
pub(crate) struct ZArchive<T: Read> {
    block_mode: bool,
    max_bits: u8,
    reader: T,
}

pub(crate) const ID: [u8; 2] = [0x1F, 0x9D];
const BLOCK_MODE: u8 = 0x80;
const BIT_MASK: u8 = 0x1f;

impl<T: Read> ZArchive<T> {
    pub(crate) fn new(mut reader: T) -> Result<Self, ArchiveError> {
        let mut header = [0; 3];
        reader.read_exact(&mut header)?;
        if header[0..2] != ID {
            return Err(ArchiveError::InvalidHeader);
        }
        let block_mode = header[2] & BLOCK_MODE != 0;
        let max_bits = header[2] & BIT_MASK;
        Ok(Self {
            block_mode,
            max_bits,
            reader,
        })
    }

    pub(crate) fn read(&mut self) -> Result<Vec<u8>, ArchiveError> {
        let mut compressed_buffer = Vec::new();
        self.reader.read_to_end(&mut compressed_buffer)?;
        let decompressed =
            lzw::Lzw::new(self.max_bits, self.block_mode).decomp(&compressed_buffer)?;
        Ok(decompressed)
    }
}
