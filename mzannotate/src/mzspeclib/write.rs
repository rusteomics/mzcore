use std::io::{self, prelude::*};

use crate::{mzspeclib::LibraryHeader, spectrum::AnnotatedSpectrum};

/// A wrapper around a writer to write out annotated spectra as mzSpecLib files.
#[derive(Debug)]
pub struct MzSpecLibTextWriter<W: Write> {
    writer: io::BufWriter<W>,
    header: LibraryHeader,
}

impl<W: Write> MzSpecLibTextWriter<W> {
    /// Create a new writer
    pub fn new(writer: io::BufWriter<W>) -> Self {
        let header: LibraryHeader = LibraryHeader::default();
        Self { writer, header }
    }

    /// Get the library header that contains all the metadata for this library
    pub const fn header(&self) -> &LibraryHeader {
        &self.header
    }

    /// Get mutable access to the library header that contains the metadata for this library
    pub const fn header_mut(&mut self) -> &mut LibraryHeader {
        &mut self.header
    }

    /// Write the header to the write stream.
    /// # Errors
    /// If writing to the underlying stream failed.
    // TODO: document in which order these have to be called, maybe look into a type state pattern to enforce this
    pub fn write_header(&mut self) -> io::Result<()> {
        self.writer.write_all(self.header.to_string().as_bytes())
    }

    /// Write a spectrum to the stream
    /// # Errors
    /// If writing to the underlying stream failed.
    pub fn write_spectrum(&mut self, spectrum: &AnnotatedSpectrum) -> io::Result<()> {
        self.writer.write_all(spectrum.to_string().as_bytes())
    }
}
