use std::io::{self, prelude::*};

use crate::mzspeclib::{LibraryHeader, LibrarySpectrum};

pub struct MzSpecLibTextWriter<W: Write> {
    writer: io::BufWriter<W>,
    header: LibraryHeader,
}

impl<W: Write> MzSpecLibTextWriter<W> {
    pub fn new(writer: io::BufWriter<W>) -> Self {
        let header: LibraryHeader = LibraryHeader::default();
        Self { writer, header }
    }

    pub fn header(&self) -> &LibraryHeader {
        &self.header
    }

    pub fn header_mut(&mut self) -> &mut LibraryHeader {
        &mut self.header
    }

    pub fn write_header(&mut self) -> io::Result<()> {
        self.writer.write_all(self.header.to_string().as_bytes())?;
        Ok(())
    }

    pub fn write_spectrum(&mut self, spectrum: &LibrarySpectrum) -> io::Result<()> {
        self.writer.write_all(spectrum.to_string().as_bytes())?;
        Ok(())
    }
}
