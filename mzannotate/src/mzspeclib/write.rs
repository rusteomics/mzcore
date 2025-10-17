use std::{
    io::{self, prelude::*},
    marker::PhantomData,
};

use mzdata::{mzpeaks::peak_set::PeakSetIter, prelude::*};

use crate::{
    mzspeclib::{Attribute, Id, LibraryHeader},
    prelude::ToMzPAF,
};

use itertools::Itertools;

/// A wrapper around a writer to write out annotated spectra as mzSpecLib files.
#[derive(Debug)]
pub struct MzSpecLibTextWriter<W: Write, S> {
    writer: W,
    header: LibraryHeader,
    state: PhantomData<S>,
}

impl<W: Write, S> MzSpecLibTextWriter<W, S> {
    /// Get the library header that contains all the metadata for this library
    pub const fn header(&self) -> &LibraryHeader {
        &self.header
    }
}

impl<W: Write> MzSpecLibTextWriter<W, Initial> {
    /// Create a new writer. It is often advisable to wrap a writer in [`io::BufWriter`] to get
    /// buffered writing.
    pub fn new(writer: W) -> Self {
        let header: LibraryHeader = LibraryHeader::default();
        Self {
            writer,
            header,
            state: PhantomData,
        }
    }

    /// Get mutable access to the library header that contains the metadata for this library
    pub const fn header_mut(&mut self) -> &mut LibraryHeader {
        &mut self.header
    }

    /// Write the header to the write stream. This has to be done before any spectra can be written
    /// and can only be done once.
    /// # Errors
    /// If writing to the underlying stream failed.
    pub fn write_header(mut self) -> io::Result<MzSpecLibTextWriter<W, HeaderWritten>> {
        self.writer.write_all(self.header.to_string().as_bytes())?;
        Ok(MzSpecLibTextWriter::<W, HeaderWritten> {
            writer: self.writer,
            header: self.header,
            state: PhantomData,
        })
    }

    // TODO: create function that takes the header and a bunch of spectra and then searches for
    // common terms+values across all spectra for more efficient encoding of the mzSpecLib file.
}

impl<W: Write> MzSpecLibTextWriter<W, HeaderWritten> {
    /// Write a spectrum to the stream. This can only be done once the header is written.
    /// # Errors
    /// If writing to the underlying stream failed.
    pub fn write_spectrum<S: MzSpecLibEncode>(&mut self, spectrum: &S) -> io::Result<()> {
        // TODO: think about uniqueness guarantees for this key
        writeln!(&mut self.writer, "<Spectrum={}>", spectrum.key())?;

        for attr in spectrum.spectrum() {
            // TODO: check for duplicates in the library header?
            writeln!(&mut self.writer, "{attr}")?;
        }
        for (id, attributes) in spectrum.analytes() {
            writeln!(&mut self.writer, "<Analyte={id}>")?;
            for attr in attributes {
                writeln!(&mut self.writer, "{attr}")?;
            }
        }
        for (id, attributes, members) in spectrum.interpretations() {
            writeln!(&mut self.writer, "<Interpretation={id}>")?;
            for attr in attributes {
                writeln!(&mut self.writer, "{attr}")?;
            }
            for (id, attributes) in members {
                writeln!(&mut self.writer, "<InterpretationMember={id}>",)?;
                for attr in attributes {
                    writeln!(&mut self.writer, "{attr}")?;
                }
            }
        }
        writeln!(&mut self.writer, "<Peaks>")?;
        let mut buffer = String::new();
        for p in spectrum.peaks() {
            write!(&mut self.writer, "{}\t{}", p.mz(), p.intensity())?;

            if p.annotations().next().is_some() {
                write!(&mut self.writer, "\t")?;
                for (i, a) in p.annotations().enumerate() {
                    if i == 0 {
                    } else {
                        write!(&mut self.writer, ",")?;
                    }
                    buffer.clear();
                    a.to_mz_paf(&mut buffer).map_err(io::Error::other)?;
                    self.writer.write_all(buffer.as_bytes())?;
                }
            }
            if p.aggregations().next().is_some() {
                write!(&mut self.writer, "\t")?;
                write!(&mut self.writer, "{}", p.aggregations().join(","))?;
            }
            writeln!(&mut self.writer)?;
        }
        Ok(())
    }

    /// Write spectra to the stream. This can only be done once the header is written.
    /// # Errors
    /// If writing to the underlying stream failed.
    pub fn write_spectra<'a, S: MzSpecLibEncode + 'a>(
        &mut self,
        spectra: impl IntoIterator<Item = &'a S>,
    ) -> io::Result<()> {
        for spectrum in spectra {
            self.write_spectrum(spectrum)?;
        }
        Ok(())
    }
}

/// The mzSpecLib file has been started but nothing has been written yet. First a header has to be
/// written before any spectra can be written.
pub struct Initial;
/// The mzSpecLib file has already been written with a header. Now spectra can be written.
pub struct HeaderWritten;

/// A single spectrum that can be encoded as an mzSpecLib file
pub trait MzSpecLibEncode {
    /// The peak type
    type P: MzSpecLibPeakEncode;
    /// The key for this spectrum
    fn key(&self) -> Id;
    /// The attributes for this spectrum
    fn spectrum(&self) -> impl IntoIterator<Item = Attribute>;
    /// The analyte attribute iterator
    type AnalyteIter: IntoIterator<Item = Attribute>;
    /// The attributes for the analytes
    fn analytes(&self) -> impl Iterator<Item = (Id, Self::AnalyteIter)>;
    /// The interpretation attribute iterator
    type InterpretationIter: IntoIterator<Item = Attribute>;
    /// The interpretation members iterator
    type InterpretationMemberIter: IntoIterator<Item = (Id, Self::InterpretationIter)>;
    /// The attributes for the interpretations
    fn interpretations(
        &self,
    ) -> impl Iterator<Item = (Id, Self::InterpretationIter, Self::InterpretationMemberIter)>;
    /// The peaks
    fn peaks(&self) -> PeakSetIter<'_, Self::P>;
}

/// A peak that can be encoded for use in an mzSpecLib file
pub trait MzSpecLibPeakEncode: CentroidLike {
    /// The annotation type, need to be able to be written as mzPAF
    type A: ToMzPAF;
    /// The annotations
    fn annotations(&self) -> impl Iterator<Item = &Self::A>;
    /// The aggregations
    fn aggregations(&self) -> impl Iterator<Item = &str>;
}
