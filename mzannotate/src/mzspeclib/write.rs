use std::{
    io::{self, prelude::*},
    marker::PhantomData,
};

use mzdata::{mzpeaks::peak_set::PeakSetIter, params::Value, prelude::*};

use crate::{
    mzspeclib::{Attribute, AttributeValue, Attributes, EntryType, Id, LibraryHeader},
    prelude::ToMzPAF,
    term,
};

use itertools::Itertools;

/// A wrapper around a writer to write out annotated spectra as mzSpecLib.txt files.
///
/// To create a valid mzSpecLib file first the header has to be written before any spectra can be
/// written. The written spectra will check the library header to not write duplicate values to the
/// file.
///
/// ```
/// use mzannotate::prelude::*;
/// // Create the file and open the writer, using a bufwriter is advised
/// let file = std::fs::File::create("../data/test.mzSpecLib.txt").unwrap();
/// let mut writer = MzSpecLibTextWriter::new(std::io::BufWriter::new(file));
/// // Set any parameters needed in the header
/// writer.header_mut().attributes[0].push(mzannotate::mzspeclib::Attribute::new(
///     term!(MS:1003188|library name),
///     mzdata::params::Value::String("Simple test library".to_string())));
/// // Write the header, note that this gives you the writer back to prevent
/// // writing the header twice or not writing the header at all.
/// let mut writer = writer.write_header().unwrap();
/// // Write the spectra, either one by one or a bigger group in one go
/// writer.write_spectrum(&AnnotatedSpectrum::default()).unwrap();
/// writer.write_spectra(&[AnnotatedSpectrum::default()]).unwrap();
/// ```
#[derive(Debug)]
pub struct MzSpecLibTextWriter<Writer: Write, State> {
    writer: Writer,
    header: LibraryHeader,
    state: PhantomData<State>,
}

impl<Writer: Write, State> MzSpecLibTextWriter<Writer, State> {
    /// Get the library header that contains all the metadata for this library
    pub const fn header(&self) -> &LibraryHeader {
        &self.header
    }
}

impl<Writer: Write> MzSpecLibTextWriter<Writer, Initial> {
    /// Create a new writer. It is often advisable to wrap a writer in [`io::BufWriter`] to get
    /// buffered writing.
    pub fn new(writer: Writer) -> Self {
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
    pub fn write_header(mut self) -> io::Result<MzSpecLibTextWriter<Writer, HeaderWritten>> {
        let version = Attribute::new(
            term!(MS:1003186|library format version),
            AttributeValue::Scalar(Value::String(self.header.format_version.clone())),
        );
        writeln!(&mut self.writer, "<mzSpecLib>\n{version}")?;
        for (id, group) in self.header.attributes.iter().enumerate() {
            for attr in group {
                if let Some(id) = id.checked_sub(1) {
                    writeln!(&mut self.writer, "[{id}]{attr}")?;
                } else {
                    writeln!(&mut self.writer, "{attr}")?;
                }
            }
        }
        for t in [
            &EntryType::Spectrum,
            &EntryType::Analyte,
            &EntryType::Interpretation,
            &EntryType::Cluster,
        ] {
            for group in self
                .header
                .attribute_classes
                .get(t)
                .iter()
                .flat_map(|s| s.iter())
            {
                writeln!(&mut self.writer, "<AttributeSet {t}={}>", group.id)?;
                for (id, group) in &group.attributes {
                    for (attr, _) in group {
                        if let Some(id) = id {
                            writeln!(&mut self.writer, "[{id}]{attr}")?;
                        } else {
                            writeln!(&mut self.writer, "{attr}")?;
                        }
                    }
                }
            }
        }
        Ok(MzSpecLibTextWriter::<Writer, HeaderWritten> {
            writer: self.writer,
            header: self.header,
            state: PhantomData,
        })
    }

    // TODO: create function that takes the header and a bunch of spectra and then searches for
    // common terms+values across all spectra for more efficient encoding of the mzSpecLib file.
}

impl<Writer: Write> MzSpecLibTextWriter<Writer, HeaderWritten> {
    /// Write a spectrum to the stream. This can only be done once the header is written.
    /// # Errors
    /// If writing to the underlying stream failed.
    pub fn write_spectrum<Spectrum: MzSpecLibEncode>(
        &mut self,
        spectrum: &Spectrum,
    ) -> io::Result<()> {
        // TODO: think about uniqueness guarantees for this key
        writeln!(&mut self.writer, "<Spectrum={}>", spectrum.key())?;

        for (id, group) in spectrum.spectrum().iter().enumerate() {
            for attr in group {
                // Check if this attribute needs to be written.
                if !self
                    .header
                    .is_already_defined(&attr, EntryType::Spectrum, &["all"])
                {
                    if let Some(id) = id.checked_sub(1) {
                        writeln!(&mut self.writer, "[{id}]{attr}")?;
                    } else {
                        writeln!(&mut self.writer, "{attr}")?;
                    }
                }
            }
        }
        for (id, attributes) in spectrum.analytes() {
            writeln!(&mut self.writer, "<Analyte={id}>")?;
            for (id, group) in attributes.iter().enumerate() {
                for attr in group {
                    if !self
                        .header
                        .is_already_defined(&attr, EntryType::Analyte, &["all"])
                    // TODO: check the protein groups additional set
                    {
                        if let Some(id) = id.checked_sub(1) {
                            writeln!(&mut self.writer, "[{id}]{attr}")?;
                        } else {
                            writeln!(&mut self.writer, "{attr}")?;
                        }
                    }
                }
            }
        }
        for (id, attributes, members) in spectrum.interpretations() {
            writeln!(&mut self.writer, "<Interpretation={id}>")?;
            for (id, group) in attributes.iter().enumerate() {
                for attr in group {
                    if !self
                        .header
                        .is_already_defined(&attr, EntryType::Interpretation, &["all"])
                    {
                        if let Some(id) = id.checked_sub(1) {
                            writeln!(&mut self.writer, "[{id}]{attr}")?;
                        } else {
                            writeln!(&mut self.writer, "{attr}")?;
                        }
                    }
                }
            }
            for (key, attributes) in members {
                writeln!(&mut self.writer, "<InterpretationMember={key}>",)?;
                for (id, group) in attributes.iter().enumerate() {
                    for attr in group {
                        if let Some(id) = id.checked_sub(1) {
                            writeln!(&mut self.writer, "[{id}]{attr}")?;
                        } else {
                            writeln!(&mut self.writer, "{attr}")?;
                        }
                    }
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
                // If there are no annotations insert an empty cell
                if p.annotations().next().is_none() {
                    write!(&mut self.writer, "\t")?;
                }
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
    pub fn write_spectra<'a, Spectrum: MzSpecLibEncode + 'a>(
        &mut self,
        spectra: impl IntoIterator<Item = &'a Spectrum>,
    ) -> io::Result<()> {
        for spectrum in spectra {
            self.write_spectrum(spectrum)?;
        }
        Ok(())
    }
}

/// The mzSpecLib file has been started but nothing has been written yet. First a header has to be
/// written before any spectra can be written.
#[allow(missing_debug_implementations, missing_copy_implementations)] // Marker ZST
pub struct Initial;
/// The mzSpecLib file has already been written with a header. Now spectra can be written.
#[allow(missing_debug_implementations, missing_copy_implementations)] // Marker ZST
pub struct HeaderWritten;

/// A single spectrum that can be encoded as an mzSpecLib file
pub trait MzSpecLibEncode {
    /// The peak type
    type P: MzSpecLibPeakEncode;
    /// The key for this spectrum
    fn key(&self) -> Id;
    /// The attributes for this spectrum
    fn spectrum(&self) -> Attributes;
    /// The attributes for the analytes
    fn analytes(&self) -> impl Iterator<Item = (Id, Attributes)>;
    /// The interpretation members iterator
    type InterpretationMemberIter: IntoIterator<Item = (Id, Attributes)>;
    /// The attributes for the interpretations
    fn interpretations(
        &self,
    ) -> impl Iterator<Item = (Id, Attributes, Self::InterpretationMemberIter)>;
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
