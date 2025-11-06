//! Loading / updating [`CVIndex`] functionality.
//! This is made as a separate module to ensure strict adherence to the internal API for `CVIndex`.
//! This ensures that the indexes of the indices and names and such cannot ever be wrongfully set
//! with code from this file.

use std::{
    io::{BufReader, BufWriter, Write},
    path::Path,
    sync::Arc,
};

use context_error::{BoxedError, Context, CreateError};

use crate::{
    CVCache, CVCompression, CVData, CVError, CVIndex, CVSource, CVVersion,
    hash_buf_reader::HashBufReader,
};

impl<CV: CVSource> CVIndex<CV> {
    /// Build this CV from the standard cache location
    fn load_from_cache() -> Result<Self, BoxedError<'static, CVError>> {
        let path = CV::default_stem().with_extension("bin");
        let file = std::fs::File::open(&path).map_err(|e| {
            BoxedError::new(
                CVError::CacheCouldNotBeOpenend,
                "CV cache file could not be openend",
                e.to_string(),
                Context::none().source(path.to_string_lossy()).to_owned(),
            )
        })?;
        let mut reader = BufReader::new(file);
        let cache: <CV::Data as CVData>::Cache =
            bincode::decode_from_std_read(&mut reader, bincode::config::standard()).map_err(
                |e| {
                    BoxedError::new(
                        CVError::CacheCouldNotBeParsed,
                        "CV cache file could not be parsed",
                        e.to_string(),
                        Context::none().source(path.to_string_lossy()).to_owned(),
                    )
                },
            )?;
        let mut result = Self::empty();
        let (version, data) = cache.deconstruct();
        result.update_skip_rebuilding_cache(data.into_iter().map(Arc::new), version);
        Ok(result)
    }

    /// Store this index in the standard cache. If the CV has
    /// [`CV::AUTOMATICALLY_WRITE_UNCOMPRESSED`] set to true this also writes the standard file.
    /// # Errors
    /// If the file could not be written to.
    fn save_to_cache(&self) -> Result<(), BoxedError<'static, CVError>> {
        self.save_to_cache_at(&CV::default_stem().with_extension("bin"))?;
        if CV::AUTOMATICALLY_WRITE_UNCOMPRESSED {
            self.save_to_file_at(
                &CV::default_stem()
                    .with_extension(CV::files().first().map_or("dat", |f| f.extension)),
            )?;
        }
        Ok(())
    }

    /// Store this index at a certain location.
    /// # Errors
    /// If the file could not be written to.
    pub fn save_to_cache_at(&self, path: &Path) -> Result<(), BoxedError<'static, CVError>> {
        let file = std::fs::File::create(path).map_err(|e| {
            BoxedError::new(
                CVError::CacheCouldNotBeOpenend,
                "CV cache file could not be openend",
                e.to_string(),
                Context::none().source(path.to_string_lossy()).to_owned(),
            )
        })?;
        let mut writer = BufWriter::new(file);
        bincode::encode_into_std_write(
            <CV::Data as CVData>::Cache::construct(
                self.version().clone(),
                self.data().map(Arc::unwrap_or_clone).collect(),
            ),
            &mut writer,
            bincode::config::standard(),
        )
        .map_err(|e| {
            BoxedError::new(
                CVError::CacheCouldNotBeMade,
                "CV cache file could not be made",
                e.to_string(),
                Context::none().source(path.to_string_lossy()).to_owned(),
            )
        })?;
        Ok(())
    }

    /// Store the uncompressed file at a certain location.
    /// # Errors
    /// If the file could not be written to.
    pub fn save_to_file_at(&self, path: &Path) -> Result<(), BoxedError<'static, CVError>> {
        let file = std::fs::File::create(path).map_err(|e| {
            BoxedError::new(
                CVError::FileCouldNotBeOpenend,
                "CV file could not be openend",
                e.to_string(),
                Context::none().source(path.to_string_lossy()).to_owned(),
            )
        })?;
        let writer = BufWriter::new(file);
        CV::write_uncompressed(writer, self.version(), self.data())
    }

    /// Initialise the data the first successful source from the list is returned:
    /// 1. Load from the binary cache.
    /// 2. Load from the default path (first compressed then uncompressed).
    /// 3. Load the static data.
    /// 4. Load an empty CV.
    ///
    /// All errors encountered in previous steps are returned.
    pub fn init() -> (Self, Vec<BoxedError<'static, CVError>>) {
        let mut errors = Vec::new();

        // Load cache
        match Self::load_from_cache() {
            Ok(data) => return (data, errors),
            Err(error) => errors.push(error),
        }

        let mut result = Self::empty();

        // If the cache failed try parsing the actual file, this could work because the cache might
        // have been built with an older version of the software with a different Data definition.
        // Inside update from path the cache is overwritten with the new info.
        match result.update_from_path([], true) {
            Ok(()) => {
                return (result, errors);
            }
            Err(error) => errors.push(error),
        }

        // Load the static data
        if let Some((version, data)) = CV::static_data() {
            result.update_skip_rebuilding_cache(data.into_iter().map(Arc::new), version);
        }

        // Fall back with empty CV
        (result, errors)
    }

    /// Initialise the data from the static data.
    ///
    /// All errors encountered in previous steps are returned.
    pub fn init_static() -> Self {
        let mut result = Self::empty();
        // Load the static data
        if let Some((version, data)) = CV::static_data() {
            result.update_skip_rebuilding_cache(data.into_iter().map(Arc::new), version);
        }
        result
    }

    /// Update the CV based on the given data, empties the CV before replacing all data with the new data.
    /// This additionally tries to save the new data to the cache.
    /// # Errors
    /// If the cache could not be saved to disk.
    pub fn update(
        &mut self,
        version: CVVersion,
        data: impl IntoIterator<Item = Arc<CV::Data>>,
    ) -> Result<(), BoxedError<'static, CVError>> {
        self.update_skip_rebuilding_cache(data, version);

        self.save_to_cache()
    }

    /// Update the CV based on the given data, empties the CV before replacing all data with the new data.
    /// This does not store the updated data to the disk, which might be useful for tests but otherwise
    /// [`Self::update`] should be used.
    /// # Errors
    /// If the cache could not be saved to disk.
    pub fn update_do_not_save_to_disk(
        &mut self,
        version: CVVersion,
        data: impl IntoIterator<Item = Arc<CV::Data>>,
    ) {
        self.update_skip_rebuilding_cache(data, version);
    }

    /// Update the CV from the default path (default name + default extension) or the given path.
    /// If no overwrite path was given this reads the default gzipped path or if that does not
    /// exist the default path. It will detect `.gz` on the overwrite path and do decompression on
    /// the fly (the same as with the default gzipped path). If the file could not be parsed a
    /// logfile will be generated at the location of the given path with the extension `.err`.
    ///
    /// This does update the cache on disk.
    ///
    /// If an overwrite path was used and the parsing was successful that file will be moved to
    /// the default path and be gzipped. This is done to keep the most current version at the
    /// default location.
    ///
    /// # Errors
    /// If the given file (or the default one if not overwritten):
    /// * Does not exist.
    /// * Could not be opened.
    /// * Could not be parsed.
    /// * (only if the filename was overwritten) could not be moved to the default location.
    pub fn update_from_path<'a>(
        &mut self,
        overwrite_path: impl IntoIterator<Item = Option<&'a Path>>,
        remove_original: bool,
    ) -> Result<(), BoxedError<'static, CVError>> {
        type Paths<'b> = (
            bool,
            std::path::PathBuf,
            std::path::PathBuf,
            Option<&'b Path>,
        );
        let stem = CV::default_stem();
        let readers: Vec<(
            HashBufReader<Box<dyn std::io::Read>, sha2::Sha256>,
            Paths<'a>,
        )> = CV::files()
            .iter()
            .zip(overwrite_path.into_iter().chain(std::iter::repeat(None)))
            .map(|(file, overwrite_path)| {
                let default_path = stem.with_extension(file.extension);
                let default_gz_path = stem.with_extension(format!("{}.gz", file.extension));
                let resolved_path = if let Some(path) = overwrite_path {
                    if !path.exists() {
                        return Err(BoxedError::new(
                            CVError::FileDoesNotExist,
                            "Given path does not exist",
                            "The given overwrite path does not exist",
                            Context::none().source(path.to_string_lossy()).to_owned(),
                        ));
                    }
                    path
                } else if default_gz_path.exists() {
                    &default_gz_path
                } else {
                    if !default_path.exists() {
                        return Err(BoxedError::new(
                            CVError::FileDoesNotExist,
                            "Default path does not exist",
                            "The default path does not exist",
                            Context::none()
                                .source(default_gz_path.to_string_lossy())
                                .to_owned(),
                        ));
                    }
                    &default_path
                };

                let compressed = resolved_path
                    .extension()
                    .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"));

                let resolved_path_string = resolved_path.to_string_lossy().to_string();
                let file = std::fs::File::open(resolved_path).map_err(|e| {
                    BoxedError::new(
                        CVError::FileCouldNotBeOpenend,
                        "CV file could not be openend",
                        e.to_string(),
                        Context::none().source(&resolved_path_string).to_owned(),
                    )
                })?;

                Ok((
                    if compressed {
                        HashBufReader::<_, sha2::Sha256>::boxed(flate2::bufread::GzDecoder::new(
                            BufReader::new(file),
                        ))
                    } else {
                        HashBufReader::<_, sha2::Sha256>::boxed(file)
                    },
                    (
                        compressed,
                        resolved_path.to_owned(),
                        default_gz_path.clone(),
                        overwrite_path,
                    ),
                ))
            })
            .collect::<Result<Vec<_>, _>>()?;

        let (readers, paths): (Vec<_>, Vec<_>) = readers.into_iter().unzip();

        let (version, data) = CV::parse(readers.into_iter()).map_err(|errors| {
            store_errors(&stem, &errors);
            BoxedError::small(
                CVError::FileCouldNotBeParsed,
                "CV file could not be parsed",
                "",
            )
            .add_underlying_errors(errors)
        })?;

        // Update the data and cache
        self.update(version, data.map(Arc::new))?;

        for (compressed, resolved_path, default_gz_path, overwrite_path) in paths {
            if let Some(path) = overwrite_path {
                let resolved_path_string = resolved_path.to_string_lossy().to_string();
                // If it all worked and this was an overwrite file store the overwrite file at the default location
                if compressed {
                    std::fs::rename(path, default_gz_path).map_err(|e| {
                        BoxedError::new(
                            CVError::FileCouldNotBeMoved,
                            "CV file could not be moved",
                            e.to_string(),
                            Context::none().source(resolved_path_string).to_owned(),
                        )
                    })?;
                } else {
                    // If the overwrite file was not compressed compress while moving to the new location
                    let source_file = std::fs::File::open(&resolved_path).map_err(|e| {
                        BoxedError::new(
                            CVError::FileCouldNotBeOpenend,
                            "CV file could not be openend",
                            e.to_string(),
                            Context::none().source(&resolved_path_string).to_owned(),
                        )
                    })?;
                    let default_gz_file = std::fs::File::create(&default_gz_path).map_err(|e| {
                        BoxedError::new(
                            CVError::FileCouldNotBeMade,
                            "CV file could not be made",
                            e.to_string(),
                            Context::none()
                                .source(default_gz_path.to_string_lossy())
                                .to_owned()
                                .to_owned(),
                        )
                    })?;
                    std::io::copy(
                        &mut BufReader::new(source_file),
                        &mut BufWriter::new(default_gz_file),
                    )
                    .map_err(|e| {
                        BoxedError::new(
                            CVError::FileCouldNotBeMoved,
                            "CV file could not be moved",
                            e.to_string(),
                            Context::none().source(&resolved_path_string).to_owned(),
                        )
                    })?;
                    if remove_original {
                        std::fs::remove_file(&resolved_path).map_err(|e| {
                            BoxedError::new(
                                CVError::FileCouldNotBeMoved,
                                "Overwrite CV file could not be deleted",
                                e.to_string(),
                                Context::none().source(&resolved_path_string).to_owned(),
                            )
                        })?;
                    }
                }
            }
        }
        Ok(())
    }

    /// Download the CV from the internet. If no overwrite URL was given it uses the default URL
    /// ([`CVSource::cv_url`]) if both are unset it errors. It downloads and compresses this file
    /// to the default location but with `.download` before the default extension. This then calls
    /// [`Self::update_from_path`] on the downloaded file.
    ///
    /// This is a blocking interface.
    ///
    /// # Errors
    ///
    /// * No URL was set with both the overwrite URL and [`CVSource::cv_url`].
    /// * The file to download to could not be made.
    /// * The file could not be downloaded or the status code of the download was not success.
    /// * The downloaded file could not be written to disk.
    ///
    /// Additionally, all errors from [`Self::update_from_path`].
    #[cfg(feature = "http")]
    pub fn update_from_url(
        &mut self,
        overwrite_urls: &[Option<&str>],
    ) -> Result<(), BoxedError<'static, CVError>> {
        let paths = CV::files()
            .iter()
            .zip(overwrite_urls.iter().chain(std::iter::repeat(&None)))
            .map(|(cv_file, overwrite_url)| {
                let url = overwrite_url.or_else(|| cv_file.url).ok_or_else(|| {
                    BoxedError::small(
                        CVError::CVUrlNotSet,
                        "Could not download CV",
                        "No URL was given or set as the default URL of this CV",
                    )
                })?;

                let download_path =
                    CV::default_stem().with_extension(format!("download.{}.gz", cv_file.extension));
                let file = std::fs::File::create(&download_path).map_err(|e| {
                    BoxedError::new(
                        CVError::FileCouldNotBeMade,
                        "CV file could not be made",
                        e.to_string(),
                        Context::none()
                            .source(download_path.to_string_lossy())
                            .to_owned(),
                    )
                })?;
                let mut writer = BufWriter::new(file);

                let url = reqwest::Url::try_from(url).map_err(|e| {
                    BoxedError::new(
                        CVError::CVUrlCouldNotBeRead,
                        "Invalid CV URL",
                        e.to_string(),
                        Context::none().source(url).to_owned(),
                    )
                })?;

                if !url.scheme().starts_with("http") {
                    return Err(BoxedError::new(
                        CVError::CVUrlCouldNotBeRead,
                        "Invalid CV URL",
                        "Only HTTP(s) files can be downloaded",
                        Context::none().source(url.to_string()).to_owned(),
                    ));
                }

                let response = reqwest::blocking::get(url.clone())
                    .map_err(|e| {
                        BoxedError::new(
                            CVError::CVUrlCouldNotBeRead,
                            "Could not download CV",
                            e.to_string(),
                            Context::none().source(url.to_string()).to_owned(),
                        )
                    })?
                    .error_for_status()
                    .map_err(|e| {
                        BoxedError::new(
                            CVError::CVUrlCouldNotBeRead,
                            "Could not download CV",
                            e.to_string(),
                            Context::none().source(url.to_string()).to_owned(),
                        )
                    })?;
                // Decompress (if needed) then compress again to gz
                match cv_file.compression {
                    CVCompression::None => {
                        let mut encoder = flate2::bufread::GzEncoder::new(
                            BufReader::new(response),
                            flate2::Compression::fast(),
                        );
                        std::io::copy(&mut encoder, &mut writer)

                        // response
                        // .copy_to(&mut writer)
                        // // .map(|_| ())
                        // .map_err(|e| e.to_string())

                        // Ok(0)
                    }
                    CVCompression::LZW => {
                        todo!()
                        // TODO: figure out LZW decompression
                        // LZW decompression did not work out
                        // let mut decoder = lzw::Decoder::new(lzw::MsbReader::new(), 8);
                        // let mut reader = BufReader::new(response);

                        // loop {
                        //     let len = {
                        //         let buf = reader.fill_buf().unwrap();
                        //         if buf.is_empty() {
                        //             break;
                        //         }
                        //         let (len, bytes) = decoder.decode_bytes(buf).unwrap();
                        //         writer.write_all(bytes).unwrap();
                        //         len
                        //     };
                        //     reader.consume(len);
                        // }

                        // Ok(())

                        // unlzw::unlzw(response, writer).unwrap();
                        // Ok(0)\
                        // let bytes: Vec<u32> = response
                        //     .bytes()
                        //     .unwrap()
                        //     .as_chunks::<4>()
                        //     .0
                        //     .iter()
                        //     .map(|bytes| u32::from_ne_bytes(*bytes))
                        //     .collect();
                        // let decompressed = lzw_compress::lzw::decompress(&bytes);
                        // writer.write_all(&decompressed).unwrap();
                        // let bytes: Vec<u32> = response
                        //     .bytes()
                        //     .unwrap()
                        //     .as_chunks::<4>()
                        //     .0
                        //     .iter()
                        //     .map(|bytes| u32::from_ne_bytes(*bytes))
                        //     .collect();
                        // let decompressed = lzw_compress::lzw::decompress(&bytes);
                        // std::io::copy(
                        //     &mut unlzw::LZWDecoder::new(BufReader::new(response)),
                        //     &mut writer,
                        // )
                        // .map_err(|e| e.to_string())
                        // let mut archive = unarc_rs::z::ZArchieve::new(response).unwrap();
                        // writer.write_all(&archive.read().unwrap()).unwrap();
                        // Ok(0)
                    }
                }
                .map_err(|e| {
                    BoxedError::new(
                        CVError::FileCouldNotBeMade,
                        "Could not download the CV file",
                        e.to_string(),
                        Context::none()
                            .source(url.to_string())
                            .lines(0, download_path.to_string_lossy())
                            .to_owned(),
                    )
                })?;
                Ok(download_path)
            })
            .collect::<Result<Vec<_>, BoxedError<'static, CVError>>>()?;

        self.update_from_path(paths.iter().map(|p| Some(Path::new(p))), true)
    }
}

/// Store these errors in a file, the path will be updated with a new extension (.err).
/// If writing failed in some way this is returned as `false`.
fn store_errors(path: &Path, errors: &[BoxedError<'_, CVError>]) -> bool {
    let err_path = path.with_extension("err");

    std::fs::File::create(&err_path).is_ok_and(|file| {
        let mut success = true;
        let mut writer = BufWriter::new(file);
        success &= writeln!(
            &mut writer,
            "Time: {}\nPath: {}",
            chrono::Local::now().to_rfc3339(),
            path.display()
        )
        .is_ok();

        for error in errors {
            success &= writeln!(&mut writer, "{error}").is_ok();
        }

        success
    })
}
