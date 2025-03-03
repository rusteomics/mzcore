#![allow(clippy::missing_panics_doc)]
use base64::Engine;
use swash::{scale::ScaleContext, CacheKey, FontRef};

use crate::{
    fragment::GlycanPosition,
    glycan::{GlycanDirection, GlycanStructure},
};
use std::{
    fmt::Write,
    io::BufWriter,
    path::{Path, PathBuf},
};

pub struct Font {
    // Full content of the font file
    data: Vec<u8>,
    // Offset to the table directory
    offset: u32,
    // Cache key
    key: CacheKey,
}

impl Font {
    pub fn from_file(path: PathBuf, index: usize) -> Option<Self> {
        // Read the full font file
        let data = std::fs::read(path).ok()?;
        // Create a temporary font reference for the first font in the file.
        // This will do some basic validation, compute the necessary offset
        // and generate a fresh cache key for us.
        let font = FontRef::from_index(&data, index)?;
        let (offset, key) = (font.offset, font.key);
        // Return our struct with the original file data and copies of the
        // offset and key from the font reference
        Some(Self { data, offset, key })
    }

    // Create the transient font reference for accessing this crate's
    // functionality.
    pub fn as_ref(&self) -> FontRef {
        // Note that you'll want to initialize the struct directly here as
        // using any of the FontRef constructors will generate a new key which,
        // while completely safe, will nullify the performance optimizations of
        // the caching mechanisms used in this crate.
        FontRef {
            data: &self.data,
            offset: self.offset,
            key: self.key,
        }
    }
}

#[test]
fn test_rendering() {
    const COLUMN_SIZE: f32 = 30.0;
    const SUGAR_SIZE: f32 = 15.0;
    const STROKE_SIZE: f32 = 1.5;

    let font = Font::from_file(
        std::fs::read_dir(
            directories::UserDirs::font_dir(
                &directories::UserDirs::new().expect("Could not find user directories"),
            )
            .unwrap_or_else(|| Path::new("C:/WINDOWS/Fonts")), // Font directory not defined for windows
        )
        .expect("Could not open font directory")
        .find(|p| {
            p.as_ref()
                .is_ok_and(|p| p.file_name().eq_ignore_ascii_case("times.ttf"))
        })
        .expect("No font files")
        .expect("Could not open font file")
        .path(),
        0,
    )
    .expect("Invalid font");

    let mut html = String::new();
    let mut footnotes = Vec::new();
    write!(&mut html, "<html lang=\"en\"><head><meta charset=\"UTF-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\"><title>Glycan render test</title></head><body>").unwrap();

    let codes = [
        ("G01670UQ", "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"),
        ("G13523IF", "Fuc(?1-?)Gal(?1-?)GalNAc(?1-"),
        ("G00613DO", "GlcN(b1-4)GlcNAc(b1-4)GlcNAc(b1-4)GlcNAc6S(?1-"),
        ("G00621IU", "Neu5Gc(a2-3/6)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)[Gal(a1-3)Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[Neu5Ac(a2-8)Neu5Ac(a2-3/6)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-2)[Neu5Ac(a2-3/6)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(?1-"),
        ("G01464QV", "Rha2,3,4Ac3(a1-2)[Xyl(b1-3)]Ara(a1-"),
        ("G04421VO", "Fruf(b2-1a)[Glc(a1-2)Glc(a1-2)Glc(a1-2)Glc(a1-2)Glc(a1-2)Glc(a1-2)]Glc"),
        ("G04458LN", "Kdn(a2-3)Gal(b1-4)ManNAc(b1-2)[Kdn(a2-3)Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[GlcNAc(b1-4)][Kdn(a2-3)Gal(b1-4)GlcNAc(b1-2)[Neu5Gc(a2-3)Gal(b1-4)GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-"),
        ("G69524KC", "Xyl(?1-?)Ara(?1-?)[Gal(?1-?)]GlcA"),
        ("G37707YH", "Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)[Gal(a1-3)Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[GlcNAc(b1-4)][Neu5Gc(a2-3/6)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-2)[Neu5Ac(a2-3/6)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(?1-"),
        ("G07370RP", "Rha(a1-3)Qui(b1-4)Rha(a1-2)Glc(b1-2)[Rha(a1-6)]Glc(b1-"),
        ("G11504PZ", "Dig3CMe(b1-3)Oli(b1-3)Oli(b1-"),
        ("G64699IM", "GlcA(b1-3)GalNAc(b1-4)4eLeg?5,7Ac2(a2-"),
        ("G14402AU", "D-Araf(b1-5)Dha(?2-3)[GalA(a1-4)GalA(a1-4)]GalA(a1-4)GalA"),
        ("G08395BZ", "Glc(b1-2a)[Ido(b1-3)]Psif"),
        ("G49642ZT", "Man(?1-?)[Man(?1-?)]Man(?1-?)[Man(?1-?)]Man(?1-?)GlcNAc(?1-?)[Fuc(?1-?)][Fuc(?1-?)]GlcNAc(?1-"),
        ("G59426OB", "Hex(?1-?)HexNAc(?1-?)HexA(?1-?)Gal(?1-?)GalNAc-ol"),
        ("G75424NV", "Hex?(?1-?)Hex?NAc(?1-?)[Hex?NAc(?1-?)]Hex?(?1-?)[Hex?(?1-?)[Hex?(?1-?)]Hex?(?1-?)][Hex?NAc(?1-?)]Hex?(?1-?)Hex?NAc(?1-?)Hex?NAc(?1-"),
        ("G36128WO", "Ido(b1-3)ManNAc(?1-3)[Ido(b1-3)L-AllNAc(b1-3)Ido(b1-4)AltNAc(b1-6)]Tal(b1-4)D-Ido(?1-"),
        ("G83422GV", "L-6dTal(a1-3)[Fuc(a1-2)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)]GlcNAc(b1-3)Gal(b1-3)[Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]GalNAc(a1-"),
        ("G09073GJ","GalNAc(?1-?)GlcA2,3NAc2(?1-?)D-FucNAc"),
        ("G00069DT","Neu(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)Glc(b1-"),
    ];

    let mut context = ScaleContext::new();
    for (index, (_, iupac)) in codes.iter().enumerate() {
        let structure = GlycanStructure::from_short_iupac(iupac, 0..iupac.len(), 0).unwrap();
        let rendered = structure
            .render(
                Some("pep".to_string()),
                COLUMN_SIZE,
                SUGAR_SIZE,
                STROKE_SIZE,
                if index % 3 == 0 {
                    GlycanDirection::LeftToRight
                } else {
                    GlycanDirection::TopDown
                },
                None,
                &[],
                [66, 66, 66],
                [255, 255, 255],
                &mut footnotes,
            )
            .unwrap();
        rendered.to_svg(&mut html).unwrap();
        let (bitmap, width) = rendered.to_bitmap(
            if index % 2 == 0 {
                zeno::Format::subpixel_bgra()
            } else {
                zeno::Format::Alpha
            },
            font.as_ref(),
            &mut context,
        );
        let mut buffer = Vec::new();
        let mut w = BufWriter::new(&mut buffer);
        let mut encoder =
            png::Encoder::new(&mut w, width as u32, (bitmap.len() / 4 / width) as u32);
        encoder.set_color(png::ColorType::Rgba);
        encoder.set_depth(png::BitDepth::Eight);
        let mut writer = encoder.write_header().unwrap();
        writer.write_image_data(&bitmap).unwrap();
        drop(writer);
        drop(w);

        write!(&mut html, "<img src=\"data:image/png;base64, ").unwrap();
        base64::engine::general_purpose::STANDARD.encode_string(&buffer, &mut html);
        write!(&mut html, "\"/>").unwrap();
    }

    for (_, iupac) in &codes {
        let structure = GlycanStructure::from_short_iupac(iupac, 0..iupac.len(), 0).unwrap();
        structure
            .render(
                Some("pep".to_string()),
                COLUMN_SIZE,
                SUGAR_SIZE,
                STROKE_SIZE,
                GlycanDirection::TopDown,
                None,
                &[],
                [66, 66, 66],
                [255, 255, 255],
                &mut footnotes,
            )
            .unwrap()
            .to_svg(&mut html)
            .unwrap();
    }

    for (index, root, breaks) in [
        (
            0,
            Some(GlycanPosition {
                inner_depth: 2,
                series_number: 2,
                branch: Vec::new(),
                attachment: None,
            }),
            Vec::new(),
        ),
        (
            0,
            Some(GlycanPosition {
                inner_depth: 2,
                series_number: 2,
                branch: Vec::new(),
                attachment: None,
            }),
            vec![GlycanPosition {
                inner_depth: 4,
                series_number: 4,
                branch: vec![1],
                attachment: None,
            }],
        ),
        (
            0,
            Some(GlycanPosition {
                inner_depth: 2,
                series_number: 2,
                branch: Vec::new(),
                attachment: None,
            }),
            vec![
                GlycanPosition {
                    inner_depth: 5,
                    series_number: 5,
                    branch: vec![0],
                    attachment: None,
                },
                GlycanPosition {
                    inner_depth: 3,
                    series_number: 3,
                    branch: vec![1],
                    attachment: None,
                },
            ],
        ),
        (
            0,
            Some(GlycanPosition {
                inner_depth: 4,
                series_number: 4,
                branch: vec![1],
                attachment: None,
            }),
            Vec::new(),
        ),
        (
            14,
            Some(GlycanPosition {
                inner_depth: 0,
                series_number: 0,
                branch: Vec::new(),
                attachment: None,
            }),
            vec![GlycanPosition {
                inner_depth: 1,
                series_number: 1,
                branch: vec![0],
                attachment: None,
            }],
        ),
        (
            14,
            Some(GlycanPosition {
                inner_depth: 1,
                series_number: 1,
                branch: vec![0],
                attachment: None,
            }),
            vec![GlycanPosition {
                inner_depth: 2,
                series_number: 2,
                branch: vec![0],
                attachment: None,
            }],
        ),
        (
            16,
            Some(GlycanPosition {
                inner_depth: 1,
                series_number: 1,
                branch: Vec::new(),
                attachment: None,
            }),
            vec![
                GlycanPosition {
                    inner_depth: 3,
                    series_number: 3,
                    branch: vec![0],
                    attachment: None,
                },
                GlycanPosition {
                    inner_depth: 4,
                    series_number: 4,
                    branch: vec![1, 1],
                    attachment: None,
                },
            ],
        ),
        (
            18,
            Some(GlycanPosition {
                inner_depth: 1,
                series_number: 1,
                branch: vec![0],
                attachment: None,
            }),
            vec![
                GlycanPosition {
                    inner_depth: 3,
                    series_number: 3,
                    branch: vec![0, 0],
                    attachment: None,
                },
                GlycanPosition {
                    inner_depth: 6,
                    series_number: 6,
                    branch: vec![0, 1],
                    attachment: None,
                },
            ],
        ),
        (
            1,
            None,
            vec![GlycanPosition {
                inner_depth: 1,
                series_number: 1,
                branch: Vec::new(),
                attachment: None,
            }],
        ),
        (
            1,
            None,
            vec![GlycanPosition {
                inner_depth: 2,
                series_number: 2,
                branch: Vec::new(),
                attachment: None,
            }],
        ),
    ] {
        let structure =
            GlycanStructure::from_short_iupac(codes[index].1, 0..codes[index].1.len(), 0).unwrap();
        let rendered = structure
            .render(
                Some("pep".to_string()),
                COLUMN_SIZE,
                SUGAR_SIZE,
                STROKE_SIZE,
                if index % 3 == 0 {
                    GlycanDirection::LeftToRight
                } else {
                    GlycanDirection::TopDown
                },
                root,
                &breaks,
                [0, 0, 0],
                [255, 255, 255],
                &mut footnotes,
            )
            .unwrap();
        rendered.to_svg(&mut html).unwrap();
        let (bitmap, width) = rendered.to_bitmap(
            if index % 2 == 0 {
                zeno::Format::subpixel_bgra()
            } else {
                zeno::Format::Alpha
            },
            font.as_ref(),
            &mut context,
        );
        let mut buffer = Vec::new();
        let mut w = BufWriter::new(&mut buffer);
        let mut encoder =
            png::Encoder::new(&mut w, width as u32, (bitmap.len() / 4 / width) as u32);
        encoder.set_color(png::ColorType::Rgba);
        encoder.set_depth(png::BitDepth::Eight);
        let mut writer = encoder.write_header().unwrap();
        writer.write_image_data(&bitmap).unwrap();
        drop(writer);
        drop(w);

        write!(&mut html, "<img src=\"data:image/png;base64, ").unwrap();
        base64::engine::general_purpose::STANDARD.encode_string(&buffer, &mut html);
        write!(&mut html, "\"/>").unwrap();
    }

    write!(&mut html, "<hr>").unwrap();
    if !footnotes.is_empty() {
        write!(&mut html, "<ol>").unwrap();
        for note in footnotes {
            write!(&mut html, "<li>{note}</li>").unwrap();
        }
        write!(&mut html, "</ol><hr>").unwrap();
    }
    for (code, _) in &codes {
        write!(
            &mut html,
            "<image src=\"https://image.glycosmos.org/snfg/png/{code}\"/>"
        )
        .unwrap();
    }
    write!(&mut html, "</body></html>").unwrap();
    std::fs::write("../rendered_glycans.html", html).unwrap();
}
