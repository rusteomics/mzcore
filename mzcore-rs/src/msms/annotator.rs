
//use itertools::Itertools;
use serde::{Deserialize, Serialize};

//use crate::ms::spectrum::SpectrumData;
use crate::msms::fragmentation::*;
use crate::msms::model::FragmentIonSeries;

/// Struct that contains the required info about a match between the exp. and theo. data.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct MatchedPeak {
    pub peak_index: usize,
    pub peak_mz: f64,
    pub peak_intensity: f64,
    pub theo_mz: f64,
    pub mz_error: f32,
    pub charge: i8,
    pub ion_type: FragmentIonSeries,
    pub frag_index: u16,  // index of m/z value in the fragmentation table (starts at 0)
    pub aa_position: u16, // AA position in amino acid sequence (starts at 1)
}

pub fn annotate_spectrum(spectrum_peaks: &[[f64;2]], frag_table: &FragmentationTable, mz_error_tol: f64) -> Vec<MatchedPeak> {

    let frag_table_n_cols = frag_table.len();
    let frag_table_n_rows = frag_table.first().unwrap().mz_values.len();

    let mut all_theo_frags = Vec::with_capacity(frag_table_n_cols * frag_table_n_rows);

    let mut frag_table_col_idx = 0 as usize;
    for frag_series_fragments in frag_table.iter() {
        let current_series = frag_series_fragments.mz_values.to_owned();

        let mut frag_table_row_idx = 0 as usize;
        for mz_value in current_series {
            //let rect = Rectangle::from_corners([mz_value - mz_error_tol, 0.0], [mz_value + mz_error_tol, 0.0]);
            //all_theo_rects.push(CustomTheoIon::new(rect, (frag_table_col_idx, frag_table_row_idx) ) );

            all_theo_frags.push((
                mz_value,
                frag_table_col_idx,
                frag_table_row_idx
            ));

            frag_table_row_idx += 1
        }

        frag_table_col_idx += 1
    }

    all_theo_frags.sort_by(|tuple_a, tuple_b| {
        tuple_a.0.partial_cmp(&tuple_b.0).unwrap()
    });

    let mut grouped_theo_frags = Vec::with_capacity(all_theo_frags.len());
    let mut theo_frags_group_buffer = Vec::new();

    let mut prev_frag_mz = 0.0;
    for theo_frag in all_theo_frags.iter() {
        let (frag_mz, _, _) = *theo_frag;

        if !theo_frags_group_buffer.is_empty() && frag_mz - prev_frag_mz > mz_error_tol {
            grouped_theo_frags.push(theo_frags_group_buffer.clone());
            theo_frags_group_buffer.clear();
        }

        theo_frags_group_buffer.push(theo_frag);

        prev_frag_mz = frag_mz;
    }

    // Append remaining fragments
    grouped_theo_frags.push(theo_frags_group_buffer);

    let mut frag_group_iter = grouped_theo_frags.iter();
    let mut cur_frag_group_opt = frag_group_iter.next();

    let mut matched_peaks: Vec<MatchedPeak> = Vec::with_capacity(all_theo_frags.len());

    let n_peaks = spectrum_peaks.len();
    let mut peak_idx = 0;
    while cur_frag_group_opt.is_some() && peak_idx < n_peaks {

        let peak = spectrum_peaks[peak_idx];
        let peak_mz = peak[0];
        let min_exp_mz = peak_mz - mz_error_tol;
        //println!("peak_mz={}",peak_mz);
        //println!("min_exp_mz={}",peak_mz);

        let mut cur_frag_group = cur_frag_group_opt.unwrap();
        let mut last_theo_mz = cur_frag_group.last().unwrap().0;
        //println!("last_theo_mz={}",last_theo_mz);

        // Check if we need to take the next theoretical fragment
        if last_theo_mz < min_exp_mz {
            let mut is_frag_above_min_mz = false;
            while !is_frag_above_min_mz && cur_frag_group_opt.is_some() {
                cur_frag_group_opt = frag_group_iter.next();
                if cur_frag_group_opt.is_some() {
                    cur_frag_group = cur_frag_group_opt.unwrap();
                    last_theo_mz = cur_frag_group.first().unwrap().0;
                    is_frag_above_min_mz = last_theo_mz >= min_exp_mz;
                }
            }
        }

        //println!("new last_theo_mz={}",last_theo_mz);

        for cur_frag in cur_frag_group {
            let (theo_mz, frag_table_col_idx, frag_table_row_idx) = **cur_frag;
            //println!("theo_mz={}",theo_mz);

            let mz_error = peak_mz - theo_mz;

            if mz_error.abs() <= mz_error_tol {
                //println!("matching fragment with m/z={}", peak_mz);

                let frag_table_col: &TheoreticalFragmentIons = frag_table.get(frag_table_col_idx).unwrap();
                let aa_position = if frag_table_col.ion_type.is_n_terminal().unwrap() {
                    (1 + frag_table_row_idx ) as u16
                } else {
                    (1 + frag_table_n_rows - frag_table_row_idx) as u16
                };

                matched_peaks.push(MatchedPeak {
                    peak_index: peak_idx,
                    peak_mz: peak_mz,
                    peak_intensity: peak[1],
                    theo_mz: theo_mz,
                    mz_error: mz_error as f32,
                    charge: frag_table_col.charge,
                    ion_type: frag_table_col.ion_type,
                    frag_index: frag_table_row_idx as u16,
                    aa_position: aa_position,
                });
            }
        }

        peak_idx += 1;
    }

    matched_peaks
}

// --- R*Tree based annotator by david-bouyssie (should be implemented in a dedicated feature)

/*
use rstar::RTree;

use rstar::primitives::Rectangle;
use rstar::PointDistance;
use rstar::primitives::GeomWithData;

pub type CustomExpPeak = GeomWithData<[f64; 2], (usize,f64)>;
type CustomTheoIon = GeomWithData<Rectangle<[f64; 2]>, (usize,usize)>;

pub enum PeakSelectionStrategy {
    NEAREST_PEAK,
    HIGHEST_PEAK,
}

fn index_peaks(spectrum_peaks: &[[f64;2]]) -> RTree<CustomExpPeak> {
    // --- Compute R*Tree indexing of experimental data --- //
    let all_custom_points: Vec<CustomExpPeak> = spectrum_peaks.iter().enumerate().map(|(peak_idx, peak)| {
        CustomExpPeak::new([peak[0], 0.0], (peak_idx,peak[1]))
    }).collect();

    RTree::bulk_load(all_custom_points)
}

fn index_spectrum(spectrum_data: &SpectrumData) -> RTree<CustomExpPeak> {
    // --- Compute R*Tree indexing of experimental data --- //
    let all_custom_points: Vec<CustomExpPeak> = spectrum_data.mz_list.iter().enumerate().map(|(peak_idx, mz)| {
        CustomExpPeak::new([*mz, 0.0], (peak_idx,spectrum_data.intensity_list[peak_idx] as f64))
    }).collect();

    RTree::bulk_load(all_custom_points)
}

fn annotate_spectrum_rtree(spectrum_peaks: &Vec<[f64;2]>, frag_table: &FragmentationTable, mz_error_tol: f64, peak_sel_strategy: PeakSelectionStrategy) -> Vec<MatchedPeak> {
    let spectrum_rtree = index_peaks(spectrum_peaks);

    annotate_indexed_spectrum(&spectrum_rtree, frag_table, mz_error_tol, peak_sel_strategy)
}

fn annotate_indexed_spectrum(indexed_spectrum: &RTree<CustomExpPeak>, frag_table: &FragmentationTable, mz_error_tol: f64, peak_sel_strategy: PeakSelectionStrategy) -> Vec<MatchedPeak> {

    let frag_table_len = frag_table.len();

    // --- Compute R*Tree indexing of theoretical data --- //
    let mut all_theo_rects = Vec::with_capacity( frag_table_len * frag_table.first().unwrap().mz_values.len());

    let mut frag_table_col_idx = 0;
    for frag_series_fragments in frag_table.clone() {
        let current_series = frag_series_fragments.mz_values;

        let mut frag_table_row_idx = 0;
        for mz_value in current_series {
            let rect = Rectangle::from_corners([mz_value - mz_error_tol, 0.0], [mz_value + mz_error_tol, 0.0]);
            all_theo_rects.push(CustomTheoIon::new(rect, (frag_table_col_idx, frag_table_row_idx) ) );

            frag_table_row_idx += 1
        }

        frag_table_col_idx += 1
    }
    let theo_tree = RTree::bulk_load(all_theo_rects);

    // --- Compute matched peaks using the intersection between the two R*Trees (experimental and theoretical) --- //
    let intersection = indexed_spectrum.intersection_candidates_with_other_tree(&theo_tree);
    /*intersection.for_each(|(exp_point, theo_line)| {
        let peak_mz = exp_point.geom()[0];
        let peak_intensity = exp_point.data;
        let theo_mz = theo_line.geom().upper()[0] + 0.02;
        let (ion_type, charge, aa_index) = theo_line.data;
        println!("{} {} {} {} {}", peak_mz, theo_mz, ion_type, charge, aa_index)
    });*/

    let matched_peaks: Vec<MatchedPeak> = intersection.map(|(exp_point, theo_line)| {
        let peak_mz = exp_point.geom()[0];
        let (peak_idx,peak_intensity) = exp_point.data;
        let (frag_table_col_idx, frag_table_row_idx) = theo_line.data;
        let frag_table_col = frag_table.get(frag_table_col_idx).unwrap();
        let theo_mz = frag_table_col.mz_values.get(frag_table_row_idx).unwrap();
        //println!("{} {} {} {} {}", peak_mz, theo_mz, ion_type, charge, aa_index)

        let mz_error = peak_mz - theo_mz;

        let aa_position = if is_ion_forward(frag_table_col.ion_type).unwrap() {
            (1 + frag_table_row_idx ) as u16
        } else {
            (1 + frag_table_len - frag_table_row_idx) as u16
        };

        MatchedPeak {
            peak_index: peak_idx,
            peak_mz: peak_mz,
            peak_intensity: peak_intensity,
            theo_mz: *theo_mz,
            mz_error: mz_error as f32,
            charge: frag_table_col.charge,
            ion_type: frag_table_col.ion_type,
            frag_index: frag_table_row_idx as u16,
            aa_position: aa_position
        }
    }).collect();

    // --- Remove matched peaks redundancy by selection the best candidate for a given theoretical m/z value --- //
    let mut non_redundant_matched_peaks = Vec::with_capacity(matched_peaks.len());

    let groups_iter = &matched_peaks.into_iter().group_by(|matched_peak| (matched_peak.theo_mz,matched_peak.ion_type));
    for (_, grouped_items_iter) in groups_iter {
        let grouped_matched_peaks: Vec<MatchedPeak> = grouped_items_iter.collect();
        let n_items = grouped_matched_peaks.len();
        //println!("{} {}", theo_ion_mz, grouped_matched_peaks.first().unwrap().peak_mz);

        // if we have more than one items in grouped_matched_peaks, it means that we have ambiguous peak matching for considered theoretical value.
        // so we need to choose one of the matched peaks according to the given selection strategy.
        if n_items > 1 {
            //println!("found one duplicate");
            let best_matched_peak = match peak_sel_strategy {
                PeakSelectionStrategy::NEAREST_PEAK => {
                    let nearest_peak_opt = grouped_matched_peaks.into_iter().min_by(|x1_peak, x2_peak| {
                        // Use partial_cmp because cmp cannot deal with floating point numbers (either f32 or f64)
                        x1_peak.mz_error.abs().partial_cmp(&x2_peak.mz_error.abs()).unwrap_or(std::cmp::Ordering::Equal)
                    });

                    nearest_peak_opt.unwrap() // safe because n_items > 1
                }
                PeakSelectionStrategy::HIGHEST_PEAK => {
                    let highest_peak_opt = grouped_matched_peaks.into_iter().max_by(|x1_peak, x2_peak| {
                        // Use partial_cmp because cmp cannot deal with floating point numbers (either f32 or f64)
                        x1_peak.peak_intensity.partial_cmp(&x2_peak.peak_intensity).unwrap_or(std::cmp::Ordering::Equal)
                    });

                    highest_peak_opt.unwrap() // safe because n_items > 1
                }
            };
            non_redundant_matched_peaks.push(best_matched_peak);
        } else if n_items == 1 {
            let single_matched_peak = grouped_matched_peaks.first().unwrap(); // safe because n_items == 1
            non_redundant_matched_peaks.push(*single_matched_peak);
        }
    }

    /*non_redundant_matched_peaks.iter().for_each(|matched_peak| {
        let peak_mz = matched_peak.peak_mz;
        let peak_intensity = matched_peak.peak_intensity;
        let theo_mz = matched_peak.theo_mz;
        let (ion_type, charge, aa_index) = (matched_peak.ion_type, matched_peak.charge, matched_peak.aa_index);
        //println!("{} {} {} {} {}", peak_mz, theo_mz, ion_type, charge, aa_index);
        println!("{}", matched_peak.mz_error);
    });*/

    non_redundant_matched_peaks
}

*/
