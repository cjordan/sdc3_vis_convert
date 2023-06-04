//! Ripped out of hyperdrive.

use std::{
    collections::BTreeSet,
    path::{Path, PathBuf},
};

use hifitime::{Duration, Epoch, TimeUnits};
use log::{debug, trace};
use marlu::{c32, rubbl_casatables, LatLngHeight, RADec, XyzGeocentric, XyzGeodetic};
use ndarray::prelude::*;
use rayon::prelude::*;
use rubbl_casatables::{Table, TableOpenMode};
use vec1::Vec1;

use super::VisRead;
use crate::{ObsContext, VisInputType};

/// Open a measurement set table read only. If `table` is `None`, then open the
/// base table.
fn read_table(ms: &Path, table: Option<&str>) -> Table {
    Table::open(
        format!("{}/{}", ms.display(), table.unwrap_or("")),
        TableOpenMode::Read,
    )
    .unwrap()
}

pub struct MsReader {
    /// Input data metadata.
    obs_context: ObsContext,

    /// The path to the measurement set on disk.
    pub ms: PathBuf,

    /// The "stride" of the data, i.e. the number of rows (baselines) before the
    /// time index changes.
    step: usize,

    /// The name of the column to be used containing visibility data in the main
    /// column.
    _data_col_name: &'static str,
}

impl MsReader {
    /// Verify and populate metadata associated with this measurement set.
    ///
    /// The measurement set is expected to be formatted in the way that
    /// cotter/Birli write measurement sets. There's a difference between a
    /// flagged antenna and an antenna which has no data. The former may be
    /// used, but its flagged status hints that maybe it shouldn't be used.
    // TODO: Handle multiple measurement sets.
    pub fn new<P: AsRef<Path>>(ms: P) -> MsReader {
        let ms = ms.as_ref();
        debug!("Using measurement set: {}", ms.display());
        assert!(ms.exists());

        let mut main_table = read_table(ms, None);
        assert!(main_table.n_rows() > 0);
        let _data_col_name = "DATA";

        // Get the tile names and XYZ positions.
        let mut antenna_table = read_table(ms, Some("ANTENNA"));
        let tile_names: Vec<String> = antenna_table.get_col_as_vec("NAME").unwrap();
        trace!("There are {} tile names", tile_names.len());
        let tile_names = Vec1::try_from_vec(tile_names).unwrap();

        let array_position = LatLngHeight {
            longitude_rad: 116.7644482_f64.to_radians(),
            latitude_rad: -26.82472208_f64.to_radians(),
            height_metres: 0.0,
        };

        let tile_xyzs: Vec<XyzGeodetic> = {
            let mut casacore_positions = Vec::with_capacity(antenna_table.n_rows() as usize);
            antenna_table
                .for_each_row(|row| {
                    let pos: Vec<f64> = row.get_cell("POSITION")?;
                    let pos_xyz = XyzGeocentric {
                        x: pos[0],
                        y: pos[1],
                        z: pos[2],
                    };
                    casacore_positions.push(pos_xyz);
                    Ok(())
                })
                .unwrap();

            let vec = XyzGeocentric::get_geocentric_vector(array_position);
            let (s_long, c_long) = array_position.longitude_rad.sin_cos();
            let xyzs = casacore_positions
                .par_iter()
                .map(|xyz| xyz.to_geodetic_inner(vec, s_long, c_long))
                .collect();

            xyzs
        };
        trace!("There are positions for {} tiles", tile_xyzs.len());
        // Not sure if this is even possible, but we'll handle it anyway.
        assert_eq!(tile_xyzs.len(), tile_names.len());
        let tile_xyzs = Vec1::try_from_vec(tile_xyzs).unwrap();
        let total_num_tiles = tile_xyzs.len();

        let num_available_tiles = total_num_tiles;
        let step = num_available_tiles * (num_available_tiles - 1) / 2;
        trace!("MS step: {}", step);

        // Work out the first and last good timesteps. This is important
        // because the observation's data may start and end with
        // visibilities that are all flagged, and (by default) we are not
        // interested in using any of those data. We work out the first and
        // last good timesteps by inspecting the flags at each timestep.
        let unflagged_timesteps: Vec<usize> = {
            // The first and last good timestep indices.
            let mut first: Option<usize> = None;
            let mut last: Option<usize> = None;

            trace!("Searching for unflagged timesteps in the MS");
            for i_step in 0..(main_table.n_rows() as usize) / step {
                trace!("Reading timestep {i_step}");
                let mut all_rows_for_step_flagged = true;
                for i_row in 0..step {
                    let vis_flags: Vec<bool> = main_table
                        .get_cell_as_vec("FLAG", (i_step * step + i_row) as u64)
                        .unwrap();
                    let all_flagged = vis_flags.into_iter().all(|f| f);
                    if !all_flagged {
                        all_rows_for_step_flagged = false;
                        if first.is_none() {
                            first = Some(i_step);
                            debug!("First good timestep: {i_step}");
                        }
                        break;
                    }
                }
                if all_rows_for_step_flagged && first.is_some() {
                    last = Some(i_step);
                    debug!("Last good timestep: {}", i_step - 1);
                    break;
                }
            }

            // Did the indices get set correctly?
            match (first, last) {
                (Some(f), Some(l)) => f..l,
                // If there weren't any flags at the end of the MS, then the
                // last timestep is fine.
                (Some(f), None) => f..main_table.n_rows() as usize / step,
                // All timesteps are flagged. The user can still use the MS,
                // but they must specify some amount of flagged timesteps.
                _ => 0..0,
            }
        }
        .collect();

        // Get the unique times in the MS.
        let utc_times: Vec<f64> = main_table.get_col_as_vec("TIME").unwrap();
        let mut utc_time_set: BTreeSet<u64> = BTreeSet::new();
        let mut timestamps = vec![];
        for utc_time in utc_times {
            let bits = utc_time.to_bits();
            if !utc_time_set.contains(&bits) {
                utc_time_set.insert(bits);

                // casacore stores the times as centroids, so no correction
                // is needed.
                let e = Epoch::from_utc_seconds(
                    // casacore stores the times as UTC seconds... but with an
                    // offset.
                    utc_time - hifitime::J1900_OFFSET * hifitime::SECONDS_PER_DAY,
                );
                // The values can be slightly off of their intended values;
                // round them to the nearest 10 milliseconds.
                timestamps.push(e.round(10.milliseconds()));
            }
        }
        let timestamps = Vec1::try_from_vec(timestamps).unwrap();
        match timestamps.as_slice() {
            // Handled above; measurement sets aren't allowed to be empty.
            [] => unreachable!(),
            [t] => debug!("Only timestep (GPS): {:.2}", t.to_gpst_seconds()),
            [t0, .., tn] => {
                debug!("First good timestep (GPS): {:.2}", t0.to_gpst_seconds());
                debug!("Last good timestep  (GPS): {:.2}", tn.to_gpst_seconds());
            }
        }

        // Get the data's time resolution. There is a possibility that the MS
        // contains only one timestep.
        let time_res = {
            // Find the minimum gap between two consecutive timestamps.
            let time_res = timestamps
                .windows(2)
                .fold(Duration::from_seconds(f64::INFINITY), |acc, ts| {
                    acc.min(ts[1] - ts[0])
                });
            debug!("Time resolution: {}s", time_res.to_seconds());
            time_res
        };

        let all_timesteps = (0..timestamps.len()).collect();
        let all_timesteps = Vec1::try_from_vec(all_timesteps).unwrap();

        // Get the frequency information.
        let mut spectral_window_table = read_table(ms, Some("SPECTRAL_WINDOW"));
        let fine_chan_freqs_f64: Vec<f64> = spectral_window_table
            .get_cell_as_vec("CHAN_FREQ", 0)
            .unwrap();
        let fine_chan_freqs = {
            let fine_chan_freqs = fine_chan_freqs_f64
                .iter()
                .map(|f| f.round() as u64)
                .collect();
            Vec1::try_from_vec(fine_chan_freqs).unwrap()
        };
        // Assume that `total_bandwidth_hz` is the total bandwidth inside the
        // measurement set, which is not necessarily the whole observation.
        let total_bandwidth_hz: f64 = spectral_window_table
            .get_cell("TOTAL_BANDWIDTH", 0)
            .unwrap();
        debug!("MS total bandwidth: {} Hz", total_bandwidth_hz);

        // Round the values in here because sometimes they have a fractional
        // component, for some reason. We're unlikely to ever have a fraction of
        // a Hz as the channel resolution.
        let freq_res = {
            let all_widths: Vec<f64> = spectral_window_table
                .get_cell_as_vec("CHAN_WIDTH", 0)
                .unwrap();
            let width = *all_widths.first().unwrap();
            // Make sure all the widths all the same.
            for w in all_widths.iter().skip(1) {
                assert!((w - width).abs() <= f64::EPSILON);
            }
            width
        };

        // Get the observation phase centre.
        let phase_centre = {
            let mut field_table = read_table(ms, Some("FIELD"));
            let phase_vec = field_table.get_cell_as_vec("PHASE_DIR", 0).unwrap();
            RADec::from_radians(phase_vec[0], phase_vec[1])
        };

        let obs_context = ObsContext {
            timestamps,
            all_timesteps,
            unflagged_timesteps,
            phase_centre,
            array_position,
            dut1: Duration::from_seconds(super::TIME_OFFSET),
            tile_names,
            tile_xyzs,
            time_res,
            freq_res,
            fine_chan_freqs,
        };

        MsReader {
            obs_context,
            ms: ms.to_path_buf(),
            step,
            _data_col_name,
        }
    }

    /// An internal method for reading visibilities. Cross- and/or
    /// auto-correlation visibilities and weights are written to the supplied
    /// arrays.
    fn read_inner(&self, mut data_bf: ArrayViewMut2<c32>, timestep: usize) {
        // When reading in a new timestep's data, these indices should be
        // multiplied by `step` to get the amount of rows to stride in the main
        // table.
        let row_range_start = timestep * self.step;
        let row_range_end = (timestep + 1) * self.step;
        let row_range = row_range_start as u64..row_range_end as u64;

        let mut main_table = read_table(&self.ms, None);
        let mut _row_index = row_range.start;
        // let mut i_bl = 0;

        for (i_bl, row) in row_range.enumerate() {
            // try get_cell_as_vec
            let ms_data: c32 = main_table.get_cell("DATA", row).unwrap();
            unsafe {
                *data_bf.uget_mut((0, i_bl)) = ms_data.conj();
            }
            // data_bf.slice_mut(s![i_bl..i_bl + 1, ..]).assign(&ms_data);
        }

        // main_table
        //     .for_each_row_in_range(row_range, |row| {
        //         // Antenna numbers are zero indexed.
        //         let ant1: i32 = row.get_cell("ANTENNA1")?;
        //         let ant2: i32 = row.get_cell("ANTENNA2")?;
        //         if ant1 == 0 && ant2 == 1 {
        //             i_bl = 0;
        //         }

        //         // The data array is arranged [frequency][instrumental_pol].
        //         let ms_data: Array2<c32> = row.get_cell(self.data_col_name)?;

        //         assert!(
        //             i_bl < data_bf.len_of(Axis(0)),
        //             "{i_bl} >= data_bf.len_of(Axis(0))"
        //         );
        //         assert!(data_bf.len_of(Axis(1)) <= ms_data.len_of(Axis(1)));
        //         data_bf.slice_mut(s![i_bl..i_bl + 1, ..]).assign(&ms_data);

        //         // Put the data and weights into the shared arrays
        //         // outside this scope. Before we can do this, we need to
        //         // remove any globally-flagged fine channels.
        //         data_bf.assign(&ms_data);

        //         row_index += 1;
        //         i_bl += 1;
        //         Ok(())
        //     })
        //     .unwrap();
    }
}

impl VisRead for MsReader {
    fn get_obs_context(&self) -> &ObsContext {
        &self.obs_context
    }

    fn get_input_data_type(&self) -> VisInputType {
        VisInputType::MeasurementSet
    }

    fn read(&self, data_bf: ArrayViewMut2<c32>, timestep: usize) {
        self.read_inner(data_bf, timestep)
    }

    fn update_file(&mut self, file: PathBuf) {
        self.ms = file;
    }
}
