//! Ripped from hyperdrive.

use std::{collections::HashSet, path::PathBuf};

use crossbeam_channel::Receiver;
use crossbeam_utils::atomic::AtomicCell;
use hifitime::{Duration, Epoch};
use indicatif::ProgressBar;
use log::{debug, trace};
use marlu::{
    c32, History, Jones, LatLngHeight, MeasurementSetWriter, ObsContext as MarluObsContext, RADec,
    UvfitsWriter, VisContext, VisWrite, XyzGeodetic,
};
use ndarray::prelude::*;
use vec1::Vec1;

use crate::averaging::Timeblock;

#[derive(Debug, Clone, Copy)]
/// All write-supported visibility formats.
pub enum VisOutputType {
    Uvfits,
    MeasurementSet,
}

/// Create the specified visibility outputs and receive visibilities to write to
/// them. This function is intended to be run concurrently with other async
/// threads, and must receive the timesteps of data as specified by `timesteps`.
///
/// # Arguments
///
/// * `outputs` - each of the output files to be written, paired with the output
///   type.
/// * `unflagged_baseline_tile_pairs` - the tile indices corresponding to
///   unflagged baselines. This includes auto-correlation "baselines" if they
///   are unflagged.
/// * `array_pos` - the position of the array that produced these visibilities.
/// * `phase_centre` - the phase centre used for the incoming visibilities.
/// * `pointing_centre` - the pointing centre used for the incoming
///   visibilities.
/// * `tile_positions` - the *un-precessed* positions of the tiles used for the
///   incoming visibilities (flagged and unflagged).
/// * `tile_names` - the names of the tiles used for the incoming visibilities
///   (flagged and unflagged).
/// * `obsid` - the MWA observation ID. If provided, it is used as the scheduled
///   start time of the observation and as an identifier. If not provided, the
///   first timestep is used as the scheduled start time and a placeholder will
///   be used for the identifier.
/// * `timestamps` - all possible timestamps that could be written out. These
///   represent the centre of the integration bin, i.e. "centroid" and not
///   "leading edge". Must be ascendingly sorted and be regularly spaced in
///   terms of `time_res`, but gaps are allowed.
/// * `timesteps` - the timesteps to be written out. These are indices into
///   `timestamps`.
/// * `time_res` - the time resolution of the incoming visibilities.
/// * `fine_chan_freqs` - all of the fine channel frequencies \[Hz\] (flagged
///   and unflagged).
/// * `freq_res` - the frequency resolution of the incoming visibilities \[Hz\].
/// * `time_average_factor` - the time average factor (i.e. average this many
///   visibilities in time before writing out).
/// * `freq_average_factor` - the frequency average factor (i.e. average this
///   many channels before writing out).
/// * `marlu_mwa_obs_context` - a tuple of [`marlu::MwaObsContext`] and a range
///   of MWA coarse channel indices. Kept optional because they're not strictly
///   needed.
/// * `rx` - the channel to receive visibilities from.
/// * `error` - a thread-safe [`bool`] to indicate if an error has occurred.
///   Receiving `true` signals that we should not continue, as another thread
///   has experienced an error.
/// * `progress_bar` - an optional progress bar to increment with writing
///   progress.
///
/// # Returns
///
/// * A neatly-formatted string reporting all of the files that got written out.
#[allow(clippy::too_many_arguments)]
pub fn write_vis<'a>(
    outputs: &'a Vec1<(PathBuf, VisOutputType)>,
    array_pos: LatLngHeight,
    phase_centre: RADec,
    tile_positions: &'a [XyzGeodetic],
    tile_names: &'a [String],
    obsid: Option<u32>,
    timestamps: &'a Vec1<Epoch>,
    timesteps: &'a Vec1<usize>,
    timeblocks: &'a Vec1<Timeblock>,
    time_res: Duration,
    dut1: Duration,
    freq_res: f64,
    fine_chan_freqs: &'a Vec1<f64>,
    unflagged_baseline_tile_pairs: &'a [(usize, usize)],
    flagged_fine_chans: &HashSet<usize>,
    time_average_factor: usize,
    freq_average_factor: usize,
    rx: Receiver<(Array2<c32>, Epoch)>,
    error: &'a AtomicCell<bool>,
    progress_bar: Option<ProgressBar>,
) {
    // Ensure our timestamps are sensible.
    for &t in timestamps {
        let diff = (t - *timestamps.first()).total_nanoseconds();
        assert!(diff % time_res.total_nanoseconds() == 0);
    }

    let start_timestamp = timestamps[*timesteps.first()];
    let vis_ctx = VisContext {
        num_sel_timesteps: timesteps.len(),
        start_timestamp,
        int_time: time_res,
        num_sel_chans: fine_chan_freqs.len(),
        start_freq_hz: *fine_chan_freqs.first(),
        freq_resolution_hz: freq_res,
        sel_baselines: unflagged_baseline_tile_pairs.to_vec(),
        avg_time: time_average_factor,
        avg_freq: freq_average_factor,
        num_vis_pols: 1,
    };

    let obs_name = obsid.map(|o| format!("{o}"));
    let sched_start_timestamp = match obsid {
        Some(gpst) => Epoch::from_gpst_seconds(f64::from(gpst)),
        None => start_timestamp,
    };
    let sched_duration = timestamps[*timesteps.last()] + time_res - sched_start_timestamp;
    let (s_lat, c_lat) = array_pos.latitude_rad.sin_cos();
    let marlu_obs_ctx = MarluObsContext {
        sched_start_timestamp,
        sched_duration,
        name: obs_name,
        phase_centre,
        pointing_centre: None,
        array_pos,
        ant_positions_enh: tile_positions
            .iter()
            .map(|xyz| xyz.to_enh_inner(s_lat, c_lat))
            .collect(),
        ant_names: tile_names.to_vec(),
        // TODO(dev): is there any value in adding this metadata via hyperdrive obs context?
        field_name: None,
        project_id: None,
        observer: None,
    };

    // Prepare history for the output vis files. It's possible that the
    // command-line call has invalid UTF-8. So use args_os and attempt to
    // convert to UTF-8 strings. If there are problems on the way, don't bother
    // trying to write the CMDLINE key.
    let cmd_line = std::env::args_os()
        .map(|a| a.into_string())
        .collect::<Result<Vec<String>, _>>()
        .map(|v| v.join(" "))
        .ok();
    let history = History {
        application: Some("mwa_hyperdrive SDC3 hacked"),
        cmd_line: cmd_line.as_deref(),
        message: None,
    };
    let mut writers = vec![];
    for (output, vis_type) in outputs {
        debug!("Setting up {} ({vis_type:?})", output.display());
        let vis_writer: Box<dyn VisWrite> = match vis_type {
            VisOutputType::Uvfits => {
                let uvfits = UvfitsWriter::from_marlu(
                    output,
                    &vis_ctx,
                    array_pos,
                    phase_centre,
                    dut1,
                    marlu_obs_ctx.name.as_deref(),
                    tile_names.to_vec(),
                    tile_positions.to_vec(),
                    false,
                    Some(&history),
                )
                .unwrap();
                Box::new(uvfits)
            }

            VisOutputType::MeasurementSet => {
                let ms = MeasurementSetWriter::new(
                    output,
                    phase_centre,
                    array_pos,
                    tile_positions.to_vec(),
                    dut1,
                    false,
                    vis_ctx.num_vis_pols,
                );
                ms.initialize(&vis_ctx, &marlu_obs_ctx, Some(&history))
                    .unwrap();
                Box::new(ms)
            }
        };
        writers.push(vis_writer);
    }

    // These arrays will contain the post-averaged values and are written out by
    // the writer when all relevant timesteps have been added.
    // [time][freq][baseline]
    let out_shape = vis_ctx.sel_dims();

    // Track a reference to the timeblock we're writing.
    let mut this_timeblock = timeblocks.first();
    // Also track the first timestamp of the tracked timeblock.
    // let mut this_start_timestamp = None;
    let mut this_average_timestamp = None;
    let mut i_timeblock = 0;
    // And the timestep into the timeblock.
    let mut this_timestep = 0;

    // Receive visibilities from another thread.
    for (i_timestep, (cross_data, timestamp)) in rx.iter().enumerate() {
        debug!(
            "Received timestep {i_timestep} (GPS {})",
            timestamp.to_gpst_seconds()
        );
        if this_average_timestamp.is_none() {
            this_average_timestamp = Some(
                timeblocks
                    .iter()
                    .find(|tb| tb.timestamps.contains(&timestamp))
                    .unwrap()
                    .median,
            );
        }

        // baseline
        assert_eq!(cross_data.len_of(Axis(1)), out_shape.2);
        // freq
        assert_eq!(
            cross_data.len_of(Axis(0)) + flagged_fine_chans.len(),
            out_shape.1
        );

        // Should we continue?
        if error.load() {
            return;
        }

        // If the next timestep doesn't belong to our tracked timeblock, write
        // out this timeblock and track the next one.
        if !this_timeblock.range.contains(&(i_timestep + 1))
            || this_timestep + 1 >= time_average_factor
        {
            debug!("Writing timeblock {i_timeblock}");
            let chunk_vis_ctx = VisContext {
                // TODO: Marlu expects "leading edge" timestamps, not centroids.
                // Fix this in Marlu.
                start_timestamp: this_average_timestamp.unwrap()
                    - time_res / 2 * time_average_factor as f64,
                num_sel_timesteps: this_timeblock.range.len(),
                ..vis_ctx.clone()
            };

            trace!("out_data dimensions: {:?}", cross_data.dim());
            trace!("this_timeblock.range: {:?}", this_timeblock.range);
            let cross_data = {
                let d = cross_data.dim();
                let v = cross_data.into_raw_vec();
                let v2 = v
                    .into_iter()
                    .map(|c| Jones::from([c, c32::default(), c32::default(), c32::default()]))
                    .collect::<Vec<_>>();
                // Array3::from_shape_vec((1, d.1, d.0), v2).unwrap()
                Array3::from_shape_vec((1, d.0, d.1), v2).unwrap()
            };

            let weights = Array3::ones(cross_data.dim());
            for vis_writer in writers.iter_mut() {
                vis_writer
                    .write_vis(cross_data.view(), weights.view(), &chunk_vis_ctx)
                    .unwrap();
                // Should we continue?
                if error.load() {
                    return;
                }
            }

            if let Some(progress_bar) = progress_bar.as_ref() {
                progress_bar.inc(1);
            }

            i_timeblock += 1;
            this_timeblock = match timeblocks.get(i_timeblock) {
                Some(t) => t,
                None => break,
            };
            this_average_timestamp = None;
            this_timestep = 0;
        } else {
            this_timestep += 1;
        }
    }

    if let Some(progress_bar) = progress_bar.as_ref() {
        progress_bar.abandon_with_message("Finished writing visibilities");
    }

    for vis_writer in writers.iter_mut() {
        vis_writer.finalise().unwrap();
    }
    debug!("Finished writing");
}
