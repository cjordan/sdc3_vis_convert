use std::{collections::HashSet, path::PathBuf, thread::scope};

use clap::{AppSettings, Parser};
use crossbeam_channel::bounded;
use crossbeam_utils::atomic::AtomicCell;
use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};
use log::{debug, info};
use ndarray::prelude::*;
use vec1::{vec1, Vec1};

use sdc3_vis_convert::{
    averaging::timesteps_to_timeblocks,
    read::{ms::MsReader, uvfits::UvfitsReader, VisRead},
    write::{write_vis, VisOutputType},
    VisInputType,
};

#[derive(Parser)]
#[clap(global_setting(AppSettings::DeriveDisplayOrder))]
#[clap(disable_help_subcommand = true)]
// #[clap(propagate_version = true)]
#[clap(infer_long_args = true)]
struct Args {
    /// The data to be flattened.
    data: Vec<PathBuf>,

    /// The flattened data to be written out.
    #[clap(short, long)]
    output: PathBuf,

    #[clap(short, long, multiple_values(true))]
    timesteps: Option<Vec<usize>>,

    /// The verbosity of the program. Increase by specifying multiple times
    /// (e.g. -vv). The default is to print only high-level information.
    #[clap(short, long, parse(from_occurrences))]
    verbosity: u8,

    /// Disable progress bars.
    #[clap(long)]
    no_progress_bars: bool,
}

fn main() {
    let mut args = Args::parse();
    args.data.sort_unstable();
    setup_logging(args.verbosity);

    // Verify all input files have the same extension.
    let input_type: VisInputType = {
        let mut input_type = None;
        for d in &args.data {
            let this_type = match d.extension().and_then(|os_str| os_str.to_str()) {
                Some("uvfits" | "uvf") => VisInputType::Uvfits,
                Some("ms") => VisInputType::MeasurementSet,
                _ => panic!("no file extension"),
            };
            if let Some(input_type) = input_type.as_ref() {
                if *input_type != this_type {
                    panic!("mismatch on input vis types");
                }
            } else {
                input_type = Some(this_type);
            }
        }
        input_type.unwrap()
    };
    info!("Input type: {input_type:?}");

    let output_type = match args.output.extension().and_then(|os_str| os_str.to_str()) {
        Some("uvfits") => VisOutputType::Uvfits,
        Some("uvf") => panic!("Refusing to write to a .uvf, use .uvfits"),
        Some("ms") => VisOutputType::MeasurementSet,
        _ => panic!("dunno what output format you want"),
    };
    info!("Output type: {output_type:?}");

    let mut reader: Box<dyn VisRead> = match input_type {
        VisInputType::MeasurementSet => Box::new(MsReader::new(&args.data[0])),

        VisInputType::Uvfits => Box::new(UvfitsReader::new(&args.data[0])),
    };
    let mut obs_context = reader.get_obs_context().clone();

    let num_stations = obs_context.tile_xyzs.len();
    let num_baselines = (num_stations * (num_stations - 1)) / 2;

    let time_average_factor = 1;
    let freq_average_factor = 1;

    // Which timesteps?
    let timesteps = match args.timesteps {
        Some(mut v) => {
            v.sort_unstable();
            Vec1::try_from_vec(v).unwrap()
        }
        None => obs_context.all_timesteps.clone(),
    };
    let timeblocks =
        timesteps_to_timeblocks(&obs_context.timestamps, time_average_factor, &timesteps);

    // Add all frequencies.
    for file in args.data.iter().skip(1) {
        let mut reader: Box<dyn VisRead> = match input_type {
            VisInputType::MeasurementSet => Box::new(MsReader::new(file)),

            VisInputType::Uvfits => Box::new(UvfitsReader::new(file)),
        };
        let file_obs_context = reader.get_obs_context();
        for freq in file_obs_context.fine_chan_freqs.iter() {
            if !obs_context.fine_chan_freqs.contains(freq) {
                obs_context.fine_chan_freqs.push(*freq);
            }
        }
    }
    obs_context.fine_chan_freqs.sort_unstable();
    let num_channels = obs_context.fine_chan_freqs.len();

    let (tx, rx) = bounded(5);
    let error = AtomicCell::new(false);
    let multi_progress = MultiProgress::with_draw_target(if args.no_progress_bars {
        ProgressDrawTarget::hidden()
    } else {
        ProgressDrawTarget::stdout()
    });
    let read_progress = multi_progress.add(
        ProgressBar::new(timeblocks.len() as _)
            .with_style(
                ProgressStyle::default_bar()
                    .template("{msg:17}: [{wide_bar:.blue}] {pos:2}/{len:2} timeblocks ({elapsed_precise}<{eta_precise})").unwrap()
                    .progress_chars("=> "),
            )
            .with_position(0)
            .with_message("Reading"),
    );
    let write_progress = multi_progress.add(
        ProgressBar::new(timeblocks.len() as _)
            .with_style(
                ProgressStyle::default_bar()
                    .template("{msg:17}: [{wide_bar:.blue}] {pos:2}/{len:2} timeblocks ({elapsed_precise}<{eta_precise})").unwrap()
                    .progress_chars("=> "),
            )
            .with_position(0)
            .with_message("Writing"),
    );
    read_progress.tick();
    write_progress.tick();

    scope(|s| {
        s.spawn(|| {
            for &timestep in &timesteps {
                debug!("Working on timestep {timestep}");
                let mut all_vis_for_a_timestep = Array2::zeros((num_channels, num_baselines));

                for (i_freq, file) in args.data.iter().enumerate() {
                    debug!("Working on {}", file.display());
                    reader.update_file(file.clone());
                    let freq_data = all_vis_for_a_timestep.slice_mut(s![i_freq..i_freq + 1, ..]);
                    reader.read(freq_data, timestep)
                }

                // Write the vis out.
                tx.send((all_vis_for_a_timestep, obs_context.timestamps[timestep]))
                    .unwrap();
                read_progress.inc(1);
            }
        });

        s.spawn(|| {
            let station_index_range = 0..num_stations;
            let baseline_pairs = {
                let mut baseline_pairs = Vec::with_capacity(num_baselines);
                for i_station1 in station_index_range.clone() {
                    for i_station2 in station_index_range.clone().skip(i_station1 + 1) {
                        baseline_pairs.push((i_station1, i_station2));
                    }
                }
                baseline_pairs
            };
            assert_eq!(baseline_pairs.len(), num_baselines);

            write_vis(
                &vec1![(args.output, output_type)],
                obs_context.array_position,
                obs_context.phase_centre,
                &obs_context.tile_xyzs,
                &obs_context.tile_names,
                None,
                &obs_context.timestamps,
                &timesteps,
                &timeblocks,
                obs_context.time_res,
                obs_context.dut1,
                obs_context.freq_res,
                &obs_context.fine_chan_freqs.mapped_ref(|i| *i as f64),
                &baseline_pairs,
                &HashSet::new(),
                time_average_factor,
                freq_average_factor,
                rx,
                &error,
                Some(write_progress),
            );
        });
    });
}

fn setup_logging(verbosity: u8) {
    let mut builder = env_logger::Builder::from_default_env();
    builder.target(env_logger::Target::Stdout);
    builder.format_target(false);
    match verbosity {
        0 => builder.filter_level(log::LevelFilter::Info),
        1 => builder.filter_level(log::LevelFilter::Debug),
        2 => builder.filter_level(log::LevelFilter::Trace),
        _ => {
            builder.filter_level(log::LevelFilter::Trace);
            builder.format(|buf, record| {
                use std::io::Write;

                // TODO: Add colours.
                let timestamp = buf.timestamp();
                let level = record.level();
                let target = record.target();
                let line = record.line().unwrap_or(0);
                let message = record.args();

                writeln!(buf, "[{timestamp} {level} {target}:{line}] {message}")
            })
        }
    };
    builder.init();
}
