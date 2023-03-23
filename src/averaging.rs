//! Ripped out of hyperdrive.

use std::ops::Range;

use hifitime::{Duration, Epoch};
use vec1::Vec1;

/// A collection of timesteps.
#[derive(Debug, Clone)]
pub struct Timeblock {
    /// The timeblock index. e.g. If all observation timesteps are being used in
    /// a single calibration timeblock, then its index is 0.
    pub index: usize,

    /// The range of indices into an *unflagged* array of visibilities.
    ///
    /// The timesteps comprising a timeblock need not be contiguous, however, we
    /// want the timestep visibilities to be contiguous. Here, `range` indicates
    /// the *unflagged* timestep indices *for this timeblock*. e.g. If timeblock
    /// 0 represents timestep 10 and timeblock 1 represents timesteps 15 and 16
    /// (and these are the only timesteps used for calibration), then timeblock
    /// 0's range is 0..1 (only one index, 0), whereas timeblock 1's range is
    /// 1..3 (two indices starting at 1).
    ///
    /// We can use a range because the timesteps belonging to a timeblock are
    /// always contiguous.
    pub range: Range<usize>,

    /// The timestamps comprising this timeblock. These are determined by the
    /// timesteps into all available timestamps.
    pub timestamps: Vec1<Epoch>,

    /// The median timestamp of the *ideal* timeblock.
    ///
    /// e.g. If we have 9 timesteps and we're averaging 3, the averaged
    /// timeblocks look like this:
    ///
    /// [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    ///
    /// But if we're only using timesteps [1, 3, 8], the timeblocks look like
    /// this.
    ///
    /// [[1, _, 3], [_, _, _], [_, 8]]
    ///
    /// In the first case, this `median` is [1, 4, 7] for each timeblock, [2, 5,
    /// 8] for the second. Note how missing timestamps don't affect it.
    pub median: Epoch,
}

pub fn timesteps_to_timeblocks(
    all_timestamps: &Vec1<Epoch>,
    time_average_factor: usize,
    timesteps_to_use: &Vec1<usize>,
) -> Vec1<Timeblock> {
    let time_res = all_timestamps
        .windows(2)
        .fold(Duration::from_seconds(f64::INFINITY), |a, t| {
            a.min(t[1] - t[0])
        });
    let timestamps_to_use = timesteps_to_use.mapped_ref(
        |&t_step|
            // TODO: Handle incorrect timestep indices.
            *all_timestamps.get(t_step).unwrap(), // Could use square brackets, but this way the unwrap is clear.
    );

    // Populate the median timestamps of all timeblocks based off of the first
    // timestamp. e.g. If there are 10 timestamps with an averaging factor of 3,
    // these are some possible situations:
    //
    // [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9]]
    //
    // [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    //
    // [[2, 3, 4], [5, 6, 7], [8, 9]]
    //
    // Each depends on the first timestep being used. For each of the respective
    // situations above, the results are:
    //
    // [(1, 0..=2), (4, 3..=5), (7, 6..=8), (10, 9..=11)]
    //
    // [(2, 1..=3), (5, 4..=6), (8, 7..=9)]
    //
    // [(3, 2..=4), (6, 5..=7), (9, 8..=10)]

    let mut timeblocks = vec![];
    let timeblock_length = Duration::from_total_nanoseconds(
        // time_average_factor as i128 * time_res.total_nanoseconds(),
        (time_average_factor - 1) as i128 * time_res.total_nanoseconds(),
    );
    let half_a_timeblock = timeblock_length / 2;
    let first_timestamp = *timestamps_to_use.first();
    let last_timestamp = *timestamps_to_use.last();
    let time_res = time_res.total_nanoseconds() as u128;
    let time_average_factor = time_average_factor as u128;
    let mut timeblock_index = 0;
    let mut timestep_index = 0;
    for i in 0.. {
        // `timeblock_start` and `timeblock_end` are not centroids but "leading
        // edge", however `timeblock_median` is a centroid.
        let timeblock_start = first_timestamp
            + Duration::from_total_nanoseconds(
                (time_res
                    .checked_mul(i)
                    .unwrap()
                    .checked_mul(time_average_factor)
                    .unwrap()) as i128,
            );
        let timeblock_end = timeblock_start + timeblock_length;
        let timeblock_median = timeblock_start + half_a_timeblock;

        if timeblock_start > last_timestamp {
            break;
        }

        let timeblock_timestamps = timestamps_to_use
            .iter()
            .filter(|ts| (timeblock_start..=timeblock_end).contains(ts))
            .copied()
            .collect::<Vec<_>>();
        if !timeblock_timestamps.is_empty() {
            let num_timeblock_timestamps = timeblock_timestamps.len();
            timeblocks.push(Timeblock {
                index: timeblock_index,
                range: timestep_index..timestep_index + num_timeblock_timestamps,
                timestamps: Vec1::try_from_vec(timeblock_timestamps).unwrap(),
                median: timeblock_median,
            });
            timeblock_index += 1;
            timestep_index += num_timeblock_timestamps;
        }
    }

    Vec1::try_from_vec(timeblocks).unwrap()
}
