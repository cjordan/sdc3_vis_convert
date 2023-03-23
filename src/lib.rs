//! Ripped out of hyperdrive.

pub mod averaging;
pub mod read;
pub mod write;

use hifitime::{Duration, Epoch};
use marlu::{LatLngHeight, RADec, XyzGeodetic};
use vec1::Vec1;

#[derive(Debug, PartialEq)]
pub enum VisInputType {
    MeasurementSet,
    Uvfits,
}

#[derive(Clone)]
pub struct ObsContext {
    /// The unique timestamps in the observation. These are stored as `hifitime`
    /// [Epoch] structs to help keep the code flexible. These include timestamps
    /// that are deemed "flagged" by the observation.
    pub timestamps: Vec1<Epoch>,

    /// The *available* timestep indices of the input data. This does not
    /// necessarily start at 0, and is not necessarily regular (e.g. a valid
    /// vector could be [1, 2, 4]).
    ///
    /// Allowing the indices to be non-regular means that we can represent input
    /// data that also isn't regular; naively reading in a dataset with 2
    /// timesteps that are separated by more than the time resolution of the
    /// data would give misleading results.
    pub all_timesteps: Vec1<usize>,

    /// The timestep indices of the input data that aren't totally flagged.
    ///
    /// This is allowed to be empty.
    pub unflagged_timesteps: Vec<usize>,

    /// The observation phase centre.
    pub phase_centre: RADec,

    /// The Earth position of the instrumental array.
    pub array_position: LatLngHeight,

    /// The difference between UT1 and UTC. If this is 0 seconds, then LSTs are
    /// wrong by up to 0.9 seconds. The code will assume that 0 seconds means
    /// that DUT1 wasn't provided and may warn the user.
    ///
    /// This is *probably* defined off of the obsid, but we don't expect DUT1 to
    /// change significantly across the course of an observation.
    pub dut1: Duration,

    /// The names of each of the tiles in the input data. This includes flagged
    /// and unavailable tiles.
    pub tile_names: Vec1<String>,

    /// The [`XyzGeodetic`] coordinates of all tiles in the array (all
    /// coordinates are specified in \[metres\]). This includes flagged and
    /// unavailable tiles.
    pub tile_xyzs: Vec1<XyzGeodetic>,

    /// The time resolution of the supplied data. This is not necessarily the
    /// native time resolution of the original observation's data, as it may
    /// have already been averaged. This is kept optional in case in the input
    /// data doesn't report the resolution and has only one timestep, and
    /// therefore no resolution.
    pub time_res: Duration,

    /// The fine-channel resolution of the supplied data \[Hz\]. This is not
    /// necessarily the fine-channel resolution of the original observation's
    /// data; this data may have applied averaging to the original observation.
    pub freq_res: f64,

    /// All of the fine-channel frequencies within the data \[Hz\]. The values
    /// reflect the frequencies at the *centre* of each channel.
    ///
    /// These are kept as ints to help some otherwise error-prone calculations
    /// using floats. By using ints, we assume there is no sub-Hz structure.
    pub fine_chan_freqs: Vec1<u64>,
}
