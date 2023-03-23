//! Ripped out of hyperdrive.

use std::{
    borrow::Cow,
    collections::{HashMap, HashSet},
    os::raw::c_char,
    path::{Path, PathBuf},
};

use fitsio::{errors::check_status as fits_check_status, hdu::FitsHdu, FitsFile};
use hifitime::{Duration, Epoch, TimeUnits};
use log::{debug, warn};
use marlu::{c32, LatLngHeight, RADec, XyzGeocentric};
use ndarray::prelude::*;
use vec1::Vec1;

use super::{
    fits::{fits_get_col, fits_get_optional_key, fits_get_required_key, fits_open, fits_open_hdu},
    VisRead,
};
use crate::{ObsContext, VisInputType};

pub struct UvfitsReader {
    /// Observation metadata.
    pub(super) obs_context: ObsContext,

    // uvfits-specific things follow.
    /// The path to the uvfits on disk.
    pub uvfits: PathBuf,

    /// The uvfits-specific metadata, like which indices contain which
    /// parameters.
    metadata: UvfitsMetadata,

    /// The "stride" of the data, i.e. the number of rows (baselines) before the
    /// time index changes.
    step: usize,
}

impl UvfitsReader {
    /// Verify and populate metadata associated with this measurement set.
    ///
    /// The measurement set is expected to be formatted in the way that
    /// cotter/Birli write measurement sets.
    pub fn new<P: AsRef<Path>>(uvfits: P) -> UvfitsReader {
        let uvfits = uvfits.as_ref();
        debug!("Using uvfits file: {}", uvfits.display());
        assert!(uvfits.exists());

        // Get the tile names, XYZ positions and antenna numbers.
        let mut uvfits_fptr = fits_open(uvfits);
        let primary_hdu = fits_open_hdu(&mut uvfits_fptr, 0);
        let antenna_table_hdu = fits_open_hdu(&mut uvfits_fptr, "AIPS AN");

        let tile_names: Vec<String> = fits_get_col(&mut uvfits_fptr, &antenna_table_hdu, "ANNAME");
        let tile_names = Vec1::try_from_vec(tile_names).unwrap();
        let total_num_tiles = tile_names.len();

        // Set up the tile map.
        let tile_nums: Vec<u32> = fits_get_col(&mut uvfits_fptr, &antenna_table_hdu, "NOSTA");
        let tile_map: HashMap<usize, usize> = tile_nums
            .into_iter()
            .zip(0..total_num_tiles)
            .map(|(a, b)| (a.try_into().expect("not larger than usize::MAX"), b))
            .collect();

        let array_position = LatLngHeight {
            longitude_rad: 116.7644482_f64.to_radians(),
            latitude_rad: -26.82472208_f64.to_radians(),
            height_metres: 0.0,
        };

        let tile_xyzs = {
            let mut tile_xyzs: Vec<XyzGeocentric> = Vec::with_capacity(total_num_tiles);
            for i in 0..total_num_tiles {
                let fits_xyz = read_cell_array(
                    &mut uvfits_fptr,
                    &antenna_table_hdu,
                    "STABXYZ",
                    i.try_into().expect("not larger than i64::MAX"),
                    3,
                );
                tile_xyzs.push(XyzGeocentric {
                    x: fits_xyz[0],
                    y: fits_xyz[1],
                    z: fits_xyz[2],
                });
            }

            // These are geocentric positions, but we need geodetic
            // ones. To convert, we need the array position. If we
            // couldn't determine it, we need to bail out here.
            let vec = XyzGeocentric::get_geocentric_vector(array_position);
            let (s_long, c_long) = array_position.longitude_rad.sin_cos();
            tile_xyzs
                .into_iter()
                .map(|geocentric| geocentric.to_geodetic_inner(vec, s_long, c_long))
                .collect()
        };
        let tile_xyzs = Vec1::try_from_vec(tile_xyzs)
            .expect("can't be empty, non-empty tile names verified above");

        let metadata = UvfitsMetadata::new(&mut uvfits_fptr, &primary_hdu);

        debug!("Number of rows in the uvfits:   {}", metadata.num_rows);
        debug!("PCOUNT:                         {}", metadata.pcount);
        debug!("Number of polarisations:        1");
        debug!("Floats per polarisation:        2");
        debug!(
            "Number of fine frequency chans: {}",
            metadata.num_fine_freq_chans
        );
        debug!("UU index:       {}", metadata.indices.u);
        debug!("VV index:       {}", metadata.indices.v);
        debug!("WW index:       {}", metadata.indices.w);
        debug!("DATE index:     {}", metadata.indices.date1);
        if let Some(d2) = metadata.indices.date2 {
            debug!("(Second) DATE index: {}", d2);
        }
        debug!("COMPLEX index:  {}", metadata.indices.complex);
        debug!("STOKES index:   {}", metadata.indices.stokes);
        debug!("FREQ index:     {}", metadata.indices.freq);
        debug!("RA index:       {}", metadata.indices.ra);
        debug!("DEC index:      {}", metadata.indices.dec);

        assert!(metadata.num_rows > 0);

        // The phase centre is described by RA and DEC if there is no SOURCE
        // table (as per the standard).
        // TODO: Check that there is no SOURCE table!
        let phase_centre = {
            let ra = fits_get_required_key(
                &mut uvfits_fptr,
                &primary_hdu,
                &format!("CRVAL{}", metadata.indices.ra),
            );
            let dec = fits_get_required_key(
                &mut uvfits_fptr,
                &primary_hdu,
                &format!("CRVAL{}", metadata.indices.dec),
            );
            RADec::from_degrees(ra, dec)
        };

        // `jd_zero` is the PZERO value of the first DATE key in the uvfits,
        // and it's supposed to "encode the Julian date at midnight of the
        // first day of the observation". That means we can round to the
        // nearest hour; doing this helps ward off float precision issues.
        let jd_zero = Epoch::from_jde_utc(0.0);

        let (all_timesteps, timestamps): (Vec<usize>, Vec<Epoch>) = metadata
            .jd_frac_timestamps
            .iter()
            .enumerate()
            .map(|(i, &jd_frac)| {
                // uvfits timestamps are in the middle of their respective
                // integration periods (centroids), so no adjustment for
                // half the integration time is needed here.
                let e = jd_zero + jd_frac;
                // Round to the nearest 10 milliseconds to avoid float
                // precision issues.
                (i, e.round(10.milliseconds()))
            })
            .unzip();
        // TODO: Determine flagging!
        let unflagged_timesteps = all_timesteps.clone();
        let all_timesteps = Vec1::try_from_vec(all_timesteps).unwrap();
        let timestamps = Vec1::try_from_vec(timestamps).unwrap();

        // Get the data's time resolution. There is a possibility that the file
        // contains only one timestep.
        let time_res = Duration::from_seconds(10.0);
        match timestamps.as_slice() {
            // Handled above; uvfits files aren't allowed to be empty.
            [] => unreachable!(),
            [t] => debug!("Only timestep (GPS): {:.2}", t.to_gpst_seconds()),
            [t0, .., tn] => {
                debug!("First good timestep (GPS): {:.2}", t0.to_gpst_seconds());
                debug!("Last good timestep  (GPS): {:.2}", tn.to_gpst_seconds());
            }
        }

        let step = metadata.num_rows / timestamps.len();

        let freq_val_str = format!("CRVAL{}", metadata.indices.freq);
        let base_freq_str: String =
            fits_get_required_key(&mut uvfits_fptr, &primary_hdu, &freq_val_str);
        let base_freq: f64 = base_freq_str.parse().unwrap();
        let base_index: isize = {
            // CRPIX might be a float. Parse it as one, then make it an int.
            let freq_val_str = format!("CRPIX{}", metadata.indices.freq);
            let f_str: String =
                fits_get_required_key(&mut uvfits_fptr, &primary_hdu, &freq_val_str);
            let f: f64 = f_str.parse().unwrap();
            f.round() as _
        };
        let freq_val_str = format!("CDELT{}", metadata.indices.freq);
        let fine_chan_width_str: String =
            fits_get_required_key(&mut uvfits_fptr, &primary_hdu, &freq_val_str);
        let freq_res: f64 = fine_chan_width_str.parse().unwrap();

        let mut fine_chan_freqs_f64 = Vec::with_capacity(metadata.num_fine_freq_chans);
        let mut fine_chan_freqs = Vec::with_capacity(metadata.num_fine_freq_chans);
        for i in 0..metadata.num_fine_freq_chans {
            let freq = (base_freq + (i as isize - base_index + 1) as f64 * freq_res).round();
            fine_chan_freqs_f64.push(freq);
            fine_chan_freqs.push(freq.round() as u64);
        }
        let fine_chan_freqs = Vec1::try_from_vec(fine_chan_freqs).unwrap();

        let obs_context = ObsContext {
            timestamps,
            all_timesteps,
            unflagged_timesteps,
            phase_centre,
            array_position,
            dut1: Duration::from_seconds(-0.1093301773071),
            tile_names,
            tile_xyzs,
            time_res,
            freq_res,
            fine_chan_freqs,
        };

        UvfitsReader {
            obs_context,
            uvfits: uvfits.to_path_buf(),
            metadata,
            step,
        }
    }

    fn read_inner(&self, mut data_fb: ArrayViewMut2<c32>, timestep: usize) {
        assert_eq!(self.metadata.num_fine_freq_chans, 1);

        let row_range_start = timestep * self.step;

        let mut uvfits = fits_open(&self.uvfits);
        fits_open_hdu(&mut uvfits, 0);
        let mut uvfits_vis = vec![0.0; self.step * 3];

        let mut status = 0;
        unsafe {
            // ffgpve = fits_read_img_flt
            fitsio_sys::ffgpve(
                uvfits.as_raw(), /* I - FITS file pointer                       */
                (1 + row_range_start)
                    .try_into()
                    .expect("not larger than i64::MAX"), /* I - group to read (1 = 1st group)           */
                1, /* I - first vector element to read (1 = 1st)  */
                uvfits_vis
                    .len()
                    .try_into()
                    .expect("not larger than i64::MAX"), /* I - number of values to read                */
                0.0,                            /* I - value for undefined pixels              */
                uvfits_vis.as_mut_ptr().cast(), /* O - array of values that are returned       */
                &mut 0,                         /* O - set to 1 if any values are null; else 0 */
                &mut status,                    /* IO - error status                           */
            );
        }
        fits_check_status(status).unwrap();

        uvfits_vis
            .chunks_exact(3)
            .zip(data_fb.iter_mut())
            .for_each(|(in_data, out_data)| {
                *out_data = c32::new(in_data[0], in_data[1]);
            });
    }
    // }
}

impl VisRead for UvfitsReader {
    fn get_obs_context(&self) -> &ObsContext {
        &self.obs_context
    }

    fn get_input_data_type(&self) -> VisInputType {
        VisInputType::Uvfits
    }

    fn read(&self, data_fb: ArrayViewMut2<c32>, timestep: usize) {
        self.read_inner(data_fb, timestep)
    }

    fn update_file(&mut self, file: PathBuf) {
        self.uvfits = file;
    }
}

struct UvfitsMetadata {
    /// The number of rows in the metafits file (hopefully equal to the number
    /// of timesteps * the number of baselines).
    num_rows: usize,

    /// The number of parameters are in each uvfits group (PCOUNT).
    pcount: usize,

    /// The... number of fine channel frequencies.
    num_fine_freq_chans: usize,

    /// The indices of various parameters (e.g. BASELINE is PTYPE4, DATE is
    /// PTYPE5, etc.)
    indices: Indices,

    /// Unique collection of JD fractions for timestamps.
    jd_frac_timestamps: Vec<Duration>,
}

impl UvfitsMetadata {
    /// Get metadata on the supplied uvfits file.
    ///
    /// This function assumes the correct HDU has already been opened (should be
    /// HDU 1, index 0).
    fn new(uvfits: &mut FitsFile, hdu: &FitsHdu) -> Self {
        let indices = Indices::new(uvfits, hdu);

        // GCOUNT tells us how many visibilities are in the file.
        let num_rows_str: String = fits_get_required_key(uvfits, hdu, "GCOUNT");
        let num_rows: usize = num_rows_str.parse().unwrap();

        // PCOUNT tells us how many parameters are in each uvfits group.
        let pcount_str: String = fits_get_required_key::<String>(uvfits, hdu, "PCOUNT");
        let pcount = pcount_str.parse::<usize>().unwrap();

        // We expect the COMPLEX index to be 2 (mandated by the standard), the
        // STOKES index to be 3, and the FREQ index to be 4. The order of these
        // indices determines the shape of the array of visibilities, and we
        // currently only support this one particular order.
        if indices.complex != 2 && indices.stokes != 3 && indices.freq != 4 {
            panic!("Bad data layout");
        }

        // NAXIS2 (COMPLEX) is how many floats are associated with a
        // polarisation. It must be either 2 or 3, as per the standard. The
        // first two floats represent the real and imag part of a complex
        // number, respectively, and the optional third is the weight. If there
        // are only 2 floats, the weight is set to 1.
        let num_floats_per_pol_str: String = fits_get_required_key(uvfits, hdu, "NAXIS2");
        let num_floats_per_pol = num_floats_per_pol_str.parse::<u8>().unwrap();
        match num_floats_per_pol {
            2 | 3 => (),
            _ => panic!("floats_per_pol wasn't 2 or 3"),
        }

        // The number of polarisations is described by the NAXIS key associated
        // with STOKES.
        let stokes_naxis_str = Cow::from(format!("NAXIS{}", indices.stokes));
        let num_pols_str: String = fits_get_required_key(uvfits, hdu, &stokes_naxis_str);
        let num_pols = num_pols_str.parse::<u8>().unwrap();

        // The pol type is described by the CRVAL key associated with STOKES.
        let stokes_crval_str = format!("CRVAL{}", indices.stokes);
        let pol_type_str: String = fits_get_required_key(uvfits, hdu, &stokes_crval_str);
        match pol_type_str.parse::<f32>() {
            Ok(pol_type) => {
                // Convert the float to an int.
                if pol_type.abs() > 127.0 {
                    panic!(
                        "STOKES {stokes_crval_str} has an unsupported value (absolute value > 127)"
                    );
                }
                let pol_type = pol_type.round() as i8;

                // We currently only support a "pol type" of -5 or -6, i.e. XX or YY.
                match (pol_type, num_pols) {
                    (-5, 1) => (),
                    _ => {
                        panic!("Bad STOKES type")
                    }
                }
            }
            Err(_) => panic!("Bad STOKES CRVAL"),
        };

        // The number of fine-frequency channels is described by the NAXIS key
        // associated with FREQ.
        let freq_naxis_str = format!("NAXIS{}", indices.freq);
        let num_fine_freq_chans_str: String = fits_get_required_key(uvfits, hdu, &freq_naxis_str);
        let num_fine_freq_chans = num_fine_freq_chans_str.parse::<usize>().unwrap();

        // Read unique group parameters (timestamps and baselines/antennas).
        let mut jd_frac_timestamp_set = HashSet::new();
        let mut jd_frac_timestamps = vec![];

        let mut group_params = Array2::zeros((num_rows, pcount));
        unsafe {
            let mut status = 0;
            // ffggpe = fits_read_grppar_flt
            fitsio_sys::ffggpe(
                uvfits.as_raw(), /* I - FITS file pointer                       */
                1,               /* I - group to read (1 = 1st group)           */
                1,               /* I - first vector element to read (1 = 1st)  */
                (pcount * num_rows)
                    .try_into()
                    .expect("not larger than i64::MAX"), /* I - number of values to read                */
                group_params.as_mut_ptr(), /* O - array of values that are returned       */
                &mut status,               /* IO - error status                           */
            );
            // Check the status.
            fits_check_status(status).unwrap();
        }

        for params in group_params.outer_iter() {
            let jd_frac = {
                let mut t = Duration::from_days(f64::from(params[usize::from(indices.date1) - 1]));
                // Use the second date, if it's there.
                if let Some(d2) = indices.date2 {
                    t += Duration::from_days(f64::from(params[usize::from(d2) - 1]));
                }
                t
            };
            if !jd_frac_timestamp_set.contains(&jd_frac.total_nanoseconds()) {
                jd_frac_timestamp_set.insert(jd_frac.total_nanoseconds());
                jd_frac_timestamps.push(jd_frac);
            }
        }

        UvfitsMetadata {
            num_rows,
            pcount,
            num_fine_freq_chans,
            indices,
            jd_frac_timestamps,
        }
    }
}

#[derive(Debug)]
struct Indices {
    /// PTYPE
    u: u8,
    /// PTYPE
    v: u8,
    /// PTYPE
    w: u8,
    /// PTYPE
    date1: u8,
    /// PTYPE
    date2: Option<u8>,
    /// CTYPE
    complex: u8,
    /// CTYPE
    stokes: u8,
    /// CTYPE
    freq: u8,
    /// CTYPE
    ra: u8,
    /// CTYPE
    dec: u8,
}

impl Indices {
    /// Find the 1-indexed indices of "PTYPE" and "CTYPE" keys we require (e.g.
    /// "UU", "VV", "WW", "RA", "DEC"). "BASELINE" will be in most uvfits files,
    /// but "ANTENNA1" and "ANTENNA2" may be used instead; exactly one of the
    /// two is ensured to be present. A second "DATE"/"_DATE" key may also be
    /// present but does not have to be.
    fn new(uvfits: &mut FitsFile, hdu: &FitsHdu) -> Self {
        // Accumulate the "PTYPE" keys.
        let mut ptypes = Vec::with_capacity(12);
        for i in 1.. {
            let ptype: Option<String> = fits_get_optional_key(uvfits, hdu, &format!("PTYPE{i}"));
            match ptype {
                Some(ptype) => ptypes.push(ptype),

                // We've found the last PTYPE.
                None => break,
            }
        }

        // We only care about UVWs, baselines and dates.
        let mut u_index = None;
        let mut v_index = None;
        let mut w_index = None;
        let mut baseline_index = None;
        let mut antenna1_index = None;
        let mut antenna2_index = None;
        let mut date1_index = None;
        let mut date2_index = None;

        for (i, key) in ptypes.into_iter().enumerate() {
            let ii = (i + 1) as u8;
            match key.as_str() {
                "UU" => {
                    if u_index.is_none() {
                        u_index = Some(ii)
                    } else {
                        warn!("Found another UU key -- only using the first");
                    }
                }
                "VV" => {
                    if v_index.is_none() {
                        v_index = Some(ii)
                    } else {
                        warn!("Found another VV key -- only using the first");
                    }
                }
                "WW" => {
                    if w_index.is_none() {
                        w_index = Some(ii)
                    } else {
                        warn!("Found another WW key -- only using the first");
                    }
                }
                "BASELINE" => {
                    if baseline_index.is_none() {
                        baseline_index = Some(ii)
                    } else {
                        warn!("Found another BASELINE key -- only using the first");
                    }
                }
                "ANTENNA1" => {
                    if antenna1_index.is_none() {
                        antenna1_index = Some(ii)
                    } else {
                        warn!("Found another ANTENNA1 key -- only using the first");
                    }
                }
                "ANTENNA2" => {
                    if antenna2_index.is_none() {
                        antenna2_index = Some(ii)
                    } else {
                        warn!("Found another ANTENNA1 key -- only using the first");
                    }
                }
                "DATE" | "_DATE" => match (date1_index, date2_index) {
                    (None, None) => date1_index = Some(ii),
                    (Some(_), None) => date2_index = Some(ii),
                    (Some(_), Some(_)) => {
                        warn!("Found more than 2 DATE/_DATE keys -- only using the first two")
                    }
                    (None, Some(_)) => unreachable!(),
                },
                _ => (),
            }
        }

        let (u, v, w, date1) = match (u_index, v_index, w_index, date1_index) {
            (Some(u), Some(v), Some(w), Some(date1)) => (u, v, w, date1),
            _ => panic!("bad indices"),
        };

        // Now find CTYPEs.
        let mut ctypes = Vec::with_capacity(12);
        for i in 2.. {
            let ctype: Option<String> = fits_get_optional_key(uvfits, hdu, &format!("CTYPE{i}"));
            match ctype {
                Some(ctype) => ctypes.push(ctype),

                // We've found the last CTYPE.
                None => break,
            }
        }

        let mut complex_index = None;
        let mut stokes_index = None;
        let mut freq_index = None;
        let mut ra_index = None;
        let mut dec_index = None;

        for (i, key) in ctypes.into_iter().enumerate() {
            let ii = (i + 2) as u8;
            match key.as_str() {
                "COMPLEX" => complex_index = Some(ii),
                "STOKES" => stokes_index = Some(ii),
                "FREQ" => freq_index = Some(ii),
                "RA" => ra_index = Some(ii),
                "DEC" => dec_index = Some(ii),
                _ => (),
            }
        }

        let (complex, stokes, freq, ra, dec) =
            match (complex_index, stokes_index, freq_index, ra_index, dec_index) {
                (Some(complex), Some(stokes), Some(freq), Some(ra), Some(dec)) => {
                    (complex, stokes, freq, ra, dec)
                }
                _ => panic!("bad indices"),
            };

        Indices {
            u,
            v,
            w,
            date1,
            date2: date2_index,
            complex,
            stokes,
            freq,
            ra,
            dec,
        }
    }
}

/// Pull out fits array-in-a-cell values; useful for e.g. STABXYZ. This function
/// assumes that the output datatype is f64, and that the fits datatype is
/// TDOUBLE, so it is not to be used generally!
fn read_cell_array(
    fits_ptr: &mut fitsio::FitsFile,
    _hdu: &fitsio::hdu::FitsHdu,
    col_name: &'static str,
    row: i64,
    n_elem: i64,
) -> Vec<f64> {
    unsafe {
        // With the column name, get the column number.
        let mut status = 0;
        let mut col_num = -1;
        let keyword = std::ffi::CString::new(col_name).expect("CString::new failed");
        // ffgcno = fits_get_colnum
        fitsio_sys::ffgcno(
            fits_ptr.as_raw(),
            0,
            keyword.as_ptr() as *mut c_char,
            &mut col_num,
            &mut status,
        );
        // Check the status.
        fits_check_status(status).unwrap();

        // Now get the specified row from that column.
        let mut array: Vec<f64> = vec![0.0; n_elem as usize];
        // ffgcv = fits_read_col
        fitsio_sys::ffgcv(
            fits_ptr.as_raw(),
            82, // TDOUBLE (fitsio.h)
            col_num,
            row + 1,
            1,
            n_elem,
            std::ptr::null_mut(),
            array.as_mut_ptr().cast(),
            &mut 0,
            &mut status,
        );
        fits_check_status(status).unwrap();

        array
    }
}
