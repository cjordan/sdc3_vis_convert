use std::path::PathBuf;

use argmin::core::{CostFunction, Error, Executor, State};
use argmin::solver::neldermead::NelderMead;
use fitsio::FitsFile;
use marlu::hifitime::{J1900_OFFSET, SECONDS_PER_DAY};
use marlu::{
    constants::VEL_C,
    hifitime::{Duration, Epoch},
    precession::precess_time,
    rubbl_casatables::{Table, TableOpenMode},
    LatLngHeight, RADec, XyzGeocentric, UVW,
};
use rayon::prelude::*;

use sdc3::fits::*;

// fn main() {
//     let n: usize = std::env::args().nth(1).unwrap().parse().unwrap();
//     let mut v = vec![XyzGeodetic::default(); n];
//     for xyz in v.iter_mut() {
//         xyz.x = 1.0 * n as f64;
//         xyz.y = 2.0 * n as f64;
//         xyz.z = 3.0 * n as f64;
//     }

//     let start = std::time::Instant::now();
//     let v2 = v
//         .iter()
//         .copied()
//         .map(|xyz| UVW::from_xyz_inner(xyz, 1.0, 2.0, 3.0, 4.0));
//     dbg!(std::time::Instant::now() - start);
//     dbg!(v2
//         .into_iter()
//         .flat_map(|uvw| [uvw.u, uvw.v, uvw.w])
//         .sum::<f64>());
// }

// fn main() {
//     let xyz1 = XyzGeocentric {
//         x: -2564976.05653377,
//         y: 5085352.7292609,
//         z: -2861040.11633764,
//     };
//     let xyz2 = XyzGeocentric {
//         x: -2564848.56192944,
//         y: 5085557.93305701,
//         z: -2860791.33394729,
//     };
//     let phase_centre = RADec::from_degrees(0.0, -30.0)
//         // .to_hadec(330.304110_f64.to_radians());
//         .to_hadec(330.70_f64.to_radians());
//     let earth_pos = LatLngHeight {
//         longitude_rad: 116.764_f64.to_radians(),
//         latitude_rad: -26.825_f64.to_radians(),
//         height_metres: 0.0,
//         // longitude_rad: 116.7644482_f64.to_radians(),
//         // latitude_rad: -26.82472208_f64.to_radians(),
//         // height_metres: 377.827,
//     };

//     let uvw1 = UVW::from_xyz(xyz1.to_geodetic(earth_pos), phase_centre);
//     let uvw2 = UVW::from_xyz(xyz2.to_geodetic(earth_pos), phase_centre);
//     dbg!(uvw1, uvw2, uvw2 - uvw1);

//     let ref_uvw = UVW {
//         u: -8.0535693e-07 * VEL_C,
//         v: 7.3264675e-07 * VEL_C,
//         w: -3.9071608e-07 * VEL_C,
//     };
//     dbg!(ref_uvw);
// }

fn main() {
    if let Err(e) = try_main() {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

const NUM_STATIONS: usize = 512;
const NUM_BASELINES: usize = (512 * 511) / 2;
const NUM_TIMESTEPS: usize = 1440;
const PCOUNT: usize = 9;

lazy_static::lazy_static! {
    static ref ARRAY_POS: LatLngHeight = LatLngHeight {
        longitude_rad: 116.764_f64.to_radians(),
        latitude_rad: -26.825_f64.to_radians(),
        height_metres: 0.0,
    };

    static ref PHASE_CENTRE: RADec = RADec {
        ra: 0.0,
        dec: -30.0_f64.to_radians(),
    };
}

pub struct GenUvw {
    xyz1: XyzGeocentric,
    xyz2: XyzGeocentric,
    ref_uvw: UVW,
    time: Epoch,
}

impl CostFunction for GenUvw {
    type Param = f64;
    type Output = f64;

    fn cost(&self, time_offset_seconds: &Self::Param) -> Result<Self::Output, Error> {
        let time_offset_seconds = Duration::from_seconds(*time_offset_seconds);
        let p = precess_time(
            ARRAY_POS.longitude_rad,
            ARRAY_POS.latitude_rad,
            *PHASE_CENTRE,
            self.time,
            time_offset_seconds,
        );

        let xyzs_gd = [
            self.xyz1.to_geodetic(*ARRAY_POS),
            self.xyz2.to_geodetic(*ARRAY_POS),
        ];
        // let xyzs_gd = p.precess_xyz(&xyzs_gd);

        let pc = PHASE_CENTRE.to_hadec(p.lmst);
        let uvw1 = UVW::from_xyz(xyzs_gd[0], pc);
        let uvw2 = UVW::from_xyz(xyzs_gd[1], pc);
        let UVW { u, v, w } = uvw2 - uvw1;
        // println!("{u} {v} {w}");
        // println!("{} {} {}", self.ref_uvw.u, self.ref_uvw.v, self.ref_uvw.w);

        let diff = (self.ref_uvw.u - u).powi(2)
            + (self.ref_uvw.v - v).powi(2)
            + (self.ref_uvw.w - w).powi(2);
        // dbg!(diff);

        Ok(diff)
    }
}

enum VisContainer {
    Uvfits {
        fptr: FitsFile,
        group_params_buffer: Vec<f32>,
    },

    Ms(PathBuf),
}

impl From<String> for VisContainer {
    fn from(value: String) -> Self {
        let value = PathBuf::from(value);
        match value.extension().and_then(|os_str| os_str.to_str()) {
            Some("uvf" | "uvfits") => VisContainer::Uvfits {
                fptr: FitsFile::open(value).unwrap(),
                group_params_buffer: vec![0.0; PCOUNT * NUM_BASELINES],
            },
            Some("ms") => VisContainer::Ms(value),
            Some(ext) => {
                eprintln!("Unrecognised extension '{ext}'");
                std::process::exit(1);
            }
            None => {
                eprintln!("No extension on input vis file");
                std::process::exit(1);
            }
        }
    }
}

impl VisContainer {
    fn get_xyzs(&mut self) -> Result<Vec<XyzGeocentric>, Box<dyn std::error::Error>> {
        match self {
            Self::Uvfits { fptr, .. } => {
                let hdu = fits_open_hdu(fptr, 0)?;
                assert_eq!(
                    PCOUNT,
                    fits_get_required_key::<usize>(fptr, &hdu, "PCOUNT")?
                );
                assert_eq!(
                    NUM_BASELINES * NUM_TIMESTEPS,
                    fits_get_required_key::<usize>(fptr, &hdu, "GCOUNT")?
                );

                let mut xyzs: Vec<XyzGeocentric> = Vec::with_capacity(NUM_STATIONS);
                let hdu = fits_open_hdu(fptr, "AIPS AN")?;
                for i in 0..NUM_STATIONS {
                    let fits_xyz = read_cell_array(
                        fptr,
                        &hdu,
                        "STABXYZ",
                        i.try_into().expect("not larger than i64::MAX"),
                        3,
                    )?;
                    xyzs.push(XyzGeocentric {
                        x: fits_xyz[0],
                        y: fits_xyz[1],
                        z: fits_xyz[2],
                    });
                }

                Ok(xyzs)
            }

            Self::Ms(pb) => {
                let mut antenna_table = Table::open(pb.join("ANTENNA"), TableOpenMode::Read)?;
                let mut casacore_positions = Vec::with_capacity(antenna_table.n_rows() as usize);
                antenna_table.for_each_row(|row| {
                    let pos: Vec<f64> = row.get_cell("POSITION")?;
                    let pos_xyz = XyzGeocentric {
                        x: pos[0],
                        y: pos[1],
                        z: pos[2],
                    };
                    casacore_positions.push(pos_xyz);
                    Ok(())
                })?;

                Ok(casacore_positions)
            }
        }
    }

    fn read_metadata_for_timestep(
        &mut self,
        i_timestep: usize,
    ) -> Result<(Vec<Epoch>, Vec<UVW>), Box<dyn std::error::Error>> {
        match self {
            VisContainer::Uvfits {
                fptr,
                group_params_buffer,
            } => {
                // Ensure we're on the data-containing HDU.
                let _ = fptr.primary_hdu()?;

                let i_row = NUM_BASELINES * i_timestep + 1;
                unsafe {
                    let mut status = 0;
                    // ffggpe = fits_read_grppar_flt
                    fitsio_sys::ffggpe(
                        fptr.as_raw(), /* I - FITS file pointer                       */
                        i_row.try_into().expect("not larger than i64::MAX"), /* I - group to read (1 = 1st group)           */
                        1, /* I - first vector element to read (1 = 1st)  */
                        group_params_buffer
                            .len()
                            .try_into()
                            .expect("not larger than i64::MAX"), /* I - number of values to read                */
                        group_params_buffer.as_mut_ptr(), /* O - array of values that are returned       */
                        &mut status, /* IO - error status                           */
                    );
                    assert_eq!(status, 0);
                }

                Ok(group_params_buffer
                    .chunks_exact(PCOUNT)
                    .map(|bl_row| {
                        let ref_uvw = UVW {
                            u: bl_row[0] as f64,
                            v: bl_row[1] as f64,
                            w: bl_row[2] as f64,
                        } * VEL_C;

                        let time = Epoch::from_jde_utc(0.0)
                            + Duration::from_days(bl_row[3] as f64)
                            + Duration::from_days(bl_row[4] as f64);

                        (time, ref_uvw)
                    })
                    .unzip())
            }

            VisContainer::Ms(pb) => {
                let mut main_table = Table::open(pb, TableOpenMode::Read)?;

                let mut times = Vec::with_capacity(NUM_BASELINES);
                let mut ref_uvws = Vec::with_capacity(NUM_BASELINES);

                let range_start = (NUM_BASELINES * i_timestep) as u64;
                let range_end = (NUM_BASELINES * (i_timestep + 1)) as u64;
                main_table.for_each_row_in_range(range_start..range_end, |row| {
                    let ref_uvw: Vec<f64> = row.get_cell("UVW")?;
                    ref_uvws.push(UVW {
                        u: ref_uvw[0],
                        v: ref_uvw[1],
                        w: ref_uvw[2],
                    });

                    let time: f64 = row.get_cell("TIME")?;
                    // casacore stores the times as UTC seconds... but with an
                    // offset.
                    times.push(Epoch::from_utc_seconds(
                        time - J1900_OFFSET * SECONDS_PER_DAY,
                    ));

                    Ok(())
                })?;

                Ok((times, ref_uvws))
            }
        }
    }
}

fn try_main() -> Result<(), Box<dyn std::error::Error>> {
    let data = std::env::args().nth(1).unwrap_or_else(|| {
        eprintln!("Supply a vis file");
        std::process::exit(1);
    });
    let i_timestep = std::env::args()
        .nth(2)
        .unwrap_or_else(|| {
            eprintln!("Supply a timestep index");
            std::process::exit(1);
        })
        .parse()?;

    let mut vis = VisContainer::from(data);
    let xyzs = vis.get_xyzs()?;

    let mut time_offsets = Vec::with_capacity(NUM_BASELINES);

    // for i_timestep in 0..NUM_TIMESTEPS {
    //     // if i_timestep != 1439 {
    //     if i_timestep > 0 {
    //         continue;
    //     }
    dbg!(i_timestep);
    let (times, ref_uvws) = vis.read_metadata_for_timestep(i_timestep)?;

    // // Serial.
    // {
    //     let mut i_station1 = 0;
    //     let mut i_station2 = 0;
    //     let mut xyz1 = xyzs[i_station1];
    //     let mut xyz2;
    //     for (i_bl, group_params) in group_params.chunks_exact(pcount).enumerate() {
    //         if i_bl % 10000 == 0 {
    //             println!("i_bl: {i_bl}");
    //         }
    //         i_station2 += 1;
    //         if i_station2 == NUM_STATIONS {
    //             i_station1 += 1;
    //             i_station2 = i_station1 + 1;
    //             xyz1 = xyzs[i_station1];
    //         }
    //         xyz2 = xyzs[i_station2];

    //         let ref_uvw = UVW {
    //             u: group_params[0] as f64,
    //             v: group_params[1] as f64,
    //             w: group_params[2] as f64,
    //         } * VEL_C;
    //         let date1 = group_params[3] as f64;
    //         let date2 = group_params[4] as f64;
    //         let time =
    //             Epoch::from_jde_utc(0.0) + Duration::from_days(date1) + Duration::from_days(date2);
    //         dbg!(ref_uvw);

    //         let offset = solve_time_offset(xyz1, xyz2, ref_uvw, time);
    //         time_offsets.push(offset);
    //     }
    // }

    // Parallel
    {
        times
            .par_iter()
            .zip(ref_uvws)
            .enumerate()
            .map(|(i_bl, (time, ref_uvw))| {
                let (i_station1, i_station2) =
                    cross_correlation_baseline_to_tiles(NUM_STATIONS, i_bl);
                let xyz1 = xyzs[i_station1];
                let xyz2 = xyzs[i_station2];

                solve_time_offset(xyz1, xyz2, ref_uvw, *time)
            })
            .collect_into_vec(&mut time_offsets);
        let av = time_offsets.iter().sum::<f64>() / time_offsets.len() as f64;
        let std = (time_offsets.iter().map(|o| (*o - av).powi(2)).sum::<f64>()
            / time_offsets.len() as f64)
            .sqrt();
        dbg!(av, std);
    }
    // }

    Ok(())
}

fn solve_time_offset(xyz1: XyzGeocentric, xyz2: XyzGeocentric, ref_uvw: UVW, time: Epoch) -> f64 {
    let cost = GenUvw {
        xyz1,
        xyz2,
        ref_uvw,
        time,
    };

    let solver = NelderMead::new(vec![-0.95, -0.96])
        .with_sd_tolerance(1e-20)
        .unwrap();

    // Run solver
    let res = Executor::new(cost, solver)
        .configure(|state| state.max_iters(100))
        // .add_observer(SlogLogger::term(), ObserverMode::Always)
        .run()
        .unwrap();

    // println!("{}", res.state.get_best_param().unwrap());
    // println!("{}", res);

    *res.state.get_best_param().unwrap()
}

fn cross_correlation_baseline_to_tiles(total_num_tiles: usize, baseline: usize) -> (usize, usize) {
    let n = (total_num_tiles - 1) as f64;
    let bl = baseline as f64;
    let tile1 = (-0.5 * (4.0 * n * (n + 1.0) - 8.0 * bl + 1.0).sqrt() + n + 0.5).floor();
    let tile2 = bl - tile1 * (n - (tile1 + 1.0) / 2.0) + 1.0;
    (tile1 as usize, tile2 as usize)
}
