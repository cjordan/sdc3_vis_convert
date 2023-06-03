// use pyo3::prelude::*;

// #[pyfunction]
// fn geocentric_pair_to_uvw(
//     earth_pos: [f64; 3],
//     phase_centre: [f64; 2],
//     offset_seconds: f64,
// ) -> [f64; 3] {
//     let time = Epoch::from_gpst_seconds(1316268778.13 + offset_seconds);
//     let lmst = (330.304110 + offset_seconds).to_radians();

//     let earth_pos = LatLngHeight {
//         longitude_rad: earth_pos[0],
//         latitude_rad: earth_pos[1],
//         height_metres: earth_pos[2],
//     };
//     let phase_centre = RADec {
//         ra: phase_centre[0],
//         dec: phase_centre[1],
//     };
//     let xyzs = [
//         XyzGeocentric {
//             x: -2564976.05653377,
//             y: 5085352.7292609,
//             z: -2861040.11633764,
//         },
//         XyzGeocentric {
//             x: -2564848.56192944,
//             y: 5085557.93305701,
//             z: -2860791.33394729,
//         },
//     ];

//     let xyzs_gd: Vec<_> = xyzs.iter().map(|xyz| xyz.to_geodetic(earth_pos)).collect();
//     let p = precess_time(
//         earth_pos.longitude_rad,
//         earth_pos.latitude_rad,
//         phase_centre,
//         time,
//         Duration::from_seconds(offset_seconds),
//     );
//     // let xyzs_gd = p.precess_xyz(&xyzs_gd);

//     let pc = phase_centre.to_hadec(lmst);
//     let uvw1 = UVW::from_xyz(xyzs_gd[0], pc);
//     let uvw2 = UVW::from_xyz(xyzs_gd[1], pc);
//     let UVW { u, v, w } = uvw2 - uvw1;

//     [u, v, w]
// }

// #[pyfunction]
// fn geocentric_pair_to_uvw(
//     earth_pos: [f64; 3],
//     phase_centre: [f64; 2],
//     date1: f64,
//     date2: f64,
//     xyz1: [f64; 3],
//     xyz2: [f64; 3],
//     offset_seconds: f64,
// ) -> [f64; 3] {
//     let earth_pos = LatLngHeight {
//         longitude_rad: earth_pos[0],
//         latitude_rad: earth_pos[1],
//         height_metres: earth_pos[2],
//     };
//     let phase_centre = RADec {
//         ra: phase_centre[0],
//         dec: phase_centre[1],
//     };

//     let time = {
//         let date = Duration::from_days(date1) + Duration::from_days(date2);
//         Epoch::from_jde_utc(0.0) + date
//     };
//     let offset_seconds = Duration::from_seconds(offset_seconds);

//     let p = precess_time(
//         earth_pos.longitude_rad,
//         earth_pos.latitude_rad,
//         phase_centre,
//         time,
//         offset_seconds,
//     );
//     let lmst = p.lmst;

//     let xyzs = [
//         XyzGeocentric {
//             x: xyz1[0],
//             y: xyz1[1],
//             z: xyz1[2],
//         },
//         XyzGeocentric {
//             x: xyz2[0],
//             y: xyz2[1],
//             z: xyz2[2],
//         },
//     ];

//     let xyzs_gd: Vec<_> = xyzs.iter().map(|xyz| xyz.to_geodetic(earth_pos)).collect();
//     // let xyzs_gd = p.precess_xyz(&xyzs_gd);

//     let pc = phase_centre.to_hadec(lmst);
//     let uvw1 = UVW::from_xyz(xyzs_gd[0], pc);
//     let uvw2 = UVW::from_xyz(xyzs_gd[1], pc);
//     let UVW { u, v, w } = uvw2 - uvw1;

//     // println!("{} {}", earth_pos, pc);
//     // dbg!(lmst.to_degrees(), offset_seconds.to_seconds());
//     // dbg!(u, v, w);

//     [u, v, w]
// }

// #[pymodule]
// fn rustuvw(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
//     m.add("__version__", env!("CARGO_PKG_VERSION"))?;
//     m.add_function(wrap_pyfunction!(geocentric_pair_to_uvw, m)?)?;

//     Ok(())
// }
