//! Ripped out of hyperdrive.

use std::{
    ffi::{CStr, CString},
    fmt::Display,
    ptr,
};

use fitsio::{hdu::*, FitsFile};

/// Open a fits file.
#[track_caller]
pub(crate) fn fits_open<P: AsRef<std::path::Path>>(file: P) -> FitsFile {
    FitsFile::open(file.as_ref()).unwrap()
}

/// Open a fits file's HDU.
#[track_caller]
pub(crate) fn fits_open_hdu<T: DescribesHdu + Display + Copy>(
    fits_fptr: &mut FitsFile,
    hdu_description: T,
) -> FitsHdu {
    fits_fptr.hdu(hdu_description).unwrap()
}

/// Given a FITS file pointer, a HDU that belongs to it, and a keyword that may
/// or may not exist, pull out the value of the keyword, parsing it into the
/// desired type.
#[track_caller]
pub(crate) fn fits_get_optional_key<T: std::str::FromStr>(
    fits_fptr: &mut FitsFile,
    hdu: &FitsHdu,
    keyword: &str,
) -> Option<T> {
    let unparsed_value: String = match hdu.read_key(fits_fptr, keyword) {
        Ok(key_value) => key_value,
        Err(e) => match &e {
            fitsio::errors::Error::Fits(fe) => match fe.status {
                202 | 204 => return None,
                _ => panic!("bad key"),
            },
            _ => panic!("bad key"),
        },
    };

    match unparsed_value.parse() {
        Ok(parsed_value) => Some(parsed_value),
        Err(_) => panic!("couldn't convert string to type"),
    }
}

/// Given a FITS file pointer, a HDU that belongs to it, and a keyword, pull out
/// the value of the keyword, parsing it into the desired type.
#[track_caller]
pub(crate) fn fits_get_required_key<T: std::str::FromStr>(
    fits_fptr: &mut FitsFile,
    hdu: &FitsHdu,
    keyword: &str,
) -> T {
    match fits_get_optional_key(fits_fptr, hdu, keyword) {
        Some(value) => value,
        None => panic!("bad key"),
    }
}

/// Get a column from a fits file's HDU.
#[track_caller]
pub(crate) fn fits_get_col<T: fitsio::tables::ReadsCol>(
    fits_fptr: &mut FitsFile,
    hdu: &FitsHdu,
    keyword: &str,
) -> Vec<T> {
    hdu.read_col(fits_fptr, keyword).unwrap()
}

/// Given a FITS file pointer, and a keyword to a long string keyword that may
/// or may not exist, pull out the long string of the keyword. This deals with
/// FITSs CONTINUE mechanism by calling a low level fits function.
#[track_caller]
pub(crate) fn fits_get_optional_key_long_string(
    fits_fptr: &mut fitsio::FitsFile,
    hdu: &FitsHdu,
    keyword: &str,
) -> Option<String> {
    // Read the long string.
    let keyword_ffi = CString::new(keyword)
        .expect("fits_get_optional_fits_long_string: CString::new() failed for keyword");
    let long_string = unsafe {
        let mut status = 0;
        let mut long_string_ptr = ptr::null_mut();
        // ffgkls = fits_read_key_longstr
        fitsio_sys::ffgkls(
            fits_fptr.as_raw(),
            keyword_ffi.as_ptr(),
            &mut long_string_ptr,
            ptr::null_mut(),
            &mut status,
        );
        // Check the call worked!
        match status {
            0 => {
                let long_string = CStr::from_ptr(long_string_ptr)
                    .to_str()
                    .expect("fits_get_optional_key_long_string: reading C string as UTF-8 failed")
                    .to_string();
                // Free the cfitsio-allocated string. The status code passed
                // isn't useful (have a look at the source if you don't believe
                // me!)
                // fffree = fits_free_memory
                fitsio_sys::fffree(long_string_ptr.cast(), &mut 0);
                Some(long_string)
            }
            202 | 204 => None,
            _ => panic!("bad long string"),
        }
    };

    long_string
}

/// Given a FITS file pointer, and a keyword to a long string keyword, pull out
/// the long string of the keyword. This deals with FITSs CONTINUE mechanism by
/// calling a low level fits function.
#[track_caller]
pub(crate) fn _fits_get_required_key_long_string(
    fits_fptr: &mut FitsFile,
    hdu: &FitsHdu,
    keyword: &str,
) -> String {
    match fits_get_optional_key_long_string(fits_fptr, hdu, keyword) {
        Some(value) => value,
        None => panic!("bad long string"),
    }
}

/// Get the size of the image on the supplied FITS file pointer and HDU.
#[track_caller]
pub(crate) fn _fits_get_image_size<'a>(
    fits_fptr: &mut FitsFile,
    hdu: &'a FitsHdu,
) -> &'a Vec<usize> {
    match &hdu.info {
        HduInfo::ImageInfo { shape, .. } => shape,
        _ => panic!("bad image"),
    }
}

/// Given a FITS file pointer and a HDU, read the associated image.
#[track_caller]
pub(crate) fn fits_get_image<T: fitsio::images::ReadImage>(
    fits_fptr: &mut FitsFile,
    hdu: &FitsHdu,
) -> T {
    match &hdu.info {
        HduInfo::ImageInfo { .. } => hdu.read_image(fits_fptr).unwrap(),
        _ => panic!("not image"),
    }
}

/// Given a FITS file pointer and a HDU, read the associated float image.
#[track_caller]
pub(crate) fn _fits_get_float_image_into_buffer(
    fits_fptr: &mut FitsFile,
    hdu: &FitsHdu,
    buffer: &mut [f32],
) {
    unsafe {
        // Get raw ptr and length to user supplied buffer
        let buffer_len = buffer.len() as i64;
        let buffer_ptr = buffer.as_mut_ptr();

        // Call the underlying cfitsio read function for floats
        let mut status = 0;
        // ffgpv = fits_read_img
        fitsio_sys::ffgpv(
            fits_fptr.as_raw(),
            fitsio_sys::TFLOAT as _,
            1,
            buffer_len,
            ptr::null_mut(),
            buffer_ptr as *mut _,
            ptr::null_mut(),
            &mut status,
        );

        // Check fits call status
        match fitsio::errors::check_status(status) {
            Ok(_) => {}
            Err(e) => panic!("{e}"),
        }
    }
}
