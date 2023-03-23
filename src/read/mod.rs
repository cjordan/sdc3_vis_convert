pub mod fits;
pub mod ms;
pub mod uvfits;

use std::path::PathBuf;

use marlu::c32;
use ndarray::ArrayViewMut2;

use super::{ObsContext, VisInputType};

pub trait VisRead: Sync + Send {
    fn get_obs_context(&self) -> &ObsContext;

    fn get_input_data_type(&self) -> VisInputType;

    fn read(&self, data_bf: ArrayViewMut2<c32>, timestep: usize);

    fn update_file(&mut self, file: PathBuf);
}
