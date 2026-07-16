pub mod weights;
pub mod fragment;
pub mod thermo_sim;

/// Fragment length distribution parameters, shared across sampling functions.
#[derive(Clone, Copy, Debug)]
pub struct FragmentLengthParams {
    pub mean: f64,
    pub min: usize,
    pub max: usize,
}

