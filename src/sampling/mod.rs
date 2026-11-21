pub mod thermo_sim;
pub mod weights;

/// Fragment length distribution parameters, shared across sampling functions.
#[derive(Clone, Copy, Debug)]
pub struct FragmentLengthParams {
    pub mean: f64,
    pub min: usize,
    pub max: usize,
}
