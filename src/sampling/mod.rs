pub mod thermo_sim;
pub mod weights;

/// Shape of the fragment length distribution.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FragmentLengthShape {
    /// Normal, truncated to `[min, max]` by redrawing.
    Normal,
    /// Log-normal, truncated the same way. Sonication and enzymatic fragmentation give
    /// right-skewed length distributions, which a Normal cannot represent: it is
    /// symmetric, so matching a long upper tail forces an equally long lower one.
    LogNormal,
}

impl FragmentLengthShape {
    pub fn from_str(s: &str) -> anyhow::Result<Self> {
        match s.to_ascii_lowercase().as_str() {
            "normal" => Ok(Self::Normal),
            "lognormal" | "log-normal" => Ok(Self::LogNormal),
            other => anyhow::bail!(
                "Unknown fragment length distribution '{}'. Use 'normal' or 'lognormal'.",
                other
            ),
        }
    }
}

/// Fragment length distribution parameters, shared across sampling functions.
#[derive(Clone, Copy, Debug)]
pub struct FragmentLengthParams {
    pub mean: f64,
    pub min: usize,
    pub max: usize,
    /// Standard deviation of the length distribution, in bases.
    ///
    /// `None` keeps the historical default of `(max - min) / 6` — the span read as a
    /// +/-3 sigma interval. That is a guess about the user's intent, not a measurement,
    /// so it is worth setting explicitly when the real distribution is known.
    pub sd: Option<f64>,
    pub shape: FragmentLengthShape,
}

impl FragmentLengthParams {
    /// Standard deviation to use, applying the historical default when unset.
    pub fn sd_or_default(&self) -> f64 {
        self.sd.unwrap_or((self.max as f64 - self.min as f64) / 6.0)
    }
}
