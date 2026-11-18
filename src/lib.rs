// Many command entry points take long parameter lists mirroring their CLI flags.
// Converting them to args structs is a separate refactor (several are already done,
// e.g. BuildProbesArgs/SimulateArgs); until then this lint is noise.
#![allow(clippy::too_many_arguments)]

pub mod sdust;
pub mod thermodynamics;
