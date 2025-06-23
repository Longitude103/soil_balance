mod boundary_cond;
mod daily_inputs;
mod hydrus;
mod root_uptake;
mod soil_hydra;

pub use hydrus::Hydrus1D;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let mut model = Hydrus1D::default();
        model.apply_bc();
        model.run();
    }
}
