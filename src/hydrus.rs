use crate::boundary_cond::{BottomBoundary, BoundaryParams};
use crate::daily_inputs::DailyInputs;
use crate::root_uptake::RootParams;
use crate::soil_hydra::SoilProfile;
use nalgebra::{DMatrix, DVector};

// Simulation parameters
struct SimParams {
    dz: f64,              // Spatial step [cm]
    dt: f64,              // Time step [day]
    z_max: f64,           // Soil column depth [cm]
    t_max: f64,           // Simulation time [day]
    root_zone_depth: f64, // Root zone depth [cm],
}

impl SimParams {
    fn new() -> Self {
        SimParams {
            dz: 10.0,              // 1 cm spatial resolution
            dt: 0.001,             // 0.001 day time step
            z_max: 150.0,          // 150 cm soil column to include groundwater
            t_max: 30.0,           // 1 day simulation
            root_zone_depth: 30.0, // 30 cm root zone depth
        }
    }
}

// Finite difference solver for Richards equation with root water uptake and variable BCs
pub struct Hydrus1D {
    soil: SoilProfile,
    root: RootParams,
    bc: BoundaryParams,
    sim: SimParams,
    inputs: DailyInputs,
    h: DVector<f64>, // Pressure head at each node [cm]
    nz: usize,       // Number of spatial nodes
}

impl Hydrus1D {
    pub fn new(soil_profile: SoilProfile) -> Self {
        let sim = SimParams::new();
        let nz = (sim.z_max / sim.dz).ceil() as usize + 1;
        // Initialize pressure head based on bottom boundary
        let h = DVector::from_fn(nz, |i, _| {
            let z = i as f64 * sim.dz;
            if z <= 60.0 { -100.0 } else { -100.0 }
        });

        let mut model = Hydrus1D {
            soil: soil_profile,
            root: RootParams::new(""),
            bc: BoundaryParams::new(),
            inputs: DailyInputs::new(),
            sim,
            h,
            nz,
        };

        // Apply initial groundwater condition if needed
        if let BottomBoundary::Groundwater(h_bot) = model.bc.bottom_boundary {
            model.h[nz - 1] = h_bot;
        }

        model
    }

    // Initialize boundary conditions
    pub fn apply_bc(&mut self) {
        if let BottomBoundary::Groundwater(h_bot) = self.bc.bottom_boundary {
            let h = self.h.as_mut_slice();
            h[self.nz - 1] = h_bot; // Apply Dirichlet for groundwater
        }
        // No action for free drainage
    }

    // Compute flux at a given depth [cm/day]
    fn compute_flux(&self, z: f64) -> f64 {
        let dz = self.sim.dz;
        // Find nodes bracketing the depth
        let i_lower = (z / dz).floor() as usize;
        let i_upper = i_lower + 1;
        if i_lower >= self.nz - 1 {
            // At bottom boundary
            match self.bc.bottom_boundary {
                BottomBoundary::Groundwater(_) => {
                    // Use backward difference for Dirichlet
                    let h_n = self.h[self.nz - 1];
                    let h_nm1 = self.h[self.nz - 2];
                    let k_n = self.soil.k(h_n, z);
                    let dh_dz = (h_n - h_nm1) / dz;
                    -k_n * (dh_dz + 1.0)
                }
                BottomBoundary::FreeDrainage => {
                    // q = -K(h_n) for free drainage
                    let h_n = self.h[self.nz - 1];
                    -self.soil.k(h_n, z)
                }
            }
        } else {
            // Interpolate pressure head and conductivity
            let z_lower = i_lower as f64 * dz;
            let z_upper = i_upper as f64 * dz;
            let frac = (z - z_lower) / (z_upper - z_lower);
            let h_lower = self.h[i_lower];
            let h_upper = self.h[i_upper];
            let k_lower = self.soil.k(h_lower, z_lower);
            let k_upper = self.soil.k(h_upper, z_upper);
            let k = k_lower + frac * (k_upper - k_lower);
            let dh_dz = (h_upper - h_lower) / dz;
            -k * (dh_dz + 1.0)
        }
    }

    // Solve one time step using implicit finite difference with nalgebra
    fn step(&mut self, time: f64) -> (f64, f64, f64, f64, f64, f64, f64) {
        let dt = self.sim.dt;
        let dz = self.sim.dz;
        let nz = self.nz;

        // Initialize tridiagonal matrix A and right-hand side b
        let mut a = DMatrix::zeros(nz, nz);
        let mut b = DVector::zeros(nz);

        // Previous time step pressure heads
        let h_old = self.h.clone();

        // Get daily inputs
        let rainfall = self.inputs.get_daily_value(time, &self.inputs.rainfall);
        let irrigation = self.inputs.get_daily_value(time, &self.inputs.irrigation);
        let evap = self.inputs.get_daily_value(time, &self.inputs.evaporation);
        let tp = self
            .inputs
            .get_daily_value(time, &self.inputs.transpiration);

        // Top boundary condition (flux or Dirichlet)
        let (top_flux, is_flux) = self.bc.top_flux(time, h_old[0], &self.inputs);

        for i in 0..nz {
            let z_i = i as f64 * dz;
            if i == 0 {
                if is_flux {
                    // Neumann BC: q = -K(h) * (dh/dz + 1) = top_flux
                    let k_i = self.soil.k(h_old[i], z_i);
                    let z_ip = (i + 1) as f64 * dz;
                    let k_ip = 0.5 * (k_i + self.soil.k(h_old[i + 1], z_ip));
                    a[(i, i)] = -k_ip / dz;
                    a[(i, i + 1)] = k_ip / dz;
                    b[i] = k_ip - top_flux;
                    // b[i] = k_ip + top_flux;
                } else {
                    // Dirichlet BC: h = h_crit
                    a[(i, i)] = 1.0;
                    b[i] = self.bc.h_crit;
                }
            } else if i == nz - 1 {
                match self.bc.bottom_boundary {
                    BottomBoundary::Groundwater(h_bot) => {
                        // Dirichlet BC for groundwater
                        a[(i, i)] = 1.0;
                        b[i] = h_bot;
                    }
                    BottomBoundary::FreeDrainage => {
                        // Free drainage: dh/dz = 0 (h_n = h_{n-1})
                        a[(i, i)] = 1.0;
                        a[(i, i - 1)] = -1.0;
                        b[i] = 0.0;
                    }
                }
            } else {
                let h_i = h_old[i];
                let h_ip1 = h_old[i + 1];
                let h_im1 = h_old[i - 1];

                // Hydraulic conductivity at interfaces
                let z_ip = (i + 1) as f64 * dz;
                let z_im = (i - 1) as f64 * dz;
                let k_ip = 0.5 * (self.soil.k(h_i, z_i) + self.soil.k(h_ip1, z_ip));
                let k_im = 0.5 * (self.soil.k(h_i, z_i) + self.soil.k(h_im1, z_im));

                // Specific moisture capacity
                let c_i = self.soil.c(h_i, z_i);

                // Root water uptake
                let root_density = self.inputs.get_root_density(time, z_i);
                let s_i = self.root.s(h_i, tp, root_density);

                // Matrix coefficients (implicit scheme)
                a[(i, i - 1)] = -k_im / dz.powi(2);
                a[(i, i)] = c_i / dt + (k_ip + k_im) / dz.powi(2);
                a[(i, i + 1)] = -k_ip / dz.powi(2);

                // Right-hand side
                b[i] = c_i * h_old[i] / dt - (k_ip - k_im) / dz - s_i;
            }
        }

        // Solve A * h_new = b
        let h_new = a.lu().solve(&b).expect("Linear system solve failed");
        self.h = h_new;
        self.apply_bc();

        // Compute fluxes
        let bottom_flux = self.compute_flux(self.sim.z_max);
        let root_zone_flux = self.compute_flux(self.sim.root_zone_depth);

        (
            top_flux,
            rainfall,
            irrigation,
            evap,
            tp,
            bottom_flux,
            root_zone_flux,
        )
    }

    // Run simulation and output results
    pub(crate) fn run(&mut self) {
        let nt = (self.sim.t_max / self.sim.dt).ceil() as usize;
        println!(
            "Time [day], Depth [cm], Pressure Head [cm], Water Content [-], Root Uptake [1/day], Root Density [-], Top Flux [cm/day], Rainfall [cm/day], Irrigation [cm/day], Evaporation [cm/day], Transpiration [cm/day], Bottom Flux [cm/day], Root Zone Flux [cm/day]"
        );
        for t in 0..=nt {
            let time = t as f64 * self.sim.dt;
            let (top_flux, rainfall, irrigation, evap, tp, bottom_flux, root_zone_flux) = if t < nt
            {
                self.step(time)
            } else {
                let flux = self.bc.top_flux(time, self.h[0], &self.inputs).0;
                let bottom_flux = self.compute_flux(self.sim.z_max);
                let root_zone_flux = self.compute_flux(self.sim.root_zone_depth);
                (
                    flux,
                    self.inputs.get_daily_value(time, &self.inputs.rainfall),
                    self.inputs.get_daily_value(time, &self.inputs.irrigation),
                    self.inputs.get_daily_value(time, &self.inputs.evaporation),
                    self.inputs
                        .get_daily_value(time, &self.inputs.transpiration),
                    bottom_flux,
                    root_zone_flux,
                )
            };
            if t % (nt / 10) == 0 || t == nt {
                for i in 0..self.nz {
                    let z = i as f64 * self.sim.dz;
                    let theta = self.soil.theta(self.h[i], z);
                    let root_density = self.inputs.get_root_density(time, z);
                    let s = self.root.s(self.h[i], tp, root_density);
                    let flux = if i == 0 { top_flux } else { 0.0 };
                    let b_flux = if i == self.nz - 1 { bottom_flux } else { 0.0 };
                    let rz_flux = if (z - self.sim.root_zone_depth).abs() < 1e-6 {
                        root_zone_flux
                    } else {
                        0.0
                    };
                    println!(
                        "{:.3}, {:.1}, {:.2}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}",
                        time,
                        z,
                        self.h[i],
                        theta,
                        s,
                        root_density,
                        flux,
                        rainfall,
                        irrigation,
                        evap,
                        tp,
                        b_flux,
                        rz_flux
                    );
                }
            }
        }
    }
}

// #[cfg(test)]
// mod tests {
//     use crate::Hydrus1D;
//     use crate::hydrus::SimParams;
//     use approx::assert_relative_eq;

//     use serde::Deserialize;
//     use std::collections::HashMap;

//     // Mock SoilParams
//     #[derive(Clone, Copy, Deserialize)]
//     struct SoilParams {
//         theta_r: f64,
//         theta_s: f64,
//         alpha: f64,
//         n: f64,
//         ks: f64,
//     }

//     impl SoilParams {
//         fn from_toml(soil_name: &str) -> Self {
//             let params = mock_soil_params();
//             params
//                 .get(soil_name)
//                 .cloned()
//                 .unwrap_or_else(|| panic!("Soil '{}' not found", soil_name))
//         }

//         fn theta(&self, h: f64) -> f64 {
//             if h >= 0.0 {
//                 self.theta_s
//             } else {
//                 let m = 1.0 - 1.0 / self.n;
//                 let ah = self.alpha * h.abs();
//                 let denom = (1.0 + ah.powf(self.n)).powf(m);
//                 self.theta_r + (self.theta_s - self.theta_r) / denom
//             }
//         }

//         fn k(&self, h: f64) -> f64 {
//             if h >= 0.0 {
//                 self.ks
//             } else {
//                 let m = 1.0 - 1.0 / self.n;
//                 let ah = self.alpha * h.abs();
//                 let se = (1.0 / (1.0 + ah.powf(self.n))).powf(m);
//                 self.ks * se.sqrt() * (1.0 - (1.0 - se.powf(1.0 / m)).powf(m)).powi(2)
//             }
//         }

//         fn c(&self, h: f64) -> f64 {
//             if h >= 0.0 {
//                 0.0
//             } else {
//                 let m = 1.0 - 1.0 / self.n;
//                 let ah = self.alpha * h.abs();
//                 let denom = (1.0 + ah.powf(self.n)).powf(m + 1.0);
//                 self.alpha * (self.theta_s - self.theta_r) * m * self.n * ah.powf(self.n - 1.0)
//                     / denom
//             }
//         }
//     }

//     // Mock SoilLayer and SoilProfile
//     #[derive(Clone)]
//     struct SoilLayer {
//         z_min: f64,
//         z_max: f64,
//         params: SoilParams,
//     }

//     #[derive(Clone)]
//     struct SoilProfile {
//         layers: Vec<SoilLayer>,
//     }

//     impl SoilProfile {
//         fn new(soil_layers: Vec<(f64, f64, &str)>) -> Self {
//             let layers = soil_layers
//                 .into_iter()
//                 .map(|(z_min, z_max, soil_name)| SoilLayer {
//                     z_min,
//                     z_max,
//                     params: SoilParams::from_toml(soil_name),
//                 })
//                 .collect();
//             SoilProfile { layers }
//         }

//         fn params_at_depth(&self, z: f64) -> &SoilParams {
//             for layer in &self.layers {
//                 if z >= layer.z_min && z <= layer.z_max {
//                     return &layer.params;
//                 }
//             }
//             &self.layers.last().expect("No soil layers defined").params
//         }

//         fn theta(&self, h: f64, z: f64) -> f64 {
//             self.params_at_depth(z).theta(h)
//         }

//         fn k(&self, h: f64, z: f64) -> f64 {
//             self.params_at_depth(z).k(h)
//         }

//         fn c(&self, h: f64, z: f64) -> f64 {
//             self.params_at_depth(z).c(h)
//         }
//     }

//     // Mock RootParams
//     #[derive(Clone, Copy)]
//     struct RootParams {
//         h1: f64,
//         h2: f64,
//         h3: f64,
//         h4: f64,
//     }

//     impl RootParams {
//         fn new(_crop: &str) -> Self {
//             RootParams {
//                 h1: -10.0,
//                 h2: -25.0,
//                 h3: -200.0,
//                 h4: -8000.0,
//             }
//         }

//         fn s(&self, h: f64, tp: f64, root_density: f64) -> f64 {
//             let alpha = if h > self.h1 || h < self.h4 {
//                 0.0
//             } else if h <= self.h2 && h >= self.h3 {
//                 1.0
//             } else if h > self.h2 && h <= self.h1 {
//                 (h - self.h1) / (self.h2 - self.h1)
//             } else {
//                 (h - self.h4) / (self.h3 - self.h4)
//             };
//             alpha * tp * root_density
//         }
//     }

//     // Mock BoundaryParams
//     #[derive(Clone)]
//     struct BoundaryParams {
//         h_crit: f64,
//         h_top: f64,
//         h_bot: f64,
//     }

//     impl BoundaryParams {
//         fn new() -> Self {
//             BoundaryParams {
//                 h_crit: -15000.0,
//                 h_top: -10.0,
//                 h_bot: 0.0,
//             }
//         }

//         fn top_flux(&self, time: f64, h_top: f64, inputs: &DailyInputs) -> (f64, bool) {
//             let rainfall = inputs.get_daily_value(time, &inputs.rainfall);
//             let irrigation = inputs.get_daily_value(time, &inputs.irrigation);
//             let total_influx = rainfall + irrigation;
//             if total_influx > 0.0 {
//                 (total_influx, true)
//             } else {
//                 let evap = inputs.get_daily_value(time, &inputs.evaporation);
//                 if h_top > self.h_crit {
//                     (-evap, true)
//                 } else {
//                     (0.0, false)
//                 }
//             }
//         }
//     }

//     // Mock DailyInputs
//     #[derive(Clone)]
//     struct DailyInputs {
//         rainfall: Vec<f64>,
//         irrigation: Vec<f64>,
//         evaporation: Vec<f64>,
//         transpiration: Vec<f64>,
//         root_dist: Vec<Vec<(f64, f64)>>,
//     }

//     impl DailyInputs {
//         fn new() -> Self {
//             DailyInputs {
//                 rainfall: vec![2.0, 0.0],
//                 irrigation: vec![1.0, 0.0],
//                 evaporation: vec![0.3, 0.5],
//                 transpiration: vec![0.5, 0.4],
//                 root_dist: vec![
//                     vec![(0.0, 0.07), (10.0, 0.02), (20.0, 0.01), (30.0, 0.0)],
//                     vec![(0.0, 0.04), (10.0, 0.03), (20.0, 0.02), (30.0, 0.0)],
//                 ],
//             }
//         }

//         fn get_daily_value(&self, time: f64, values: &Vec<f64>) -> f64 {
//             let day = time.floor() as usize;
//             if day < values.len() {
//                 values[day]
//             } else {
//                 *values.last().unwrap_or(&0.0)
//             }
//         }

//         fn get_root_density(&self, time: f64, z: f64) -> f64 {
//             let day = time.floor() as usize;
//             let default_dist = vec![(0.0, 0.0)];
//             let profile = if day < self.root_dist.len() {
//                 &self.root_dist[day]
//             } else {
//                 self.root_dist.last().unwrap_or(&default_dist)
//             };

//             for i in 0..profile.len() - 1 {
//                 let (z1, d1) = profile[i];
//                 let (z2, d2) = profile[i + 1];
//                 if z >= z1 && z <= z2 {
//                     return d1 + (d2 - d1) * (z - z1) / (z2 - z1);
//                 }
//             }
//             profile.last().map(|&(_, d)| d).unwrap_or(0.0)
//         }
//     }

//     // Mock soil parameters
//     fn mock_soil_params() -> HashMap<&'static str, SoilParams> {
//         let mut params = HashMap::new();
//         params.insert(
//             "loam",
//             SoilParams {
//                 theta_r: 0.078,
//                 theta_s: 0.43,
//                 alpha: 0.036,
//                 n: 1.56,
//                 ks: 24.96,
//             },
//         );
//         params.insert(
//             "loamy-sand",
//             SoilParams {
//                 theta_r: 0.057,
//                 theta_s: 0.41,
//                 alpha: 0.124,
//                 n: 2.28,
//                 ks: 350.2,
//             },
//         );
//         params.insert(
//             "sand",
//             SoilParams {
//                 theta_r: 0.045,
//                 theta_s: 0.43,
//                 alpha: 0.145,
//                 n: 2.68,
//                 ks: 712.8,
//             },
//         );
//         params
//     }

//     // Test SimParams
//     #[test]
//     fn test_sim_params_new() {
//         let sim = SimParams::new();
//         assert_eq!(sim.dz, 10.0);
//         assert_eq!(sim.dt, 0.1);
//         assert_eq!(sim.z_max, 150.0);
//         assert_eq!(sim.t_max, 30.0);
//         assert_eq!(sim.root_zone_depth, 30.0);
//     }

//     // Test Hydrus1D methods
//     #[test]
//     fn test_hydrus1d_new() {
//         let model = Hydrus1D::new();
//         let sim = model.sim;
//         let nz = ((sim.z_max / sim.dz).ceil() as usize) + 1;
//         assert_eq!(model.nz, nz);
//         assert_eq!(model.h.len(), nz);
//         assert_eq!(model.soil.layers.len(), 3);
//         assert_eq!(model.soil.layers[0].z_min, 0.0);
//         assert_eq!(model.soil.layers[0].z_max, 60.0);
//         assert_eq!(model.soil.layers[1].z_min, 60.0);
//         assert_eq!(model.soil.layers[1].z_max, 80.0);
//         assert_eq!(model.soil.layers[2].z_min, 80.0);
//         assert_eq!(model.soil.layers[2].z_max, 150.0);
//         // Check initial pressure heads
//         assert_eq!(model.h[nz - 1], 0.0); // Groundwater at z = 150 cm
//         assert_eq!(model.h[0], -33.0); // Loam at z = 0 cm
//         assert_eq!(model.h[6], -33.0); // Loam at z = 60 cm
//         assert_eq!(model.h[7], -15.0); // Loamy-sand at z = 70 cm
//         assert_eq!(model.h[10], -15.0); // Sand at z = 100 cm
//     }

//     #[test]
//     fn test_apply_bc() {
//         let mut model = Hydrus1D::new();
//         model.h[model.nz - 1] = -100.0; // Set different value
//         model.apply_bc();
//         assert_eq!(model.h[model.nz - 1], 0.0); // Should reset to h_bot = 0.0
//     }

//     #[test]
//     fn test_compute_flux() {
//         let model = Hydrus1D::new();
//         // Test at bottom boundary (z = 150 cm)
//         let bottom_flux = model.compute_flux(model.sim.z_max);
//         let h_n = model.h[model.nz - 1];
//         let h_nm1 = model.h[model.nz - 2];
//         let dz = model.sim.dz;
//         let k_n = model.soil.k(h_n, 150.0);
//         let dh_dz = (h_n - h_nm1) / dz;
//         let expected_flux = -k_n * (dh_dz + 1.0);
//         assert_relative_eq!(bottom_flux, expected_flux, epsilon = 1e-6);

//         // Test at root zone (z = 30 cm)
//         let root_zone_flux = model.compute_flux(model.sim.root_zone_depth);
//         let i_lower = (30.0 / dz).floor() as usize;
//         let i_upper = i_lower + 1;
//         let z_lower = i_lower as f64 * dz;
//         let z_upper = i_upper as f64 * dz;
//         let frac = (30.0 - z_lower) / (z_upper - z_lower);
//         let h_lower = model.h[i_lower];
//         let h_upper = model.h[i_upper];
//         let k_lower = model.soil.k(h_lower, z_lower);
//         let k_upper = model.soil.k(h_upper, z_upper);
//         let k = k_lower + frac * (k_upper - k_lower);
//         let dh_dz = (h_upper - h_lower) / dz;
//         let expected_root_flux = -k * (dh_dz + 1.0);
//         assert_relative_eq!(root_zone_flux, expected_root_flux, epsilon = 1e-6);
//     }

//     #[test]
//     fn test_step() {
//         let mut model = Hydrus1D::new();
//         let time = 0.0;
//         let (top_flux, rainfall, irrigation, evap, tp, _bottom_flux, _root_zone_flux) =
//             model.step(time);
//         // Check input values
//         assert_eq!(rainfall, 1.0); // From DailyInputs
//         assert_eq!(irrigation, 0.0);
//         assert_eq!(evap, 0.3);
//         assert_eq!(tp, 0.0);
//         // Check top flux (rainfall + irrigation since total_influx > 0)
//         assert_eq!(top_flux, 0.7);
//         // Check that pressure heads are updated
//         assert!(model.h.iter().all(|&h| h.is_finite()));
//         // Check bottom boundary condition
//         assert_relative_eq!(model.h[model.nz - 1], -100.0);
//     }

//     // #[test]
//     // fn test_run() {
//     //     let mut model = Hydrus1D::new();
//     //     // Redirect stdout to capture output
//     //     let stdout = std::io::stdout();
//     //     let mut output = Vec::new();
//     //     {
//     //         let mut buffer = std::io::BufWriter::new(&mut output);
//     //         std::io::set_output_capture(Some(&mut buffer));
//     //         model.run();
//     //         std::io::set_output_capture(None);
//     //     }
//     //     let output_str = String::from_utf8(output).expect("Invalid UTF-8");
//     //     // Check header
//     //     assert!(output_str.contains("Time [day], Depth [cm], Pressure Head [cm], Water Content [-], Root Uptake [1/day], Root Density [-], Top Flux [cm/day], Rainfall [cm/day], Irrigation [cm/day], Evaporation [cm/day], Transpiration [cm/day], Bottom Flux [cm/day], Root Zone Flux [cm/day]"));
//     //     // Check that output contains data rows
//     //     assert!(output_str.lines().count() > 1);
//     // }
// }
