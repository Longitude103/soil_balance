use crate::boundary_cond::BoundaryParams;
use crate::daily_inputs::DailyInputs;
use crate::root_uptake::RootParams;
use crate::soil_hydra::SoilProfile;
use nalgebra::{DMatrix, DVector};

// Simulation parameters
struct SimParams {
    dz: f64,    // Spatial step [cm]
    dt: f64,    // Time step [day]
    z_max: f64, // Soil column depth [cm]
    t_max: f64, // Simulation time [day]
    root_zone_depth: f64, // Root zone depth [cm]
}

impl SimParams {
    fn new() -> Self {
        SimParams {
            dz: 10.0,      // 1 cm spatial resolution
            dt: 0.1,    // 0.001 day time step
            z_max: 150.0,  // 150 cm soil column to include groundwater
            t_max: 30.0,   // 1 day simulation
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

impl Default for Hydrus1D {
    fn default() -> Self {
        Self::new()
    }
}

impl Hydrus1D {
    pub fn new() -> Self {
        let sim = SimParams::new();
        let nz = (sim.z_max / sim.dz).ceil() as usize + 1;
        // let h_init = -50.0; // Initial pressure head [cm]

        let mut h = DVector::from_fn(nz, |i, _| {
            let z = i as f64 * sim.dz;
            if z >= 150.0 { 0.0 } else if z <= 60.0 { -33.0 } else { -15.0 }
        });

        // Set h = 0 at z = 150 cm (bottom node)
        h[nz - 1] = 0.0;

        // Define soil layers: loam (0–60 cm), sand (60–150 cm)
        let soil_layers = vec![
            (0.0, 60.0, "loam"),
            (60.0, 80.0, "loamy-sand"),
            (80.0, 150.0, "sand"),
        ];
        let soil_profile = SoilProfile::new(soil_layers);

        Hydrus1D {
            soil: soil_profile,
            root: RootParams::new(""),
            bc: BoundaryParams::new(),
            inputs: DailyInputs::new(),
            sim,
            h,
            nz,
        }
    }

    // Initialize boundary conditions
    pub fn apply_bc(&mut self) {
        let h = self.h.as_mut_slice();
        h[self.nz - 1] = self.bc.h_bot;
    }

    // Compute flux at a given depth [cm/day]
    fn compute_flux(&self, z: f64) -> f64 {
        let dz = self.sim.dz;
        // Find nodes bracketing the depth
        let i_lower = (z / dz).floor() as usize;
        let i_upper = i_lower + 1;
        if i_lower >= self.nz - 1 {
            // At bottom boundary, use backward difference
            let h_n = self.h[self.nz - 1];
            let h_nm1 = self.h[self.nz - 2];
            let k_n = self.soil.k(h_n, z);
            let dh_dz = (h_n - h_nm1) / dz;
            -k_n * (dh_dz + 1.0)
        } else {
            // Interpolate pressure head and conductivity
            let z_lower = i_lower as f64 * dz;
            let z_upper = i_upper as f64 * dz;
            let frac = (z - z_lower) / (z_upper - z_lower);
            let h_lower = self.h[i_lower];
            let h_upper = self.h[i_upper];
            let k_lower = self.soil.k(h_lower, z_lower);
            let k_upper = self.soil.k(h_upper, z_upper);
            // let h = h_lower + frac * (h_upper - h_lower);
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
        let tp = self.inputs.get_daily_value(time, &self.inputs.transpiration);

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
                } else {
                    // Dirichlet BC: h = h_crit
                    a[(i, i)] = 1.0;
                    b[i] = self.bc.h_crit;
                }
            } else if i == nz - 1 {
                // Dirichlet BC at bottom (groundwater)
                a[(i, i)] = 1.0;
                b[i] = self.bc.h_bot; // h = 0
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

        (top_flux, rainfall, irrigation, evap, tp, bottom_flux, root_zone_flux)
    }

    // Run simulation and output results
    pub(crate) fn run(&mut self) {
        let nt = (self.sim.t_max / self.sim.dt).ceil() as usize;
        println!("Time [day], Depth [cm], Pressure Head [cm], Water Content [-], Root Uptake [1/day], Root Density [-], Top Flux [cm/day], Rainfall [cm/day], Irrigation [cm/day], Evaporation [cm/day], Transpiration [cm/day], Bottom Flux [cm/day], Root Zone Flux [cm/day]");
        for t in 0..=nt {
            let time = t as f64 * self.sim.dt;
            let (top_flux, rainfall, irrigation, evap, tp, bottom_flux, root_zone_flux) = if t < nt {
                self.step(time)
            } else {
                let flux = self.bc.top_flux(time, self.h[0], &self.inputs).0;
                let bottom_flux = self.compute_flux(self.sim.z_max);
                let root_zone_flux = self.compute_flux(self.sim.root_zone_depth);
                (flux,
                 self.inputs.get_daily_value(time, &self.inputs.rainfall),
                 self.inputs.get_daily_value(time, &self.inputs.irrigation),
                 self.inputs.get_daily_value(time, &self.inputs.evaporation),
                 self.inputs.get_daily_value(time, &self.inputs.transpiration),
                 bottom_flux,
                 root_zone_flux)
            };
            if t % (nt / 10) == 0 || t == nt {
                for i in 0..self.nz {
                    let z = i as f64 * self.sim.dz;
                    let theta = self.soil.theta(self.h[i], z);
                    let root_density = self.inputs.get_root_density(time, z);
                    let s = self.root.s(self.h[i], tp, root_density);
                    let flux = if i == 0 { top_flux } else { 0.0 };
                    let b_flux = if i == self.nz - 1 { bottom_flux } else { 0.0 };
                    let rz_flux = if (z - self.sim.root_zone_depth).abs() < 1e-6 { root_zone_flux } else { 0.0 };
                    println!("{:.3}, {:.1}, {:.2}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}",
                             time, z, self.h[i], theta, s, root_density, flux, rainfall, irrigation, evap, tp, b_flux, rz_flux);
                }
            }
        }
    }
}
