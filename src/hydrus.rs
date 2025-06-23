use crate::boundary_cond::BoundaryParams;
use crate::daily_inputs::DailyInputs;
use crate::root_uptake::RootParams;
use crate::soil_hydra::SoilParams;
use nalgebra::{DMatrix, DVector};

// Simulation parameters
struct SimParams {
    dz: f64,    // Spatial step [cm]
    dt: f64,    // Time step [day]
    z_max: f64, // Soil column depth [cm]
    t_max: f64, // Simulation time [day]
}

impl SimParams {
    fn new() -> Self {
        SimParams {
            dz: 1.0,      // 1 cm spatial resolution
            dt: 0.001,    // 0.001 day time step
            z_max: 100.0, // 100 cm soil column
            t_max: 1.0,   // 1 day simulation
        }
    }
}

// Finite difference solver for Richards equation with root water uptake and variable BCs
pub struct Hydrus1D {
    soil: SoilParams,
    root: RootParams,
    bc: BoundaryParams,
    sim: SimParams,
    inputs: DailyInputs,
    h: DVector<f64>, // Pressure head at each node [cm]
    nz: usize,       // Number of spatial nodes
}

impl Hydrus1D {
    pub fn new() -> Self {
        let sim = SimParams::new();
        let nz = (sim.z_max / sim.dz).ceil() as usize + 1;
        let h_init = -50.0; // Initial pressure head [cm]
        let h = DVector::from_element(nz, h_init);
        Hydrus1D {
            soil: SoilParams::new(),
            root: RootParams::new(),
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

    // Solve one time step using implicit finite difference with nalgebra
    fn step(&mut self, time: f64) -> (f64, f64, f64, f64, f64) {
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
            if i == 0 {
                if is_flux {
                    // Neumann BC: q = -K(h) * (dh/dz + 1) = top_flux
                    let k_i = self.soil.k(h_old[i]);
                    let k_ip = 0.5 * (k_i + self.soil.k(h_old[i + 1]));
                    a[(i, i)] = -k_ip / dz;
                    a[(i, i + 1)] = k_ip / dz;
                    b[i] = k_ip - top_flux;
                } else {
                    // Dirichlet BC: h = h_crit
                    a[(i, i)] = 1.0;
                    b[i] = self.bc.h_crit;
                }
            } else if i == nz - 1 {
                // Dirichlet BC at bottom
                a[(i, i)] = 1.0;
                b[i] = self.bc.h_bot;
            } else {
                let h_i = h_old[i];
                let h_ip1 = h_old[i + 1];
                let h_im1 = h_old[i - 1];

                // Hydraulic conductivity at interfaces
                let k_ip = 0.5 * (self.soil.k(h_i) + self.soil.k(h_ip1));
                let k_im = 0.5 * (self.soil.k(h_i) + self.soil.k(h_im1));

                // Specific moisture capacity
                let c_i = self.soil.c(h_i);

                // Root water uptake
                let z = i as f64 * dz;
                let root_density = self.inputs.get_root_density(time, z);
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

        (top_flux, rainfall, irrigation, evap, tp)
    }

    // Run simulation and output results
    pub fn run(&mut self) {
        let nt = (self.sim.t_max / self.sim.dt).ceil() as usize;
        println!(
            "Time [day], Depth [cm], Pressure Head [cm], Water Content [-], Root Uptake [1/day], Root Density [-], Top Flux [cm/day], Rainfall [cm/day], Irrigation [cm/day], Evaporation [cm/day], Transpiration [cm/day]"
        );
        for t in 0..=nt {
            let time = t as f64 * self.sim.dt;
            let (top_flux, rainfall, irrigation, evap, tp) = if t < nt {
                self.step(time)
            } else {
                let flux = self.bc.top_flux(time, self.h[0], &self.inputs).0;
                (
                    flux,
                    self.inputs.get_daily_value(time, &self.inputs.rainfall),
                    self.inputs.get_daily_value(time, &self.inputs.irrigation),
                    self.inputs.get_daily_value(time, &self.inputs.evaporation),
                    self.inputs
                        .get_daily_value(time, &self.inputs.transpiration),
                )
            };
            if t % (nt / 10) == 0 || t == nt {
                for i in 0..self.nz {
                    let z = i as f64 * self.sim.dz;
                    let theta = self.soil.theta(self.h[i]);
                    let root_density = self.inputs.get_root_density(time, z);
                    let s = self.root.s(self.h[i], tp, root_density);
                    let flux = if i == 0 { top_flux } else { 0.0 };
                    println!(
                        "{:.3}, {:.1}, {:.2}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}, {:.3}",
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
                        tp
                    );
                }
            }
        }
    }
}
