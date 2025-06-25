use crate::daily_inputs::DailyInputs;


// Bottom boundary condition types
#[derive(Clone)]
pub(crate) enum BottomBoundary {
    Groundwater(f64), // Dirichlet: fixed pressure head (cm)
    FreeDrainage,     // Neumann: q = -K(h), dh/dz = 0
}

// Boundary condition parameters
#[derive(Clone)]
pub struct BoundaryParams {
    pub(crate) h_crit: f64,          // Critical pressure head for evaporation [cm]
    h_top: f64,           // Fallback top pressure head [cm]
    pub(crate) bottom_boundary: BottomBoundary, // Bottom boundary condition
}

impl BoundaryParams {
    pub fn new() -> Self {
        BoundaryParams {
            h_crit: -15000.0, // Critical pressure head for evaporation
            h_top: -10.0,     // Fallback top pressure head
            bottom_boundary: BottomBoundary::FreeDrainage, // Default to free drainage
            // For groundwater, use: BottomBoundary::Groundwater(0.0)
        }
    }

    // Top boundary flux [cm/day] (positive = rainfall + irrigation, negative = evaporation)
    pub(crate) fn top_flux(&self, time: f64, h_top: f64, inputs: &DailyInputs) -> (f64, bool) {
        let rainfall = inputs.get_daily_value(time, &inputs.rainfall);
        let irrigation = inputs.get_daily_value(time, &inputs.irrigation);
        let evap = inputs.get_daily_value(time, &inputs.evaporation);
        let total_influx = rainfall + irrigation - evap;
        if total_influx > 0.0 {
            (total_influx, true) // Rainfall and/or irrigation as flux
        } else if h_top > self.h_crit {
            (-evap, true) // Evaporation as flux
        } else {
            (0.0, false) // Switch to Dirichlet (h = h_crit)
        }
    }
}
