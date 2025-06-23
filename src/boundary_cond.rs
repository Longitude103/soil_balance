use crate::daily_inputs::DailyInputs;

// Boundary condition parameters
#[derive(Clone)]
pub struct BoundaryParams {
    pub h_crit: f64, // Critical pressure head for evaporation [cm]
    // h_top: f64,  // Fallback top pressure head [cm]
    pub h_bot: f64, // Bottom pressure head [cm]
}

impl BoundaryParams {
    pub fn new() -> Self {
        BoundaryParams {
            h_crit: -15000.0, // Critical pressure head for evaporation
            // h_top: -10.0,     // Fallback top pressure head
            h_bot: -100.0, // Constant bottom pressure head
        }
    }

    // Top boundary flux [cm/day] (positive = rainfall + irrigation, negative = evaporation)
    pub fn top_flux(&self, time: f64, h_top: f64, inputs: &DailyInputs) -> (f64, bool) {
        let rainfall = inputs.get_daily_value(time, &inputs.rainfall);
        let irrigation = inputs.get_daily_value(time, &inputs.irrigation);
        let total_influx = rainfall + irrigation;
        if total_influx > 0.0 {
            (total_influx, true) // Rainfall and/or irrigation as flux
        } else {
            let evap = inputs.get_daily_value(time, &inputs.evaporation);
            if h_top > self.h_crit {
                (-evap, true) // Evaporation as flux
            } else {
                (0.0, false) // Switch to Dirichlet (h = h_crit)
            }
        }
    }
}
