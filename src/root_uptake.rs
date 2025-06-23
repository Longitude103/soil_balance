// Root water uptake parameters (Feddes model)
#[derive(Clone, Copy)]
pub struct RootParams {
    h1: f64, // Pressure head below which uptake starts [cm]
    h2: f64, // Optimal uptake range start [cm]
    h3: f64, // Optimal uptake range end [cm]
    h4: f64, // Pressure head below which uptake stops [cm]
}

impl RootParams {
    pub fn new() -> Self {
        RootParams {
            h1: -10.0, // Feddes parameters for generic crop
            h2: -25.0,
            h3: -200.0,
            h4: -8000.0,
        }
    }

    // Root water uptake sink term S(h, z, t) [1/day]
    pub fn s(&self, h: f64, tp: f64, root_density: f64) -> f64 {
        let alpha = if h > self.h1 || h < self.h4 {
            0.0
        } else if h <= self.h2 && h >= self.h3 {
            1.0
        } else if h > self.h2 && h <= self.h1 {
            (h - self.h1) / (self.h2 - self.h1)
        } else {
            (h - self.h4) / (self.h3 - self.h4)
        };
        alpha * tp * root_density // Root density scales the transpiration
    }
}
