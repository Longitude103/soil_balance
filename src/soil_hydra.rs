// Soil hydraulic parameters (van Genuchten-Mualem model)
#[derive(Clone, Copy)]
pub struct SoilParams {
    theta_r: f64, // Residual water content [cm³/cm³]
    theta_s: f64, // Saturated water content [cm³/cm³]
    alpha: f64,   // van Genuchten parameter [1/cm]
    n: f64,       // van Genuchten parameter [-]
    ks: f64,      // Saturated hydraulic conductivity [cm/day]
}

impl SoilParams {
    pub fn new() -> Self {
        SoilParams {
            theta_r: 0.045, // Example values for sandy loam
            theta_s: 0.43,
            alpha: 0.145,
            n: 1.89,
            ks: 104.8,
        }
    }

    // Water content as a function of pressure head (h, cm)
    pub fn theta(&self, h: f64) -> f64 {
        if h >= 0.0 {
            self.theta_s
        } else {
            let m = 1.0 - 1.0 / self.n;
            let ah = self.alpha * h.abs();
            let denom = (1.0 + ah.powf(self.n)).powf(m);
            self.theta_r + (self.theta_s - self.theta_r) / denom
        }
    }

    // Hydraulic conductivity as a function of pressure head (h, cm)
    pub fn k(&self, h: f64) -> f64 {
        if h >= 0.0 {
            self.ks
        } else {
            let m = 1.0 - 1.0 / self.n;
            let ah = self.alpha * h.abs();
            let se = (1.0 / (1.0 + ah.powf(self.n))).powf(m); // Effective saturation
            self.ks * se.sqrt() * (1.0 - (1.0 - se.powf(1.0 / m)).powf(m)).powi(2)
        }
    }

    // Specific moisture capacity, C(h) = d(theta)/dh [1/cm]
    pub fn c(&self, h: f64) -> f64 {
        if h >= 0.0 {
            0.0
        } else {
            let m = 1.0 - 1.0 / self.n;
            let ah = self.alpha * h.abs();
            let denom = (1.0 + ah.powf(self.n)).powf(m + 1.0);
            self.alpha * (self.theta_s - self.theta_r) * m * self.n * ah.powf(self.n - 1.0) / denom
        }
    }
}
