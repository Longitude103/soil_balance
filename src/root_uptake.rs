use std::fs;

// Root water uptake parameters (Feddes model)
#[derive(Clone, Copy)]
pub struct RootParams {
    h1: f64, // Pressure head below which uptake starts [cm]
    h2: f64, // Optimal uptake range start [cm]
    h3: f64, // Optimal uptake range end [cm]
    h4: f64, // Pressure head below which uptake stops [cm]
}

impl RootParams {
    pub(crate) fn new(crop: &str) -> Self {
        let toml_str = fs::read_to_string("feddes_parameters.toml").expect("Failed to read TOML file");
        let params: toml::Value = toml::from_str(&toml_str).expect("Failed to parse TOML");

        let mut crop_name = crop.to_lowercase();
        if crop_name.is_empty() {
            crop_name = "grass".to_string();
        }

        let mut crop_params = params
            .get(crop_name)
            .and_then(|v| v.as_table());

        if crop_params.is_none() {
            crop_params = params.get("grass").and_then(|v| v.as_table());
        }

        let crop_params = crop_params.expect("Missing crop parameters");

        let h1 = crop_params.get("h1").and_then(|v| v.as_float()).expect("Missing h1");
        let h2 = crop_params.get("h2").and_then(|v| v.as_float()).expect("Missing h2");
        let h3 = crop_params.get("h3").and_then(|v| v.as_float()).expect("Missing h3");
        let h4 = crop_params.get("h4").and_then(|v| v.as_float()).expect("Missing h4");

        RootParams { h1, h2, h3, h4 }
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
