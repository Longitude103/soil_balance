use std::fs;

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
    pub fn new(theta_r: f64, theta_s: f64, alpha: f64, n: f64, ks: f64) -> Self {
        SoilParams {
            theta_r,
            theta_s,
            alpha,
            n,
            ks,
        }
    }

    pub(crate) fn from_toml(soil_name: &str) -> Self {
        let toml_str = fs::read_to_string("soil_parameters.toml")
            .expect("Failed to read soil_parameters.toml");
        let params: toml::Value = toml::from_str(&toml_str).expect("Failed to parse TOML");

        let mut soil_name = soil_name.to_lowercase();
        if soil_name.is_empty() {
            soil_name = "loam".to_string();
        }

        let mut soil_params = params.get(soil_name).and_then(|v| v.as_table());

        if soil_params.is_none() {
            soil_params = params.get("loam").and_then(|v| v.as_table());
        }

        let soil_params = soil_params.expect("Invalid soils");

        let theta_r = soil_params
            .get("theta_r")
            .and_then(|v| v.as_float())
            .expect("Missing theta_r");
        let theta_s = soil_params
            .get("theta_s")
            .and_then(|v| v.as_float())
            .expect("Missing theta_s");
        let alpha = soil_params
            .get("alpha")
            .and_then(|v| v.as_float())
            .expect("Missing alpha");
        let n = soil_params
            .get("n")
            .and_then(|v| v.as_float())
            .expect("Missing n");
        let ks = soil_params
            .get("ks")
            .and_then(|v| v.as_float())
            .expect("Missing ks");

        SoilParams {
            theta_r,
            theta_s,
            alpha,
            n,
            ks,
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

// Soil layer definition
#[derive(Clone)]
pub struct SoilLayer {
    pub(crate) z_min: f64,    // Start depth of layer [cm]
    pub(crate) z_max: f64,    // End depth of layer [cm]
    pub params: SoilParams,   // Soil parameters
    pub horizon_name: String, // Horizon name like "H1", "A", or "B"
}

impl SoilLayer {
    // This function creates a new soil layer from a vector of soil parameters
    pub fn new(z_min: f64, z_max: f64, params: SoilParams, horizon_name: String) -> Self {
        SoilLayer {
            z_min,
            z_max,
            params,
            horizon_name,
        }
    }
}

// Soil profile managing multiple layers
#[derive(Clone)]
pub struct SoilProfile {
    pub layers: Vec<SoilLayer>,
    pub hydro_group: String, // Hydrological group like "A", "B", "C"
    pub soil_name: String,   // name of soils like "Dailey loamy sand, 3 to 9 percent slopes"
}

impl SoilProfile {
    // This function creates a new soil profile from a vector of soil layers
    pub fn new() -> Self {
        SoilProfile {
            layers: Vec::new(),
            hydro_group: String::new(),
            soil_name: String::new(),
        }
    }

    // this function creates a new soil profile from a vector of soil layers that are derived from the default TOML file of soil_parameters
    pub fn new_from_toml(
        soil_layers: Vec<(f64, f64, &str, &str)>,
        hydro_group: String,
        soil_name: String,
    ) -> Self {
        let layers = soil_layers
            .into_iter()
            .map(|(z_min, z_max, soil_name, horizon_name)| SoilLayer {
                z_min,
                z_max,
                params: SoilParams::from_toml(soil_name),
                horizon_name: horizon_name.to_string(),
            })
            .collect();
        SoilProfile {
            layers,
            hydro_group,
            soil_name,
        }
    }

    // Get soil parameters at a given depth
    fn params_at_depth(&self, z: f64) -> &SoilParams {
        for layer in &self.layers {
            if z >= layer.z_min && z <= layer.z_max {
                return &layer.params;
            }
        }
        // Fallback to last layer if z exceeds all layers
        &self.layers.last().expect("No soil layers defined").params
    }

    // Hydraulic functions
    pub(crate) fn theta(&self, h: f64, z: f64) -> f64 {
        self.params_at_depth(z).theta(h)
    }

    pub(crate) fn k(&self, h: f64, z: f64) -> f64 {
        self.params_at_depth(z).k(h)
    }

    pub(crate) fn c(&self, h: f64, z: f64) -> f64 {
        self.params_at_depth(z).c(h)
    }
}
