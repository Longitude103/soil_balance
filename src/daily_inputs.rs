// Daily input data for rainfall, irrigation, evaporation, transpiration, and root distribution
#[derive(Clone)]
pub struct DailyInputs {
    pub rainfall: Vec<f64>,              // Daily rainfall rates [cm/day]
    pub irrigation: Vec<f64>,            // Daily irrigation rates [cm/day]
    pub evaporation: Vec<f64>,           // Daily potential evaporation rates [cm/day]
    pub transpiration: Vec<f64>,         // Daily potential transpiration rates [cm/day]
    pub root_dist: Vec<Vec<(f64, f64)>>, // Daily root distribution profiles [(depth [cm], density [-])]
}

impl DailyInputs {
    pub fn new() -> Self {
        // Example: 2 days of data
        DailyInputs {
            // create a vector with 30 values of 1
            rainfall: vec![1.0; 30], // cm/day
            irrigation: vec![0.0; 30], // cm/day
            evaporation: vec![0.3; 30], // cm/day
            transpiration: vec![0.0; 30], // cm/day

            // (cm to depth, % root density) tuples for each day
            root_dist: vec![vec![(10.0, 0.00), (20.0, 0.000), (30.0, 0.0)]; 30],

            // rainfall: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0],      // 2 cm/day on day 0, 0 cm/day on day 1
            // irrigation: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],    // 1 cm/day on day 0, 0 cm/day on day 1
            // evaporation: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.5],   // 0.3 cm/day on day 0, 0.5 cm/day on day 1
            // transpiration: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.4], // 0.5 cm/day on day 0, 0.4 cm/day on day 1
            // root_dist: vec![
            //     // Day 0: No Roots (70% in top 10 cm)
            //     vec![(0.0, 0.00), (10.0, 0.00), (20.0, 0.00), (30.0, 0.0)],
            //     // Day 1: No Roots (70% in top 10 cm)
            //     vec![(0.0, 0.00), (10.0, 0.00), (20.0, 0.00), (30.0, 0.0)],
            //     // Day 2: No Roots (70% in top 10 cm)
            //     vec![(0.0, 0.00), (10.0, 0.00), (20.0, 0.00), (30.0, 0.0)],
            //     // Day 3: No Roots (70% in top 10 cm)
            //     vec![(0.0, 0.00), (10.0, 0.00), (20.0, 0.00), (30.0, 0.0)],
            //     // Day 4: No Roots (70% in top 10 cm)
            //     vec![(0.0, 0.00), (10.0, 0.00), (20.0, 0.00), (30.0, 0.0)],
            //     // Day 5: No Roots (70% in top 10 cm)
            //     vec![(0.0, 0.00), (10.0, 0.00), (20.0, 0.00), (30.0, 0.0)],
            //     // Day 6: No Roots (70% in top 10 cm)
            //     vec![(0.0, 0.00), (10.0, 0.00), (20.0, 0.00), (30.0, 0.0)],
            //     // Day 7: No Roots (70% in top 10 cm)
            //     vec![(0.0, 0.00), (10.0, 0.00), (20.0, 0.00), (30.0, 0.0)],
            //     // Day 8: Surface-heavy distribution (70% in top 10 cm)
            //     vec![(0.0, 0.07), (10.0, 0.02), (20.0, 0.01), (30.0, 0.0)],
            //     // Day 9: More uniform distribution
            //     vec![(0.0, 0.04), (10.0, 0.03), (20.0, 0.02), (30.0, 0.0)],
            // ],
        }
    }

    // Get value for the current day (floor of time in days)
    pub(crate) fn get_daily_value(&self, time: f64, values: &Vec<f64>) -> f64 {
        let day = time.floor() as usize;
        if day < values.len() {
            values[day]
        } else {
            *values.last().unwrap_or(&0.0) // Use last value if time exceeds input length
        }
    }

    // Get root density for the current day and interpolate at depth z
    pub(crate) fn get_root_density(&self, time: f64, z: f64) -> f64 {
        let day = time.floor() as usize;
        let r_dist = vec![(0.0, 0.0)];
        let profile = if day < self.root_dist.len() {
            &self.root_dist[day]
        } else {
            self.root_dist.last().unwrap_or(&r_dist)
        };

        // Linear interpolation between depth points
        for i in 0..profile.len() - 1 {
            let (z1, d1) = profile[i];
            let (z2, d2) = profile[i + 1];
            if z >= z1 && z <= z2 {
                return d1 + (d2 - d1) * (z - z1) / (z2 - z1);
            }
        }
        profile.last().map(|&(_, d)| d).unwrap_or(0.0)
    }
}
