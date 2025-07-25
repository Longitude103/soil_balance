/**
Module for calculating NRCS runoff curve number and daily runoff.

This module provides functions to estimate the NRCS (formerly SCS) Curve Number (CN)
based on soil group and land use, and to calculate runoff depth using the CN method.
Units are in millimeters for metric consistency with hydrological models.

Calculates the NRCS Curve Number (CN) for average antecedent moisture condition (AMC II).

This is a simplified lookup based on common agricultural land uses and hydrologic soil groups.
For more comprehensive tables, consider expanding this function or using external data.

# Arguments
* `soil_group` - Hydrologic soil group ("A", "B", "C", "D").
* `land_use` - Land use description (e.g., "row crops good", "pasture fair").

# Returns
The curve number as f64.

# Panics
If the soil group or land use is not recognized.
*/
pub fn calculate_curve_number(soil_group: &str, land_use: &str) -> Result<f64, String> {
    match (
        soil_group.to_uppercase().as_str(),
        land_use.to_lowercase().as_str(),
    ) {
        // Row crops, straight row, good condition
        ("A", "row crops good") => Ok(67.0),
        ("B", "row crops good") => Ok(78.0),
        ("C", "row crops good") => Ok(85.0),
        ("D", "row crops good") => Ok(89.0),

        // Row crops, straight row, fair condition
        ("A", "row crops fair") => Ok(72.0),
        ("B", "row crops fair") => Ok(83.0),
        ("C", "row crops fair") => Ok(90.0),
        ("D", "row crops fair") => Ok(94.0),

        // Row crops, straight row, poor condition
        ("A", "row crops poor") => Ok(77.0),
        ("B", "row crops poor") => Ok(88.0),
        ("C", "row crops poor") => Ok(95.0),
        ("D", "row crops poor") => Ok(99.0),

        // Row crops, straight row, excellent condition
        ("A", "row crops excellent") => Ok(97.0),
        ("B", "row crops excellent") => Ok(108.0),
        ("C", "row crops excellent") => Ok(115.0),
        ("D", "row crops excellent") => Ok(119.0),

        // Pasture, fair condition
        ("A", "pasture fair") => Ok(49.0),
        ("B", "pasture fair") => Ok(69.0),
        ("C", "pasture fair") => Ok(79.0),
        ("D", "pasture fair") => Ok(84.0),

        // Pasture, poor condition
        ("A", "pasture poor") => Ok(34.0),
        ("B", "pasture poor") => Ok(54.0),
        ("C", "pasture poor") => Ok(64.0),
        ("D", "pasture poor") => Ok(69.0),

        // Pasture, excellent condition
        ("A", "pasture excellent") => Ok(29.0),
        ("B", "pasture excellent") => Ok(49.0),
        ("C", "pasture excellent") => Ok(59.0),
        ("D", "pasture excellent") => Ok(64.0),

        // Pasture, good condition
        ("A", "pasture good") => Ok(39.0),
        ("B", "pasture good") => Ok(59.0),
        ("C", "pasture good") => Ok(69.0),
        ("D", "pasture good") => Ok(74.0),

        // Add more land uses as needed
        _ => Err(format!(
            "Unknown combination of soil group '{}' and land use '{}'",
            soil_group, land_use
        )),
    }
}

/**
Calculates the daily runoff depth using the NRCS Curve Number method.

This uses the metric version of the formula, with depths in millimeters.

# Arguments
- `water_application` - Total water input (e.g., precipitation + irrigation) in mm for the day.
- `cn` - Curve number (from `calculate_curve_number` or provided).

# Returns
The estimated runoff depth in mm.
*/
pub fn calculate_runoff(water_application: f64, cn: f64) -> f64 {
    if water_application <= 0.0 {
        return 0.0;
    }

    // Maximum potential retention (S) in mm
    let s = (25400.0 / cn) - 254.0;

    // Initial abstraction (Ia = 0.2 * S)
    let ia = 0.2 * s;

    if water_application <= ia {
        return 0.0;
    }

    // Runoff depth Q = (P - Ia)^2 / (P - Ia + S)
    (water_application - ia).powi(2) / (water_application - ia + s)
}
