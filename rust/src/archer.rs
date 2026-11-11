//! Archer reference overlays (NCRP 147 Table A.4) and radionuclide fits.

type Row = (f64, (f64, f64, f64));

pub const ARCHER_LEAD: [Row; 11] = [
    (50.0, (8.801, 27.28, 0.2957)),
    (60.0, (6.951, 24.89, 0.4198)),
    (70.0, (5.369, 23.49, 0.5881)),
    (80.0, (4.040, 21.69, 0.7187)),
    (90.0, (3.067, 18.83, 0.7726)),
    (100.0, (2.500, 15.28, 0.7557)),
    (110.0, (2.296, 11.70, 0.6827)),
    (120.0, (2.246, 8.950, 0.5873)),
    (125.0, (2.219, 7.923, 0.5386)),
    (140.0, (2.009, 5.916, 0.4018)),
    (150.0, (1.757, 5.177, 0.3156)),
];
pub const ARCHER_CONCRETE: [Row; 11] = [
    (50.0, (0.09, 0.17, 0.232)),
    (60.0, (0.063, 0.169, 0.273)),
    (70.0, (0.051, 0.169, 0.385)),
    (80.0, (0.046, 0.154, 0.493)),
    (90.0, (0.042, 0.113, 0.469)),
    (100.0, (0.039, 0.085, 0.427)),
    (110.0, (0.037, 0.074, 0.475)),
    (120.0, (0.036, 0.071, 0.607)),
    (125.0, (0.035, 0.071, 0.697)),
    (140.0, (0.03345, 0.07476, 1.047)),
    (150.0, (0.03243, 0.08599, 1.467)),
];
pub const ARCHER_GYPSUM: [Row; 11] = [
    (50.0, (0.038, 0.087, 0.51)),
    (60.0, (0.030, 0.080, 0.617)),
    (70.0, (0.023, 0.072, 0.729)),
    (80.0, (0.019, 0.061, 0.810)),
    (90.0, (0.016, 0.050, 0.859)),
    (100.0, (0.015, 0.042, 0.894)),
    (110.0, (0.013, 0.035, 0.924)),
    (120.0, (0.012, 0.030, 0.957)),
    (125.0, (0.012, 0.029, 0.968)),
    (140.0, (0.01088, 0.02436, 0.9964)),
    (150.0, (0.0103, 0.02198, 1.013)),
];

pub fn table(material: &str) -> Option<&'static [Row]> {
    match material {
        "lead" => Some(&ARCHER_LEAD),
        "concrete" => Some(&ARCHER_CONCRETE),
        "gypsum" => Some(&ARCHER_GYPSUM),
        _ => None,
    }
}

pub fn archer_transmission(x_mm: f64, al: f64, be: f64, ga: f64) -> f64 {
    ((1.0 + be / al) * (al * ga * x_mm).exp() - be / al).powf(-1.0 / ga)
}

/// Archer transmission for overlay; None if the material has no table.
pub fn archer_reference(x_mm: f64, kvp: f64, material: &str) -> Option<f64> {
    let tbl = table(material)?;
    if x_mm <= 0.0 {
        return Some(1.0);
    }
    let kv = kvp.max(tbl[0].0).min(tbl[tbl.len() - 1].0);
    let mut a = tbl[tbl.len() - 1].1;
    for i in 0..tbl.len() - 1 {
        let (lo, hi) = (tbl[i].0, tbl[i + 1].0);
        if lo <= kv && kv <= hi {
            let f = (kv - lo) / (hi - lo);
            let (p, q) = (tbl[i].1, tbl[i + 1].1);
            a = (p.0 + f * (q.0 - p.0), p.1 + f * (q.1 - p.1), p.2 + f * (q.2 - p.2));
            break;
        }
    }
    Some(archer_transmission(x_mm, a.0, a.1, a.2))
}

/// Archer alpha at the tabulated kVp nearest the (clamped) beam kVp.
pub fn archer_alpha_near(kvp: f64, material: &str) -> Option<f64> {
    let tbl = table(material)?;
    let kv = kvp.max(tbl[0].0).min(tbl[tbl.len() - 1].0);
    let mut best = tbl[0];
    for r in tbl {
        if (r.0 - kv).abs() < (best.0 - kv).abs() {
            best = *r;
        }
    }
    Some(best.1 .0)
}

/// (name, energies, lead, concrete) with (alpha, beta, gamma, unit, basis)
pub type Fit = Option<(f64, f64, f64, &'static str, &'static str)>;
pub const RADIONUCLIDES: [(&str, &str, Fit, Fit); 10] = [
    ("F-18", "511 keV annihilation", Some((1.543, -0.4408, 2.136, "cm", "verified")), Some((0.1539, -0.1161, 2.0752, "cm", "verified"))),
    ("Ga-67", "93 / 185 / 300 / 393 keV", Some((0.04412, 0.7326, 0.3869, "mm", "physics")), Some((0.1865, -0.1529, 1.683, "cm", "physics"))),
    ("Ga-68", "511 keV annihilation", Some((1.41, -0.402, 4.47, "cm", "physics")), Some((0.134, -0.094, 3.2, "cm", "physics"))),
    ("Ge-68/Ga-68 gen", "511 keV annihilation", Some((1.41, -0.402, 4.47, "cm", "physics")), Some((0.134, -0.094, 3.2, "cm", "physics"))),
    ("I-123", "159 keV + 529 keV at 1.4 %", Some((0.1322, 14.85, 1.488, "mm", "inferred")), Some((0.1866, 1253.0, 74.93, "cm", "inferred"))),
    ("I-131", "364 / 637 / 723 keV", Some((-0.009, 2.439, 0.117, "cm", "verified")), Some((0.152, -0.099, 2.634, "cm", "physics"))),
    ("In-111", "171 / 245 keV", Some((0.7094, 30.63, 5.979, "mm", "inferred")), Some((0.2339, -0.1265, 0.4618, "cm", "inferred"))),
    ("Lu-177", "113 / 208 keV", Some((0.3855, 1.071, 0.2822, "mm", "verified")), Some((0.2477, -0.1173, 0.4404, "cm", "verified"))),
    ("Tc-99m", "140 keV", Some((2.558, 1.010, 4.344, "mm", "verified")), Some((0.3102, -0.1729, 0.3622, "cm", "verified"))),
    ("Tl-201", "69-81 keV x-rays + 167 keV", Some((1.67, 1.358, 0.3692, "mm", "inferred")), None),
];

fn nuclide_raw(name: &str, material: &str) -> Fit {
    let e = RADIONUCLIDES.iter().find(|r| r.0 == name)?;
    match material {
        "lead" => e.2,
        "concrete" => e.3,
        _ => None,
    }
}

/// (alpha, beta, gamma) in mm^-1, or None if unusable/unavailable.
pub fn nuclide_params(name: &str, material: &str) -> Option<(f64, f64, f64)> {
    let (mut a, mut b, g, unit, _basis) = nuclide_raw(name, material)?;
    if unit == "cm" {
        a /= 10.0;
        b /= 10.0;
    }
    if a == 0.0 || (1.0 + b / a).abs() < 1e-12 {
        return None;
    }
    Some((a, b, g))
}

/// (basis, unit) for display.
pub fn nuclide_info(name: &str, material: &str) -> Option<(&'static str, &'static str)> {
    nuclide_raw(name, material).map(|v| (v.4, v.3))
}

pub fn archer_reference_nuclide(x_mm: f64, name: &str, material: &str) -> Option<f64> {
    let (al, be, ga) = nuclide_params(name, material)?;
    if x_mm <= 0.0 {
        return Some(1.0);
    }
    Some(archer_transmission(x_mm, al, be, ga))
}
