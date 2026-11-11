//! SECTION 2 - spectrum generation, filtration, beam quality, RQA and
//! radioisotope sources (port of `Spectrum`, `spectrum_kramers`,
//! `spectrum_spekpy`, `solve_filtration_from_hvl`, `rqa_*`, isotopes).

use crate::num::{argmin, gradient, pairwise_sum, trapz};
use crate::xs::{air_muen, air_muen1, CrossSections, Kind};

#[derive(Clone, Debug)]
#[allow(non_snake_case)]
pub struct Spectrum {
    pub e: Vec<f64>,
    pub n: Vec<f64>,
    pub kvp: f64,
    pub model: String,
    /// (material, mm, density)
    pub filters: Vec<(String, f64, f64)>,
    pub absolute: bool,
    pub kerma_uGy_per_mAs: Option<f64>,
}

#[allow(non_snake_case)]
impl Spectrum {
    pub fn new(e: Vec<f64>, n: Vec<f64>, kvp: f64, model: &str) -> Spectrum {
        Spectrum { e, n, kvp, model: model.to_string(), filters: vec![], absolute: false, kerma_uGy_per_mAs: None }
    }

    pub fn fluence(&self) -> f64 {
        trapz(&self.n, &self.e)
    }

    pub fn energy_fluence(&self) -> f64 {
        let ne: Vec<f64> = self.n.iter().zip(&self.e).map(|(n, e)| n * e).collect();
        trapz(&ne, &self.e)
    }

    pub fn mean_energy(&self) -> f64 {
        let f = self.fluence();
        if f > 0.0 {
            self.energy_fluence() / f
        } else {
            0.0
        }
    }

    /// Beer-Lambert per bin. NOT Monte Carlo - no scatter.
    pub fn filtered(&self, xs: &CrossSections, material: &str, thickness_mm: f64, density: Option<f64>) -> Result<Spectrum, String> {
        if thickness_mm <= 0.0 {
            return Ok(self.clone());
        }
        let mu = xs.mu(material, &self.e, Kind::Total, density)?;
        let mut out = self.clone();
        out.n = self.n.iter().zip(&mu).map(|(n, m)| n * (-m * thickness_mm / 10.0).exp()).collect();
        let rho = match density {
            Some(d) if d != 0.0 => d,
            _ => xs.density(material)?,
        };
        out.filters.push((material.to_string(), thickness_mm, rho));
        if let Some(k) = self.kerma_uGy_per_mAs {
            let k0 = self.kerma_integral();
            let k1 = out.kerma_integral();
            out.kerma_uGy_per_mAs = Some(if k0 > 0.0 { k * k1 / k0 } else { 0.0 });
        }
        Ok(out)
    }

    /// SUM Phi.E.(mu_en/rho)_air in the spectrum's own units.
    pub fn kerma_integral(&self) -> f64 {
        let de = gradient(&self.e);
        let v: Vec<f64> = (0..self.e.len()).map(|i| self.n[i] * de[i] * self.e[i] * air_muen1(self.e[i])).collect();
        pairwise_sum(&v)
    }

    /// Unshielded primary air kerma, mGy per mA.min at 1 m (Archer 1994 fit).
    pub fn archer_K_poly(kvp: f64) -> f64 {
        1.222170186 - 0.056644073 * kvp + 0.001227127 * kvp.powi(2) - 3.13579e-6 * kvp.powi(3)
    }

    /// Thickness of `material` halving the beam (narrow beam), mm.
    pub fn hvl(&self, xs: &CrossSections, material: &str, target: f64, density: Option<f64>, kerma: bool) -> Result<f64, String> {
        let mu = xs.mu(material, &self.e, Kind::Total, density)?;
        let de = gradient(&self.e);
        let mue = air_muen(&self.e);
        let w: Vec<f64> = (0..self.e.len())
            .map(|i| if kerma { self.n[i] * self.e[i] * mue[i] * de[i] } else { self.n[i] * self.e[i] * de[i] })
            .collect();
        let base = pairwise_sum(&w);
        if base <= 0.0 {
            return Ok(0.0);
        }
        let mut tmp = vec![0.0; w.len()];
        let mut frac = |t_mm: f64| -> f64 {
            for i in 0..w.len() {
                tmp[i] = w[i] * (-mu[i] * t_mm / 10.0).exp();
            }
            pairwise_sum(&tmp) / base
        };
        let (mut lo, mut hi) = (0.0f64, 1.0f64);
        for _ in 0..80 {
            if frac(hi) <= target {
                break;
            }
            hi *= 2.0;
            if hi > 1e5 {
                return Ok(f64::INFINITY);
            }
        }
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if frac(mid) > target {
                lo = mid;
            } else {
                hi = mid;
            }
            if hi - lo < 1e-7 {
                break;
            }
        }
        Ok(0.5 * (lo + hi))
    }

    pub fn hvl_al(&self, xs: &CrossSections) -> Result<f64, String> {
        self.hvl(xs, "aluminium", 0.5, None, true)
    }

    /// Monoenergetic energy with the same air-kerma HVL, keV.
    pub fn effective_energy(&self, xs: &CrossSections, material: &str) -> Result<f64, String> {
        let h_mm = self.hvl(xs, material, 0.5, None, true)?;
        if !h_mm.is_finite() || h_mm <= 0.0 {
            return Ok(0.0);
        }
        let target_mu = 2f64.ln() / (h_mm / 10.0);
        let a = xs.emin.max(1.0).log10();
        let b = xs.emax.min(1000.0).log10();
        let num = 4000usize;
        let step = (b - a) / (num as f64 - 1.0);
        let grid: Vec<f64> = (0..num).map(|i| if i == num - 1 { 10f64.powf(b) } else { 10f64.powf(i as f64 * step + a) }).collect();
        let mu = xs.mu(material, &grid, Kind::Total, None)?;
        let d: Vec<f64> = mu.iter().map(|m| (m - target_mu).abs()).collect();
        Ok(grid[argmin(&d)])
    }

    pub fn homogeneity(&self, xs: &CrossSections, material: &str) -> Result<f64, String> {
        let h1 = self.hvl(xs, material, 0.5, None, true)?;
        let h2 = self.hvl(xs, material, 0.25, None, true)?;
        Ok(if (h2 - h1) > 0.0 { h1 / (h2 - h1) } else { 0.0 })
    }
}

// -- tungsten K lines -------------------------------------------------------
pub const W_K_LINES: [(f64, f64); 6] =
    [(57.981, 0.2877), (59.318, 0.5006), (66.950, 0.0564), (67.244, 0.1088), (67.653, 0.0025), (69.035, 0.0431)];
pub const W_K_EDGE: f64 = 69.525;

/// Built-in fallback: Kramers' law + target self-filtration + W K lines.
pub fn spectrum_kramers(xs: &CrossSections, kvp: f64, filters: &[(String, f64)], target_angle_deg: f64) -> Result<Spectrum, String> {
    if !(kvp > 1.0) || !kvp.is_finite() {
        return Err("Tube potential must be greater than 1 kVp.".into());
    }
    let n_e = kvp.ceil() as usize; // len(np.arange(1, kvp+1, 1))
    let e: Vec<f64> = (0..n_e).map(|i| 1.0 + i as f64).collect();
    let mut n: Vec<f64> = e.iter().map(|&x| if x < kvp { (kvp / x.max(1e-6) - 1.0).max(0.0) } else { 0.0 }).collect();
    let mu_w = xs.mu("tungsten", &e, Kind::Total, None)?;
    let depth_cm = 2.0e-4 / target_angle_deg.to_radians().sin().max(1e-3);
    for i in 0..n.len() {
        n[i] *= (-mu_w[i] * depth_cm).exp();
    }
    if kvp > W_K_EDGE {
        let brem_area = trapz(&n, &e);
        let drive = (kvp - W_K_EDGE) / kvp.max(1.0);
        let k_total = 0.55 * brem_area * drive;
        for &(e_line, inten) in W_K_LINES.iter() {
            if e_line < kvp {
                let d: Vec<f64> = e.iter().map(|x| (x - e_line).abs()).collect();
                let j = argmin(&d);
                n[j] += k_total * inten;
            }
        }
    }
    let mut s = Spectrum::new(e, n, kvp, "Kramers (approximate)");
    for (m, mm) in filters {
        s = s.filtered(xs, m, *mm, None)?;
    }
    Ok(s)
}

// -- SpekPy, from the pre-computed table ------------------------------------
pub struct SpekTable {
    pub kv_min: u32,
    pub kv_max: u32,
    pub th_min: u32,
    pub th_max: u32,
    /// entries[th - th_min][kv - kv_min] = (kerma, e_first, N)
    pub entries: Vec<Vec<(f64, f64, Vec<f32>)>>,
}

static SPEK_BIN: &[u8] = include_bytes!("../data/spekpy_table.bin");

thread_local! {
    static SPEK: Result<SpekTable, String> = parse_spek(SPEK_BIN);
}

fn parse_spek(z: &[u8]) -> Result<SpekTable, String> {
    let raw = miniz_oxide::inflate::decompress_to_vec_zlib(z).map_err(|e| format!("SpekPy table: {:?}", e))?;
    if raw.len() < 20 || &raw[0..4] != b"SPK1" {
        return Err("SpekPy table: bad header".into());
    }
    let u = |o: usize| u32::from_le_bytes(raw[o..o + 4].try_into().unwrap());
    let f = |o: usize| f64::from_le_bytes(raw[o..o + 8].try_into().unwrap());
    let (kv_min, kv_max, th_min, th_max) = (u(4), u(8), u(12), u(16));
    let mut p = 20;
    let mut entries = vec![];
    for _ in th_min..=th_max {
        let mut row = vec![];
        for _ in kv_min..=kv_max {
            let k = f(p);
            let nb = u(p + 8) as usize;
            let e0 = f(p + 12);
            p += 20;
            let n: Vec<f32> = (0..nb).map(|i| f32::from_le_bytes(raw[p + 4 * i..p + 4 * i + 4].try_into().unwrap())).collect();
            p += 4 * nb;
            row.push((k, e0, n));
        }
        entries.push(row);
    }
    Ok(SpekTable { kv_min, kv_max, th_min, th_max, entries })
}

pub fn spek_range() -> Option<(u32, u32, u32, u32)> {
    SPEK.with(|t| t.as_ref().ok().map(|t| (t.kv_min, t.kv_max, t.th_min, t.th_max)))
}

/// SpekPy's unfiltered spectrum (dk = 1 keV, 1 mAs) and its air kerma, read
/// from the table. kVp is rounded to the nearest whole kV (the table holds
/// SpekPy's exact output at whole kV); the anode angle is interpolated
/// linearly between whole degrees.
pub fn spectrum_spekpy(kvp: f64, anode_angle: f64) -> Result<Spectrum, String> {
    SPEK.with(|t| {
        let t = t.as_ref().map_err(|e| e.clone())?;
        let kv = kvp.round();
        if kv < t.kv_min as f64 || kv > t.kv_max as f64 || !kv.is_finite() {
            return Err(format!(
                "The SpekPy table covers {}-{} kVp. Use the built-in model outside that range.",
                t.kv_min, t.kv_max
            ));
        }
        if !(anode_angle >= t.th_min as f64 && anode_angle <= t.th_max as f64) {
            return Err(format!("The SpekPy table covers anode angles {}-{} deg.", t.th_min, t.th_max));
        }
        let ki = (kv as u32 - t.kv_min) as usize;
        let a0 = anode_angle.floor().min(t.th_max as f64);
        let f = anode_angle - a0;
        let i0 = (a0 as u32 - t.th_min) as usize;
        let i1 = (i0 + 1).min(t.entries.len() - 1);
        let (k0, e0, n0) = &t.entries[i0][ki];
        let (k1, _, n1) = &t.entries[i1][ki];
        let n: Vec<f64> = n0.iter().zip(n1).map(|(a, b)| (1.0 - f) * (*a as f64) + f * (*b as f64)).collect();
        let e: Vec<f64> = (0..n.len()).map(|i| e0 + i as f64).collect();
        let k = (1.0 - f) * k0 + f * k1;
        let mut s = Spectrum::new(e, n, kv, "SpekPy");
        s.absolute = true;
        s.kerma_uGy_per_mAs = Some(k);
        Ok(s)
    })
}

// SpekPy's own filter attenuation (mu per mm) and air mu_en, for the
// materials app.py maps to SpekPy names - see tools/build_spekpy_table.py.
pub struct SpekMu {
    pub e0: f64,
    pub kerma_const: f64,
    pub muen: Vec<f64>,
    pub mats: Vec<(String, Vec<f64>)>,
}

static SPEKMU_BIN: &[u8] = include_bytes!("../data/spekpy_mu.bin");

thread_local! {
    static SPEKMU: Result<SpekMu, String> = parse_spekmu(SPEKMU_BIN);
}

fn parse_spekmu(z: &[u8]) -> Result<SpekMu, String> {
    let raw = miniz_oxide::inflate::decompress_to_vec_zlib(z).map_err(|e| format!("SpekPy mu table: {:?}", e))?;
    if raw.len() < 24 || &raw[0..4] != b"SPM1" {
        return Err("SpekPy mu table: bad header".into());
    }
    let u = |o: usize| u32::from_le_bytes(raw[o..o + 4].try_into().unwrap()) as usize;
    let f = |o: usize| f64::from_le_bytes(raw[o..o + 8].try_into().unwrap());
    let nb = u(4);
    let (e0, c) = (f(8), f(16));
    let mut p = 24;
    let muen: Vec<f64> = (0..nb).map(|i| f(p + 8 * i)).collect();
    p += 8 * nb;
    let nm = u(p);
    p += 4;
    let mut mats = vec![];
    for _ in 0..nm {
        let l = raw[p] as usize;
        let name = String::from_utf8_lossy(&raw[p + 1..p + 1 + l]).to_string();
        p += 1 + l;
        mats.push((name, (0..nb).map(|i| f(p + 8 * i)).collect()));
        p += 8 * nb;
    }
    Ok(SpekMu { e0, kerma_const: c, muen, mats })
}

/// app.py's `spectrum_spekpy(kvp, angle, filters)`: SpekPy applies the
/// filters with its own attenuation data and reports its own air kerma.
/// Materials app.py maps to SpekPy names (Al, Cu, Be, water, air, tissue,
/// Pb, W) use SpekPy's coefficients; any other material - which app.py
/// would silently drop - is applied with the Elam data instead.
pub fn spectrum_spekpy_filtered(xs: &CrossSections, kvp: f64, anode_angle: f64, filters: &[(String, f64)]) -> Result<Spectrum, String> {
    let mut s = spectrum_spekpy(kvp, anode_angle)?;
    if filters.is_empty() {
        return Ok(s);
    }
    SPEKMU.with(|t| -> Result<(), String> {
        let t = t.as_ref().map_err(|e| e.clone())?;
        let idx = |e: f64| ((e - t.e0).round() as usize).min(t.muen.len() - 1);
        for (m, mm) in filters {
            if *mm <= 0.0 {
                continue;
            }
            match t.mats.iter().find(|(k, _)| k == m) {
                Some((_, mu)) => {
                    for i in 0..s.n.len() {
                        s.n[i] *= (-mu[idx(s.e[i])] * mm).exp();
                    }
                    s.filters.push((m.clone(), *mm, xs.density(m).unwrap_or(0.0)));
                }
                None => {
                    let k = s.kerma_uGy_per_mAs;
                    s = s.filtered(xs, m, *mm, None)?;
                    s.kerma_uGy_per_mAs = k;
                }
            }
        }
        let v: Vec<f64> = (0..s.n.len()).map(|i| s.n[i] * s.e[i] * t.muen[idx(s.e[i])]).collect();
        s.kerma_uGy_per_mAs = Some(t.kerma_const * pairwise_sum(&v));
        Ok(())
    })?;
    Ok(s)
}

pub fn make_spectrum(xs: &CrossSections, kvp: f64, filters: &[(String, f64)], anode_angle: f64, use_spek: bool) -> Result<Spectrum, String> {
    if use_spek {
        spectrum_spekpy_filtered(xs, kvp, anode_angle, filters)
    } else {
        spectrum_kramers(xs, kvp, filters, anode_angle)
    }
}

/// Total Al filtration reproducing a measured HVL at this kVp.
/// Returns (Some(mm), achieved) or (None, achievable).
pub fn solve_filtration_from_hvl(
    xs: &CrossSections,
    kvp: f64,
    measured_hvl_mm: f64,
    anode_angle: f64,
    use_spek: bool,
    max_mm: f64,
) -> Result<(Option<f64>, f64), String> {
    if measured_hvl_mm <= 0.0 {
        return Err("Measured HVL must be greater than zero.".into());
    }
    let hvl_at = |mm: f64| -> Result<f64, String> {
        make_spectrum(xs, kvp, &[("aluminium".to_string(), mm)], anode_angle, use_spek)?.hvl_al(xs)
    };
    let lo_hvl = hvl_at(0.0)?;
    if lo_hvl > measured_hvl_mm {
        return Ok((None, lo_hvl));
    }
    let (mut lo, mut hi) = (0.0, 1.0);
    while hvl_at(hi)? < measured_hvl_mm {
        hi *= 2.0;
        if hi > max_mm {
            return Ok((None, hvl_at(max_mm)?));
        }
    }
    for _ in 0..50 {
        let mid = 0.5 * (lo + hi);
        if hvl_at(mid)? < measured_hvl_mm {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-4 {
            break;
        }
    }
    let mm = 0.5 * (lo + hi);
    Ok((Some(mm), hvl_at(mm)?))
}

/// Air kerma in uGy per mAs at 1 m, or None without an absolute scale.
/// A user calibration (scale factor) takes precedence over SpekPy's own.
pub fn spectrum_output(sp: &Spectrum, cal_scale: Option<f64>) -> Option<f64> {
    match cal_scale {
        Some(s) => Some(sp.kerma_integral() * s),
        None => sp.kerma_uGy_per_mAs,
    }
}

// -- IEC 61267 RQA ----------------------------------------------------------
/// (name, kV, nominal added Al mm, target HVL mm Al, confidence)
pub const RQA_ED2: [(&str, f64, f64, f64, &str); 9] = [
    ("RQA2", 40.0, 4.0, 2.2, "verified"),
    ("RQA3", 50.0, 10.0, 3.8, "verified"),
    ("RQA4", 60.0, 16.0, 5.4, "unverified"),
    ("RQA5", 70.0, 21.0, 6.8, "verified"),
    ("RQA6", 80.0, 26.0, 8.2, "unverified"),
    ("RQA7", 90.0, 30.0, 9.2, "unverified"),
    ("RQA8", 100.0, 34.0, 10.1, "unverified"),
    ("RQA9", 120.0, 40.0, 11.6, "verified"),
    ("RQA10", 150.0, 45.0, 13.3, "unverified"),
];
pub const RQA_EDITION: &str = "IEC 61267 Ed.2 (2005)";
/// IEC 62220-1 Table 3, SNR_in^2 (mm^-2.uGy^-1)
pub const RQA_SNR_IN2: [(&str, f64); 4] = [("RQA3", 20673.0), ("RQA5", 29653.0), ("RQA7", 32490.0), ("RQA9", 31007.0)];

pub fn rqa_lookup(q: &str) -> Result<(f64, f64, f64, &'static str), String> {
    RQA_ED2
        .iter()
        .find(|r| r.0 == q)
        .map(|r| (r.1, r.2, r.3, r.4))
        .ok_or_else(|| format!("{} is not a defined RQA quality.", q))
}

pub fn snr_in2(q: &str) -> Option<f64> {
    RQA_SNR_IN2.iter().find(|r| r.0 == q).map(|r| r.1)
}

/// Photon fluence per unit air kerma (mm^-2.uGy^-1), or None if not absolute.
pub fn rqa_fluence_per_kerma(sp: &Spectrum, cal_scale: Option<f64>) -> Option<f64> {
    let k = spectrum_output(sp, cal_scale)?;
    if k == 0.0 {
        return None;
    }
    Some(sp.fluence() / 100.0 / k)
}

fn solve_al(xs: &CrossSections, base: &Spectrum, target: f64, max_mm: f64) -> Result<(Option<f64>, f64), String> {
    let hvl_after = |mm: f64| -> Result<f64, String> { base.filtered(xs, "aluminium", mm, None)?.hvl_al(xs) };
    let h0 = hvl_after(0.0)?;
    if h0 > target {
        return Ok((None, h0));
    }
    let (mut lo, mut hi) = (0.0, 1.0);
    while hvl_after(hi)? < target {
        hi *= 2.0;
        if hi > max_mm {
            return Ok((None, hvl_after(max_mm)?));
        }
    }
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if hvl_after(mid)? < target {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-4 {
            break;
        }
    }
    let mm = 0.5 * (lo + hi);
    Ok((Some(mm), hvl_after(mm)?))
}

/// Added Al to reach a quality's target HVL. (mm or None, achieved, target)
pub fn rqa_solve_filtration(xs: &CrossSections, base: &Spectrum, quality: &str, max_mm: f64) -> Result<(Option<f64>, f64, f64), String> {
    let (_kv, _nom, target, _c) = rqa_lookup(quality)?;
    let (mm, ach) = solve_al(xs, base, target, max_mm)?;
    Ok((mm, ach, target))
}

/// As above, on top of a fixed copper thickness.
pub fn rqa_solve_filtration_cu(
    xs: &CrossSections,
    base: &Spectrum,
    quality: &str,
    cu_mm: f64,
    max_mm: f64,
) -> Result<(Option<f64>, f64, f64), String> {
    let (_kv, _nom, target, _c) = rqa_lookup(quality)?;
    let base_cu = base.filtered(xs, "copper", cu_mm, None)?;
    let (mm, ach) = solve_al(xs, &base_cu, target, max_mm)?;
    Ok((mm, ach, target))
}

// -- radioisotope sources ---------------------------------------------------
pub const ISOTOPE_LINES: [(&str, &[(f64, f64)]); 15] = [
    ("Am-241", &[(26.34, 0.024), (59.541, 0.359)]),
    ("Ba-133", &[(53.2, 0.022), (81.0, 0.367), (276.4, 0.0716), (302.9, 0.183), (356.0, 0.619), (383.8, 0.089)]),
    ("Co-57", &[(122.06, 0.856), (136.47, 0.1068)]),
    ("Co-60", &[(1173.2, 0.9985), (1332.5, 0.9998)]),
    ("Cs-137", &[(661.657, 0.851)]),
    ("F-18", &[(511.0, 1.936)]),
    ("Ga-67", &[(93.3, 0.39), (184.6, 0.214), (300.2, 0.168), (393.5, 0.046)]),
    ("I-123", &[(159.0, 0.833), (529.0, 0.014)]),
    ("I-125", &[(27.4, 0.716), (31.0, 0.196), (35.5, 0.067)]),
    ("I-131", &[(80.2, 0.026), (284.3, 0.061), (364.5, 0.815), (636.9, 0.072), (722.9, 0.0179)]),
    ("In-111", &[(171.3, 0.907), (245.4, 0.941)]),
    ("Ir-192", &[(206.0, 0.032), (296.0, 0.283), (308.5, 0.293), (316.5, 0.83), (468.1, 0.477), (604.4, 0.0823), (612.5, 0.0534)]),
    ("Se-75", &[(136.0, 0.586), (264.7, 0.589), (279.5, 0.25), (400.7, 0.114)]),
    ("Tc-99m", &[(140.5, 0.89)]),
    ("Xe-133", &[(81.0, 0.37)]),
];

pub fn isotope_lines(name: &str) -> Result<&'static [(f64, f64)], String> {
    ISOTOPE_LINES.iter().find(|r| r.0 == name).map(|r| r.1).ok_or_else(|| format!("{} is not in the isotope line list.", name))
}

pub fn isotope_lines_beyond_xsdata(xs: &CrossSections, name: &str) -> Result<Vec<f64>, String> {
    Ok(isotope_lines(name)?.iter().filter(|(e, _)| *e > xs.emax).map(|(e, _)| *e).collect())
}

pub fn make_isotope_spectrum(name: &str) -> Result<Spectrum, String> {
    let lines = isotope_lines(name)?;
    let top = lines.iter().map(|(e, _)| *e).fold(f64::MIN, f64::max).ceil() as i64 + 10;
    let e: Vec<f64> = (0..top).map(|i| 1.0 + i as f64).collect();
    let mut n = vec![0.0; e.len()];
    for &(el, inten) in lines {
        let d: Vec<f64> = e.iter().map(|x| (x - el).abs()).collect();
        n[argmin(&d)] += inten;
    }
    Ok(Spectrum::new(e, n, top as f64, &format!("isotope: {}", name)))
}
