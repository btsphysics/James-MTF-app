//! SECTION 1 - cross-section data (port of `CrossSections`) and air kerma.

use crate::npz;
use crate::num::{interp1, pairwise_sum};
use std::collections::BTreeMap;

#[derive(Clone)]
pub struct Fluor {
    pub kedge: f64,
    pub kyield: f64,
    pub klines: Vec<f64>,
    pub kcum: Vec<f64>,
}

pub struct Material {
    pub photo: Vec<f64>,
    pub coh: Vec<f64>,
    pub incoh: Vec<f64>,
    /// log(max(y, 1e-300)) cached for the log-log interpolation
    lphoto: Vec<f64>,
    lcoh: Vec<f64>,
    lincoh: Vec<f64>,
    pub density: f64,
    pub ffcum: Option<Vec<f64>>,
    pub fluor: Option<Fluor>,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Kind {
    Total,
    Photo,
    Coh,
    Incoh,
}

impl Kind {
    pub fn parse(s: &str) -> Kind {
        match s {
            "photo" => Kind::Photo,
            "coh" => Kind::Coh,
            "incoh" => Kind::Incoh,
            _ => Kind::Total,
        }
    }
}

pub struct CrossSections {
    pub e: Vec<f64>,
    pub log_e: Vec<f64>,
    pub emin: f64,
    pub emax: f64,
    pub ff_x: Option<Vec<f64>>,
    pub mats: BTreeMap<String, Material>,
}

fn logs(y: &[f64]) -> Vec<f64> {
    y.iter().map(|&v| v.max(1e-300).ln()).collect()
}

impl CrossSections {
    pub fn from_npz(buf: &[u8]) -> Result<CrossSections, String> {
        let z = npz::load_npz(buf)?;
        let e = z.get("energy_keV").ok_or("xsdata.npz has no energy_keV array")?.clone();
        if e.len() < 2 {
            return Err("energy_keV is too short".into());
        }
        let log_e: Vec<f64> = e.iter().map(|v| v.ln()).collect();
        let mut names: Vec<String> = z.keys().filter_map(|k| k.split_once("__").map(|(m, _)| m.to_string())).collect();
        names.sort();
        names.dedup();
        let mut mats = BTreeMap::new();
        for m in names {
            let g = |f: &str| z.get(&format!("{}__{}", m, f));
            let (photo, coh, incoh, dens) = match (g("photo"), g("coh"), g("incoh"), g("density")) {
                (Some(a), Some(b), Some(c), Some(d)) => (a.clone(), b.clone(), c.clone(), d[0]),
                _ => continue,
            };
            let fluor = match (g("kedge"), g("kyield"), g("klines"), g("kcum")) {
                (Some(a), Some(b), Some(c), Some(d)) => Some(Fluor { kedge: a[0], kyield: b[0], klines: c.clone(), kcum: d.clone() }),
                _ => None,
            };
            mats.insert(
                m.clone(),
                Material {
                    lphoto: logs(&photo),
                    lcoh: logs(&coh),
                    lincoh: logs(&incoh),
                    photo,
                    coh,
                    incoh,
                    density: dens,
                    ffcum: g("ffcum").cloned(),
                    fluor,
                },
            );
        }
        let emin = e[0];
        let emax = *e.last().unwrap();
        Ok(CrossSections { ff_x: z.get("formfactor_x").cloned(), e, log_e, emin, emax, mats })
    }

    pub fn mat(&self, m: &str) -> Result<&Material, String> {
        self.mats.get(m).ok_or_else(|| format!("Material '{}' is not in xsdata.npz.", m))
    }

    pub fn density(&self, m: &str) -> Result<f64, String> {
        Ok(self.mat(m)?.density)
    }

    #[inline]
    fn interp_one(&self, ly: &[f64], e_kev: f64) -> f64 {
        let eq = e_kev.max(self.emin).min(self.emax);
        interp1(eq.ln(), &self.log_e, ly).exp()
    }

    /// Mass attenuation coefficient, cm2/g, at one energy.
    pub fn mass_atten1(&self, m: &Material, e_kev: f64, kind: Kind) -> f64 {
        match kind {
            Kind::Total => self.interp_one(&m.lphoto, e_kev) + self.interp_one(&m.lcoh, e_kev) + self.interp_one(&m.lincoh, e_kev),
            Kind::Photo => self.interp_one(&m.lphoto, e_kev),
            Kind::Coh => self.interp_one(&m.lcoh, e_kev),
            Kind::Incoh => self.interp_one(&m.lincoh, e_kev),
        }
    }

    /// Linear attenuation coefficient, cm^-1, over an energy array.
    pub fn mu(&self, material: &str, e_kev: &[f64], kind: Kind, density: Option<f64>) -> Result<Vec<f64>, String> {
        let m = self.mat(material)?;
        let rho = density.unwrap_or(m.density);
        Ok(e_kev.iter().map(|&e| self.mass_atten1(m, e, kind) * rho).collect())
    }

    pub fn mu1(&self, material: &str, e_kev: f64, kind: Kind, density: Option<f64>) -> Result<f64, String> {
        let m = self.mat(material)?;
        let rho = density.unwrap_or(m.density);
        Ok(self.mass_atten1(m, e_kev, kind) * rho)
    }
}

// -- NIST Hubbell & Seltzer, dry air mu_en/rho -----------------------------
pub const AIR_MUEN_KEV: [f64; 27] = [
    1.0, 1.5, 2.0, 3.0, 3.2029, 3.20291, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 150.0, 200.0,
    300.0, 400.0, 500.0, 600.0, 800.0, 1000.0,
];
pub const AIR_MUEN_CM2G: [f64; 27] = [
    3.599e3, 1.188e3, 5.262e2, 1.614e2, 1.330e2, 1.460e2, 7.636e1, 3.931e1, 2.270e1, 9.446, 4.742, 1.334, 5.389e-1, 1.537e-1,
    6.833e-2, 4.098e-2, 3.041e-2, 2.407e-2, 2.325e-2, 2.496e-2, 2.672e-2, 2.872e-2, 2.949e-2, 2.966e-2, 2.953e-2, 2.882e-2,
    2.789e-2,
];

pub const KEV_TO_J: f64 = 1.602176634e-16;
pub const KERMA_FACTOR: f64 = KEV_TO_J * 1000.0;

thread_local! {
    static AIR_LOG: (Vec<f64>, Vec<f64>) = (
        AIR_MUEN_KEV.iter().map(|v| v.ln()).collect(),
        AIR_MUEN_CM2G.iter().map(|v| v.ln()).collect(),
    );
}

/// Mass energy-absorption coefficient of dry air, cm2/g.
pub fn air_muen1(e_kev: f64) -> f64 {
    let e = e_kev.max(AIR_MUEN_KEV[0]).min(AIR_MUEN_KEV[26]);
    AIR_LOG.with(|(le, lm)| interp1(e.ln(), le, lm).exp())
}

pub fn air_muen(e: &[f64]) -> Vec<f64> {
    e.iter().map(|&v| air_muen1(v)).collect()
}

/// Air kerma in uGy from photons/cm2 per bin.
pub fn air_kerma(e: &[f64], phi: &[f64]) -> f64 {
    let v: Vec<f64> = e.iter().zip(phi).map(|(&e, &p)| p * e * air_muen1(e)).collect();
    pairwise_sum(&v) * KERMA_FACTOR * 1e6
}
