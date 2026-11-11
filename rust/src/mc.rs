//! SECTION 3 - Monte Carlo photon transport through a 1D slab.
//! Line-for-line port of `_mc_slab`, `_mc_batch` and `run_monte_carlo`.
//!
//! The random generator is MT19937 seeded exactly like numpy's legacy
//! `np.random.seed(int)` with `np.random.random()` = 53-bit doubles, so a
//! run consumes the same random stream as the pure-numpy path of app.py.

use crate::num::{searchsorted_left, trapz};
use crate::spectrum::Spectrum;
use crate::xs::{air_muen1, CrossSections, Kind};

pub const ELECTRON_MASS_KEV: f64 = 510.998950;

pub struct Mt19937 {
    mt: [u32; 624],
    idx: usize,
}

impl Mt19937 {
    pub fn new(seed: u32) -> Mt19937 {
        let mut mt = [0u32; 624];
        let mut s = seed;
        for (pos, slot) in mt.iter_mut().enumerate() {
            *slot = s;
            s = 1812433253u32.wrapping_mul(s ^ (s >> 30)).wrapping_add(pos as u32 + 1);
        }
        Mt19937 { mt, idx: 624 }
    }

    fn gen(&mut self) {
        const UPPER: u32 = 0x8000_0000;
        const LOWER: u32 = 0x7fff_ffff;
        const MATRIX_A: u32 = 0x9908_b0df;
        for i in 0..624 {
            let y = (self.mt[i] & UPPER) | (self.mt[(i + 1) % 624] & LOWER);
            let mut v = self.mt[(i + 397) % 624] ^ (y >> 1);
            if y & 1 != 0 {
                v ^= MATRIX_A;
            }
            self.mt[i] = v;
        }
        self.idx = 0;
    }

    #[inline]
    pub fn next_u32(&mut self) -> u32 {
        if self.idx >= 624 {
            self.gen();
        }
        let mut y = self.mt[self.idx];
        self.idx += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }

    /// numpy random_sample(): uniform on [0, 1) with 53-bit resolution.
    #[inline]
    pub fn random(&mut self) -> f64 {
        let a = (self.next_u32() >> 5) as f64;
        let b = (self.next_u32() >> 6) as f64;
        (a * 67108864.0 + b) / 9007199254740992.0
    }
}

/// Log-log interpolation of a cross-section array at energy x (keV).
#[inline]
fn interp_log(x: f64, log_e: &[f64], y: &[f64]) -> f64 {
    let lx = x.ln();
    let n = log_e.len();
    if lx <= log_e[0] {
        return y[0];
    }
    if lx >= log_e[n - 1] {
        return y[n - 1];
    }
    let (mut lo, mut hi) = (0usize, n - 1);
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if log_e[mid] <= lx {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let (y0, y1) = (y[lo], y[hi]);
    let f = (lx - log_e[lo]) / (log_e[hi] - log_e[lo]);
    if y0 <= 0.0 || y1 <= 0.0 {
        return y0 + f * (y1 - y0);
    }
    (y0.ln() + f * (y1.ln() - y0.ln())).exp()
}

/// Klein-Nishina by direct rejection. Returns (E_out, cos_theta).
#[inline]
fn sample_compton(e: f64, rng: &mut Mt19937) -> (f64, f64) {
    let a = e / ELECTRON_MASS_KEV;
    for _ in 0..500 {
        let ct = 2.0 * rng.random() - 1.0;
        let r = 1.0 / (1.0 + a * (1.0 - ct));
        let f = r * r * (r + 1.0 / r - (1.0 - ct * ct));
        if 2.0 * rng.random() <= f {
            return (e * r, ct);
        }
    }
    (e, 1.0)
}

#[inline]
fn lin_at(v: f64, xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len();
    if v <= xs[0] {
        return ys[0];
    }
    if v >= xs[n - 1] {
        return ys[n - 1];
    }
    let (mut lo, mut hi) = (0usize, n - 1);
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if xs[mid] <= v {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let f = (v - xs[lo]) / (xs[hi] - xs[lo]);
    ys[lo] + f * (ys[hi] - ys[lo])
}

#[inline]
fn invert(v: f64, ys: &[f64], xs: &[f64]) -> f64 {
    let n = ys.len();
    if v <= ys[0] {
        return xs[0];
    }
    if v >= ys[n - 1] {
        return xs[n - 1];
    }
    let (mut lo, mut hi) = (0usize, n - 1);
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if ys[mid] <= v {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let d = ys[hi] - ys[lo];
    let f = if d <= 0.0 { 0.0 } else { (v - ys[lo]) / d };
    xs[lo] + f * (xs[hi] - xs[lo])
}

/// Coherent scattering angle from the real atomic form factor.
#[inline]
fn sample_rayleigh(e: f64, ff_x: &[f64], ff_cum: &[f64], rng: &mut Mt19937) -> f64 {
    let lam = 12.398 / e;
    let c_max = lin_at(1.0 / lam, ff_x, ff_cum);
    if c_max <= 0.0 {
        return 1.0;
    }
    for _ in 0..200 {
        let x = invert(rng.random() * c_max, ff_cum, ff_x);
        let mut st2 = x * lam;
        if st2 > 1.0 {
            st2 = 1.0;
        }
        let ct = 1.0 - 2.0 * st2 * st2;
        if rng.random() <= 0.5 * (1.0 + ct * ct) {
            return ct;
        }
    }
    1.0
}

#[inline]
fn new_direction(u: f64, cos_t: f64, rng: &mut Mt19937) -> f64 {
    let st = (1.0 - cos_t * cos_t).max(0.0).sqrt();
    let su = (1.0 - u * u).max(0.0).sqrt();
    let phi = 2.0 * std::f64::consts::PI * rng.random();
    u * cos_t + su * st * phi.cos()
}

/// Everything the transport kernel needs, prepared once per run.
pub struct Kernel {
    pub d_cm: f64,
    pub log_e: Vec<f64>,
    pub mu_ph: Vec<f64>,
    pub mu_co: Vec<f64>,
    pub mu_in: Vec<f64>,
    pub resp: Vec<f64>,
    pub u_min: f64,
    pub e_min: f64,
    pub ff_x: Vec<f64>,
    pub ff_cum: Vec<f64>,
    pub k_edge: f64,
    pub k_yield: f64,
    pub k_lines: Vec<f64>,
    pub k_cum: Vec<f64>,
    pub w_cut: f64,
    pub max_scatter: usize,
    pub edges: Vec<f64>,
}

impl Kernel {
    /// One photon history: (transmitted weight, transmitted response, backscatter, collisions)
    fn slab(&self, e0: f64, hist: &mut [f64], rng: &mut Mt19937) -> (f64, f64, f64, u64) {
        let (le, ph_a, co_a, in_a, resp) = (&self.log_e, &self.mu_ph, &self.mu_co, &self.mu_in, &self.resp);
        let d_cm = self.d_cm;
        let mut e = e0;
        let mut z = 0.0f64;
        let mut u = 1.0f64;
        let mut w = 1.0f64;
        let (mut t_w, mut t_e, mut b_w) = (0.0f64, 0.0f64, 0.0f64);
        let mut ncoll = 0u64;
        let edges = &self.edges;
        let nb = edges.len() - 1;

        for _ in 0..self.max_scatter {
            let ph = interp_log(e, le, ph_a);
            let co = interp_log(e, le, co_a);
            let inc = interp_log(e, le, in_a);
            let tot = ph + co + inc;
            if tot <= 0.0 {
                if u >= self.u_min && e >= self.e_min {
                    t_w += w;
                    t_e += w * interp_log(e, le, resp);
                }
                break;
            }
            let (dist, forward) = if u > 1e-9 {
                ((d_cm - z) / u, true)
            } else if u < -1e-9 {
                (z / (-u), false)
            } else {
                (1e30, true)
            };
            let p_esc = if dist < 1e29 { (-tot * dist).exp() } else { 0.0 };
            if forward {
                if u >= self.u_min && e >= self.e_min {
                    let c = w * p_esc;
                    t_w += c;
                    t_e += c * interp_log(e, le, resp);
                    if edges[0] <= e && e < edges[nb] {
                        let (mut lo2, mut hi2) = (0usize, nb);
                        while hi2 - lo2 > 1 {
                            let md = (lo2 + hi2) / 2;
                            if edges[md] <= e {
                                lo2 = md;
                            } else {
                                hi2 = md;
                            }
                        }
                        hist[lo2] += c;
                    }
                }
            } else {
                b_w += w * p_esc;
            }

            w *= 1.0 - p_esc; // forced collision
            if w < self.w_cut {
                if rng.random() > 0.25 {
                    break;
                }
                w *= 4.0;
            }

            let r = rng.random();
            let mut s = -(1.0 - r * (1.0 - p_esc)).ln() / tot;
            if s > dist {
                s = dist;
            }
            z += s * u;
            if z < 0.0 {
                z = 0.0;
            } else if z > d_cm {
                z = d_cm;
            }
            ncoll += 1;

            let mut p_fl = 0.0;
            if e > self.k_edge && self.k_yield > 0.0 {
                p_fl = ph * self.k_yield;
            }
            w *= (co + inc + p_fl) / tot;
            if w < self.w_cut {
                if rng.random() > 0.25 {
                    break;
                }
                w *= 4.0;
            }

            let chan = co + inc + p_fl;
            if chan <= 0.0 {
                break;
            }
            let r2 = rng.random() * chan;
            let cos_t;
            if r2 < p_fl {
                let xi = rng.random();
                let mut k = 0usize;
                while k < self.k_cum.len() - 1 && self.k_cum[k] < xi {
                    k += 1;
                }
                e = self.k_lines[k];
                u = 2.0 * rng.random() - 1.0;
                continue;
            } else if r2 < p_fl + co {
                cos_t = sample_rayleigh(e, &self.ff_x, &self.ff_cum, rng);
            } else {
                let (e2, ct) = sample_compton(e, rng);
                e = e2;
                cos_t = ct;
                if e < 1.0 {
                    break;
                }
            }
            u = new_direction(u, cos_t, rng);
        }
        (t_w, t_e, b_w, ncoll)
    }
}

#[derive(Clone, Copy, PartialEq)]
pub enum Score {
    Kerma,
    Fluence,
}

pub struct McParams {
    pub material: String,
    pub thickness_mm: f64,
    pub n_hist: u64,
    pub density: Option<f64>,
    pub seed: u32,
    pub batch: u64,
    pub score: Score,
    pub acceptance_deg: f64,
    pub e_threshold: f64,
    pub fluorescence: bool,
}

/// A run in progress, advanced batch by batch so a UI can show progress.
pub struct McRun {
    pub p: McParams,
    k: Kernel,
    rng: Mt19937,
    e: Vec<f64>,
    n: Vec<f64>,
    cdf: Vec<f64>,
    hist: Vec<f64>,
    tot_w: f64,
    tot_w2: f64,
    tot_e: f64,
    tot_b: f64,
    tot_c: u64,
    pub done: u64,
}

pub struct McResult {
    pub b_mc: f64,
    pub b_narrow: f64,
    pub buildup: f64,
    pub rel_std_err: f64,
    pub backscatter: f64,
    pub mean_collisions: f64,
    pub hist_e: Vec<f64>,
    pub hist_w: Vec<f64>,
    pub hist_norm: Vec<f64>,
}

impl McRun {
    pub fn new(xs: &CrossSections, sp: &Spectrum, p: McParams) -> Result<McRun, String> {
        let rng = Mt19937::new(p.seed);
        let m = xs.mat(&p.material)?;
        let rho = p.density.unwrap_or(m.density);
        let d_cm = p.thickness_mm / 10.0;
        let u_min = p.acceptance_deg.max(0.01).min(90.0).to_radians().cos();
        let ff_x = xs.ff_x.clone();
        let ff_cum = m.ffcum.clone();
        let (ff_x, ff_cum) = match (ff_x, ff_cum) {
            (Some(a), Some(b)) => (a, b),
            _ => {
                return Err(format!(
                    "xsdata.npz has no form factor data for {}. Replace it with the version that came with this build.",
                    p.material
                ))
            }
        };
        let fl = m.fluor.clone().ok_or_else(|| format!("xsdata.npz has no K-edge data for {}.", p.material))?;
        let (k_edge, k_yield) = if p.fluorescence { (fl.kedge, fl.kyield) } else { (1e9, 0.0) };

        let e = sp.e.clone();
        let n: Vec<f64> = sp.n.iter().map(|v| v.max(0.0)).collect();
        if crate::num::sum(&n) <= 0.0 {
            return Err("Spectrum is empty.".into());
        }
        let grid: Vec<f64> = xs.log_e.iter().map(|v| v.exp()).collect();
        let resp: Vec<f64> = match p.score {
            Score::Kerma => grid.iter().map(|&g| g * air_muen1(g)).collect(),
            Score::Fluence => grid.clone(),
        };
        let de = crate::num::gradient(&e);
        let wgt: Vec<f64> = n.iter().zip(&de).map(|(a, b)| a * b).collect();
        if crate::num::sum(&wgt) <= 0.0 {
            return Err("Spectrum has no photons.".into());
        }
        let mut cdf = Vec::with_capacity(wgt.len());
        let mut acc = 0.0;
        for w in &wgt {
            acc += w;
            cdf.push(acc);
        }
        let last = *cdf.last().unwrap();
        for c in cdf.iter_mut() {
            *c /= last;
        }
        let deg = if e.len() > 1 {
            let d: Vec<f64> = (0..e.len() - 1).map(|i| e[i + 1] - e[i]).collect();
            crate::num::mean(&d)
        } else {
            1.0
        };
        let mut edges: Vec<f64> = e.iter().map(|v| v - 0.5 * deg).collect();
        edges.push(e[e.len() - 1] + 0.5 * deg);
        let hist = vec![0.0; edges.len() - 1];
        let k = Kernel {
            d_cm,
            log_e: xs.log_e.clone(),
            mu_ph: m.photo.iter().map(|v| v * rho).collect(),
            mu_co: m.coh.iter().map(|v| v * rho).collect(),
            mu_in: m.incoh.iter().map(|v| v * rho).collect(),
            resp,
            u_min,
            e_min: p.e_threshold,
            ff_x,
            ff_cum,
            k_edge,
            k_yield,
            k_lines: fl.klines.clone(),
            k_cum: fl.kcum.clone(),
            w_cut: 1e-8,
            max_scatter: 400,
            edges,
        };
        Ok(McRun { p, k, rng, e, n, cdf, hist, tot_w: 0.0, tot_w2: 0.0, tot_e: 0.0, tot_b: 0.0, tot_c: 0, done: 0 })
    }

    pub fn finished(&self) -> bool {
        self.done >= self.p.n_hist
    }

    /// Run one batch (as app.py does: draw m energies, then track them).
    pub fn step(&mut self) {
        if self.finished() {
            return;
        }
        let m = self.p.batch.min(self.p.n_hist - self.done) as usize;
        let ne = self.e.len();
        let mut e0 = Vec::with_capacity(m);
        for _ in 0..m {
            let u = self.rng.random();
            let i = searchsorted_left(&self.cdf, u).min(ne - 1);
            e0.push(self.e[i]);
        }
        let (mut sw, mut sw2, mut se, mut sb) = (0.0, 0.0, 0.0, 0.0);
        let mut nc = 0u64;
        for &en in &e0 {
            let (_tw, te, bw, ncoll) = self.k.slab(en, &mut self.hist, &mut self.rng);
            sw += te;
            sw2 += te * te;
            se += interp_log(en, &self.k.log_e, &self.k.resp);
            sb += bw;
            nc += ncoll;
        }
        self.tot_w += sw;
        self.tot_w2 += sw2;
        self.tot_e += se;
        self.tot_b += sb;
        self.tot_c += nc;
        self.done += m as u64;
    }

    pub fn result(&self, xs: &CrossSections) -> Result<McResult, String> {
        let nh = self.p.n_hist as f64;
        let b_mc = if self.tot_e > 0.0 { self.tot_w / self.tot_e } else { 0.0 };
        let mean_r = self.tot_w / nh;
        let var = (self.tot_w2 / nh - mean_r * mean_r).max(0.0);
        let rse = if mean_r > 0.0 { (var / nh).sqrt() / mean_r } else { f64::INFINITY };

        let mu = xs.mu(&self.p.material, &self.e, Kind::Total, self.p.density)?;
        let resp_in: Vec<f64> = self.e.iter().map(|&e| if self.p.score == Score::Kerma { e * air_muen1(e) } else { e }).collect();
        let y0: Vec<f64> = (0..self.e.len()).map(|i| self.n[i] * resp_in[i]).collect();
        let y1: Vec<f64> = (0..self.e.len()).map(|i| self.n[i] * resp_in[i] * (-mu[i] * self.k.d_cm).exp()).collect();
        let ef0 = trapz(&y0, &self.e);
        let ef1 = trapz(&y1, &self.e);
        let b_narrow = if ef0 > 0.0 { ef1 / ef0 } else { 0.0 };
        let edges = &self.k.edges;
        let centres: Vec<f64> = (0..edges.len() - 1).map(|i| 0.5 * (edges[i] + edges[i + 1])).collect();
        let hs = crate::num::sum(&self.hist);
        let hist_norm = if hs > 0.0 { self.hist.iter().map(|v| v / hs).collect() } else { self.hist.clone() };
        Ok(McResult {
            b_mc,
            b_narrow,
            buildup: if b_narrow > 0.0 { b_mc / b_narrow } else { f64::NAN },
            rel_std_err: rse,
            backscatter: self.tot_b / nh,
            mean_collisions: self.tot_c as f64 / nh,
            hist_e: centres,
            hist_w: self.hist.clone(),
            hist_norm,
        })
    }
}

/// Convenience: whole run in one call.
pub fn run_monte_carlo(xs: &CrossSections, sp: &Spectrum, p: McParams) -> Result<McResult, String> {
    let mut r = McRun::new(xs, sp, p)?;
    while !r.finished() {
        r.step();
    }
    r.result(xs)
}

