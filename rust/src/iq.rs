//! SECTION 6 - image quality: slanted-edge MTF, flat-field NNPS and DQE
//! (port of `slanted_edge_mtf`, `_detrend_2d`, `flat_field_nnps`,
//! `nnps_profile`, `compute_dqe`).

use crate::num::{arange, gradient, hanning, interp, mean, percentiles, searchsorted_right};
use rustfft::{num_complex::Complex, FftPlanner};

#[derive(Clone)]
pub struct Image {
    pub w: usize,
    pub h: usize,
    pub data: Vec<f32>,
}

impl Image {
    #[inline]
    pub fn at(&self, r: usize, c: usize) -> f32 {
        self.data[r * self.w + c]
    }
    /// rows r0..r1, cols c0..c1 (Python slice semantics, clamped)
    pub fn crop(&self, r0: i64, r1: i64, c0: i64, c1: i64) -> Image {
        let cl = |v: i64, n: usize| v.max(0).min(n as i64) as usize;
        let (r0, r1, c0, c1) = (cl(r0, self.h), cl(r1, self.h), cl(c0, self.w), cl(c1, self.w));
        let (h, w) = (r1.saturating_sub(r0), c1.saturating_sub(c0));
        let mut data = Vec::with_capacity(w * h);
        for r in r0..r0 + h {
            data.extend_from_slice(&self.data[r * self.w + c0..r * self.w + c0 + w]);
        }
        Image { w, h, data }
    }
    pub fn transpose(&self) -> Image {
        let mut data = vec![0f32; self.data.len()];
        for r in 0..self.h {
            for c in 0..self.w {
                data[c * self.h + r] = self.data[r * self.w + c];
            }
        }
        Image { w: self.h, h: self.w, data }
    }
}

pub struct MtfResult {
    pub freq: Vec<f64>,
    pub mtf: Vec<f64>,
    pub angle_deg: f64,
    pub esf: Vec<f64>,
    pub lsf: Vec<f64>,
    pub bin_mm: f64,
}

pub fn rfft_abs(x: &[f64]) -> Vec<f64> {
    let n = x.len();
    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(n);
    let mut buf: Vec<Complex<f64>> = x.iter().map(|&v| Complex::new(v, 0.0)).collect();
    fft.process(&mut buf);
    buf[..n / 2 + 1].iter().map(|c| c.norm()).collect()
}

/// IEC-style slanted-edge MTF from a single ROI containing one edge.
/// `horizontal` = the edge separates top/bottom (scan columns).
pub fn slanted_edge_mtf(roi: &Image, px_mm: f64, oversample: usize, horizontal: bool) -> Result<MtfResult, String> {
    let work = if horizontal { roi.transpose() } else { roi.clone() };
    let (n_rows, n_cols) = (work.h, work.w);
    if n_rows < 2 || n_cols < 2 {
        return Err("ROI is too small.".into());
    }
    let all: Vec<f64> = work.data.iter().map(|&v| v as f64).collect();
    let p = percentiles(&all, &[2.0, 98.0]);
    let (lo, hi) = (p[0], p[1]);
    if hi <= lo {
        return Err("No contrast found in this ROI - check the edge is actually inside the selected region.".into());
    }
    let mid = 0.5 * (lo + hi);

    let mut cr: Vec<(f64, f64)> = vec![];
    for i in 0..n_rows {
        let line = &all[i * n_cols..(i + 1) * n_cols];
        let idx: Vec<usize> = (0..n_cols - 1).filter(|&j| (line[j] > mid) != (line[j + 1] > mid)).collect();
        if idx.is_empty() {
            continue;
        }
        let j = idx[idx.len() / 2];
        let (y0, y1) = (line[j], line[j + 1]);
        if y1 == y0 {
            continue;
        }
        cr.push((i as f64, j as f64 + (mid - y0) / (y1 - y0)));
    }
    let frac_ok = cr.len() as f64 / n_rows as f64;
    if frac_ok < 0.5 {
        let (dir, wrong) = if horizontal { ("columns", "vertical") } else { ("rows", "horizontal") };
        return Err(format!(
            "Only {:.0}% of {} found a clean crossing (need at least 50%). Contrast range in this ROI: {:.0} to {:.0}. \
             Most likely causes: (1) the ROI doesn't actually straddle an edge - it may sit entirely on one side; \
             (2) it includes a corner of the aperture, where the edge direction changes; (3) this ROI actually has a {} \
             edge in it, which belongs in the OTHER MTF direction, not this one.",
            100.0 * frac_ok,
            dir,
            lo,
            hi,
            wrong
        ));
    }
    // least-squares line col = a*row + b (np.polyfit deg 1)
    let nn = cr.len() as f64;
    let mx = cr.iter().map(|c| c.0).sum::<f64>() / nn;
    let my = cr.iter().map(|c| c.1).sum::<f64>() / nn;
    let sxx: f64 = cr.iter().map(|c| (c.0 - mx) * (c.0 - mx)).sum();
    let sxy: f64 = cr.iter().map(|c| (c.0 - mx) * (c.1 - my)).sum();
    if sxx <= 0.0 {
        return Err("Too few edge crossings to fit the edge.".into());
    }
    let a = sxy / sxx;
    let b = my - a * mx;
    let angle_deg = a.atan().to_degrees();

    let mut dist = Vec::with_capacity(n_rows * n_cols);
    for r in 0..n_rows {
        for c in 0..n_cols {
            dist.push(c as f64 - (a * r as f64 + b));
        }
    }
    let bin_w = 1.0 / oversample as f64;
    let q = percentiles(&dist, &[1.0, 99.0]);
    let bins = arange(q[0], q[1], bin_w);
    if bins.len() < 2 {
        return Err("Too few valid ESF bins - ROI may be too small.".into());
    }
    let nb = bins.len() - 1;
    let mut sums = vec![0.0f64; nb];
    let mut cnt = vec![0usize; nb];
    for (k, &d) in dist.iter().enumerate() {
        let bi = searchsorted_right(&bins, d); // np.digitize
        if bi >= 1 && bi <= nb {
            sums[bi - 1] += all[k];
            cnt[bi - 1] += 1;
        }
    }
    let esf: Vec<f64> = (0..nb).filter(|&k| cnt[k] > 0).map(|k| sums[k] / cnt[k] as f64).collect();
    if esf.len() < 16 {
        return Err("Too few valid ESF bins - ROI may be too small.".into());
    }
    let lsf = gradient(&esf);
    let win = hanning(lsf.len());
    let lsf_w: Vec<f64> = lsf.iter().zip(&win).map(|(a, b)| a * b).collect();
    let mut mtf = rfft_abs(&lsf_w);
    if mtf[0] != 0.0 {
        let m0 = mtf[0];
        for v in mtf.iter_mut() {
            *v /= m0;
        }
    }
    let n = lsf_w.len();
    let d = bin_w * px_mm;
    let freq: Vec<f64> = (0..n / 2 + 1).map(|k| k as f64 / (n as f64 * d)).collect();
    Ok(MtfResult { freq, mtf, angle_deg, esf, lsf, bin_mm: d })
}

/// Subtract a 2nd-order 2D polynomial least-squares fit.
/// Columns [1, y, y^2, x, xy, x^2] as in app.py; coordinates are rescaled to
/// [-1, 1] before solving, which spans the same space (identical fit) but
/// keeps the normal equations well conditioned.
fn detrend2(block: &[f64], n: usize) -> Vec<f64> {
    let s = |i: usize| if n > 1 { 2.0 * i as f64 / (n as f64 - 1.0) - 1.0 } else { 0.0 };
    let mut ata = [[0.0f64; 6]; 6];
    let mut atz = [0.0f64; 6];
    for yi in 0..n {
        let y = s(yi);
        for xi in 0..n {
            let x = s(xi);
            let t = [1.0, y, y * y, x, x * y, x * x];
            let z = block[yi * n + xi];
            for i in 0..6 {
                atz[i] += t[i] * z;
                for j in i..6 {
                    ata[i][j] += t[i] * t[j];
                }
            }
        }
    }
    for i in 0..6 {
        for j in 0..i {
            ata[i][j] = ata[j][i];
        }
    }
    let c = solve6(ata, atz);
    let mut out = vec![0.0; n * n];
    for yi in 0..n {
        let y = s(yi);
        for xi in 0..n {
            let x = s(xi);
            let fit = c[0] + c[1] * y + c[2] * y * y + c[3] * x + c[4] * x * y + c[5] * x * x;
            out[yi * n + xi] = block[yi * n + xi] - fit;
        }
    }
    out
}

fn solve6(mut a: [[f64; 6]; 6], mut b: [f64; 6]) -> [f64; 6] {
    for k in 0..6 {
        let mut p = k;
        for i in k + 1..6 {
            if a[i][k].abs() > a[p][k].abs() {
                p = i;
            }
        }
        a.swap(k, p);
        b.swap(k, p);
        let d = a[k][k];
        if d == 0.0 {
            continue;
        }
        for i in k + 1..6 {
            let f = a[i][k] / d;
            for j in k..6 {
                a[i][j] -= f * a[k][j];
            }
            b[i] -= f * b[k];
        }
    }
    let mut x = [0.0; 6];
    for k in (0..6).rev() {
        let mut s = b[k];
        for j in k + 1..6 {
            s -= a[k][j] * x[j];
        }
        x[k] = if a[k][k] != 0.0 { s / a[k][k] } else { 0.0 };
    }
    x
}

pub struct NnpsResult {
    pub nnps2d: Vec<f64>,
    pub n: usize,
    pub freq: Vec<f64>,
    pub mean_signal: f64,
    pub count: usize,
}

/// fftshift(fftfreq(n, d))
pub fn shifted_freq(n: usize, d: f64) -> Vec<f64> {
    let mut f: Vec<f64> = (0..n).map(|i| if i < (n + 1) / 2 { i as f64 } else { i as f64 - n as f64 } / (d * n as f64)).collect();
    f.rotate_left((n + 1) / 2);
    f
}

/// 2D NNPS: overlapping detrended, windowed, FFT'd ROIs, averaged and
/// normalised by the squared mean signal.
pub fn flat_field_nnps(
    img: &Image,
    px_mm: f64,
    roi: usize,
    step: usize,
    mask_below: f64,
    mut progress: impl FnMut(f64),
) -> Result<NnpsResult, String> {
    let (ny, nx) = (img.h, img.w);
    let wh = hanning(roi);
    let mut win = vec![0.0; roi * roi];
    for i in 0..roi {
        for j in 0..roi {
            win[i * roi + j] = wh[i] * wh[j];
        }
    }
    let wsq: Vec<f64> = win.iter().map(|v| v * v).collect();
    let win_mean_sq = mean(&wsq);
    let mut positions = vec![];
    if ny > roi && nx > roi {
        let mut y0 = 0;
        while y0 < ny - roi {
            let mut x0 = 0;
            while x0 < nx - roi {
                positions.push((y0, x0));
                x0 += step;
            }
            y0 += step;
        }
    }
    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(roi);
    let mut accum = vec![0.0f64; roi * roi];
    let mut count = 0usize;
    let mut valid_means = vec![];
    let mut block = vec![0.0f64; roi * roi];
    let mut buf = vec![Complex::new(0.0, 0.0); roi * roi];
    let mut col = vec![Complex::new(0.0, 0.0); roi];
    let npos = positions.len();
    for (i, &(y0, x0)) in positions.iter().enumerate() {
        let mut masked = false;
        for r in 0..roi {
            let row = &img.data[(y0 + r) * nx + x0..(y0 + r) * nx + x0 + roi];
            for (c, &v) in row.iter().enumerate() {
                if (v as f64) <= mask_below {
                    masked = true;
                }
                block[r * roi + c] = v as f64;
            }
        }
        if masked {
            continue;
        }
        valid_means.push(mean(&block));
        let d = detrend2(&block, roi);
        for k in 0..roi * roi {
            buf[k] = Complex::new(d[k] * win[k], 0.0);
        }
        for r in 0..roi {
            fft.process(&mut buf[r * roi..(r + 1) * roi]);
        }
        for c in 0..roi {
            for r in 0..roi {
                col[r] = buf[r * roi + c];
            }
            fft.process(&mut col);
            for r in 0..roi {
                buf[r * roi + c] = col[r];
            }
        }
        let scale = px_mm * px_mm / (roi * roi) as f64 / win_mean_sq;
        let sh = (roi + 1) / 2; // fftshift
        for r in 0..roi {
            let rs = (r + roi - sh) % roi; // destination of row r after shift
            for c in 0..roi {
                let cs = (c + roi - sh) % roi;
                accum[rs * roi + cs] += buf[r * roi + c].norm_sqr() * scale;
            }
        }
        count += 1;
        if i % 20 == 0 {
            progress(i as f64 / npos as f64);
        }
    }
    if count == 0 {
        return Err("No usable ROIs - flat field may be smaller than expected, or mostly masked-out background.".into());
    }
    let ms = mean(&valid_means);
    let nnps2d: Vec<f64> = accum.iter().map(|v| (v / count as f64) / (ms * ms)).collect();
    Ok(NnpsResult { nnps2d, n: roi, freq: shifted_freq(roi, px_mm), mean_signal: ms, count })
}

/// 1D profile from a band of rows (h) or columns (v) about the centre.
pub fn nnps_profile(nnps2d: &[f64], n: usize, horizontal: bool, band: usize) -> Vec<f64> {
    let c = n / 2;
    let (a, b) = (c.saturating_sub(band), (c + band).min(n));
    let cnt = (b - a) as f64;
    (0..n)
        .map(|k| {
            let mut s = 0.0;
            for j in a..b {
                s += if horizontal { nnps2d[j * n + k] } else { nnps2d[k * n + j] };
            }
            s / cnt
        })
        .collect()
}

/// DQE(f) = MTF^2 / (q K NNPS) on the positive-frequency NNPS grid.
pub fn compute_dqe(freq_mtf: &[f64], mtf: &[f64], freq_nnps: &[f64], nnps: &[f64], q: f64, k_ugy: f64) -> (Vec<f64>, Vec<f64>) {
    let mut fn_ = vec![];
    let mut nn = vec![];
    for (f, v) in freq_nnps.iter().zip(nnps) {
        if *f >= 0.0 {
            fn_.push(*f);
            nn.push(*v);
        }
    }
    let fmax = crate::num::max(freq_mtf);
    let fc: Vec<f64> = fn_.iter().cloned().filter(|f| *f <= fmax).collect();
    let mi = interp(&fc, freq_mtf, mtf);
    let ni = interp(&fc, &fn_, &nn);
    let dqe = mi.iter().zip(&ni).map(|(m, n)| if *n > 0.0 { m * m / (q * k_ugy * n) } else { f64::NAN }).collect();
    (fc, dqe)
}
