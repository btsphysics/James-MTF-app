//! NumPy-compatible numeric helpers, written so results track numpy closely
//! (pairwise summation, np.interp's exact formula, np.gradient, etc.).

/// numpy's pairwise summation (used by np.sum / np.mean on contiguous f64).
pub fn pairwise_sum(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut res = 0.0;
        for &x in a {
            res += x;
        }
        res
    } else if n <= 128 {
        let mut r = [0.0f64; 8];
        r.copy_from_slice(&a[0..8]);
        let mut i = 8;
        while i < n - (n % 8) {
            for j in 0..8 {
                r[j] += a[i + j];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        while i < n {
            res += a[i];
            i += 1;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        pairwise_sum(&a[..n2]) + pairwise_sum(&a[n2..])
    }
}

pub fn sum(a: &[f64]) -> f64 {
    pairwise_sum(a)
}

pub fn mean(a: &[f64]) -> f64 {
    if a.is_empty() {
        f64::NAN
    } else {
        pairwise_sum(a) / a.len() as f64
    }
}

/// np.trapz(y, x)
pub fn trapz(y: &[f64], x: &[f64]) -> f64 {
    if y.len() < 2 {
        return 0.0;
    }
    let v: Vec<f64> = (0..y.len() - 1).map(|i| (x[i + 1] - x[i]) * (y[i + 1] + y[i]) / 2.0).collect();
    pairwise_sum(&v)
}

/// np.gradient(f) with unit spacing.
pub fn gradient(f: &[f64]) -> Vec<f64> {
    let n = f.len();
    if n < 2 {
        return vec![0.0; n];
    }
    let mut g = vec![0.0; n];
    g[0] = f[1] - f[0];
    g[n - 1] = f[n - 1] - f[n - 2];
    for i in 1..n - 1 {
        g[i] = (f[i + 1] - f[i - 1]) / 2.0;
    }
    g
}

/// Index j with xp[j] <= x < xp[j+1], for xp[0] <= x < xp[n-1].
#[inline]
fn bsearch(xp: &[f64], x: f64) -> usize {
    let mut lo = 0usize;
    let mut hi = xp.len() - 1;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if xp[mid] <= x {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

/// np.interp(x, xp, fp) for a single value.
#[inline]
pub fn interp1(x: f64, xp: &[f64], fp: &[f64]) -> f64 {
    let n = xp.len();
    if n == 0 {
        return f64::NAN;
    }
    if n == 1 {
        return fp[0];
    }
    if x.is_nan() {
        return x;
    }
    if x < xp[0] {
        return fp[0];
    }
    if x > xp[n - 1] {
        return fp[n - 1];
    }
    if x == xp[n - 1] {
        return fp[n - 1];
    }
    let j = bsearch(xp, x);
    let slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]);
    let mut r = slope * (x - xp[j]) + fp[j];
    if r.is_nan() {
        r = slope * (x - xp[j + 1]) + fp[j + 1];
        if r.is_nan() && fp[j] == fp[j + 1] {
            r = fp[j];
        }
    }
    r
}

pub fn interp(x: &[f64], xp: &[f64], fp: &[f64]) -> Vec<f64> {
    x.iter().map(|&v| interp1(v, xp, fp)).collect()
}

/// np.percentile(a, q) with the default 'linear' method.
pub fn percentile_sorted(sorted: &[f64], q: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    let idx = q / 100.0 * (n as f64 - 1.0);
    let lo = idx.floor() as usize;
    let hi = (lo + 1).min(n - 1);
    let t = idx - lo as f64;
    let a = sorted[lo];
    let b = sorted[hi];
    let d = b - a;
    if t >= 0.5 {
        b - d * (1.0 - t)
    } else {
        a + d * t
    }
}

pub fn percentiles(a: &[f64], qs: &[f64]) -> Vec<f64> {
    let mut s: Vec<f64> = a.to_vec();
    s.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    qs.iter().map(|&q| percentile_sorted(&s, q)).collect()
}

/// np.hanning(M)
pub fn hanning(m: usize) -> Vec<f64> {
    if m == 0 {
        return vec![];
    }
    if m == 1 {
        return vec![1.0];
    }
    (0..m)
        .map(|n| 0.5 - 0.5 * (2.0 * std::f64::consts::PI * n as f64 / (m as f64 - 1.0)).cos())
        .collect()
}

/// np.arange(start, stop, step) for floats.
pub fn arange(start: f64, stop: f64, step: f64) -> Vec<f64> {
    let n = ((stop - start) / step).ceil();
    let n = if n.is_finite() && n > 0.0 { n as usize } else { 0 };
    (0..n).map(|i| start + i as f64 * step).collect()
}

/// np.searchsorted(a, v, side='right')
pub fn searchsorted_right(a: &[f64], v: f64) -> usize {
    let (mut lo, mut hi) = (0usize, a.len());
    while lo < hi {
        let mid = (lo + hi) / 2;
        if a[mid] <= v {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

/// np.searchsorted(a, v, side='left')
pub fn searchsorted_left(a: &[f64], v: f64) -> usize {
    let (mut lo, mut hi) = (0usize, a.len());
    while lo < hi {
        let mid = (lo + hi) / 2;
        if a[mid] < v {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

/// np.argmin
pub fn argmin(a: &[f64]) -> usize {
    let mut best = 0;
    for i in 1..a.len() {
        if a[i] < a[best] || (a[best].is_nan() && !a[i].is_nan()) {
            best = i;
        }
    }
    best
}

pub fn max(a: &[f64]) -> f64 {
    a.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
}
