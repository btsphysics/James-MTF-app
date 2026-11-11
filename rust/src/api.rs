//! JSON command interface used by the browser front end (and the native
//! CLI used for validation). Each command takes a JSON object and returns a
//! JSON object; large binary outputs (image views) go to a side buffer.

use crate::archer;
use crate::dicom::{read_dicom, Dicom};
use crate::iq;
use crate::mc::{McParams, McRun, Score};
use crate::spectrum::{self as spm, Spectrum};
use crate::xs::{CrossSections, Kind};
use serde_json::{json, Value};
use std::collections::HashMap;

pub struct State {
    pub xs: Option<CrossSections>,
    pub images: HashMap<String, Dicom>,
    pub mc: Option<McRun>,
    pub bin: Vec<u8>,
    pub progress: Box<dyn FnMut(f64)>,
}

impl State {
    pub fn new(progress: Box<dyn FnMut(f64)>) -> State {
        State { xs: None, images: HashMap::new(), mc: None, bin: vec![], progress }
    }
}

/// f64 -> JSON, mapping NaN to null and infinities to "inf"/"-inf".
pub fn num(v: f64) -> Value {
    if v.is_finite() {
        json!(v)
    } else if v.is_nan() {
        Value::Null
    } else if v > 0.0 {
        json!("inf")
    } else {
        json!("-inf")
    }
}
pub fn arr(v: &[f64]) -> Value {
    Value::Array(v.iter().map(|&x| num(x)).collect())
}
fn opt(v: Option<f64>) -> Value {
    v.map(num).unwrap_or(Value::Null)
}

fn f(req: &Value, k: &str) -> Result<f64, String> {
    req.get(k).and_then(|v| v.as_f64()).ok_or_else(|| format!("missing number '{}'", k))
}
fn fo(req: &Value, k: &str) -> Option<f64> {
    req.get(k).and_then(|v| v.as_f64())
}
fn s<'a>(req: &'a Value, k: &str) -> Result<&'a str, String> {
    req.get(k).and_then(|v| v.as_str()).ok_or_else(|| format!("missing string '{}'", k))
}
fn b(req: &Value, k: &str) -> bool {
    req.get(k).and_then(|v| v.as_bool()).unwrap_or(false)
}
fn farr(req: &Value, k: &str) -> Result<Vec<f64>, String> {
    let a = req.get(k).and_then(|v| v.as_array()).ok_or_else(|| format!("missing array '{}'", k))?;
    Ok(a.iter().map(|v| v.as_f64().unwrap_or(f64::NAN)).collect())
}
fn filters(req: &Value) -> Vec<(String, f64)> {
    req.get("filters")
        .and_then(|v| v.as_array())
        .map(|a| {
            a.iter()
                .filter_map(|p| {
                    let p = p.as_array()?;
                    Some((p.first()?.as_str()?.to_string(), p.get(1)?.as_f64()?))
                })
                .filter(|(_, mm)| *mm > 0.0)
                .collect()
        })
        .unwrap_or_default()
}

fn filters_json(sp: &Spectrum) -> Value {
    Value::Array(sp.filters.iter().map(|(m, t, r)| json!([m, t, r])).collect())
}

fn spec_json(sp: &Spectrum) -> Value {
    json!({"E": arr(&sp.e), "N": arr(&sp.n), "kvp": sp.kvp, "model": sp.model, "filters": filters_json(sp),
           "absolute": sp.absolute, "kerma_uGy_per_mAs": opt(sp.kerma_uGy_per_mAs)})
}

fn spec_from(req: &Value, key: &str) -> Result<Spectrum, String> {
    let o = req.get(key).ok_or_else(|| format!("missing spectrum '{}'", key))?;
    let mut sp = Spectrum::new(farr(o, "E")?, farr(o, "N")?, f(o, "kvp")?, o.get("model").and_then(|v| v.as_str()).unwrap_or(""));
    sp.kerma_uGy_per_mAs = fo(o, "kerma_uGy_per_mAs");
    sp.absolute = b(o, "absolute");
    Ok(sp)
}

/// unfiltered and filtered spectrum from the usual tube inputs
fn tube(xs: &CrossSections, req: &Value) -> Result<(Spectrum, Spectrum), String> {
    let kvp = f(req, "kvp")?;
    let ang = fo(req, "angle").unwrap_or(12.0);
    let use_spek = req.get("engine").and_then(|v| v.as_str()).unwrap_or("spekpy") == "spekpy";
    let unf = if use_spek { spm::spectrum_spekpy(kvp, ang)? } else { spm::spectrum_kramers(xs, kvp, &[], ang)? };
    let mut sp = unf.clone();
    for (m, mm) in filters(req) {
        sp = sp.filtered(xs, &m, mm, None)?;
    }
    Ok((unf, sp))
}

pub fn dispatch(st: &mut State, req: &Value, blob: Option<Vec<u8>>) -> Result<Value, String> {
    let cmd = s(req, "cmd")?;
    st.bin.clear();
    if cmd == "load_xs" {
        let data = blob.ok_or("load_xs needs the file bytes")?;
        let xs = CrossSections::from_npz(&data)?;
        let mats: Vec<Value> = xs
            .mats
            .iter()
            .map(|(k, m)| json!({"key": k, "density": m.density, "ff": m.ffcum.is_some(), "fl": m.fluor.is_some()}))
            .collect();
        let r = json!({"materials": mats, "emin": xs.emin, "emax": xs.emax, "n_energy": xs.e.len(),
                       "spek_range": spm::spek_range().map(|t| json!([t.0, t.1, t.2, t.3]))});
        st.xs = Some(xs);
        return Ok(r);
    }
    if cmd == "dicom_load" {
        let data = blob.ok_or("dicom_load needs the file bytes")?;
        let slot = s(req, "slot")?.to_string();
        let d = read_dicom(&data)?;
        let mut mn = f32::INFINITY;
        let mut mx = f32::NEG_INFINITY;
        for &v in &d.image.data {
            if v < mn {
                mn = v;
            }
            if v > mx {
                mx = v;
            }
        }
        // robust window from a subsample
        let stride = (d.image.data.len() / 200_000).max(1);
        let sub: Vec<f64> = d.image.data.iter().step_by(stride).map(|&v| v as f64).collect();
        let p = crate::num::percentiles(&sub, &[1.0, 99.0]);
        let r = json!({"w": d.image.w, "h": d.image.h, "px_mm": d.px_mm, "sign": d.sign, "photometric": d.photometric,
                       "vmin": mn, "vmax": mx, "p1": p[0], "p99": p[1],
                       "info": d.info.iter().map(|(k, v)| json!([k, v])).collect::<Vec<_>>()});
        st.images.insert(slot, d);
        return Ok(r);
    }
    if cmd == "image_view" {
        let d = st.images.get(s(req, "slot")?).ok_or("image not loaded")?;
        let (x0, y0, x1, y1) = (f(req, "x0")?, f(req, "y0")?, f(req, "x1")?, f(req, "y1")?);
        let (ow, oh) = (f(req, "out_w")? as usize, f(req, "out_h")? as usize);
        let (vmin, vmax) = (f(req, "vmin")? as f32, f(req, "vmax")? as f32);
        let inv = b(req, "invert");
        let span = if vmax > vmin { vmax - vmin } else { 1.0 };
        let img = &d.image;
        let mut out = vec![0u8; ow * oh];
        for j in 0..oh {
            let y = y0 + (j as f64 + 0.5) * (y1 - y0) / oh as f64;
            for i in 0..ow {
                let x = x0 + (i as f64 + 0.5) * (x1 - x0) / ow as f64;
                let v = if x >= 0.0 && y >= 0.0 && (x as usize) < img.w && (y as usize) < img.h {
                    let t = ((img.at(y as usize, x as usize) - vmin) / span).clamp(0.0, 1.0);
                    let t = if inv { 1.0 - t } else { t };
                    (t * 255.0 + 0.5) as u8
                } else {
                    40
                };
                out[j * ow + i] = v;
            }
        }
        st.bin = out;
        return Ok(json!({"bin": st.bin.len()}));
    }
    if cmd == "mtf" {
        let d = st.images.get(s(req, "slot")?).ok_or("Load the edge image first.")?;
        let roi = farr(req, "roi")?;
        let crop = d.image.crop(roi[0] as i64, roi[1] as i64, roi[2] as i64, roi[3] as i64);
        let ov = fo(req, "oversample").unwrap_or(8.0) as usize;
        let r = iq::slanted_edge_mtf(&crop, d.px_mm, ov, b(req, "horizontal"))?;
        return Ok(json!({"freq": arr(&r.freq), "mtf": arr(&r.mtf), "angle_deg": r.angle_deg, "esf": arr(&r.esf),
                         "lsf": arr(&r.lsf), "bin_mm": r.bin_mm, "px_mm": d.px_mm}));
    }
    if cmd == "nnps" {
        let d = st.images.get(s(req, "slot")?).ok_or("Load the flat-field image first.")?;
        let roi = fo(req, "roi_size").unwrap_or(256.0) as usize;
        let step = fo(req, "step").unwrap_or(128.0) as usize;
        let mask = fo(req, "mask_below").unwrap_or(1.0);
        let prog = &mut st.progress;
        let r = iq::flat_field_nnps(&d.image, d.px_mm, roi, step, mask, |x| prog(x))?;
        let ph = iq::nnps_profile(&r.nnps2d, r.n, true, 3);
        let pv = iq::nnps_profile(&r.nnps2d, r.n, false, 3);
        return Ok(json!({"freq": arr(&r.freq), "profH": arr(&ph), "profV": arr(&pv), "mean_sig": r.mean_signal,
                         "n_used": r.count, "px_mm": d.px_mm}));
    }
    if cmd == "dqe" {
        let (q, k) = (f(req, "q")?, f(req, "K")?);
        let (fh, dh) = iq::compute_dqe(&farr(req, "fmh")?, &farr(req, "mtfh")?, &farr(req, "fH")?, &farr(req, "profH")?, q, k);
        let (fv, dv) = iq::compute_dqe(&farr(req, "fmv")?, &farr(req, "mtfv")?, &farr(req, "fV")?, &farr(req, "profV")?, q, k);
        return Ok(json!({"fDh": arr(&fh), "dqe_h": arr(&dh), "fDv": arr(&fv), "dqe_v": arr(&dv)}));
    }
    if cmd == "tables" {
        return Ok(json!({
            "rqa": spm::RQA_ED2.iter().map(|r| json!([r.0, r.1, r.2, r.3, r.4])).collect::<Vec<_>>(),
            "rqa_edition": spm::RQA_EDITION,
            "snr_in2": spm::RQA_SNR_IN2.iter().map(|r| json!([r.0, r.1])).collect::<Vec<_>>(),
            "isotopes": spm::ISOTOPE_LINES.iter().map(|r| json!([r.0, r.1.iter().map(|l| json!([l.0, l.1])).collect::<Vec<_>>()])).collect::<Vec<_>>(),
            "spek_range": spm::spek_range().map(|t| json!([t.0, t.1, t.2, t.3])),
        }));
    }

    let xs = st.xs.as_ref().ok_or("xsdata.npz has not been loaded.")?;
    match cmd {
        "spectrum" => {
            let (unf, sp) = tube(xs, req)?;
            let hm = req.get("hvl_mat").and_then(|v| v.as_str()).unwrap_or("aluminium");
            let cal = fo(req, "cal_scale");
            let h1 = sp.hvl(xs, hm, 0.5, None, true)?;
            let h2 = sp.hvl(xs, hm, 0.25, None, true)?;
            Ok(json!({
                "unf": spec_json(&unf), "sp": spec_json(&sp),
                "mean": sp.mean_energy(), "eff": sp.effective_energy(xs, hm)?,
                "hvl1": num(h1), "hvl2": num(h2 - h1), "homog": sp.homogeneity(xs, hm)?,
                "fluence": sp.fluence(), "efluence": sp.energy_fluence(),
                "unf_fluence": unf.fluence(), "unf_mean": unf.mean_energy(),
                "k0": unf.kerma_integral(), "k1": sp.kerma_integral(),
                "output": opt(spm::spectrum_output(&sp, cal)),
                "archer_K": spm::Spectrum::archer_K_poly(sp.kvp),
            }))
        }
        "solve_hvl" => {
            let use_spek = s(req, "engine")? == "spekpy";
            let (mm, ach) = spm::solve_filtration_from_hvl(xs, f(req, "kvp")?, f(req, "hvl")?, fo(req, "angle").unwrap_or(12.0), use_spek, 40.0)?;
            Ok(json!({"mm": opt(mm), "achieved": num(ach)}))
        }
        "set_calibration" => {
            let measured = f(req, "measured")?;
            let dist = f(req, "distance")?;
            if measured <= 0.0 {
                return Err("Measured output must be greater than zero.".into());
            }
            if dist <= 0.0 {
                return Err("Distance must be greater than zero.".into());
            }
            let use_spek = s(req, "engine")? == "spekpy";
            let sp = spm::make_spectrum(xs, f(req, "kvp")?, &filters(req), fo(req, "angle").unwrap_or(12.0), use_spek)?;
            let at_1m = measured * dist * dist;
            let integral = sp.kerma_integral();
            if integral <= 0.0 {
                return Err("Spectrum carries no kerma - cannot calibrate.".into());
            }
            let spek_ratio = sp.kerma_uGy_per_mAs.filter(|k| *k != 0.0).map(|k| at_1m / k);
            Ok(json!({"scale": at_1m / integral, "spek_ratio": opt(spek_ratio), "kvp": sp.kvp, "model": sp.model}))
        }
        "sweep" => {
            let ang = fo(req, "angle").unwrap_or(12.0);
            let use_spek = s(req, "engine")? == "spekpy";
            let cal = fo(req, "cal_scale");
            let fl = filters(req);
            let mut rows = vec![];
            for kv in (50..155).step_by(5) {
                let kv = kv as f64;
                let base = if use_spek { spm::spectrum_spekpy(kv, ang)? } else { spm::spectrum_kramers(xs, kv, &[], ang)? };
                let mut s2 = base;
                for (m, mm) in &fl {
                    s2 = s2.filtered(xs, m, *mm, None)?;
                }
                let k = spm::spectrum_output(&s2, cal);
                rows.push(json!({"kVp": kv, "K_app": opt(k.map(|v| v * 0.06)), "K_archer": Spectrum::archer_K_poly(kv),
                                 "HVL": num(s2.hvl_al(xs)?), "mean": s2.mean_energy()}));
            }
            Ok(json!({"rows": rows}))
        }
        "rqa_identify" => {
            let (_unf, sp) = tube(xs, req)?;
            let cal = fo(req, "cal_scale");
            let hvl = sp.hvl_al(xs)?;
            let mut rows: Vec<(f64, Value)> = spm::RQA_ED2
                .iter()
                .map(|r| {
                    let dev = 100.0 * (hvl - r.3) / r.3;
                    (dev, json!({"quality": r.0, "kV": r.1, "nominal_added_Al_mm": r.2, "target_HVL_mm": r.3,
                                 "your_HVL_mm": (hvl * 1000.0).round() / 1000.0, "deviation_pc": (dev * 100.0).round() / 100.0,
                                 "confidence": r.4}))
                })
                .collect();
            rows.sort_by(|a, b| {
                let da = a.1["deviation_pc"].as_f64().unwrap_or(0.0).abs();
                let db = b.1["deviation_pc"].as_f64().unwrap_or(0.0).abs();
                da.partial_cmp(&db).unwrap()
            });
            let q = spm::rqa_fluence_per_kerma(&sp, cal);
            let mut frows = vec![];
            if let Some(qv) = q {
                let mut fr: Vec<Value> = spm::RQA_SNR_IN2
                    .iter()
                    .map(|(n, t)| {
                        json!({"quality": n, "target_SNRin2": t, "your_SNRin2": (qv * 10.0).round() / 10.0,
                               "deviation_pc": (100.0 * (qv - t) / t * 100.0).round() / 100.0})
                    })
                    .collect();
                fr.sort_by(|a, b| {
                    a["deviation_pc"].as_f64().unwrap().abs().partial_cmp(&b["deviation_pc"].as_f64().unwrap().abs()).unwrap()
                });
                frows = fr;
            }
            let _ = &mut rows;
            Ok(json!({"rows": rows.into_iter().map(|r| r.1).collect::<Vec<_>>(), "hvl": num(hvl), "mean": sp.mean_energy(),
                      "kvp": sp.kvp, "filters": filters_json(&sp), "q": opt(q), "frows": frows}))
        }
        "rqa_solve" => {
            let (_unf, sp) = tube(xs, req)?;
            let q = s(req, "quality")?;
            let (kv, nominal, _target, conf) = spm::rqa_lookup(q)?;
            let (mm, ach, tgt) = spm::rqa_solve_filtration(xs, &sp, q, 80.0)?;
            let (mmc, achc, _) = spm::rqa_solve_filtration_cu(xs, &sp, q, 0.5, 80.0)?;
            let hard = match mm {
                Some(m) => {
                    let h = sp.filtered(xs, "aluminium", m, None)?;
                    json!({"E": arr(&h.e), "N": arr(&h.n), "mean": h.mean_energy()})
                }
                None => Value::Null,
            };
            let hard_cu = match mmc {
                Some(m) => {
                    let h = sp.filtered(xs, "copper", 0.5, None)?.filtered(xs, "aluminium", m, None)?;
                    json!({"E": arr(&h.e), "N": arr(&h.n), "mean": h.mean_energy()})
                }
                None => Value::Null,
            };
            Ok(json!({"kv": kv, "nominal": nominal, "target": tgt, "conf": conf, "mm": opt(mm), "achieved": num(ach),
                      "mm_cu": opt(mmc), "achieved_cu": num(achc), "start_hvl": num(sp.hvl_al(xs)?),
                      "sp": spec_json(&sp), "hard": hard, "hard_cu": hard_cu}))
        }
        "isotope" => {
            let name = s(req, "name")?;
            let sp = spm::make_isotope_spectrum(name)?;
            let beyond = spm::isotope_lines_beyond_xsdata(xs, name)?;
            let mut v = spec_json(&sp);
            v["mean"] = json!(sp.mean_energy());
            v["beyond"] = arr(&beyond);
            Ok(v)
        }
        "spec_info" => {
            let sp = spec_from(req, "sp")?;
            Ok(json!({"mean": sp.mean_energy()}))
        }
        "mu_curves" => {
            let mat = s(req, "mat")?;
            let dens = fo(req, "density");
            let (e0, e1, n) = (f(req, "emin")?, f(req, "emax")?, f(req, "n")? as usize);
            let e: Vec<f64> = (0..n).map(|i| if n > 1 { e0 + (e1 - e0) * i as f64 / (n as f64 - 1.0) } else { e0 }).collect();
            Ok(json!({"E": arr(&e), "total": arr(&xs.mu(mat, &e, Kind::Total, dens)?), "photo": arr(&xs.mu(mat, &e, Kind::Photo, dens)?),
                      "coh": arr(&xs.mu(mat, &e, Kind::Coh, dens)?), "incoh": arr(&xs.mu(mat, &e, Kind::Incoh, dens)?)}))
        }
        "density" => Ok(json!({"density": xs.density(s(req, "mat")?)?})),
        "archer" => {
            let mat = s(req, "mat")?;
            let xs_mm = farr(req, "x")?;
            let iso = req.get("iso").and_then(|v| v.as_str());
            let (bvals, has_fit, info, alpha, mu_top): (Vec<Option<f64>>, bool, Option<(&str, &str)>, Option<f64>, Option<f64>) = match iso {
                Some(name) => {
                    let p = archer::nuclide_params(name, mat);
                    let mean_e = fo(req, "mean_energy").unwrap_or(0.0);
                    (
                        xs_mm.iter().map(|&x| archer::archer_reference_nuclide(x, name, mat)).collect(),
                        p.is_some(),
                        archer::nuclide_info(name, mat),
                        p.map(|v| v.0),
                        if p.is_some() { Some(xs.mu1(mat, mean_e, Kind::Total, None)? / 10.0) } else { None },
                    )
                }
                None => {
                    let kvp = f(req, "kvp")?;
                    let a = archer::archer_alpha_near(kvp, mat);
                    (
                        xs_mm.iter().map(|&x| archer::archer_reference(x, kvp, mat)).collect(),
                        archer::table(mat).is_some(),
                        None,
                        a,
                        if a.is_some() { Some(xs.mu1(mat, kvp, Kind::Total, None)? / 10.0) } else { None },
                    )
                }
            };
            Ok(json!({"B": bvals.iter().map(|v| opt(*v)).collect::<Vec<_>>(), "has_fit": has_fit,
                      "basis": info.map(|i| i.0), "unit": info.map(|i| i.1), "alpha": opt(alpha), "mu_top": opt(mu_top)}))
        }
        "mc_start" => {
            let sp = spec_from(req, "sp")?;
            let p = McParams {
                material: s(req, "material")?.to_string(),
                thickness_mm: f(req, "thickness")?,
                n_hist: f(req, "n_hist")? as u64,
                density: fo(req, "density").filter(|d| *d > 0.0),
                seed: fo(req, "seed").unwrap_or(1234.0) as u32,
                batch: fo(req, "batch").unwrap_or(20000.0) as u64,
                score: if req.get("score").and_then(|v| v.as_str()) == Some("fluence") { Score::Fluence } else { Score::Kerma },
                acceptance_deg: fo(req, "acceptance").unwrap_or(90.0),
                e_threshold: fo(req, "e_threshold").unwrap_or(0.0),
                fluorescence: req.get("fluorescence").and_then(|v| v.as_bool()).unwrap_or(true),
            };
            st.mc = Some(McRun::new(xs, &sp, p)?);
            Ok(json!({"ok": true}))
        }
        "mc_step" => {
            let r = st.mc.as_mut().ok_or("no Monte Carlo run in progress")?;
            let n = fo(req, "batches").unwrap_or(1.0).max(1.0) as usize;
            for _ in 0..n {
                r.step();
            }
            Ok(json!({"done": r.done, "n_hist": r.p.n_hist, "finished": r.finished()}))
        }
        "mc_result" => {
            let r = st.mc.as_ref().ok_or("no Monte Carlo run in progress")?;
            let res = r.result(xs)?;
            Ok(json!({"B_mc": num(res.b_mc), "B_narrow": num(res.b_narrow), "buildup": num(res.buildup),
                      "rel_std_err": num(res.rel_std_err), "backscatter": num(res.backscatter),
                      "mean_collisions": num(res.mean_collisions), "hist_E": arr(&res.hist_e), "hist_w": arr(&res.hist_w),
                      "hist_norm": arr(&res.hist_norm), "thickness_mm": r.p.thickness_mm, "n_hist": r.p.n_hist}))
        }
        _ => Err(format!("unknown command '{}'", cmd)),
    }
}
