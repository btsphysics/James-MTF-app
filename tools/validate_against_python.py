#!/usr/bin/env python3
"""
Cross-check the Rust port against the original app.py functions.

    cargo build --release --manifest-path rust/Cargo.toml
    python3 tools/validate_against_python.py [workdir]

Needs numpy, spekpy and pydicom for the full set.
"""
import json, os, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pyref import ref

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
CLI = os.path.join(ROOT, "rust", "target", "release", "xray_cli")
WORK = sys.argv[1] if len(sys.argv) > 1 else "/tmp/xray_validate"
os.makedirs(WORK, exist_ok=True)

reqs = []
def q(**kw):
    reqs.append(kw)
    return len(reqs) - 1

q(cmd="load_xs", file=os.path.join(ROOT, "xsdata.npz"))
results = []
fails = 0

def check(name, py, rs, tol):
    global fails
    py = np.atleast_1d(np.asarray(py, dtype=float))
    rs = np.atleast_1d(np.asarray([np.nan if v is None else (np.inf if v == "inf" else v) for v in np.atleast_1d(rs)], dtype=float))
    if py.shape != rs.shape:
        print("FAIL %-48s shape %s vs %s" % (name, py.shape, rs.shape)); fails += 1; return
    both_nan = np.isnan(py) & np.isnan(rs)
    den = np.maximum(np.abs(py), 1e-300)
    rel = np.where(both_nan, 0.0, np.abs(py - rs) / den)
    worst = float(np.nanmax(rel)) if rel.size else 0.0
    ok = worst <= tol
    fails += (not ok)
    print("%s %-48s max rel diff %.2e (tol %.0e)" % ("ok  " if ok else "FAIL", name, worst, tol))

# ---------------- spectra ---------------------------------------------------
cases = [("builtin", 80, 12, [("aluminium", 2.5)], "aluminium"),
         ("builtin", 125, 16, [("aluminium", 2.5), ("copper", 0.1)], "copper"),
         ("builtin", 72.5, 10, [("aluminium", 3.0), ("pmma", 20)], "aluminium"),
         ("spekpy", 80, 12, [("aluminium", 2.5)], "aluminium"),
         ("spekpy", 120, 14, [("aluminium", 3.0), ("copper", 0.2)], "aluminium")]
spec_ids = [q(cmd="spectrum", engine=e, kvp=k, angle=a, filters=[list(f) for f in fl], hvl_mat=hm)
            for e, k, a, fl, hm in cases]
solve_ids = [q(cmd="solve_hvl", engine="builtin", kvp=80, hvl=3.5, angle=12),
             q(cmd="solve_hvl", engine="spekpy", kvp=70, hvl=2.8, angle=12)]
cal_id = q(cmd="set_calibration", engine="spekpy", kvp=80, angle=12, filters=[["aluminium", 3.1], ["copper", 0.1]],
           measured=40.0, distance=1.2)
sweep_id = q(cmd="sweep", engine="spekpy", angle=12, filters=[["aluminium", 2.5]])
rqa_id = q(cmd="rqa_identify", engine="spekpy", kvp=70, angle=12, filters=[["aluminium", 23.5]])
rqs_id = q(cmd="rqa_solve", engine="spekpy", kvp=70, angle=12, filters=[["aluminium", 2.5]], quality="RQA5")

# ---------------- Monte Carlo (same seed, same random stream) ---------------
mc_cases = [("builtin", 100, [("aluminium", 2.5)], "lead", 1.0, 1500, "kerma", 90.0, 0.0),
            ("builtin", 125, [("aluminium", 2.5)], "concrete", 40.0, 1500, "fluence", 30.0, 10.0),
            ("iso:Cs-137", None, [], "lead", 10.0, 1500, "kerma", 90.0, 0.0)]
mc_ids = []
for eng, kvp, fl, mat, t, nh, sc, acc, eth in mc_cases:
    if eng.startswith("iso:"):
        s = ref.make_isotope_spectrum(eng[4:])
    else:
        s = ref.spectrum_kramers(kvp, fl)
    spj = {"E": s.E.tolist(), "N": s.N.tolist(), "kvp": s.kvp, "model": s.model}
    q(cmd="mc_start", sp=spj, material=mat, thickness=t, n_hist=nh, seed=1234, batch=20000,
      score=sc, acceptance=acc, e_threshold=eth)
    q(cmd="mc_step", batches=1000)
    mc_ids.append((q(cmd="mc_result"), s, mat, t, nh, sc, acc, eth))

arch_id = q(cmd="archer", mat="lead", x=[0.25, 1.0, 2.0, 3.0], kvp=87.5)
archn_id = q(cmd="archer", mat="lead", x=[1.0, 5.0], iso="Tc-99m", mean_energy=140.5)

# ---------------- DICOM / MTF / NNPS ---------------------------------------
dicom_ok = ref.HAVE_PYDICOM
if dicom_ok:
    import pydicom
    from pydicom.dataset import Dataset, FileMetaDataset
    from pydicom.uid import ExplicitVRLittleEndian, ImplicitVRLittleEndian, ExplicitVRBigEndian, generate_uid

    def write(path, arr, ts, px=0.139, slope=None):
        meta = FileMetaDataset()
        meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.1.1"
        meta.MediaStorageSOPInstanceUID = generate_uid()
        meta.TransferSyntaxUID = ts
        ds = Dataset(); ds.file_meta = meta
        ds.SOPClassUID = meta.MediaStorageSOPClassUID; ds.SOPInstanceUID = meta.MediaStorageSOPInstanceUID
        ds.Modality = "DX"; ds.Rows, ds.Columns = arr.shape
        ds.SamplesPerPixel = 1; ds.PhotometricInterpretation = "MONOCHROME2"
        ds.BitsAllocated = 16; ds.BitsStored = 16; ds.HighBit = 15; ds.PixelRepresentation = 0
        ds.ImagerPixelSpacing = [px, px]
        if slope:
            ds.RescaleSlope = slope; ds.RescaleIntercept = 3
        seq = Dataset(); seq.PixelSpacing = [9.9, 9.9]          # nested - must be ignored
        ds.ReferencedImageSequence = [seq]
        ds.PixelData = arr.astype(np.uint16).tobytes()
        if ts == ExplicitVRBigEndian:
            ds.PixelData = arr.astype(">u2").tobytes()
        ds.save_as(path, enforce_file_format=True)

    rng = np.random.default_rng(5)
    yy, xx = np.mgrid[0:600, 0:500]
    ang = np.radians(2.3)
    d = (xx - 250) * np.cos(ang) - (yy - 300) * np.sin(ang)       # near-vertical edge
    from math import erf
    verf = np.vectorize(erf)
    edge = 1000 + 2500 * 0.5 * (1 + verf(d / 1.6)) + rng.normal(0, 8, d.shape)
    edge = np.clip(edge, 0, 65535)
    p_edge = os.path.join(WORK, "edge.dcm"); write(p_edge, edge, ExplicitVRLittleEndian)
    p_edge_i = os.path.join(WORK, "edge_implicit.dcm"); write(p_edge_i, edge.T.copy(), ImplicitVRLittleEndian, slope=1.5)
    p_edge_b = os.path.join(WORK, "edge_be.dcm")
    try:
        write(p_edge_b, edge, ExplicitVRBigEndian)
    except Exception as e:
        p_edge_b = None
    p_edge_r = os.path.join(WORK, "edge_rle.dcm")
    try:
        from pydicom.uid import RLELossless
        dsr = pydicom.dcmread(p_edge); dsr.compress(RLELossless); dsr.save_as(p_edge_r)
    except Exception as e:
        print("(RLE test skipped: %s)" % e); p_edge_r = None
    flat = 2000 + 40 * (xx[:560] / 500.0) + rng.normal(0, 20, (560, 500))
    flat[:, :30] = 0                                               # masked border
    flat = np.clip(flat, 0, 65535)
    p_flat = os.path.join(WORK, "flat.dcm"); write(p_flat, flat, ExplicitVRLittleEndian, px=0.15)
    roi_v = [100, 500, 200, 300]                                   # r0 r1 c0 c1
    q(cmd="dicom_load", slot="edge", file=p_edge)
    mtfv_id = q(cmd="mtf", slot="edge", roi=roi_v, horizontal=False)
    q(cmd="dicom_load", slot="edgei", file=p_edge_i)
    mtfh_id = q(cmd="mtf", slot="edgei", roi=[200, 300, 100, 500], horizontal=True)
    if p_edge_b:
        q(cmd="dicom_load", slot="edgeb", file=p_edge_b)
        mtfb_id = q(cmd="mtf", slot="edgeb", roi=roi_v, horizontal=False)
    if p_edge_r:
        q(cmd="dicom_load", slot="edger", file=p_edge_r)
        mtfr_id = q(cmd="mtf", slot="edger", roi=roi_v, horizontal=False)
    q(cmd="dicom_load", slot="flat", file=p_flat)
    nnps_id = q(cmd="nnps", slot="flat")

# ---------------- run the Rust side ----------------------------------------
inp = "\n".join(json.dumps(r) for r in reqs) + "\n"
out = subprocess.run([CLI], input=inp, capture_output=True, text=True, check=True).stdout.splitlines()
R = [json.loads(l) for l in out]
for i, r in enumerate(R):
    if "err" in r:
        print("Rust error on request", i, reqs[i]["cmd"], r["err"]); fails += 1
R = [r.get("ok") for r in R]

# ---------------- compare ---------------------------------------------------
print("\n-- spectra --")
for (eng, kvp, ang, fl, hm), i in zip(cases, spec_ids):
    r = R[i]
    unf = ref.spectrum_spekpy(kvp, ang, ()) if eng == "spekpy" else ref.spectrum_kramers(kvp, (), target_angle_deg=ang)
    s = unf
    for m, mm in fl:
        s = s.filtered(m, mm)
    tag = "%s %s kVp" % (eng, kvp)
    tol = 1e-6 if eng == "spekpy" else 1e-11      # spekpy table is stored as float32
    check(tag + " N(E)", s.N, r["sp"]["N"], tol)
    check(tag + " mean energy", s.mean_energy, r["mean"], tol)
    check(tag + " HVL1 in " + hm, s.hvl(hm, 0.5), r["hvl1"], 1e-6)
    check(tag + " HVL2", s.hvl(hm, 0.25) - s.hvl(hm, 0.5), r["hvl2"], 1e-5)
    check(tag + " effective energy", s.effective_energy(hm), r["eff"], 1e-9)
    check(tag + " homogeneity", s.homogeneity(hm), r["homog"], 1e-5)
    check(tag + " kerma integral", s.kerma_integral(), r["k1"], tol)
    if eng == "spekpy":
        check(tag + " SpekPy air kerma uGy/mAs", s.kerma_uGy_per_mAs, r["output"], 1e-6)

print("\n-- filtration solve / sweep / RQA --")
mm, ach = ref.solve_filtration_from_hvl(80, 3.5, 12, prefer_spekpy=False)
check("solve HVL builtin: total Al", mm, R[solve_ids[0]]["mm"], 1e-9)
mm, ach = ref.solve_filtration_from_hvl(70, 2.8, 12, prefer_spekpy=True)
check("solve HVL SpekPy: total Al", mm, R[solve_ids[1]]["mm"], 1e-6)
msp = ref.make_spectrum(80, [("aluminium", 3.1), ("copper", 0.1)], 12, True)
ref.set_calibration(msp, 40.0, 1.2)
check("calibration scale factor", ref.CALIBRATION["scale"], R[cal_id]["scale"], 1e-6)
check("calibration correction to SpekPy", ref.CALIBRATION["spek_ratio"], R[cal_id]["spek_ratio"], 1e-6)
ref.clear_calibration()
rows = R[sweep_id]["rows"]
pyK, pyH = [], []
for kv in range(50, 155, 5):
    s2 = ref.spectrum_spekpy(kv, 12, ()).filtered("aluminium", 2.5)
    pyK.append(ref.spectrum_output(s2) * 0.06); pyH.append(s2.hvl("aluminium"))
check("sweep K (mGy/mA.min)", pyK, [r["K_app"] for r in rows], 1e-6)
check("sweep HVL", pyH, [r["HVL"] for r in rows], 1e-6)
s70 = ref.spectrum_spekpy(70, 12, ()).filtered("aluminium", 23.5)
q_py, fr = ref.rqa_fluence_match(s70)
check("RQA identify: Q/K_air", q_py, R[rqa_id]["q"], 1e-6)
check("RQA identify: HVL", s70.hvl("aluminium"), R[rqa_id]["hvl"], 1e-6)
print("     nearest quality: python %s / rust %s" % (ref.rqa_match(s70)[0]["quality"], R[rqa_id]["rows"][0]["quality"]))
b70 = ref.spectrum_spekpy(70, 12, ()).filtered("aluminium", 2.5)
mm, a, t = ref.rqa_solve_filtration(b70, "RQA5")
mmc, ac, _ = ref.rqa_solve_filtration_cu(b70, "RQA5")
check("RQA5 solve: added Al", mm, R[rqs_id]["mm"], 1e-6)
check("RQA5 solve: Al on 0.5 Cu", mmc, R[rqs_id]["mm_cu"], 1e-6)

print("\n-- Monte Carlo (numpy MT19937 stream reproduced) --")
for i, s, mat, t, nh, sc, acc, eth in mc_ids:
    pr = ref.run_monte_carlo(s, mat, t, n_hist=nh, seed=1234, score=sc, acceptance_deg=acc, e_threshold=eth)
    r = R[i]
    tag = "%s %s mm %s" % (s.model[:14], t, mat)
    check(tag + " B_mc", pr["B_mc"], r["B_mc"], 1e-9)
    check(tag + " B_narrow", pr["B_narrow"], r["B_narrow"], 1e-12)
    check(tag + " rel std err", pr["rel_std_err"], r["rel_std_err"], 1e-8)
    check(tag + " backscatter", pr["backscatter"], r["backscatter"], 1e-9)
    check(tag + " mean collisions", pr["mean_collisions"], r["mean_collisions"], 1e-12)
    check(tag + " transmitted spectrum", pr["hist_w"], r["hist_w"], 1e-8)

print("\n-- Archer --")
check("Archer lead 87.5 kVp", [ref.archer_reference(x, 87.5, "lead") for x in [0.25, 1, 2, 3]], R[arch_id]["B"], 1e-12)
check("Archer Tc-99m lead", [ref.archer_reference_nuclide(x, "Tc-99m", "lead") for x in [1, 5]], R[archn_id]["B"], 1e-12)

if dicom_ok:
    print("\n-- DICOM / MTF / NNPS --")
    arr, px, _ = ref.load_dicom_pixels(p_edge)
    r0, r1, c0, c1 = roi_v
    pm = ref.slanted_edge_mtf(arr[r0:r1, c0:c1], px, edge_axis="vertical")
    rm = R[mtfv_id]
    check("vertical MTF edge angle", pm["angle_deg"], rm["angle_deg"], 1e-5)
    check("vertical MTF ESF", pm["esf"], rm["esf"], 1e-5)
    n = min(len(pm["mtf"]), len(rm["mtf"]))
    check("vertical MTF curve (to Nyquist)", pm["mtf"][: n // 8], rm["mtf"][: n // 8], 1e-4)
    check("vertical MTF freq", pm["freq"], rm["freq"], 1e-12)
    arr, px, _ = ref.load_dicom_pixels(p_edge_i)
    pm = ref.slanted_edge_mtf(arr[200:300, 100:500], px, edge_axis="horizontal")
    rm = R[mtfh_id]
    check("horizontal MTF (implicit VR, rescaled) angle", pm["angle_deg"], rm["angle_deg"], 1e-5)
    n = min(len(pm["mtf"]), len(rm["mtf"]))
    check("horizontal MTF curve", pm["mtf"][: n // 8], rm["mtf"][: n // 8], 1e-4)
    if p_edge_b:
        check("big-endian file gives same MTF", R[mtfv_id]["mtf"], R[mtfb_id]["mtf"], 0)
    if p_edge_r:
        check("RLE-compressed file gives same MTF", R[mtfv_id]["mtf"], R[mtfr_id]["mtf"], 0)
    arr, px, _ = ref.load_dicom_pixels(p_flat)
    nn, fr, ms, cnt = ref.flat_field_nnps(arr, px)
    fH, pH = ref.nnps_profile(nn, fr, "h"); fV, pV = ref.nnps_profile(nn, fr, "v")
    rn = R[nnps_id]
    print("     ROIs used: python %d / rust %d" % (cnt, rn["n_used"]))
    check("NNPS mean signal", ms, rn["mean_sig"], 1e-6)
    check("NNPS freq", fr, rn["freq"], 1e-12)
    check("NNPS horizontal profile", pH, rn["profH"], 1e-5)
    check("NNPS vertical profile", pV, rn["profV"], 1e-5)
    fDh, dh = ref.compute_dqe(R[mtfv_id]["freq"] and np.array(R[mtfv_id]["freq"]), np.array(R[mtfv_id]["mtf"]), fH, pH, 29653.0, 2.5)
    rq = subprocess.run([CLI], input=json.dumps(dict(cmd="dqe", q=29653.0, K=2.5, fmh=R[mtfv_id]["freq"], mtfh=R[mtfv_id]["mtf"],
                         fH=fH.tolist(), profH=pH.tolist(), fmv=R[mtfv_id]["freq"], mtfv=R[mtfv_id]["mtf"], fV=fV.tolist(),
                         profV=pV.tolist())) + "\n", capture_output=True, text=True).stdout
    check("DQE", dh, json.loads(rq)["ok"]["dqe_h"], 1e-12)

print("\n%s" % ("ALL CHECKS PASSED" if fails == 0 else "%d CHECK(S) FAILED" % fails))
sys.exit(1 if fails else 0)
