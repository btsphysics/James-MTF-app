#!/usr/bin/env python3
"""
Pre-compute SpekPy spectra for the browser build.

SpekPy (~200 MB of physics tables) cannot run inside a web page, so the
Rust/WebAssembly port carries SpekPy's *output* instead: the unfiltered
tungsten spectrum (dk = 1 keV, mas = 1, SpekPy's default reference distance)
and SpekPy's own air-kerma figure, for every integer kVp 20-150 and every
integer anode angle 5-25 deg. Filtration is then applied in the app by
Beer-Lambert, exactly as app.py's run() does.

app.py's calibration path (make_spectrum with SpekPy) lets SpekPy apply the
filters itself, with SpekPy's own attenuation data, and reads SpekPy's own
air kerma. To reproduce that, a second file carries SpekPy's mu*t per mm
for the materials app.py maps to SpekPy names, SpekPy's air mu_en/rho, and
the kerma constant, all on SpekPy's 1 keV bin centres (1.5 .. 149.5 keV):

Output 2: rust/data/spekpy_mu.bin  (zlib-compressed, little-endian)
    u32 magic 'SPM1', u32 n_bins, f64 E_first, f64 kerma_const
    f64 muen_air[n_bins]
    u32 n_materials, then per material: u8 name_len, name (app key),
        f64 mu_per_mm[n_bins]

Pass --mu-only to rebuild just the second file.

Output: rust/data/spekpy_table.bin  (zlib-compressed, little-endian)
    u32 magic 'SPK1', u32 kv_min, u32 kv_max, u32 th_min, u32 th_max
    then for each angle (th_min..th_max), for each kVp (kv_min..kv_max):
        f64 kerma_uGy_per_mAs
        u32 n_bins
        f64 E_first_centre
        f32 N[n_bins]            photons / keV / cm2 / mAs
"""
import os, struct, sys, zlib
import numpy as np
import spekpy as sp

KV_MIN, KV_MAX = 20, 150
TH_MIN, TH_MAX = 5, 25

here = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(here, "..", "rust", "data")
os.makedirs(data_dir, exist_ok=True)

# ---- SpekPy filtration / kerma data (app.py's name_map) -------------------
NAME_MAP = {"aluminium": "Al", "copper": "Cu", "beryllium": "Be",
            "water": "Water", "air": "Air", "tissue": "Water",
            "lead": "Pb", "tungsten": "W"}
s = sp.Spek(kvp=KV_MAX, th=12, dk=1.0, mas=1.0)
k = s.model.k
E, N = s.get_spectrum(edges=False)
muen = s.muen_air_data.get_muen_over_rho_air(k)
C = float(np.asarray(s.get_kerma()).ravel()[0]) / float(np.sum(N * E * muen))
mu_out = bytearray(struct.pack("<4sIdd", b"SPM1", len(k), float(k[0]), C))
mu_out += np.asarray(muen, dtype="<f8").tobytes()
mu_out += struct.pack("<I", len(NAME_MAP))
for key, name in NAME_MAP.items():
    mu_out += struct.pack("<B", len(key)) + key.encode()
    mu_out += np.asarray(s.mu_data.get_mu_t(name, k, 1.0), dtype="<f8").tobytes()
with open(os.path.join(data_dir, "spekpy_mu.bin"), "wb") as f:
    f.write(zlib.compress(bytes(mu_out), 9))
print("wrote spekpy_mu.bin", file=sys.stderr)
if "--mu-only" in sys.argv:
    sys.exit(0)

out = bytearray(struct.pack("<4sIIII", b"SPK1", KV_MIN, KV_MAX, TH_MIN, TH_MAX))
for th in range(TH_MIN, TH_MAX + 1):
    for kv in range(KV_MIN, KV_MAX + 1):
        s = sp.Spek(kvp=kv, th=th, dk=1.0, mas=1.0)
        E, N = s.get_spectrum(edges=False)
        k = float(np.asarray(s.get_kerma()).ravel()[0])
        out += struct.pack("<dId", k, len(E), float(E[0]))
        out += np.asarray(N, dtype="<f4").tobytes()
    print("angle", th, "done", flush=True)

dst = os.path.join(data_dir, "spekpy_table.bin")
with open(dst, "wb") as f:
    f.write(zlib.compress(bytes(out), 9))
print("wrote", dst, len(out), "raw bytes", file=sys.stderr)
