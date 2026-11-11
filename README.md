# X-ray Spectrum and Monte Carlo Transmission

`app.py` is the original Tk desktop app, which ran on an embedded Python with
numpy, SpekPy, numba and pydicom. This repository adds a **Rust rewrite** of
it that runs in a web browser.

## Using it

Open **`app.html`** in Chrome, Edge or Firefox by double-clicking it. You
don't need to install anything, run a server, or have a network connection.
The page is a single 4.9 MB file. It contains:

* the Rust code compiled to WebAssembly,
* `xsdata.npz`, the cross-section data, built in. You can load a different
  copy with **Load xsdata.npz** (top right),
* SpekPy's spectra, pre-computed (see below).

It has the same four tabs as the Python app: **Spectrum**, **Monte Carlo**,
**MTF / NNPS / DQE** and **Notes**. Everything is calculated on your own
computer and no image or data leaves the machine.

Differences from the Tk app that you will notice:

* **Monte Carlo is parallel.** Each thickness point runs on its own Web
  Worker. 1.4 million histories take about 2 s on a 4-core laptop.
* **Edge-image viewer.** Left-drag draws an ROI. The mouse wheel zooms in at
  full resolution. Right-drag or Shift+drag pans. "Contrast" switches
  between the full range and 1-99 %.
* **Export** saves CSV files. The Python app could also save `.xlsx`.
* **Plots** have a "Save PNG" button and show the coordinates under the
  mouse.

## Layout

| Path | What it is |
|---|---|
| `app.py` | Original Python/Tk application (unchanged) |
| `app.html` | **Built browser app** |
| `rust/` | The Rust crate: all physics and image analysis |
| `rust/src/xs.rs` | Cross sections (log-log interpolation), NIST air mu_en, air kerma |
| `rust/src/spectrum.rs` | Spectra (SpekPy table, built-in Kramers), filtration, HVL, effective energy, calibration, RQA, isotopes |
| `rust/src/mc.rs` | Monte Carlo slab transport (forced collisions, Klein-Nishina, form-factor Rayleigh, K fluorescence) |
| `rust/src/archer.rs` | Archer / NCRP 147 and radionuclide overlay fits |
| `rust/src/iq.rs` | Slanted-edge MTF, flat-field NNPS, DQE |
| `rust/src/dicom.rs` | DICOM reader (replaces pydicom) |
| `rust/src/npz.rs` | `.npz` / `.npy` reader |
| `rust/src/api.rs` | JSON command interface used by the page |
| `rust/data/` | Pre-computed SpekPy tables |
| `web/` | Front-end source (HTML/CSS/JS: UI, canvas plotting, worker pool) |
| `tools/build_html.py` | Builds `app.html` |
| `tools/build_spekpy_table.py` | Regenerates `rust/data/` from SpekPy |
| `tools/validate_against_python.py` | Checks the Rust port against the original `app.py` functions |

## Rebuilding

```bash
rustup target add wasm32-unknown-unknown
```

```bash
python3 tools/build_html.py
```

This builds the WebAssembly module and writes `app.html`. To check the port
against the Python original, which needs numpy, spekpy and pydicom:

```bash
cargo build --release --manifest-path rust/Cargo.toml
```

```bash
python3 tools/validate_against_python.py
```

## How faithful it is

`tools/validate_against_python.py` runs the functions in `app.py` itself,
with the GUI removed, side by side with the Rust code on the same inputs.
Current results:

* **Built-in spectra, HVL 1/2, effective energy, homogeneity, kerma,
  filtration-from-HVL solve, RQA identify/solve:** identical to rounding
  (relative difference about 1e-14).
* **SpekPy spectra, SpekPy air kerma, kVp sweep, calibration scale factor
  and "correction to SpekPy":** match to about 1e-8. The SpekPy table is
  stored as float32.
* **Monte Carlo:** bit-identical. The Rust code uses numpy's MT19937
  generator with the same `np.random.seed(1234)` seeding, so it uses exactly
  the same random numbers as `app.py`'s numpy path. B_mc, the standard
  error, backscatter, collision counts and the transmitted spectrum all
  match exactly, for tube spectra, isotopes, kerma and fluence scoring, and
  with a narrowed acceptance angle and an energy threshold.
  * When numba was installed, `app.py`'s kernel used numba's own unseeded
    generator. Those results are only statistically equivalent, not
    identical.
* **MTF, NNPS, DQE:** match to about 1e-5 or better. pydicom works in
  float32, while the Rust code works in float64.
* **DICOM:** these files all decode to the same pixels:
  * implicit VR little endian,
  * explicit VR little endian,
  * explicit VR big endian,
  * RLE Lossless,
  * rescale slope/intercept applied,
  * a decoy PixelSpacing nested inside a sequence (correctly ignored).

### SpekPy

SpekPy cannot run in a browser: its physics tables are about 200 MB. The
page therefore carries SpekPy 2.5's own output instead:

* the unfiltered spectrum (dk = 1 keV, 1 mAs) and SpekPy's air kerma,
* for every whole kVp from **20 to 150** and every whole anode angle from
  **5 to 25 degrees**.

How that differs from running SpekPy live:

* A kVp that isn't a whole number is run at the nearest whole kV, and the
  report says so.
* An anode angle that isn't a whole number is interpolated between whole
  degrees.
* In some places `app.py` let SpekPy do the filtering: when solving
  filtration from a measured HVL, and for the spectrum the output
  calibration is taken on. There, SpekPy's own attenuation data and air
  kerma are reproduced for Al, Cu, Be, water, tissue, air, Pb and W.
  * For other materials, `app.py` silently skipped the filter in that path.
    The Rust code applies them using the Elam data instead.

### DICOM

The reader handles:

* uncompressed images (implicit or explicit VR, little or big endian),
* deflated files,
* RLE Lossless,
* JPEG baseline and JPEG lossless.

JPEG-LS and JPEG 2000 files are refused with a message asking for an
uncompressed export. Multi-frame files use the first frame.
