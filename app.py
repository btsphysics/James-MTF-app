#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
X-RAY SPECTRUM AND MONTE CARLO TRANSMISSION
===========================================
Companion to the NCRP 147 Shielding Calculator. Where that app evaluates the
Archer empirical model, this one works from first principles:

    1. GENERATE a tube spectrum          (semi-analytical, not Monte Carlo)
    2. FILTER it through materials       (Beer-Lambert per energy bin)
    3. TRANSPORT it through a barrier    (Monte Carlo - the only MC step)

WHY THE SPLIT MATTERS
    Steps 1 and 2 are instant and need no Monte Carlo. SpekCalc and SpekPy
    are semi-analytical: Poludniowski used Monte Carlo once, upstream, to
    derive the electron energy/depth distribution in the tungsten target,
    then fitted an analytical bremsstrahlung expression to it. Filtering a
    spectrum is just N'(E) = N(E) exp(-mu(E) x) bin by bin.

    Step 2 alone gives NARROW-BEAM transmission, which always over-attenuates
    because it discards scattered photons that do reach the far side of the
    barrier. Archer parameters are BROAD-BEAM - they include that buildup.
    Monte Carlo is what bridges the two, so plotting narrow-beam, Monte Carlo
    and Archer together shows the buildup factor directly and gives an
    independent check on the Archer fits.

REQUIRES: numpy, scipy, matplotlib, pandas       (as the main app)
OPTIONAL: spekpy   - better spectra; a rough built-in model is used without it
          numba    - 50-100x faster Monte Carlo; falls back to numpy

DATA: xsdata.npz must sit beside this file. It holds photoelectric, coherent
      and incoherent mass attenuation coefficients for 15 curated materials
      (shielding/filtration compounds and a few pure elements) plus every
      element Z=1-92, all from 1 to 800 keV, derived from the Elam tables
      via xraydb. Verified to 0.011 % against xraydb and to published NIST
      values (water 0.1707 cm2/g at 100 keV). The 15 curated materials have
      absorption edges explicitly resolved on the energy grid; the added
      elements share that same grid without extra points at their own
      edges, so expect a little more smoothing right at an added element's
      K-edge than at lead's or tungsten's. K-fluorescence yield for the
      added elements is a jump-ratio approximation, reproducing lead and
      tungsten's own curated values to about 3 % rather than exactly.

SCORING NOTE: transmission here is an ENERGY-FLUENCE ratio. A true air-kerma
      ratio needs mass energy-absorption coefficients, which are not in the
      Elam tables, so they are deliberately not faked. The difference is
      second order for a ratio but it is not zero - see the Notes tab.
"""

import os
import threading
import queue
import time
import sys
import math
import datetime
import traceback

_MISSING = []
try:
    import numpy as np
except ImportError:
    _MISSING.append("numpy")
try:
    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                   NavigationToolbar2Tk)
    from matplotlib.widgets import RectangleSelector
    import matplotlib.patches as mpatches
except ImportError:
    _MISSING.append("matplotlib")
try:
    import pandas as pd
except ImportError:
    _MISSING.append("pandas")

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

if _MISSING:
    msg = ("This application needs:\n\n    " + "\n    ".join(_MISSING)
           + "\n\npip install " + " ".join(_MISSING))
    try:
        _r = tk.Tk(); _r.withdraw(); messagebox.showerror("Missing packages", msg)
    except Exception:
        print(msg)
    sys.exit(1)

# -- optional accelerator ---------------------------------------------------
try:
    from numba import njit
    HAVE_NUMBA = True
except ImportError:
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        def deco(fn):
            return fn
        return deco

# -- optional spectrum engine ----------------------------------------------
try:
    import spekpy as sp
    HAVE_SPEKPY = True
except Exception:
    HAVE_SPEKPY = False

# -- optional DICOM reader ---------------------------------------------------
try:
    import pydicom
    HAVE_PYDICOM = True
except Exception:
    HAVE_PYDICOM = False

# numpy 2.0 renamed trapz to trapezoid and removed the old name.
# Bind once so the app runs on both 1.x and 2.x.
_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")

APP_TITLE = "X-ray Spectrum and Monte Carlo Transmission"
APP_VERSION = "1.0"

# =============================================================================
#  SECTION 1 -- CROSS SECTION DATA
# =============================================================================

_HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() \
    else os.getcwd()


class CrossSections:
    """Log-log interpolated photo / coherent / incoherent coefficients."""

    def __init__(self, path=None):
        path = path or os.path.join(_HERE, "xsdata.npz")
        if not os.path.exists(path):
            raise FileNotFoundError(
                "xsdata.npz not found next to this script (%s)." % path)
        z = np.load(path)
        self.E = z["energy_keV"]
        self.logE = np.log(self.E)
        self.materials = sorted({k.split("__")[0] for k in z.files
                                 if "__" in k})
        self._d = {}
        for m in self.materials:
            self._d[m] = {
                "photo": z[m + "__photo"],
                "coh": z[m + "__coh"],
                "incoh": z[m + "__incoh"],
                "density": float(z[m + "__density"][0]),
            }
        self.Emin, self.Emax = float(self.E[0]), float(self.E[-1])
        self.ff_x = z["formfactor_x"] if "formfactor_x" in z.files else None
        self._ff = {m: z[m + "__ffcum"] for m in self.materials
                    if (m + "__ffcum") in z.files}
        self.has_form_factors = bool(self._ff)
        self._fl = {}
        for mm_ in self.materials:
            if (mm_ + "__kedge") in z.files:
                self._fl[mm_] = (z[mm_ + "__kedge"], z[mm_ + "__kyield"],
                                 z[mm_ + "__klines"], z[mm_ + "__kcum"])
        self.has_fluorescence = bool(self._fl)

    def fluor(self, material):
        """(K edge keV, net yield, line energies, cumulative intensities).

        Note the .item() calls. numpy deprecated float() on a one-element
        array in 1.25 and it now raises, so a stored scalar must be
        extracted explicitly rather than coerced.
        """
        e, y, ln, cu = self._fl[material]
        return (float(np.asarray(e).ravel()[0]),
                float(np.asarray(y).ravel()[0]),
                np.asarray(ln, dtype=np.float64),
                np.asarray(cu, dtype=np.float64))

    def form_factor(self, material):
        """Cumulative of x.F(x)^2 for coherent scattering, or None."""
        return self._ff.get(material)

    def density(self, material):
        return self._d[material]["density"]

    def _interp(self, material, kind, E_keV):
        y = self._d[material][kind]
        Eq = np.clip(np.asarray(E_keV, dtype=float), self.Emin, self.Emax)
        # log-log is the right space: these are near power laws between edges
        with np.errstate(divide="ignore", invalid="ignore"):
            out = np.exp(np.interp(np.log(Eq), self.logE,
                                   np.log(np.maximum(y, 1e-300))))
        return out

    def mass_atten(self, material, E_keV, kind="total"):
        """Mass attenuation coefficient, cm2/g."""
        if kind == "total":
            return (self._interp(material, "photo", E_keV)
                    + self._interp(material, "coh", E_keV)
                    + self._interp(material, "incoh", E_keV))
        return self._interp(material, kind, E_keV)

    def mu(self, material, E_keV, kind="total", density=None):
        """Linear attenuation coefficient, cm^-1."""
        rho = self.density(material) if density is None else density
        return self.mass_atten(material, E_keV, kind) * rho

    def arrays(self, material, density=None):
        """(logE, mu_photo, mu_coh, mu_incoh) in cm^-1, for the MC kernel."""
        rho = self.density(material) if density is None else density
        d = self._d[material]
        return (self.logE,
                d["photo"] * rho, d["coh"] * rho, d["incoh"] * rho)


# -- Mass energy-absorption coefficients for dry air -----------------------
# Hubbell & Seltzer, NIST tables of X-ray mass attenuation and mass
# energy-absorption coefficients (NISTIR 5632), "Air, Dry (Near Sea Level)".
# Taken verbatim from physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/air.html
# Needed to convert a photon spectrum into air kerma. Energies in keV,
# coefficients in cm2/g. The argon K edge at 3.203 keV is bracketed.
_AIR_MUEN_KEV = np.array([
    1.0, 1.5, 2.0, 3.0, 3.2029, 3.20291, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0,
    20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 150.0, 200.0, 300.0, 400.0,
    500.0, 600.0, 800.0, 1000.0])
_AIR_MUEN_CM2G = np.array([
    3.599e3, 1.188e3, 5.262e2, 1.614e2, 1.330e2, 1.460e2, 7.636e1, 3.931e1,
    2.270e1, 9.446, 4.742, 1.334, 5.389e-1, 1.537e-1, 6.833e-2, 4.098e-2,
    3.041e-2, 2.407e-2, 2.325e-2, 2.496e-2, 2.672e-2, 2.872e-2, 2.949e-2,
    2.966e-2, 2.953e-2, 2.882e-2, 2.789e-2])
_AIR_LOGE = np.log(_AIR_MUEN_KEV)
_AIR_LOGMU = np.log(_AIR_MUEN_CM2G)

# keV -> J, and cm2/g -> the factor turning fluence x energy x mu_en into Gy
_KEV_TO_J = 1.602176634e-16
_KERMA_FACTOR = _KEV_TO_J * 1000.0      # x1000 converts per gram to per kg


# -- IEC 61267 RQA beam qualities ------------------------------------------
# RQA beams are RQR beams hardened with added aluminium, representing the
# beam that has passed through a patient. They are DEFINED BY THEIR HVL -
# the added filtration is nominal guidance, and the procedure in the
# standard is to adjust filtration until the target HVL is met.
#
# Entries marked "unverified" could not be confirmed from an accessible
# published source and should be checked against the standard before use.
#   (kV, nominal added Al mm, target HVL mm Al, confidence)
RQA_ED2 = {                       # IEC 61267:2005, 2nd edition
    "RQA2":  (40,  4.0,  2.2,  "verified"),
    "RQA3":  (50, 10.0,  3.8,  "verified"),
    "RQA4":  (60, 16.0,  5.4,  "unverified"),
    "RQA5":  (70, 21.0,  6.8,  "verified"),
    "RQA6":  (80, 26.0,  8.2,  "unverified"),
    "RQA7":  (90, 30.0,  9.2,  "unverified"),
    "RQA8":  (100, 34.0, 10.1, "unverified"),
    "RQA9":  (120, 40.0, 11.6, "verified"),
    "RQA10": (150, 45.0, 13.3, "unverified"),
}
RQA_EDITION = "IEC 61267 Ed.2 (2005)"

# -- IEC 62220-1 Table 3: mandatory SNR_in^2 ---------------------------------
# SNR_in^2, in mm^-2.uGy^-1, is numerically the incident photon fluence per
# unit air kerma (Q/K_air): for an ideal photon-counting detector SNR^2
# equals the mean photon count, so this is a fluence figure wearing a DQE
# name. It is independent of HVL - two spectra can share an RQA's target
# HVL while differing in Q/K_air, because IEC 61267 constrains the beam by
# HVL and kV, not by the fluence shape.
#
# Only these four qualities carry a MANDATORY value in IEC 62220-1 Table 3.
# Values are exactly as supplied from the user's copy of the standard, not
# re-derived, and no confidence flag is attached because there was nothing
# to independently check them against.
RQA_SNR_IN2 = {
    "RQA3": 20673.0,
    "RQA5": 29653.0,
    "RQA7": 32490.0,
    "RQA9": 31007.0,
}


def rqa_fluence_per_kerma(spectrum):
    """
    Photon fluence per unit air kerma, mm^-2.uGy^-1 - what IEC 62220-1 calls
    SNR_in^2.

    Needs an ABSOLUTE spectrum (Calibrate, or SpekPy's own normalisation).
    This is a ratio of two absolute quantities - fluence and kerma - not
    two shape-only ones like HVL, so a relative-only spectrum has no
    meaningful answer here. None is returned rather than a number that
    looks real but is not comparable to the standard.
    """
    k = spectrum_output(spectrum)
    if not k:
        return None
    return spectrum.fluence / 100.0 / k     # photons/cm2 -> /mm2, per uGy/mAs


def rqa_fluence_match(spectrum):
    """
    Compare the spectrum's own Q/K_air against the IEC 62220-1 Table 3
    mandatory values. Returns (your_value_or_None, rows_sorted_by_closeness).
    """
    q = rqa_fluence_per_kerma(spectrum)
    if q is None:
        return None, []
    out = []
    for name, target in RQA_SNR_IN2.items():
        out.append({"quality": name, "target_SNRin2": target,
                    "your_SNRin2": round(q, 1),
                    "deviation_pc": round(100 * (q - target) / target, 2)})
    out.sort(key=lambda d: abs(d["deviation_pc"]))
    return q, out


def rqa_match(spectrum):
    """
    Which RQA quality does this spectrum most resemble?
    Compares measured HVL against each target. Returns a list of dicts
    sorted by how close the HVL is.
    """
    tbl = RQA_ED2
    hvl = spectrum.hvl("aluminium")
    out = []
    for name, (kv, add, target, conf) in sorted(tbl.items(),
                                                key=lambda kv_: kv_[1][0]):
        out.append({"quality": name, "kV": kv, "nominal_added_Al_mm": add,
                    "target_HVL_mm": target, "your_HVL_mm": round(hvl, 3),
                    "deviation_pc": round(100 * (hvl - target) / target, 2),
                    "confidence": conf})
    out.sort(key=lambda d: abs(d["deviation_pc"]))
    return out


def rqa_solve_filtration(base_spectrum, quality, max_mm=80.0):
    """
    How much aluminium must be added to `base_spectrum` to reach the target
    HVL of `quality`? This mirrors the procedure in the standard, which
    adjusts filtration until the HVL target is met rather than trusting the
    nominal figure.

    Returns (added_mm, achieved_HVL, target_HVL) or (None, ...) if the
    target cannot be reached.
    """
    if quality not in RQA_ED2:
        raise ValueError("%s is not a defined RQA quality." % quality)
    kv, nominal, target, _conf = RQA_ED2[quality]

    def hvl_after(mm):
        return base_spectrum.filtered("aluminium", mm).hvl("aluminium")

    if hvl_after(0.0) > target:
        return None, hvl_after(0.0), target      # already too hard
    lo, hi = 0.0, 1.0
    while hvl_after(hi) < target:
        hi *= 2.0
        if hi > max_mm:
            return None, hvl_after(max_mm), target
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if hvl_after(mid) < target:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-4:
            break
    mm = 0.5 * (lo + hi)
    return mm, hvl_after(mm), target


def rqa_solve_filtration_cu(base_spectrum, quality, cu_mm=0.5, max_mm=80.0):
    """
    Same idea as rqa_solve_filtration, but for the IEC-recognised
    alternative recipe: a fixed thickness of copper (0.5 mm by default -
    IEC 61267's own alternative for RQA5, 0.5 mm Cu + 2 mm Al) plus
    whatever aluminium on top is needed to reach the target HVL.

    A different filter composition landing on the same HVL can still
    produce a measurably different spectrum, and a different photon
    fluence per unit air kerma - see the fluence comparison in "Identify
    beam" - so this is offered ALONGSIDE the pure-Al answer, not as a
    like-for-like replacement.

    Returns (added_al_mm, achieved_HVL, target_HVL) or (None, ...) if the
    target cannot be reached with this copper thickness.
    """
    if quality not in RQA_ED2:
        raise ValueError("%s is not a defined RQA quality." % quality)
    kv, nominal, target, _conf = RQA_ED2[quality]
    base_cu = base_spectrum.filtered("copper", cu_mm)

    def hvl_after(mm):
        return base_cu.filtered("aluminium", mm).hvl("aluminium")

    if hvl_after(0.0) > target:
        return None, hvl_after(0.0), target      # 0.5 mm Cu alone overshoots
    lo, hi = 0.0, 1.0
    while hvl_after(hi) < target:
        hi *= 2.0
        if hi > max_mm:
            return None, hvl_after(max_mm), target
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if hvl_after(mid) < target:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-4:
            break
    mm = 0.5 * (lo + hi)
    return mm, hvl_after(mm), target


# -- Output calibration ----------------------------------------------------
# A measurement at ONE operating point fixes the absolute scale of the model
# everywhere. The model supplies the SHAPE - how output varies with kVp and
# filtration - and the measurement supplies the one number the shape cannot
# provide. So the stored quantity is a SCALE FACTOR, not the measurement.
#
#     scale = measured_uGy_per_mAs_at_1m / model_kerma_integral
#
# Every later spectrum then reports  kerma_integral x scale.
CALIBRATION = {
    "scale": None,        # uGy per mAs at 1 m, per unit of kerma_integral
    "kvp": None,          # where it was taken
    "filters": None,
    "measured": None,
    "distance": None,
    "engine": None,
    "spek_ratio": None,   # measured / SpekPy's own value, if available
}


def clear_calibration():
    for k in CALIBRATION:
        CALIBRATION[k] = None


def set_calibration(spectrum, measured_uGy_per_mAs, distance_m):
    """
    Fix the absolute scale from one measurement. Returns the scale factor.
    """
    if measured_uGy_per_mAs <= 0:
        raise ValueError("Measured output must be greater than zero.")
    if distance_m <= 0:
        raise ValueError("Distance must be greater than zero.")
    at_1m = measured_uGy_per_mAs * distance_m ** 2
    integral = spectrum.kerma_integral()
    if integral <= 0:
        raise ValueError("Spectrum carries no kerma - cannot calibrate.")
    CALIBRATION["scale"] = at_1m / integral
    CALIBRATION["kvp"] = spectrum.kvp
    CALIBRATION["filters"] = list(spectrum.filters)
    CALIBRATION["measured"] = measured_uGy_per_mAs
    CALIBRATION["distance"] = distance_m
    CALIBRATION["engine"] = spectrum.model
    CALIBRATION["spek_ratio"] = (
        at_1m / spectrum.kerma_uGy_per_mAs
        if spectrum.kerma_uGy_per_mAs else None)
    return CALIBRATION["scale"]


def spectrum_output(spectrum):
    """
    Air kerma in uGy per mAs at 1 m, or None if there is no absolute scale.
    A user calibration takes precedence over SpekPy's own normalisation.
    """
    if CALIBRATION["scale"] is not None:
        return spectrum.kerma_integral() * CALIBRATION["scale"]
    return spectrum.kerma_uGy_per_mAs


def output_basis():
    if CALIBRATION["scale"] is not None:
        return ("calibrated: %s uGy/mAs at %s m, %s kVp"
                % (CALIBRATION["measured"], CALIBRATION["distance"],
                   CALIBRATION["kvp"]))
    return "SpekPy absolute scale"


def air_muen(E_keV):
    """Mass energy-absorption coefficient of dry air, cm2/g."""
    E = np.clip(np.asarray(E_keV, dtype=float), _AIR_MUEN_KEV[0],
                _AIR_MUEN_KEV[-1])
    return np.exp(np.interp(np.log(E), _AIR_LOGE, _AIR_LOGMU))


def air_kerma(E_keV, fluence_per_cm2):
    """
    Air kerma in uGy from a photon fluence spectrum.

        K = SUM  Phi(E) . E . (mu_en/rho)_air(E)

    fluence_per_cm2 must be photons per cm2 in each bin (not per keV).
    """
    E = np.asarray(E_keV, dtype=float)
    phi = np.asarray(fluence_per_cm2, dtype=float)
    gy = float(np.sum(phi * E * air_muen(E)) * _KERMA_FACTOR)
    return gy * 1e6                      # Gy -> uGy


XS = None          # populated in main()

MATERIAL_LABEL = {
    "lead": "Lead", "concrete": "Concrete", "gypsum": "Gypsum",
    "steel": "Steel", "glass": "Glass", "borosilicate": "Borosilicate",
    "leadglass": "Lead glass",
    "aluminium": "Aluminium", "copper": "Copper", "water": "Water",
    "pmma": "PMMA", "tungsten": "Tungsten", "beryllium": "Beryllium",
    "air": "Air", "tissue": "Tissue",
}

# -- pure elements, Z=1-92, added via xsdata.npz --------------------------
# Same source and method as the curated materials above (Elam tables via
# xraydb) - see build_xsdata.py. K-fluorescence yield for these is a
# jump-ratio approximation (P_K x omega_K) rather than whatever compilation
# gave the lead/tungsten values above, so it reproduces those two to about
# 3%, not exactly - a second-order effect on an already-optional channel
# (the "include fluorescence" toggle), not the main photon transport.
ELEMENT_LABEL = {
    "hydrogen": "Hydrogen", "helium": "Helium", "lithium": "Lithium",
    "boron": "Boron", "carbon": "Carbon", "nitrogen": "Nitrogen",
    "oxygen": "Oxygen", "fluorine": "Fluorine", "neon": "Neon",
    "sodium": "Sodium", "magnesium": "Magnesium", "silicon": "Silicon",
    "phosphorus": "Phosphorus", "sulfur": "Sulfur", "chlorine": "Chlorine",
    "argon": "Argon", "potassium": "Potassium", "calcium": "Calcium",
    "scandium": "Scandium", "titanium": "Titanium", "vanadium": "Vanadium",
    "chromium": "Chromium", "manganese": "Manganese", "iron": "Iron",
    "cobalt": "Cobalt", "nickel": "Nickel", "zinc": "Zinc",
    "gallium": "Gallium", "germanium": "Germanium", "arsenic": "Arsenic",
    "selenium": "Selenium", "bromine": "Bromine", "krypton": "Krypton",
    "rubidium": "Rubidium", "strontium": "Strontium", "yttrium": "Yttrium",
    "zirconium": "Zirconium", "niobium": "Niobium", "molybdenum": "Molybdenum",
    "technetium": "Technetium", "ruthenium": "Ruthenium", "rhodium": "Rhodium",
    "palladium": "Palladium", "silver": "Silver", "cadmium": "Cadmium",
    "indium": "Indium", "tin": "Tin", "antimony": "Antimony",
    "tellurium": "Tellurium", "iodine": "Iodine", "xenon": "Xenon",
    "cesium": "Caesium", "barium": "Barium", "lanthanum": "Lanthanum",
    "cerium": "Cerium", "praseodymium": "Praseodymium",
    "neodymium": "Neodymium", "promethium": "Promethium",
    "samarium": "Samarium", "europium": "Europium", "gadolinium": "Gadolinium",
    "terbium": "Terbium", "dysprosium": "Dysprosium", "holmium": "Holmium",
    "erbium": "Erbium", "thulium": "Thulium", "ytterbium": "Ytterbium",
    "lutetium": "Lutetium", "hafnium": "Hafnium", "tantalum": "Tantalum",
    "rhenium": "Rhenium", "osmium": "Osmium", "iridium": "Iridium",
    "platinum": "Platinum", "gold": "Gold", "mercury": "Mercury",
    "thallium": "Thallium", "bismuth": "Bismuth", "polonium": "Polonium",
    "astatine": "Astatine", "radon": "Radon", "francium": "Francium",
    "radium": "Radium", "actinium": "Actinium", "thorium": "Thorium",
    "protactinium": "Protactinium", "uranium": "Uranium",
}
MATERIAL_LABEL.update(ELEMENT_LABEL)


# =============================================================================
#  SECTION 2 -- SPECTRUM GENERATION
# =============================================================================

# Tungsten K lines: energies and relative intensities from the Elam tables.
W_K_LINES = [(57.981, 0.2877), (59.318, 0.5006), (66.950, 0.0564),
             (67.244, 0.1088), (67.653, 0.0025), (69.035, 0.0431)]
W_K_EDGE = 69.525


class Spectrum:
    """Photon fluence per keV on a 1 keV grid."""

    def __init__(self, E, N, kvp, model, filters=None, absolute=False,
                 kerma_uGy_per_mAs=None):
        self.E = np.asarray(E, dtype=float)          # bin centres, keV
        self.N = np.asarray(N, dtype=float)          # photons per keV
        self.kvp = kvp
        self.model = model
        self.filters = list(filters or [])           # [(material, mm, rho)]
        # absolute=True means N carries a real fluence scale, so air kerma
        # can be quoted in uGy. The built-in model is arbitrary units.
        self.absolute = absolute
        self.kerma_uGy_per_mAs = kerma_uGy_per_mAs

    # -- derived quantities -------------------------------------------
    def copy(self):
        return Spectrum(self.E.copy(), self.N.copy(), self.kvp, self.model,
                        list(self.filters), self.absolute,
                        self.kerma_uGy_per_mAs)

    @property
    def fluence(self):
        return float(_trapz(self.N, self.E))

    @property
    def energy_fluence(self):
        return float(_trapz(self.N * self.E, self.E))

    @property
    def mean_energy(self):
        f = self.fluence
        return self.energy_fluence / f if f > 0 else 0.0

    def filtered(self, material, thickness_mm, density=None):
        """Beer-Lambert per bin. This is NOT Monte Carlo - no scatter."""
        if thickness_mm <= 0:
            return self.copy()
        mu = XS.mu(material, self.E, "total", density)     # cm^-1
        out = self.copy()
        out.N = self.N * np.exp(-mu * thickness_mm / 10.0)
        out.filters.append((material, thickness_mm,
                            density if density else XS.density(material)))
        if self.kerma_uGy_per_mAs is not None:
            k0 = self.kerma_integral()
            k1 = out.kerma_integral()
            out.kerma_uGy_per_mAs = (self.kerma_uGy_per_mAs * k1 / k0
                                     if k0 > 0 else 0.0)
        return out

    def kerma_integral(self):
        """SUM Phi.E.(mu_en/rho)_air over the spectrum, in the spectrum's
        own units. A ratio of two of these is unit-free."""
        dE = np.gradient(self.E)
        return float(np.sum(self.N * dE * self.E * air_muen(self.E)))

    def output_uGy_per_mAs(self, scale=None):
        """
        Air kerma at the reference distance, uGy per mAs.
        `scale` overrides with a measured value at the unfiltered point.
        Returns None if the spectrum has no absolute scale.
        """
        if scale is not None:
            return scale
        return self.kerma_uGy_per_mAs

    @staticmethod
    def archer_K_poly(kvp):
        """
        Unshielded primary air kerma per workload, mGy per mA.min at 1 m.

        Fitted directly to the measured data in Archer, Fewell, Conway and
        Quinn (1994), Med Phys 21, 1499-507. NCRP 147 reproduces the same
        coefficients but with one sign wrong, so this is the corrected form
        and the 1994 paper is the primary source, not the report.

        Valid over 50-150 kVp. A cubic has no obligation to behave outside
        its fitted range and does not - do not use it below 50 kVp.
        """
        return (1.222170186 - 0.056644073 * kvp
                + 0.001227127 * kvp ** 2 - 3.13579e-6 * kvp ** 3)

    def hvl(self, material="aluminium", target=0.5, density=None,
            quantity="kerma"):
        """
        Thickness of `material` halving the beam, in mm. Narrow beam.

        quantity = "kerma"   AIR KERMA halving - the conventional definition
                             and what an ion chamber or a RaySafe reports
                   "fluence" energy-fluence halving

        These are NOT interchangeable. Weighting by air kerma brings in the
        air mass energy-absorption coefficient, which varies steeply at low
        energy, so the two definitions disagree by tens of per cent across
        the diagnostic range. Anything compared against a measurement must
        use the air-kerma definition.
        """
        mu = XS.mu(material, self.E, "total", density)
        dE = np.gradient(self.E)
        if quantity == "kerma":
            w = self.N * self.E * air_muen(self.E) * dE
        else:
            w = self.N * self.E * dE
        base = float(np.sum(w))
        if base <= 0:
            return 0.0

        def frac(t_mm):
            return float(np.sum(w * np.exp(-mu * t_mm / 10.0))) / base

        lo, hi = 0.0, 1.0
        for _ in range(80):
            if frac(hi) <= target:
                break
            hi *= 2.0
            if hi > 1e5:
                return float("inf")
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if frac(mid) > target:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-7:
                break
        return 0.5 * (lo + hi)

    def effective_energy(self, material="aluminium"):
        """Monoenergetic energy with the same air-kerma HVL, keV."""
        h_mm = self.hvl(material)
        if not np.isfinite(h_mm) or h_mm <= 0:
            return 0.0
        target_mu = math.log(2.0) / (h_mm / 10.0)          # cm^-1
        grid = np.logspace(math.log10(max(XS.Emin, 1.0)),
                           math.log10(min(XS.Emax, 1000.0)), 4000)
        mu = XS.mu(material, grid, "total")
        # mu falls with E over the diagnostic range, so invert on a
        # monotonic decreasing branch
        idx = np.argmin(np.abs(mu - target_mu))
        return float(grid[idx])

    def homogeneity(self, material="aluminium"):
        h1 = self.hvl(material, 0.5)
        h2 = self.hvl(material, 0.25)
        return h1 / (h2 - h1) if (h2 - h1) > 0 else 0.0


def spectrum_spekpy(kvp, anode_angle=12.0, filters=(), mas=1.0, ripple=0.0):
    """SpekPy spectrum. filters = [(material_name, mm)]."""
    s = sp.Spek(kvp=kvp, th=anode_angle, dk=1.0, mas=mas)
    if ripple:
        try:
            s = sp.Spek(kvp=kvp, th=anode_angle, dk=1.0, mas=mas,
                        physics="spekpy-v2", ripple=ripple)
        except Exception:
            pass
    name_map = {"aluminium": "Al", "copper": "Cu", "beryllium": "Be",
                "water": "Water", "air": "Air", "tissue": "Water",
                "lead": "Pb", "tungsten": "W"}
    for mat, mm in filters:
        try:
            s.filter(name_map.get(mat, mat), mm)
        except Exception:
            pass
    E, N = s.get_spectrum(edges=False)
    k = None
    try:
        # SpekPy reports air kerma directly at its reference distance, which
        # avoids having to guess at its fluence normalisation convention
        k = float(np.asarray(s.get_kerma()).ravel()[0])
    except Exception:
        k = None
    return Spectrum(np.asarray(E), np.asarray(N), kvp, "SpekPy",
                    absolute=(k is not None), kerma_uGy_per_mAs=k)


def spectrum_kramers(kvp, filters=(), target="tungsten",
                     target_angle_deg=12.0):
    """
    Fallback bremsstrahlung model: Kramers' law with target self-filtration
    and tungsten K lines at their tabulated energies and intensities.

    This is ROUGH. It reproduces the shape and the K lines but not the
    accuracy of SpekPy or SpekCalc. Install spekpy for anything quantitative.
    """
    E = np.arange(1.0, float(kvp) + 1.0, 1.0)
    # Kramers: photon number per unit energy ~ Z (E0/E - 1)
    N = np.where(E < kvp, (kvp / np.maximum(E, 1e-6) - 1.0), 0.0)
    N = np.maximum(N, 0.0)

    # crude target self-absorption: photons born at mean depth, escaping at
    # the anode angle, attenuated in tungsten
    mu_w = XS.mu("tungsten", E, "total")
    depth_cm = 2.0e-4 / max(math.sin(math.radians(target_angle_deg)), 1e-3)
    N *= np.exp(-mu_w * depth_cm)

    # characteristic K lines, only if the tube can excite them
    if kvp > W_K_EDGE:
        brem_area = float(_trapz(N, E))
        drive = (kvp - W_K_EDGE) / max(kvp, 1.0)
        k_total = 0.55 * brem_area * drive
        for e_line, inten in W_K_LINES:
            if e_line < kvp:
                j = int(np.argmin(np.abs(E - e_line)))
                N[j] += k_total * inten          # one 1 keV bin
    s = Spectrum(E, N, kvp, "Kramers (approximate)")
    for mat, mm in filters:
        s = s.filtered(mat, mm)
    return s


def solve_filtration_from_hvl(kvp, measured_hvl_mm, anode_angle=12.0,
                              prefer_spekpy=True, max_mm=40.0):
    """
    Find the total aluminium-equivalent filtration that reproduces a
    MEASURED HVL at this kVp.

    This is a far stronger determination of a beam than a label. It captures
    inherent filtration, added filtration, anode ageing and anything else
    hardening the beam, all in one number, because it works from what the
    beam actually does rather than from what the paperwork claims.

    Returns (total_mm_Al, achieved_HVL) or (None, achievable_HVL) if the
    measured value cannot be reached - which usually means the spectrum
    model is already harder than the real tube.
    """
    if measured_hvl_mm <= 0:
        raise ValueError("Measured HVL must be greater than zero.")

    def hvl_at(mm):
        return make_spectrum(kvp, [("aluminium", mm)], anode_angle,
                             prefer_spekpy).hvl("aluminium")

    lo_hvl = hvl_at(0.0)
    if lo_hvl > measured_hvl_mm:
        return None, lo_hvl              # model already too hard
    lo, hi = 0.0, 1.0
    while hvl_at(hi) < measured_hvl_mm:
        hi *= 2.0
        if hi > max_mm:
            return None, hvl_at(max_mm)
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if hvl_at(mid) < measured_hvl_mm:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-4:
            break
    mm = 0.5 * (lo + hi)
    return mm, hvl_at(mm)


def make_spectrum(kvp, filters=(), anode_angle=12.0, prefer_spekpy=True):
    if prefer_spekpy and HAVE_SPEKPY:
        try:
            return spectrum_spekpy(kvp, anode_angle, filters)
        except Exception:
            pass
    return spectrum_kramers(kvp, filters, target_angle_deg=anode_angle)


# =============================================================================
#  SECTION 2B -- RADIOISOTOPE SOURCES
#
#  An alternative to a tube spectrum for the Monte Carlo tab: a sealed or
#  unsealed source's photon lines, for shielding calculations around
#  nuclear medicine, brachytherapy or industrial radiography sources
#  rather than an X-ray tube.
# =============================================================================

# (energy_keV, intensity per decay) for each line with intensity >= ~1%.
# Energies and intensities are widely published values (NNDC/ENSDF-derived,
# cross-checked against Ir-192 patent/reference literature during this
# session) but have NOT been individually re-verified line by line against
# ENSDF here the way the RQA HVL table was. Treat these as good for a
# shielding cross-check, not as a substitute for NNDC/ICRP-107 if the
# result needs to stand up to compliance review.
#
# F-18 is a positron emitter: intensity 1.936 reflects TWO 511 keV photons
# per decay (annihilation), not a single-line branching ratio like the others.
ISOTOPE_LINES = {
    "Am-241": [(26.34, 0.024), (59.541, 0.359)],
    "Ba-133": [(53.2, 0.022), (81.0, 0.367), (276.4, 0.0716), (302.9, 0.183),
              (356.0, 0.619), (383.8, 0.089)],
    "Co-57": [(122.06, 0.856), (136.47, 0.1068)],
    "Co-60": [(1173.2, 0.9985), (1332.5, 0.9998)],
    "Cs-137": [(661.657, 0.851)],
    "F-18": [(511.0, 1.936)],
    "Ga-67": [(93.3, 0.39), (184.6, 0.214), (300.2, 0.168), (393.5, 0.046)],
    "I-123": [(159.0, 0.833), (529.0, 0.014)],
    "I-125": [(27.4, 0.716), (31.0, 0.196), (35.5, 0.067)],
    "I-131": [(80.2, 0.026), (284.3, 0.061), (364.5, 0.815), (636.9, 0.072),
             (722.9, 0.0179)],
    "In-111": [(171.3, 0.907), (245.4, 0.941)],
    "Ir-192": [(206.0, 0.032), (296.0, 0.283), (308.5, 0.293), (316.5, 0.83),
              (468.1, 0.477), (604.4, 0.0823), (612.5, 0.0534)],
    "Se-75": [(136.0, 0.586), (264.7, 0.589), (279.5, 0.25), (400.7, 0.114)],
    "Tc-99m": [(140.5, 0.89)],
    "Xe-133": [(81.0, 0.37)],
}


def isotope_lines_beyond_xsdata(name):
    """Lines above XS.Emax (800 keV, the Elam table's own validity limit).
    xsdata.npz's cross-section arrays clip silently above this - a photon
    at these energies is transported using the frozen 800 keV attenuation
    coefficient, which is measurably wrong (real mu keeps falling with
    energy in the Compton-dominated regime). Co-60's two lines both sit
    here. Flag rather than silently clip."""
    return [e for e, _ in ISOTOPE_LINES[name] if e > XS.Emax]


def make_isotope_spectrum(name):
    """
    Build a Spectrum-compatible source from a radioisotope's photon lines.

    Each line becomes a single 1 keV histogram bin holding its full
    relative intensity - exactly the "line in a bin" convention already
    used for tungsten's characteristic K lines in spectrum_kramers().
    Every downstream calculation (fluence, kerma, HVL, Beer-Lambert
    filtering, the Monte Carlo sampler itself) already treats a spectrum
    as a histogram and reads bin centres directly, so a handful of
    isolated one-line bins is a valid degenerate case of the same
    machinery - nothing else needed to change to support it.

    Not absolute: these are relative intensities per decay, not a real
    fluence at a distance (that needs source activity and geometry, which
    this tool does not model). The Monte Carlo buildup factor and
    transmission ratio are unaffected, since both are ratios.

    "kvp" is set to the top of the constructed grid, purely so the
    existing energy-labelled plot axes and reports have something
    sensible to show - it is not a tube potential.
    """
    if name not in ISOTOPE_LINES:
        raise ValueError("%s is not in the isotope line list." % name)
    lines = ISOTOPE_LINES[name]
    top = int(math.ceil(max(e for e, _ in lines))) + 10
    E = np.arange(1.0, top + 1.0, 1.0)
    N = np.zeros_like(E)
    for e_line, inten in lines:
        j = int(np.argmin(np.abs(E - e_line)))
        N[j] += inten
    return Spectrum(E, N, float(top), "isotope: %s" % name, absolute=False)


# =============================================================================
#  SECTION 3 -- MONTE CARLO PHOTON TRANSPORT (1D SLAB)
#
#  Geometry: an infinite slab of thickness d, photons entering the front face
#  travelling along +z. Transmission is scored at the back face.
#
#  Variance reduction, which is what makes B = 1e-6 reachable at all:
#    * forced collisions  - at every step the probability of escaping the
#      slab without another interaction is scored analytically, then the
#      photon is forced to interact inside, carrying weight (1 - p_escape).
#      Nothing is ever "lost" to a lucky escape, so the estimator has far
#      lower variance than counting photons that happen to get through.
#    * implicit absorption - weight is multiplied by the scattering
#      probability rather than terminating on photoelectric absorption.
#    * Russian roulette   - kills negligible weights without bias.
#
#  Physics:
#    photoelectric  absorbed, with K-fluorescence tracked (toggle) for
#                   materials that have K-edge data - see the Notes tab
#    incoherent     Klein-Nishina, sampled by direct rejection,
#                   free-electron approximation (no binding correction)
#    coherent       real atomic form factor, sampled from the tabulated
#                   cumulative of x.F(x)^2 - see the Notes tab
# =============================================================================

ELECTRON_MASS_KEV = 510.998950


@njit(cache=True)
def _interp_log(x, logE, y):
    """Log-log interpolation of a cross-section array at energy x (keV)."""
    lx = math.log(x)
    n = logE.shape[0]
    if lx <= logE[0]:
        return y[0]
    if lx >= logE[n - 1]:
        return y[n - 1]
    lo = 0
    hi = n - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if logE[mid] <= lx:
            lo = mid
        else:
            hi = mid
    y0 = y[lo]
    y1 = y[hi]
    if y0 <= 0.0 or y1 <= 0.0:
        f = (lx - logE[lo]) / (logE[hi] - logE[lo])
        return y0 + f * (y1 - y0)
    f = (lx - logE[lo]) / (logE[hi] - logE[lo])
    return math.exp(math.log(y0) + f * (math.log(y1) - math.log(y0)))


@njit(cache=True)
def _sample_compton(E, seed_state):
    """
    Sample the Klein-Nishina distribution by direct rejection.

    The differential cross section in cos(theta), dropping constants, is

        f(ct) = r^2 (r + 1/r - sin^2 theta),   r = E'/E = 1/(1 + a(1 - ct))

    Since r <= 1 always, f <= 2, so uniform sampling of ct on [-1, 1] with
    acceptance probability f/2 is exact. Efficiency is about 65 per cent at
    low energy and degrades gently.

    This replaces a hand-written Kahn rejection method that was wrong: it
    over-produced forward scatter, giving a mean cos(theta) of 0.41 against
    the analytic 0.23 at 300 keV, and chi-squared per degree of freedom in
    the hundreds. Too-forward scattering keeps photons in the beam and
    inflates transmission, so the error propagated straight into the
    barrier results.

    Returns (E_out, cos_theta). Free-electron approximation, no binding.
    """
    a = E / ELECTRON_MASS_KEV
    for _ in range(500):
        ct = 2.0 * np.random.random() - 1.0
        r = 1.0 / (1.0 + a * (1.0 - ct))
        f = r * r * (r + 1.0 / r - (1.0 - ct * ct))
        if 2.0 * np.random.random() <= f:
            return E * r, ct
    return E, 1.0


@njit(cache=True)
def _lin_at(v, xs, ys):
    n = xs.shape[0]
    if v <= xs[0]:
        return ys[0]
    if v >= xs[n - 1]:
        return ys[n - 1]
    lo = 0
    hi = n - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= v:
            lo = mid
        else:
            hi = mid
    f = (v - xs[lo]) / (xs[hi] - xs[lo])
    return ys[lo] + f * (ys[hi] - ys[lo])


@njit(cache=True)
def _invert(v, ys, xs):
    n = ys.shape[0]
    if v <= ys[0]:
        return xs[0]
    if v >= ys[n - 1]:
        return xs[n - 1]
    lo = 0
    hi = n - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ys[mid] <= v:
            lo = mid
        else:
            hi = mid
    d = ys[hi] - ys[lo]
    f = 0.0 if d <= 0.0 else (v - ys[lo]) / d
    return xs[lo] + f * (xs[hi] - xs[lo])


@njit(cache=True)
def _sample_rayleigh(E, ff_x, ff_cum):
    """
    Coherent scattering angle from the real atomic form factor.

        dsigma/dOmega ~ (1 + cos^2 theta)/2 . F(x)^2
        x = sin(theta/2)/lambda,  lambda[A] = 12.398/E[keV]

    In x the shape is energy independent - only the upper limit 1/lambda
    moves - so sample x from the tabulated cumulative of x F(x)^2, convert
    to an angle, and apply the polarisation factor by rejection.

    This replaces a pair of crude bounds (zero deflection on one side,
    Thomson with no form factor on the other) with the actual distribution.
    Validated against exact integration to better than 0.0021 in mean
    cos(theta) across four materials and 40-120 keV.
    """
    lam = 12.398 / E
    c_max = _lin_at(1.0 / lam, ff_x, ff_cum)
    if c_max <= 0.0:
        return 1.0
    for _ in range(200):
        x = _invert(np.random.random() * c_max, ff_cum, ff_x)
        st2 = x * lam
        if st2 > 1.0:
            st2 = 1.0
        ct = 1.0 - 2.0 * st2 * st2
        if np.random.random() <= 0.5 * (1.0 + ct * ct):
            return ct
    return 1.0


@njit(cache=True)
def _sample_thomson():
    """Thomson angular distribution, dsigma/dOmega ~ (1 + cos^2 theta)."""
    for _ in range(200):
        ct = 2.0 * np.random.random() - 1.0
        if np.random.random() <= 0.5 * (1.0 + ct * ct):
            return ct
    return 1.0


@njit(cache=True)
def _new_direction(u, cos_t):
    """
    New direction cosine relative to +z after scattering by cos_t from a
    direction with cosine u, with the azimuth sampled uniformly.
    """
    st = math.sqrt(max(0.0, 1.0 - cos_t * cos_t))
    su = math.sqrt(max(0.0, 1.0 - u * u))
    phi = 2.0 * math.pi * np.random.random()
    return u * cos_t + su * st * math.cos(phi)


@njit(cache=True)
def _mc_slab(E0, d_cm, logE, mu_ph, mu_co, mu_in, resp, u_min, e_min,
             ff_x, ff_cum, k_edge, k_yield, k_lines, k_cum,
             w_cut, max_scatter, out_E_edges, hist):
    """
    Track one photon history.

    `resp` is the detector response on the same energy grid: E for an
    energy-fluence score, or E x (mu_en/rho)_air for an air-kerma score.

    `u_min` is the cosine of the detector's acceptance half-angle. A real
    broad-beam measurement uses a finite detector at a finite distance, so
    photons emerging at a wide angle miss it. Scoring every forward photon
    regardless of angle counts scatter no experiment ever sees, and inflates
    transmission. u_min = 0 collects the whole forward hemisphere.

    Returns (transmitted_weight, transmitted_response, backscatter, n_coll).
    """
    E = E0
    z = 0.0
    u = 1.0
    w = 1.0
    t_w = 0.0
    t_e = 0.0
    b_w = 0.0
    ncoll = 0

    for _ in range(max_scatter):
        ph = _interp_log(E, logE, mu_ph)
        co = _interp_log(E, logE, mu_co)
        inc = _interp_log(E, logE, mu_in)
        tot = ph + co + inc
        if tot <= 0.0:
            if u >= u_min and E >= e_min:
                t_w += w
                t_e += w * _interp_log(E, logE, resp)
            break

        # distance to whichever face the photon is heading for
        if u > 1e-9:
            dist = (d_cm - z) / u
            forward = True
        elif u < -1e-9:
            dist = z / (-u)
            forward = False
        else:
            # travelling parallel to the faces: it will interact eventually
            dist = 1e30
            forward = True

        p_esc = math.exp(-tot * dist) if dist < 1e29 else 0.0
        if forward:
            # a real detector sees neither wide angles nor very soft photons
            if u >= u_min and E >= e_min:
                c = w * p_esc
                t_w += c
                t_e += c * _interp_log(E, logE, resp)
                # bin by the energy the photon actually has when it leaves,
                # which is the transmitted spectrum
                nb = out_E_edges.shape[0] - 1
                if out_E_edges[0] <= E < out_E_edges[nb]:
                    lo2 = 0
                    hi2 = nb
                    while hi2 - lo2 > 1:
                        md = (lo2 + hi2) // 2
                        if out_E_edges[md] <= E:
                            lo2 = md
                        else:
                            hi2 = md
                    hist[lo2] += c
        else:
            b_w += w * p_esc

        w *= (1.0 - p_esc)                      # forced collision
        if w < w_cut:
            if np.random.random() > 0.25:       # Russian roulette
                break
            w *= 4.0

        # sample the collision point within the remaining path
        r = np.random.random()
        s = -math.log(1.0 - r * (1.0 - p_esc)) / tot
        if s > dist:
            s = dist
        z += s * u
        if z < 0.0:
            z = 0.0
        elif z > d_cm:
            z = d_cm
        ncoll += 1

        # Implicit absorption. Photoelectric absorption above the K edge is
        # not the end of the story: the atom re-emits a K fluorescence
        # photon at 73-85 keV for lead, which sits in the window where lead
        # is most transparent. So that channel survives too, with the net
        # yield (K-shell share of photoelectric x fluorescence yield).
        p_fl = 0.0
        if E > k_edge and k_yield > 0.0:
            p_fl = ph * k_yield
        w *= (co + inc + p_fl) / tot
        if w < w_cut:
            if np.random.random() > 0.25:
                break
            w *= 4.0

        # which scatter?
        chan = co + inc + p_fl
        if chan <= 0.0:
            break                               # pure absorber, nothing left
        r2 = np.random.random() * chan
        if r2 < p_fl:
            # K fluorescence: new energy from the line spectrum, emitted
            # isotropically, so the direction is completely reset
            xi = np.random.random()
            k = 0
            while k < k_cum.shape[0] - 1 and k_cum[k] < xi:
                k += 1
            E = k_lines[k]
            u = 2.0 * np.random.random() - 1.0
            continue
        elif r2 < p_fl + co:
            cos_t = _sample_rayleigh(E, ff_x, ff_cum)
        else:
            E, cos_t = _sample_compton(E, 0)
            if E < 1.0:
                break
        u = _new_direction(u, cos_t)

    return t_w, t_e, b_w, ncoll


@njit(cache=True)
def _mc_batch(E0_arr, d_cm, logE, mu_ph, mu_co, mu_in, resp, u_min, e_min,
              ff_x, ff_cum, k_edge, k_yield, k_lines, k_cum,
              w_cut, max_scatter, out_E_edges, hist):
    """
    Run a batch. Returns (sum_resp, sum_resp2, sum_incident_resp, back, coll)
    so the caller can form an unbiased ratio and its standard error.
    """
    n = E0_arr.shape[0]
    sw = 0.0
    sw2 = 0.0
    se = 0.0
    sb = 0.0
    nc = 0
    nb = out_E_edges.shape[0] - 1
    for i in range(n):
        tw, te, bw, ncoll = _mc_slab(E0_arr[i], d_cm, logE, mu_ph, mu_co,
                                     mu_in, resp, u_min, e_min, ff_x, ff_cum,
                                     k_edge, k_yield, k_lines, k_cum,
                                     w_cut, max_scatter, out_E_edges, hist)
        # te is the response-weighted transmitted signal for this history;
        # se accumulates the incident response so the ratio is unbiased
        sw += te
        sw2 += te * te
        se += _interp_log(E0_arr[i], logE, resp)
        sb += bw
        nc += ncoll
    return sw, sw2, se, sb, nc


def run_monte_carlo(spectrum, material, thickness_mm, n_hist=200000,
                    density=None, seed=1234, progress=None, batch=20000,
                    score="kerma", acceptance_deg=90.0, e_threshold=0.0,
                    fluorescence=True):
    """
    Broad-beam transmission of `spectrum` through `thickness_mm` of material.

    score = "kerma"   air-kerma ratio, the quantity shielding uses
            "fluence" energy-fluence ratio

    acceptance_deg is the detector's acceptance half-angle. 90 degrees
    collects the entire forward hemisphere, which is what an infinite plane
    detector would do and is NOT what a real measurement does. Narrowing it
    discards wide-angle scatter that would have missed a real detector.

    Returns a dict with narrow-beam and Monte Carlo transmission, the
    buildup factor, a relative standard error, and the transmitted spectrum.
    """
    np.random.seed(seed)
    logE, mu_ph, mu_co, mu_in = XS.arrays(material, density)
    d_cm = thickness_mm / 10.0
    u_min = math.cos(math.radians(min(max(acceptance_deg, 0.01), 90.0)))
    ff_x = XS.ff_x
    ff_cum = XS.form_factor(material)
    k_edge, k_yield, k_lines, k_cum = XS.fluor(material)
    if not fluorescence:
        k_edge, k_yield = 1e9, 0.0
    if ff_x is None or ff_cum is None:
        raise ValueError(
            "xsdata.npz has no form factor data for %s. Replace it with the "
            "version that came with this build." % material)

    E = spectrum.E
    N = np.maximum(spectrum.N, 0.0)
    if N.sum() <= 0:
        raise ValueError("Spectrum is empty.")

    # detector response on the cross-section grid
    grid = np.exp(logE)
    if score == "kerma":
        resp = grid * air_muen(grid)
    else:
        resp = grid.copy()

    # Sample incident energies from the PHOTON fluence distribution; the
    # response weighting is applied when scoring, which keeps the ratio
    # unbiased for either choice of score.
    #
    # The spectrum is a HISTOGRAM, so pick a bin and take its centre.
    # Interpolating a cumulative sum instead (np.interp on cumsum) silently
    # shifts every sampled photon down by half a bin, which biases the
    # result low by a few per cent per mean free path.
    dE = np.gradient(E)
    wgt = N * dE
    if wgt.sum() <= 0:
        raise ValueError("Spectrum has no photons.")
    cdf = np.cumsum(wgt)
    cdf = cdf / cdf[-1]

    # Bin the transmitted spectrum on the INCIDENT energy grid. Uncollided
    # photons keep their exact incident energy, so a finer output grid
    # leaves them landing in every other bin and produces a comb artifact.
    dEg = float(np.mean(np.diff(E))) if E.size > 1 else 1.0
    edges = np.concatenate([E - 0.5 * dEg, [E[-1] + 0.5 * dEg]])
    hist = np.zeros(len(edges) - 1)

    tot_w = tot_w2 = tot_e = tot_b = 0.0
    tot_c = 0
    done = 0
    while done < n_hist:
        m = min(batch, n_hist - done)
        u = np.random.random(m)
        E0 = E[np.searchsorted(cdf, u, side="left").clip(0, len(E) - 1)]
        sw, sw2, se, sb, nc = _mc_batch(
            E0.astype(np.float64), d_cm, logE, mu_ph, mu_co, mu_in, resp,
            u_min, float(e_threshold), ff_x, ff_cum,
            k_edge, k_yield, k_lines, k_cum, 1e-8, 400, edges, hist)
        tot_w += sw; tot_w2 += sw2; tot_e += se; tot_b += sb; tot_c += nc
        done += m
        if progress:
            progress(done / n_hist)

    # ratio of transmitted response to incident response
    B_mc = tot_w / tot_e if tot_e > 0 else 0.0
    mean_r = tot_w / n_hist
    var = max(tot_w2 / n_hist - mean_r ** 2, 0.0)
    rse = (math.sqrt(var / n_hist) / mean_r) if mean_r > 0 else float("inf")

    mu = XS.mu(material, E, "total", density)
    resp_in = E * air_muen(E) if score == "kerma" else E
    ef0 = float(_trapz(N * resp_in, E))
    ef1 = float(_trapz(N * resp_in * np.exp(-mu * d_cm), E))
    B_narrow = ef1 / ef0 if ef0 > 0 else 0.0

    centres = 0.5 * (edges[:-1] + edges[1:])
    return {
        "B_mc": B_mc,
        "B_narrow": B_narrow,
        "buildup": (B_mc / B_narrow) if B_narrow > 0 else float("nan"),
        "rel_std_err": rse,
        "backscatter": tot_b / n_hist,
        "mean_collisions": tot_c / n_hist,
        "hist_E": centres,
        "hist_w": hist,
        "hist_norm": (hist / hist.sum()) if hist.sum() > 0 else hist,
        "n_hist": n_hist,
        "thickness_mm": thickness_mm,
        "material": material,
        "score": score,
        "acceptance_deg": acceptance_deg,
        "e_threshold": e_threshold,
    }


# -- Archer reference, for overlay only ------------------------------------
# Primary NCRP 147 Table A.4 parameters, carried here purely so the Monte
# Carlo result can be plotted against the empirical model. The shielding
# calculator remains the authority for design work.
ARCHER_REF = {
    "lead": {50: (8.801, 27.28, 0.2957), 60: (6.951, 24.89, 0.4198),
             70: (5.369, 23.49, 0.5881), 80: (4.040, 21.69, 0.7187),
             90: (3.067, 18.83, 0.7726), 100: (2.500, 15.28, 0.7557),
             110: (2.296, 11.70, 0.6827), 120: (2.246, 8.950, 0.5873),
             125: (2.219, 7.923, 0.5386), 140: (2.009, 5.916, 0.4018),
             150: (1.757, 5.177, 0.3156)},
    "concrete": {50: (0.09, 0.17, 0.232), 60: (0.063, 0.169, 0.273),
                 70: (0.051, 0.169, 0.385), 80: (0.046, 0.154, 0.493),
                 90: (0.042, 0.113, 0.469), 100: (0.039, 0.085, 0.427),
                 110: (0.037, 0.074, 0.475), 120: (0.036, 0.071, 0.607),
                 125: (0.035, 0.071, 0.697), 140: (0.03345, 0.07476, 1.047),
                 150: (0.03243, 0.08599, 1.467)},
    "gypsum": {50: (0.038, 0.087, 0.51), 60: (0.030, 0.080, 0.617),
               70: (0.023, 0.072, 0.729), 80: (0.019, 0.061, 0.810),
               90: (0.016, 0.050, 0.859), 100: (0.015, 0.042, 0.894),
               110: (0.013, 0.035, 0.924), 120: (0.012, 0.030, 0.957),
               125: (0.012, 0.029, 0.968), 140: (0.01088, 0.02436, 0.9964),
               150: (0.0103, 0.02198, 1.013)},
}


def _archer_transmission(x_mm, al, be, ga):
    """The Archer three-parameter transmission equation itself, shared by
    both the kVp-indexed table below and the radionuclide table."""
    return ((1 + be / al) * math.exp(al * ga * x_mm) - be / al) ** (-1.0 / ga)


def archer_reference(x_mm, kvp, material):
    """Archer transmission for overlay; None if the material has no table."""
    tbl = ARCHER_REF.get(material)
    if not tbl or x_mm <= 0:
        return None if not tbl else 1.0
    ks = sorted(tbl)
    kv = min(max(kvp, ks[0]), ks[-1])
    for i in range(len(ks) - 1):
        lo, hi = ks[i], ks[i + 1]
        if lo <= kv <= hi:
            f = (kv - lo) / (hi - lo)
            a = tuple(tbl[lo][j] + f * (tbl[hi][j] - tbl[lo][j])
                      for j in range(3))
            break
    else:
        a = tbl[ks[-1]]
    al, be, ga = a
    return _archer_transmission(x_mm, al, be, ga)


# -- Radionuclide Archer parameters ------------------------------------------
# The same Archer three-parameter model, but fitted directly to a
# radionuclide's emission spectrum rather than an x-ray tube's kVp - so
# there is no kVp to interpolate across, just one fixed (alpha, beta,
# gamma) per nuclide/material. Lead and concrete only.
#
# Ported as supplied, including its own unit tags and basis flags:
#   verified - matches a published mm^-1 fit, or the implied TVL matches
#              a published TVL
#   physics  - only one unit reading gives a physically possible TVL
#              (e.g. the 511 keV annihilation sources must shield like a
#              known ~15 mm Pb HVL source)
#   inferred - both readings were arguable; the one consistent with the
#              nuclide's photon energies was taken - WORTH RE-CHECKING
#              before relying on it for a real shielding decision.
#
# Only F-18, Ga-67, I-123, I-131, In-111 and Tc-99m overlap with the
# isotopes this app can actually run a Monte Carlo source for (see
# ISOTOPE_LINES above) - Ga-68, the Ge-68/Ga-68 generator, Lu-177 and
# Tl-201 are carried here for completeness but have no gamma-line data
# wired in, so no MC run (and therefore no comparison curve) is available
# for them yet.
#
# Entry format: (alpha, beta, gamma, unit, basis)
RADIONUCLIDES = {
    "F-18":     {"energies": "511 keV annihilation",
                "lead":     (1.543, -0.4408, 2.136, "cm", "verified"),
                "concrete": (0.1539, -0.1161, 2.0752, "cm", "verified")},
    "Ga-67":    {"energies": "93 / 185 / 300 / 393 keV",
                "lead":     (0.04412, 0.7326, 0.3869, "mm", "physics"),
                "concrete": (0.1865, -0.1529, 1.683, "cm", "physics")},
    "Ga-68":    {"energies": "511 keV annihilation",
                "lead":     (1.41, -0.402, 4.47, "cm", "physics"),
                "concrete": (0.134, -0.094, 3.2, "cm", "physics")},
    "Ge-68/Ga-68 gen": {"energies": "511 keV annihilation",
                "lead":     (1.41, -0.402, 4.47, "cm", "physics"),
                "concrete": (0.134, -0.094, 3.2, "cm", "physics")},
    "I-123":    {"energies": "159 keV + 529 keV at 1.4 %",
                "lead":     (0.1322, 14.85, 1.488, "mm", "inferred"),
                "concrete": (0.1866, 1253.0, 74.93, "cm", "inferred")},
    "I-131":    {"energies": "364 / 637 / 723 keV",
                "lead":     (-0.009, 2.439, 0.117, "cm", "verified"),
                "concrete": (0.152, -0.099, 2.634, "cm", "physics")},
    "In-111":   {"energies": "171 / 245 keV",
                "lead":     (0.7094, 30.63, 5.979, "mm", "inferred"),
                "concrete": (0.2339, -0.1265, 0.4618, "cm", "inferred")},
    "Lu-177":   {"energies": "113 / 208 keV",
                "lead":     (0.3855, 1.071, 0.2822, "mm", "verified"),
                "concrete": (0.2477, -0.1173, 0.4404, "cm", "verified")},
    "Tc-99m":   {"energies": "140 keV",
                "lead":     (2.558, 1.010, 4.344, "mm", "verified"),
                "concrete": (0.3102, -0.1729, 0.3622, "cm", "verified")},
    "Tl-201":   {"energies": "69-81 keV x-rays + 167 keV",
                "lead":     (1.67, 1.358, 0.3692, "mm", "inferred"),
                # alpha = -beta exactly, so 1 + beta/alpha = 0 and the
                # Archer equation divides by zero - flagged unusable in
                # the source sheet, nothing to load for concrete.
                "concrete": None},
}
NUCLIDE_MATERIALS = ["lead", "concrete"]


def nuclide_params(name, material):
    """(alpha, beta, gamma) in mm^-1, or None if unusable/unavailable."""
    entry = RADIONUCLIDES.get(name)
    if not entry:
        return None
    vals = entry.get(material)
    if not vals:
        return None
    a, b, g, unit, _basis = vals
    if unit == "cm":                     # per cm -> per mm
        a, b = a / 10.0, b / 10.0
    if a == 0 or abs(1.0 + b / a) < 1e-12:
        return None                      # degenerate fit
    return (a, b, g)


def nuclide_info(name, material):
    """(basis, unit) for display, or (None, None)."""
    vals = RADIONUCLIDES.get(name, {}).get(material)
    return (vals[4], vals[3]) if vals else (None, None)


def archer_reference_nuclide(x_mm, name, material):
    """
    Archer transmission for a radionuclide source; None if this
    nuclide/material combination has no usable fit. No kVp involved - a
    nuclide has one fixed spectrum, so there is nothing to interpolate.
    """
    p = nuclide_params(name, material)
    if p is None or x_mm <= 0:
        return None if p is None else 1.0
    al, be, ga = p
    return _archer_transmission(x_mm, al, be, ga)


# =============================================================================
#  SECTION 4 -- GUI
# =============================================================================

BG, CARD, ACCENT, OK, WARN = "#f4f6f8", "#ffffff", "#1f5f8b", "#1a7f4b", "#b3541e"
MONO = ("Consolas", 10)


def fmt(v, sig=4):
    if v is None:
        return "n/a"
    if isinstance(v, str):
        return v
    try:
        v = float(v)
    except (TypeError, ValueError):
        return str(v)
    if v != v:
        return "n/a"
    if math.isinf(v):
        return "inf"
    if v == 0:
        return "0"
    a = abs(v)
    if a >= 1e5 or a < 1e-3:
        return "%.*e" % (sig - 1, v)
    if a >= 100:
        return "%.1f" % v
    if a >= 10:
        return "%.2f" % v
    if a >= 1:
        return "%.3f" % v
    return "%.*f" % (sig + 1, v)


class LabelledEntry(ttk.Frame):
    def __init__(self, parent, label, default="", unit="", width=12, tip=""):
        super().__init__(parent)
        self.var = tk.StringVar(value=str(default))
        ttk.Label(self, text=label, width=24, anchor="w").pack(side="left")
        self.entry = ttk.Entry(self, textvariable=self.var, width=width,
                               justify="right")
        self.entry.pack(side="left")
        if unit:
            ttk.Label(self, text=" " + unit, width=9, anchor="w",
                      foreground="#666").pack(side="left")

    def get(self, default=None):
        t = self.var.get().strip()
        if not t:
            if default is None:
                raise ValueError("A required value is blank.")
            return default
        try:
            return float(t)
        except ValueError:
            raise ValueError("'%s' is not a number." % t)

    def set(self, v):
        self.var.set(str(v))


class LabelledCombo(ttk.Frame):
    def __init__(self, parent, label, values, default=None, width=15):
        super().__init__(parent)
        self.var = tk.StringVar(value=default if default else values[0])
        ttk.Label(self, text=label, width=24, anchor="w").pack(side="left")
        self.combo = ttk.Combobox(self, textvariable=self.var, values=values,
                                  width=width, state="readonly")
        self.combo.pack(side="left")

    def get(self):
        return self.var.get()


class PlotPane(ttk.Frame):
    def __init__(self, parent, figsize=(7.4, 4.4)):
        super().__init__(parent)
        self.fig = Figure(figsize=figsize, dpi=100, facecolor=CARD)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        try:
            tb = NavigationToolbar2Tk(self.canvas, self, pack_toolbar=False)
            tb.pack(side="bottom", fill="x")
        except TypeError:
            holder = ttk.Frame(self); holder.pack(side="bottom", fill="x")
            tb = NavigationToolbar2Tk(self.canvas, holder)
        tb.update()
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=True)

    def reset(self):
        self.fig.clf()
        return self.fig

    def done(self):
        try:
            self.fig.tight_layout()
        except Exception:
            pass
        self.canvas.draw_idle()


class ResultPane(ttk.Frame):
    def __init__(self, parent, title="Results"):
        super().__init__(parent)
        self._df = None
        bar = ttk.Frame(self); bar.pack(fill="x", pady=(0, 4))
        ttk.Label(bar, text=title, font=("Segoe UI", 10, "bold")).pack(side="left")
        ttk.Button(bar, text="Copy", width=8, command=self.copy).pack(side="right")
        self.btn = ttk.Button(bar, text="Export", width=9, command=self.export,
                              state="disabled")
        self.btn.pack(side="right", padx=4)
        wrap = ttk.Frame(self); wrap.pack(fill="both", expand=True)
        self.text = tk.Text(wrap, wrap="none", font=MONO, height=13,
                            background=CARD, relief="solid", borderwidth=1,
                            padx=10, pady=8)
        sb = ttk.Scrollbar(wrap, orient="vertical", command=self.text.yview)
        self.text.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y")
        self.text.pack(side="left", fill="both", expand=True)
        t = self.text
        t.tag_configure("h", font=("Consolas", 10, "bold"), foreground=ACCENT,
                        spacing1=8, spacing3=3)
        t.tag_configure("val", font=("Consolas", 10, "bold"))
        t.tag_configure("big", font=("Consolas", 13, "bold"), foreground=OK)
        t.tag_configure("note", foreground="#777", font=("Consolas", 9))
        t.tag_configure("warn", foreground=WARN, font=("Consolas", 10, "bold"))
        t.configure(state="disabled")

    def clear(self):
        self.text.configure(state="normal"); self.text.delete("1.0", "end")
        self.text.configure(state="disabled")
        self._df = None; self.btn.configure(state="disabled")

    def _ins(self, s, tag=None):
        self.text.configure(state="normal"); self.text.insert("end", s, tag)
        self.text.configure(state="disabled")

    def head(self, s):
        self._ins("\n" + s + "\n" + "-" * max(len(s), 48) + "\n", "h")

    def row(self, k, v, u="", w=30):
        self._ins("  " + str(k).ljust(w)); self._ins(str(v), "val")
        self._ins((" " + u if u else "") + "\n")

    def big(self, k, v, u=""):
        self._ins("  " + k.ljust(30)); self._ins(str(v) + (" " + u if u else ""), "big")
        self._ins("\n")

    def note(self, s):
        self._ins("  " + s + "\n", "note")

    def warn(self, s):
        self._ins("  ! " + s + "\n", "warn")

    def frame(self, df):
        self._df = df; self.btn.configure(state="normal")
        for line in df.to_string(index=False).splitlines():
            self._ins("  " + line + "\n")

    def top(self):
        self.text.see("1.0")

    def copy(self):
        self.clipboard_clear()
        self.clipboard_append(self.text.get("1.0", "end").rstrip())

    def export(self):
        if self._df is None:
            return
        p = filedialog.asksaveasfilename(defaultextension=".csv",
                                         initialfile="montecarlo.csv",
                                         filetypes=[("CSV", "*.csv"),
                                                    ("Excel", "*.xlsx")])
        if not p:
            return
        try:
            self._df.to_excel(p, index=False) if p.endswith(".xlsx") \
                else self._df.to_csv(p, index=False)
        except Exception as e:
            messagebox.showerror("Export failed", str(e))


def mat_key(label):
    """Map a display label back to its data key, e.g. PMMA -> pmma."""
    for k, v in MATERIAL_LABEL.items():
        if v == label:
            return k
    return str(label).lower()


def mat_labels(keys):
    return [MATERIAL_LABEL.get(k, k) for k in keys]


_ELEMENTS_BY_LABEL = sorted(ELEMENT_LABEL, key=lambda k: ELEMENT_LABEL[k])

# Both lists now share the same base set of curated materials - there is
# no physics reason PMMA/tissue/water/air/beryllium/borosilicate can't be
# used as a Monte Carlo barrier material even though the original design
# only offered them as beam filters. tungsten already had cross-section
# data (used for the target calc) but was never exposed as a selectable
# filter/barrier material either - fixed here too.
_BASE_MATERIALS = ["aluminium", "copper", "pmma", "beryllium", "borosilicate",
                  "glass", "water", "tissue", "air", "lead", "leadglass",
                  "concrete", "gypsum", "steel", "tungsten"]
FILTER_MATERIALS = _BASE_MATERIALS + _ELEMENTS_BY_LABEL
BARRIER_MATERIALS = _BASE_MATERIALS + _ELEMENTS_BY_LABEL


class SpectrumTab(ttk.Frame):
    title = "  Spectrum  "

    def __init__(self, parent):
        super().__init__(parent, padding=10)
        self._shape_set = False
        ttk.Label(self, text=(
            "Generate a tube spectrum and filter it. "
            "Filtering is Beer-Lambert per energy bin, which gives "
            "narrow-beam attenuation."),
            foreground="#555", wraplength=1200, justify="left").pack(
                fill="x", pady=(0, 8))
        body = ttk.Frame(self); body.pack(fill="both", expand=True)
        left = ttk.Frame(body); left.pack(side="left", fill="y", padx=(0, 12))

        s = ttk.LabelFrame(left, text="Tube", padding=8); s.pack(fill="x", pady=(0, 8))
        self.kvp = LabelledEntry(s, "Tube potential", 125, "kVp"); self.kvp.pack(fill="x", pady=2)
        self.angle = LabelledEntry(s, "Anode angle", 12, "deg"); self.angle.pack(fill="x", pady=2)
        self.engine = LabelledCombo(s, "Spectrum model",
                                    (["SpekPy"] if HAVE_SPEKPY else []) + ["Built-in (approx)"],
                                    "SpekPy" if HAVE_SPEKPY else "Built-in (approx)",
                                    width=18)
        self.engine.pack(fill="x", pady=2)

        s2 = ttk.LabelFrame(left, text="Total filtration (inherent + added)",
                            padding=8)
        s2.pack(fill="x", pady=(0, 8))
        ttk.Label(s2, text="2.5 mm Al is the regulatory minimum total\n"
                           "above 70 kV; typical general radiography is\n"
                           "2.5-3.5 mm. Tube labels usually quote total,\n"
                           "not inherent - inherent alone is 0.5-1.5 mm.",
                  foreground="#777", font=("Segoe UI", 8),
                  justify="left").pack(anchor="w", pady=(0, 6))
        self.filters = []
        for i, (mat, mm) in enumerate([("aluminium", 2.5), ("copper", 0.0),
                                       ("water", 0.0)]):
            row = ttk.Frame(s2); row.pack(fill="x", pady=1)
            c = LabelledCombo(row, "Filter %d" % (i + 1),
                              mat_labels(FILTER_MATERIALS),
                              MATERIAL_LABEL.get(mat, mat), width=12)
            c.pack(side="left")
            e = LabelledEntry(row, "", mm, "mm", width=7); e.pack(side="left")
            self.filters.append((c, e))

        s3 = ttk.LabelFrame(left, text="Calibration", padding=8)
        s3.pack(fill="x", pady=(0, 8))
        ttk.Label(s3, text="Measured on the tube. Enter either or\n"
                           "both, then press Calibrate once.",
                  foreground="#777", font=("Segoe UI", 8),
                  justify="left").pack(anchor="w", pady=(0, 6))
        self.hvl_mat = LabelledCombo(s3, "Measure HVL in",
                                     mat_labels(FILTER_MATERIALS),
                                     "Aluminium", width=12)
        self.hvl_mat.pack(fill="x", pady=2)
        self.hvl_meas = LabelledEntry(s3, "Measured HVL", "", "mm")
        self.hvl_meas.pack(fill="x", pady=2)
        self.meas_kvp = LabelledEntry(s3, "measured at", "", "kVp",
                                      tip="Blank uses the tube potential "
                                          "above.")
        self.meas_kvp.pack(fill="x", pady=2)
        self.meas = LabelledEntry(s3, "Measured output", "", "uGy/mAs")
        self.meas.pack(fill="x", pady=2)
        self.meas_d = LabelledEntry(s3, "at distance", 1.0, "m")
        self.meas_d.pack(fill="x", pady=2)
        brow = ttk.Frame(s3); brow.pack(fill="x", pady=(7, 0))
        ttk.Button(brow, text="Calibrate", command=self.calibrate).pack(
            side="left", expand=True, fill="x", padx=(0, 3))
        ttk.Button(brow, text="Clear", command=self.uncalibrate).pack(
            side="left", expand=True, fill="x", padx=(3, 0))
        self.callbl = ttk.Label(s3, text="not calibrated", foreground="#777",
                                font=("Segoe UI", 8), justify="left",
                                wraplength=210)
        self.callbl.pack(anchor="w", pady=(6, 0))
        ttk.Label(s3, text="HVL sets the shape by solving filtration.\n"
                           "Output sets the scale. Both are properties of\n"
                           "the tube, so once solved at the kVp\n"
                           "measured at, the calibration applies at every other kVp.",
                  foreground="#777", font=("Segoe UI", 8),
                  justify="left").pack(anchor="w", pady=(4, 0))




        ttk.Button(left, text="Generate", style="Accent.TButton",
                   command=self.run).pack(fill="x", ipady=4, pady=(8, 0))
        ttk.Button(left, text="Sweep kVp vs Archer 1994",
                   command=self.sweep).pack(fill="x", ipady=3, pady=(4, 0))

        s5 = ttk.LabelFrame(left, text="RQA beam qualities", padding=8)
        s5.pack(fill="x", pady=(10, 0))
        self.quality = LabelledCombo(s5, "Quality",
                                     sorted(RQA_ED2.keys(),
                                            key=lambda q: RQA_ED2[q][0]),
                                     "RQA5", width=10)
        self.quality.pack(fill="x", pady=2)
        qrow = ttk.Frame(s5); qrow.pack(fill="x", pady=(6, 0))
        ttk.Button(qrow, text="Identify beam",
                   command=self.rqa_identify).pack(side="left", expand=True,
                                                   fill="x", padx=(0, 3))
        ttk.Button(qrow, text="Solve for target",
                   command=self.rqa_solve).pack(side="left", expand=True,
                                                fill="x", padx=(3, 0))
        ttk.Label(s5, text="IEC 61267 Ed.2 (2005). Qualities are\n"
                           "defined by HVL, not by the nominal\n"
                           "added filtration.",
                  foreground="#777", font=("Segoe UI", 8),
                  justify="left").pack(anchor="w", pady=(5, 0))

        pane = ttk.PanedWindow(body, orient="vertical"); pane.pack(side="left", fill="both", expand=True)
        self.results = ResultPane(pane); self.plot = PlotPane(pane)
        pane.add(self.results, weight=2); pane.add(self.plot, weight=3)
        self.spectrum = None

    def current_filters(self):
        out = []
        for c, e in self.filters:
            try:
                mm = e.get(0.0)
            except ValueError:
                mm = 0.0
            if mm > 0:
                out.append((mat_key(c.get()), mm))
        return out


    def _cal_label(self):
        if CALIBRATION["scale"] is None and not self._shape_set:
            self.callbl.configure(text="not calibrated", foreground="#777")
            return
        bits = []
        if self._shape_set:
            bits.append("shape")
        if CALIBRATION["scale"] is not None:
            bits.append("scale")
        self.callbl.configure(
            text="%s set at %s kVp" % (" and ".join(bits), fmt(self.kvp.get(0))),
            foreground=OK)

    def calibrate(self):
        """
        One action, two jobs.

        HVL fixes the SHAPE: filtration is solved so the model reproduces the
        measured HVL. Output fixes the SCALE: a factor is derived so the
        model reports real uGy.

        Shape is solved FIRST, because the scale factor is derived against
        whichever spectrum is loaded - deriving it before correcting the
        filtration would attach it to the wrong beam. Doing both in one
        action means that ordering cannot be got wrong.
        """
        r = self.results
        try:
            kvp_show = self.kvp.get()
            # HVL and output are measured at a particular tube potential.
            # Filtration and the scale factor are properties of the TUBE, so
            # they must be derived at the kVp they were measured at, then
            # applied at whatever potential is being simulated.
            kvp = self.meas_kvp.get(0.0) or kvp_show
            ang = self.angle.get(12.0)
            hvl_t = self.hvl_meas.get(0.0)
            out_t = self.meas.get(0.0)
            dist = self.meas_d.get(1.0)
            use_spek = self.engine.get().startswith("SpekPy") and HAVE_SPEKPY

            if hvl_t <= 0 and out_t <= 0:
                raise ValueError("Enter a measured HVL, a measured output, "
                                 "or both.")

            report = []          # collected now, written AFTER the rebuild,
                                 # because run() clears the results pane

            # ---- 1. shape -------------------------------------------
            shape_changed = False
            if hvl_t > 0:
                mm, ach = solve_filtration_from_hvl(kvp, hvl_t, ang, use_spek)
                if mm is None:
                    report.append(("warn",
                        "Cannot reach an HVL of %s mm. With no filtration at "
                        "all the model already gives %s mm, so it is making a "
                        "harder beam than your tube." % (fmt(hvl_t), fmt(ach))))
                    report.append(("note",
                        "On the built-in model, install spekpy. On SpekPy, "
                        "check the anode angle."))
                    self._shape_set = False
                else:
                    self.filters[0][0].var.set(MATERIAL_LABEL["aluminium"])
                    self.filters[0][1].set(round(mm, 3))
                    for _c, e in self.filters[1:]:
                        e.set(0.0)
                    self._shape_set = True
                    shape_changed = True
                    report.append(("row", ("Measured HVL", fmt(hvl_t), "mm")))
                    report.append(("big", ("Total filtration", fmt(mm),
                                           "mm Al equivalent")))
                    report.append(("row", ("Reproduced HVL", fmt(ach), "mm")))
                    report.append(("row", ("Deviation",
                                   "%+.3f %%" % (100 * (ach - hvl_t) / hvl_t),
                                   "")))
            else:
                report.append(("note",
                    "No HVL given - filtration left as entered above, so the "
                    "spectral shape is only as good as that figure."))

            # a re-solved shape invalidates any scale taken on the old beam
            if shape_changed and out_t <= 0 and CALIBRATION["scale"] is not None:
                clear_calibration()
                report.append(("warn",
                    "The filtration changed, so the previous scale factor "
                    "referred to a different beam and has been cleared. "
                    "Re-enter your measured output to restore it."))

            # derive the scale against the spectrum AS MEASURED, which may
            # not be the one currently on display
            meas_sp = make_spectrum(kvp, self.current_filters(), ang,
                                    use_spek)
            self.run()           # rebuild the displayed spectrum

            # ---- 2. scale -------------------------------------------
            if out_t > 0:
                scale = set_calibration(meas_sp, out_t, dist)
                self.run()
                report.append(("row", ("Measured output",
                               "%s uGy/mAs at %s m" % (fmt(out_t), fmt(dist)),
                               "")))
                report.append(("row", ("Scale factor", fmt(scale, 6), "")))
                if CALIBRATION["spek_ratio"]:
                    report.append(("row", ("Correction to SpekPy",
                                   "x %s" % fmt(CALIBRATION["spek_ratio"], 4),
                                   "")))
                    report.append(("note",
                        "A factor near 1 means SpekPy already agreed with "
                        "your measurement."))
            else:
                report.append(("note",
                    "No output given - beam quality is calibrated but dose "
                    "cannot be quoted in uGy."))
            self._cal_label()

            # ---- write the report, after the rebuild ------------------
            r.head("Calibration from measurements at %s kVp" % fmt(kvp))
            if abs(kvp - kvp_show) > 1e-9:
                r.note("Tube potential on display is %s kVp. Filtration and "
                       "scale were derived at %s kVp and carried across, "
                       "which is what makes them tube properties rather than "
                       "settings." % (fmt(kvp_show), fmt(kvp)))
            for kind, payload in report:
                if kind == "row":
                    r.row(payload[0], payload[1], payload[2])
                elif kind == "big":
                    r.big(payload[0], payload[1], payload[2])
                elif kind == "note":
                    r.note(payload)
                else:
                    r.warn(payload)

            r.head("What was set")
            r.row("Shape", "from measured HVL" if self._shape_set
                  else "from the filtration entered above")
            r.row("Scale", "from measured output"
                  if CALIBRATION["scale"] is not None else "not set")
            r.head("How to check it")
            r.note("Change the kVp and press Generate. Both the HVL and the "
                   "output should MOVE. If the output stays put, the scale "
                   "is not being applied.")
            r.note("Add filtration and press Generate. HVL should rise and "
                   "output should fall.")
            r.note("Press Sweep to see the whole kVp range on your scale.")
            r.top()
        except ValueError as e:
            r.clear(); r.head("Calibration error"); r.warn(str(e))

    def uncalibrate(self):
        clear_calibration()
        self._shape_set = False
        self._cal_label()
        self.run()

    def _current_spectrum(self):
        if self.spectrum is None:
            self.run()
        if self.spectrum is None:
            raise ValueError("Generate a spectrum first.")
        return self.spectrum

    def rqa_identify(self):
        r = self.results; r.clear()
        try:
            sp = self._current_spectrum()
            ed = RQA_EDITION
            rows = rqa_match(sp)
            hvl = sp.hvl("aluminium")
            r.head("Your beam against %s" % ed)
            r.row("Tube potential", fmt(sp.kvp), "kVp")
            r.row("Filtration", ", ".join("%s %s mm"
                                          % (MATERIAL_LABEL.get(m, m), fmt(t))
                                          for m, t, _ in sp.filters) or "none")
            r.big("Measured HVL", fmt(hvl), "mm Al")
            r.row("Mean energy", fmt(sp.mean_energy), "keV")
            best = rows[0]
            r.head("Closest match")
            if abs(best["deviation_pc"]) < 5:
                r.big("Nearest quality", best["quality"])
                r.row("Target HVL", fmt(best["target_HVL_mm"]), "mm Al")
                r.row("Deviation", "%+.2f %%" % best["deviation_pc"])
            else:
                r.row("Nearest quality", "%s (%.1f %% off)"
                      % (best["quality"], best["deviation_pc"]))
                r.warn("Your beam is not close to any RQA quality. RQA beams "
                       "are heavily hardened - RQA5 needs about 21 mm of "
                       "added aluminium on top of a normal tube.")
            if best["confidence"] != "verified":
                r.warn("The %s entry could not be confirmed from an "
                       "accessible published source. Check it against the "
                       "standard before relying on it." % best["quality"])
            if best["kV"] != sp.kvp:
                r.note("Note the kV differs: %s is specified at %d kV and "
                       "your spectrum is at %s kVp. HVL alone does not make "
                       "a beam an RQA quality - the tube potential must "
                       "match too." % (best["quality"], best["kV"],
                                       fmt(sp.kvp)))
            df = pd.DataFrame(rows)
            r.head("All qualities, sorted by closeness")
            r.frame(df)
            r._df = df

            # ---- fluence / air-kerma (IEC 62220-1 Table 3) -------------
            r.head("Photon fluence per air kerma (IEC 62220-1 Table 3)")
            q, frows = rqa_fluence_match(sp)
            if q is None:
                r.note("Not available - this spectrum has no absolute "
                       "scale. Press Calibrate with a measured output, or "
                       "use SpekPy's own normalisation, to get a real "
                       "photons/mm2/uGy figure. HVL alone does not fix "
                       "this: two beams can share an RQA's target HVL "
                       "while differing in Q/K_air, since the standard "
                       "constrains the beam by HVL and kV, not by fluence.")
            else:
                r.big("Your Q/K_air (SNR_in^2)", fmt(q, 1), "mm^-2.uGy^-1")
                fbest = frows[0]
                r.row("Closest mandatory quality", "%s (%+.2f %%)"
                      % (fbest["quality"], fbest["deviation_pc"]))
                r.note("Only RQA3/5/7/9 carry a mandatory value in Table 3 "
                       "- this is not compared against RQA2/4/6/8/10.")
                fdf = pd.DataFrame(frows)
                r.frame(fdf)
            r.top()

            fig = self.plot.reset(); ax = fig.add_subplot(111)
            names = [d["quality"] for d in sorted(rows, key=lambda d: d["kV"])]
            tgts = [d["target_HVL_mm"] for d in sorted(rows,
                                                       key=lambda d: d["kV"])]
            kvs = [d["kV"] for d in sorted(rows, key=lambda d: d["kV"])]
            ax.plot(kvs, tgts, "s-", color=WARN, lw=1.6, ms=6,
                    label="%s targets" % ed)
            for k, t, n in zip(kvs, tgts, names):
                ax.annotate(n, (k, t), fontsize=7, xytext=(3, -10),
                            textcoords="offset points", color=WARN)
            ax.plot([sp.kvp], [hvl], "o", color=ACCENT, ms=11, zorder=5,
                    label="your beam")
            ax.grid(alpha=0.25); ax.set_axisbelow(True)
            ax.set_xlabel("tube potential (kV)", fontsize=9)
            ax.set_ylabel("HVL (mm Al)", fontsize=9)
            ax.set_title("Beam quality against the RQA series", fontsize=10,
                         color=ACCENT, fontweight="bold")
            ax.legend(fontsize=8); ax.tick_params(labelsize=8)
            self.plot.done()
        except ValueError as e:
            r.head("Input error"); r.warn(str(e))

    def rqa_solve(self):
        r = self.results; r.clear()
        try:
            sp = self._current_spectrum()
            ed = RQA_EDITION
            q = self.quality.get()
            kv, nominal, target, conf = RQA_ED2[q]
            mm, achieved, tgt = rqa_solve_filtration(sp, q)
            mm_al_cu, achieved_cu, _ = rqa_solve_filtration_cu(sp, q)

            r.head("Reaching %s from your beam" % q)
            r.row("Standard", ed)
            r.row("Specified at", "%d kV" % kv)
            r.row("Target HVL", fmt(tgt), "mm Al")
            r.row("Nominal added Al", fmt(nominal), "mm")
            r.row("Your starting HVL", fmt(sp.hvl("aluminium")), "mm Al")
            if mm is None:
                r.warn("Pure Al: cannot reach the target. Your beam is "
                       "either already harder than the target, or 80 mm "
                       "of aluminium is not enough. Check the tube "
                       "potential matches.")
            else:
                r.big("Add (pure Al)", fmt(mm), "mm Al")
                r.row("Achieved HVL", fmt(achieved), "mm Al")
                r.row("Deviation from target",
                      "%+.3f %%" % (100 * (achieved - tgt) / tgt))
                r.row("Versus nominal", "%+.2f mm" % (mm - nominal))
                r.note("The nominal figure assumes the standard's reference "
                       "tube. Yours differs, which is exactly why the "
                       "standard specifies HVL rather than filtration.")
            r.head("Alternative: 0.5 mm Cu + Al")
            if mm_al_cu is None:
                r.warn("0.5 mm Cu alone already exceeds the target HVL - "
                       "no amount of added Al can bring it back down. This "
                       "Cu thickness cannot reach %s from your beam." % q)
            else:
                r.big("Add (on top of 0.5 mm Cu)", fmt(mm_al_cu), "mm Al")
                r.row("Achieved HVL", fmt(achieved_cu), "mm Al")
                r.row("Deviation from target",
                      "%+.3f %%" % (100 * (achieved_cu - tgt) / tgt))
                r.note("Same target HVL, different filter composition. A "
                       "beam that matches by HVL alone is not guaranteed "
                       "to match in photon fluence per unit air kerma - "
                       "see the fluence comparison in Identify beam. This "
                       "is offered as a genuinely different route to the "
                       "same HVL, not a more accurate version of the "
                       "pure-Al answer above.")
            if sp.kvp != kv:
                r.warn("Your spectrum is at %s kVp but %s is specified at "
                       "%d kV. Matching HVL alone does not reproduce the "
                       "quality." % (fmt(sp.kvp), q, kv))
            if conf != "verified":
                r.warn("This entry is unverified - confirm against the "
                       "standard before use.")
            r.top()

            fig = self.plot.reset(); ax = fig.add_subplot(111)
            ax.plot(sp.E, sp.N / max(float(np.max(sp.N)), 1e-30),
                    color="#999", lw=1.4, label="your beam")
            if mm is not None:
                hard = sp.filtered("aluminium", mm)
                ax.plot(hard.E, hard.N / max(float(np.max(hard.N)), 1e-30),
                        color=ACCENT, lw=2.0,
                        label="%s  (+%s mm Al)" % (q, fmt(mm)))
                ax.axvline(hard.mean_energy, color=OK, ls="--", lw=1.0,
                           label="mean %.1f keV" % hard.mean_energy)
            if mm_al_cu is not None:
                hard_cu = sp.filtered("copper", 0.5).filtered("aluminium",
                                                              mm_al_cu)
                ax.plot(hard_cu.E, hard_cu.N
                        / max(float(np.max(hard_cu.N)), 1e-30),
                        color="#c1666b", lw=2.0, ls="--",
                        label="%s  (0.5 mm Cu + %s mm Al)"
                        % (q, fmt(mm_al_cu)))
            ax.set_xlim(0, sp.kvp * 1.02); ax.set_ylim(0, 1.15)
            ax.grid(alpha=0.25); ax.set_axisbelow(True)
            ax.set_xlabel("photon energy (keV)", fontsize=9)
            ax.set_ylabel("relative fluence per keV", fontsize=9)
            ax.set_title("Hardening to %s" % q, fontsize=10, color=ACCENT,
                         fontweight="bold")
            ax.legend(fontsize=8); ax.tick_params(labelsize=8)
            self.plot.done()
        except ValueError as e:
            r.head("Input error"); r.warn(str(e))

    def draw_spectrum(self, unf, sp, kvp):
        """
        Plot the beam as it actually is.

        Showing the unfiltered spectrum alongside it does not work: an
        unfiltered spectrum rises steeply toward low energy, so its peak is
        both far to the left and far larger, and putting the two on one
        scale flattens the real beam into nothing. Rescaling the unfiltered
        curve to match would misrepresent it, and a log axis buries the
        shape under the low-energy tail.

        The unfiltered spectrum is a fiction in any case - no tube emits it.
        What the filtration did is reported as numbers in the results pane,
        where it can be stated precisely instead of implied by a curve.
        """
        fig = self.plot.reset()
        ax = fig.add_subplot(111)
        peak = float(np.max(sp.N)) if sp.N.size else 1.0

        ax.fill_between(sp.E, sp.N, color=ACCENT, alpha=0.22)
        ax.plot(sp.E, sp.N, color=ACCENT, lw=1.9)
        ax.axvline(sp.mean_energy, color=OK, ls="--", lw=1.1)
        ax.annotate("mean %.1f keV" % sp.mean_energy,
                    (sp.mean_energy, peak * 0.94), fontsize=8, color=OK,
                    xytext=(5, 0), textcoords="offset points")

        filt = ", ".join("%s %s mm" % (MATERIAL_LABEL.get(mm_, mm_), fmt(t))
                         for mm_, t, _ in sp.filters) or "no filtration"
        ax.set_xlim(0, kvp)
        ax.set_ylim(0, peak * 1.15 if peak > 0 else 1)
        ax.grid(alpha=0.25); ax.set_axisbelow(True)
        ax.set_xlabel("photon energy (keV)", fontsize=9)
        ax.set_ylabel("fluence per keV", fontsize=9)
        ax.set_title("%s kVp, %s" % (fmt(kvp), filt), fontsize=10,
                     color=ACCENT, fontweight="bold")
        ax.tick_params(labelsize=8)
        self.plot.done()

    def sweep(self):
        """
        Hold filtration fixed and sweep kVp, comparing output against the
        fit to Archer 1994 and its tabulated anchor values.

        These are ONE source, not two - the curve is the fit to the data,
        so they are plotted as a single reference rather than as if they
        were independent checks.

        This is the test that means something. Matching output at ONE kVp by
        adjusting filtration is guaranteed - one free parameter against one
        target. Reproducing the SHAPE across kVp with nothing further
        adjusted is not guaranteed, and the ratio K(125)/K(100) is
        independent of absolute normalisation, so it isolates filtration
        from anode angle and generator waveform.
        """
        r = self.results; r.clear()
        try:
            ang = self.angle.get(12.0)
            fl = self.current_filters()
            use_spek = self.engine.get().startswith("SpekPy") and HAVE_SPEKPY
            kvs = list(range(50, 155, 5))
            rows = []
            for kv in kvs:
                base = (spectrum_spekpy(kv, ang, ()) if use_spek
                        else spectrum_kramers(kv, (), target_angle_deg=ang))
                s2 = base
                for mat, mm in fl:
                    s2 = s2.filtered(mat, mm)
                k = spectrum_output(s2)
                rows.append({
                    "kVp": kv,
                    "K_app_mGy_per_mAmin": (k * 0.06) if k is not None else None,
                    "K_Archer_mGy_per_mAmin": Spectrum.archer_K_poly(kv),
                    "HVL_Al_mm": s2.hvl("aluminium"),
                    "mean_keV": s2.mean_energy})
            df = pd.DataFrame(rows)
            have_abs = df["K_app_mGy_per_mAmin"].notna().all()

            r.head("Output against kVp, filtration held fixed")
            r.row("Filtration", ", ".join("%s %s mm" % (MATERIAL_LABEL.get(m, m),
                                                        fmt(t))
                                          for m, t in fl) or "none")
            r.row("Anode angle", fmt(ang), "deg")
            if not have_abs:
                r.warn("No absolute scale, so only HVL and mean energy are "
                       "meaningful here. Install spekpy to compare output.")
            else:
                df["ratio"] = (df["K_app_mGy_per_mAmin"]
                               / df["K_Archer_mGy_per_mAmin"])
                k100 = float(df.loc[df.kVp == 100, "K_app_mGy_per_mAmin"].iloc[0])
                k125 = float(df.loc[df.kVp == 125, "K_app_mGy_per_mAmin"].iloc[0])
                r.head("The test that actually discriminates")
                r.row("app K(125)/K(100)", fmt(k125 / k100, 4))
                r.row("Archer tabulated 7.170/4.722", fmt(7.170 / 4.722, 4))
                r.row("Archer fitted ratio",
                      fmt(Spectrum.archer_K_poly(125)
                          / Spectrum.archer_K_poly(100), 4))
                r.note("This ratio is independent of absolute normalisation, "
                       "so it is insensitive to anode angle and generator "
                       "waveform. Matching it is evidence about FILTRATION. "
                       "Matching a single output value is not - that is one "
                       "free parameter fitted to one number.")
                r.row("Spread of app/Archer ratio",
                      "%s to %s" % (fmt(df["ratio"].min(), 4),
                                    fmt(df["ratio"].max(), 4)))
                r.note("A flat ratio near 1 across the whole range means the "
                       "filtration is right. A ratio that trends with kVp "
                       "means it is compensating for something else.")
                r.note("The reference is a fit to the measured data in Archer "
                       "1994. NCRP 147 carries the same coefficients with one "
                       "sign wrong, so the paper is the source used here.")

            r.head("Sweep")
            disp = df.copy()
            for c in disp.columns:
                if c != "kVp":
                    disp[c] = disp[c].map(lambda v: fmt(v, 4)
                                          if v is not None else "-")
            r.frame(disp)
            r._df = df
            r.head("Archer 1994 tabulated anchor values")
            r.row("100 kVp", "4.722", "mGy per mA.min at 1 m")
            r.row("125 kVp", "7.170", "mGy per mA.min at 1 m")
            r.note("These sit alongside the lead parameters a=2.500, "
                   "b=15.28, g=0.7557 - the same pair your shielding "
                   "calculator uses at 100 kVp. The curve above is the "
                   "fit to this data, so the two are one source, not two.")
            r.top()

            fig = self.plot.reset()
            ax = fig.add_subplot(121)
            if have_abs:
                ax.plot(df.kVp, df["K_app_mGy_per_mAmin"], "o-", color=ACCENT,
                        lw=1.8, ms=4, label="this app")
            ax.plot(df.kVp, df["K_Archer_mGy_per_mAmin"], "s--", color=WARN,
                    lw=1.5, ms=4, label="Archer 1994 (fitted)")
            ax.grid(alpha=0.25); ax.set_axisbelow(True)
            ax.set_xlabel("kVp", fontsize=9)
            ax.set_ylabel("mGy per mA.min at 1 m", fontsize=9)
            ax.set_title("Output", fontsize=10, color=ACCENT,
                         fontweight="bold")
            ax.legend(fontsize=7); ax.tick_params(labelsize=8)

            ax2 = fig.add_subplot(122)
            if have_abs:
                ax2.plot(df.kVp, df["ratio"], "o-", color=ACCENT, lw=1.8, ms=4)
                ax2.axhline(1.0, color="#c33", ls="--", lw=1.0)
                ax2.set_ylabel("app / Archer", fontsize=9)
                ax2.set_title("Ratio - flat is good", fontsize=10,
                              color=ACCENT, fontweight="bold")
            else:
                ax2.plot(df.kVp, df["HVL_Al_mm"], "o-", color=OK, lw=1.8, ms=4)
                ax2.set_ylabel("HVL in Al (mm)", fontsize=9)
                ax2.set_title("Beam quality", fontsize=10, color=ACCENT,
                              fontweight="bold")
            ax2.grid(alpha=0.25); ax2.set_axisbelow(True)
            ax2.set_xlabel("kVp", fontsize=9); ax2.tick_params(labelsize=8)
            self.plot.done()
        except ValueError as e:
            r.head("Input error"); r.warn(str(e))

    def run(self):
        r = self.results; r.clear()
        try:
            kvp = self.kvp.get(); ang = self.angle.get(12.0)
            fl = self.current_filters()
            use_spek = self.engine.get().startswith("SpekPy") and HAVE_SPEKPY
            unf = (spectrum_spekpy(kvp, ang, ()) if use_spek
                   else spectrum_kramers(kvp, (), target_angle_deg=ang))
            sp = unf
            for mat, mm in fl:
                sp = sp.filtered(mat, mm)
            self.spectrum = sp
            hm = mat_key(self.hvl_mat.get())

            r.head("Spectrum at %s kVp" % fmt(kvp))
            r.row("Model", sp.model)
            if not use_spek and HAVE_SPEKPY:
                r.note("SpekPy is installed - select it for better accuracy.")
            if not HAVE_SPEKPY:
                r.warn("SpekPy not installed. The built-in model reproduces "
                       "the shape and the tungsten K lines but is NOT "
                       "quantitative. pip install spekpy.")
            r.row("Filtration", ", ".join("%s %s mm" % (MATERIAL_LABEL.get(m, m), fmt(t))
                                          for m, t in fl) or "none")
            r.blank_line = None
            r.head("Beam quality")
            r.big("Mean energy", fmt(sp.mean_energy), "keV")
            r.row("Effective energy", fmt(sp.effective_energy(hm)), "keV")
            r.row("HVL 1 in %s" % MATERIAL_LABEL.get(hm, hm),
                  fmt(sp.hvl(hm, 0.5)), "mm")
            r.row("HVL 2 (0.5 -> 0.25)",
                  fmt(sp.hvl(hm, 0.25) - sp.hvl(hm, 0.5)), "mm")
            r.row("Homogeneity coefficient", fmt(sp.homogeneity(hm)))
            r.note("HVL is air-kerma defined, matching what an ion chamber "
                   "or a RaySafe measures.")
            r.row("Photon fluence (rel.)", fmt(sp.fluence))
            r.row("Energy fluence (rel.)", fmt(sp.energy_fluence))
            if unf.fluence > 0 and fl:
                r.head("What the filtration removed")
                r.row("Photons", "%.1f %%" % (100 * (1 - sp.fluence / unf.fluence)))
                k0, k1 = unf.kerma_integral(), sp.kerma_integral()
                if k0 > 0:
                    r.row("Air kerma", "%.1f %%" % (100 * (1 - k1 / k0)))
                r.row("Mean energy raised",
                      "%s -> %s keV" % (fmt(unf.mean_energy),
                                        fmt(sp.mean_energy)))
                r.note("The plot shows the beam as it is. An unfiltered "
                       "spectrum peaks far to the left and far higher, so "
                       "drawing both on one scale hides the real beam - and "
                       "no tube emits an unfiltered spectrum anyway.")

            # ---- air kerma output ----------------------------------
            r.head("Output at 1 m")
            k_uGy_mAs = spectrum_output(sp)

            if k_uGy_mAs is None:
                r.warn("No absolute scale. The built-in spectrum model is in "
                       "arbitrary units, so dose cannot be quoted. Either "
                       "install spekpy, or measure your tube output once and "
                       "press Calibrate - one measurement fixes the scale at "
                       "every kVp and filtration.")
                r.note("Beam quality figures above are unaffected - HVL, "
                       "mean and effective energy are all ratios and need "
                       "no absolute scale.")
            else:
                k_mGy_mAmin = k_uGy_mAs * 60.0 / 1000.0
                r.row("Basis", output_basis())
                if CALIBRATION["scale"] is not None:
                    if CALIBRATION["kvp"] != sp.kvp:
                        r.note("Calibrated at %s kVp, now showing %s kVp - "
                               "the model is carrying the scale across, "
                               "which is the point of calibrating."
                               % (fmt(CALIBRATION["kvp"]), fmt(sp.kvp)))
                    if CALIBRATION["spek_ratio"]:
                        r.row("Correction to SpekPy",
                              "x %s" % fmt(CALIBRATION["spek_ratio"], 4))
                        r.note("A factor near 1 means SpekPy already agreed "
                               "with your measurement.")
                r.big("Air kerma", fmt(k_uGy_mAs), "uGy per mAs at 1 m")
                r.big("Same as", fmt(k_mGy_mAmin), "mGy per mA.min at 1 m")
                ncrp = Spectrum.archer_K_poly(sp.kvp)
                r.row("Archer 1994", fmt(ncrp), "mGy per mA.min")
                if ncrp > 0:
                    diff = 100.0 * (k_mGy_mAmin - ncrp) / ncrp
                    r.row("Difference", "%+.1f %%" % diff)
                r.note("Archer measured a particular beam quality, so "
                       "agreement depends on matching its filtration. A "
                       "difference of tens of per cent is normal and does "
                       "not imply either is wrong.")
                r.note("K = SUM Phi(E) . E . (mu_en/rho)_air(E), with air "
                       "coefficients from the NIST Hubbell and Seltzer "
                       "tables. No Monte Carlo is involved - this is an "
                       "integral over the spectrum.")

            df = pd.DataFrame({"E_keV": sp.E,
                               "unfiltered": unf.N,
                               "filtered": sp.N})
            r.head("Spectrum data")
            r.frame(df.iloc[::max(1, len(df) // 25)])
            r.note("Full spectrum available via Export.")
            r._df = df
            r.top()

            self.draw_spectrum(unf, sp, kvp)
        except ValueError as e:
            r.head("Input error"); r.warn(str(e))


class MonteCarloTab(ttk.Frame):
    title = "  Monte Carlo  "

    figsize = (9.4, 4.2)

    def __init__(self, parent, spectrum_tab):
        super().__init__(parent, padding=10)
        self.spectrum_tab = spectrum_tab
        ttk.Label(self, text=(
            "Monte Carlo photon transport through a slab. "),
            foreground="#555", wraplength=1200, justify="left").pack(
                fill="x", pady=(0, 8))
        body = ttk.Frame(self); body.pack(fill="both", expand=True)
        left = ttk.Frame(body); left.pack(side="left", fill="y", padx=(0, 12))

        s0 = ttk.LabelFrame(left, text="Source", padding=8)
        s0.pack(fill="x", pady=(0, 8))
        self.source = LabelledCombo(s0, "Source", ["X-ray spectrum", "Isotope"],
                                    "X-ray spectrum", width=16)
        self.source.pack(fill="x", pady=2)
        self.isotope_row = ttk.Frame(s0)
        self.isotope = LabelledCombo(self.isotope_row, "Isotope",
                                     sorted(ISOTOPE_LINES), "Cs-137", width=16)
        self.isotope.pack(fill="x", pady=2)
        self.source.combo.bind("<<ComboboxSelected>>", self._source_changed)
        self.isotope.combo.bind("<<ComboboxSelected>>",
                               self._update_archer_availability)

        s = ttk.LabelFrame(left, text="Barrier", padding=8); s.pack(fill="x", pady=(0, 8))
        self.material = LabelledCombo(s, "Material",
                                      mat_labels(BARRIER_MATERIALS),
                                      "Lead", width=12)
        self.material.pack(fill="x", pady=2)
        self.material.combo.bind("<<ComboboxSelected>>",
                                 self._update_archer_availability)
        self.tmin = LabelledEntry(s, "Thickness from", 0.25, "mm"); self.tmin.pack(fill="x", pady=2)
        self.tmax = LabelledEntry(s, "Thickness to", 3.0, "mm"); self.tmax.pack(fill="x", pady=2)
        self.nsteps = LabelledEntry(s, "Number of steps", 7, ""); self.nsteps.pack(fill="x", pady=2)
        self.density = LabelledEntry(s, "Density override", "", "g/cm3"); self.density.pack(fill="x", pady=2)

        s2 = ttk.LabelFrame(left, text="Simulation", padding=8); s2.pack(fill="x", pady=(0, 8))
        self.nhist = LabelledEntry(s2, "Histories per point",
                                   1000000 if HAVE_NUMBA else 20000, "")
        self.nhist.pack(fill="x", pady=2)
        ttk.Label(s2, text="Coherent scatter uses real atomic form factors,\n"
                           "so it is computed, not bracketed.",
                  foreground="#777", font=("Segoe UI", 8),
                  justify="left").pack(anchor="w", pady=(2, 4))
        self.score = LabelledCombo(s2, "Score", ["Air kerma", "Energy fluence"],
                                   "Air kerma", width=24)
        self.score.pack(fill="x", pady=2)
        self.accept = LabelledEntry(s2, "Detector acceptance angle", 90, "deg")
        self.accept.pack(fill="x", pady=2)
        self.ethr = LabelledEntry(s2, "Detector energy threshold", 0, "keV")
        self.ethr.pack(fill="x", pady=2)
        self.show_mu = tk.BooleanVar(value=True)
        ttk.Checkbutton(s2, text="overlay attenuation coefficients",
                        variable=self.show_mu).pack(anchor="w", pady=(4, 0))
        self.show_archer = tk.BooleanVar(value=True)
        self.archer_chk = ttk.Checkbutton(s2, text="overlay Archer reference",
                                          variable=self.show_archer)
        self.archer_chk.pack(anchor="w", pady=(4, 0))

        if not HAVE_NUMBA:
            ttk.Label(left, text="numba not installed - the transport loop\n"
                                 "runs in pure Python and is roughly 50x\n"
                                 "slower. pip install numba.",
                      foreground=WARN, font=("Segoe UI", 8),
                      justify="left").pack(fill="x", pady=(0, 6))

        self.btn = ttk.Button(left, text="Run", style="Accent.TButton",
                              command=self.run)
        self.btn.pack(fill="x", ipady=4)
        self.prog = ttk.Progressbar(left, mode="determinate", maximum=100)
        self.prog.pack(fill="x", pady=(6, 0))

        sl = ttk.LabelFrame(left, text="Spectrum shown", padding=8)
        sl.pack(fill="x", pady=(10, 0))
        # tk.Scale rather than ttk.Scale: it takes a resolution, so the
        # handle snaps to the simulated thicknesses instead of sliding
        # through fractional positions that round to the same one.
        self.slider_var = tk.IntVar(value=0)
        self.slider = tk.Scale(sl, from_=0, to=0, orient="horizontal",
                               resolution=1, showvalue=False,
                               variable=self.slider_var,
                               command=self._slider_moved,
                               background=BG, highlightthickness=0,
                               troughcolor="#dfe3e8", borderwidth=0,
                               sliderlength=18, width=12)
        self.slider.pack(fill="x")
        self.slider_lbl = ttk.Label(sl, text="run first", foreground="#777",
                                    font=("Segoe UI", 8))
        self.slider_lbl.pack(anchor="w", pady=(4, 0))
        self._state = None

        pane = ttk.PanedWindow(body, orient="vertical"); pane.pack(side="left", fill="both", expand=True)
        self.results = ResultPane(pane)
        self.plot = PlotPane(pane, self.figsize)
        pane.add(self.results, weight=2); pane.add(self.plot, weight=3)

        self._source_changed()

    def _source_changed(self, _=None):
        """Show/hide the isotope picker; Archer availability is handled
        separately since it depends on both source AND material."""
        is_iso = self.source.get() == "Isotope"
        if is_iso:
            self.isotope_row.pack(fill="x", pady=(4, 0))
        else:
            self.isotope_row.pack_forget()
        self._update_archer_availability()

    def _update_archer_availability(self, _=None):
        """
        X-ray spectrum source: Archer's kVp-indexed NCRP 147 table is
        always available.

        Isotope source: that kVp table means nothing for a monoenergetic
        gamma line, so it is never used here. Instead, a SEPARATE table of
        Archer fits made directly against radionuclide emission spectra
        (RADIONUCLIDES, lead/concrete only) is checked for the current
        isotope + material combination. Only enable the overlay where a
        fit actually exists, rather than leaving it on to silently show
        nothing or (worse) the wrong table.
        """
        is_iso = self.source.get() == "Isotope"
        if not is_iso:
            self.archer_chk.configure(state="normal")
            return
        mat = mat_key(self.material.get())
        has_fit = nuclide_params(self.isotope.get(), mat) is not None
        self.archer_chk.configure(state="normal" if has_fit else "disabled")
        self.show_archer.set(has_fit)

    def _slider_moved(self, _=None):
        st = getattr(self, "_state", None)
        if not st or not st["specs"]:
            return
        try:
            i = int(round(float(self.slider_var.get())))
        except (TypeError, ValueError):
            return
        self._draw(max(0, min(i, len(st["specs"]) - 1)))

    def _draw(self, idx):
        """
        Redraw both panels for one of the thicknesses already simulated.
        Every step of the run stores its transmitted spectrum, so moving the
        slider is a redraw, not a re-run.
        """
        st = self._state
        df, sp, mat, dens, specs = (st["df"], st["sp"], st["mat"],
                                    st["dens"], st["specs"])
        idx = max(0, min(idx, len(specs) - 1)) if specs else 0
        sel = specs[idx] if specs else None
        if sel:
            self.slider_lbl.configure(
                text="%s mm of %s" % (fmt(sel[0]),
                                      MATERIAL_LABEL.get(mat, mat)),
                foreground=OK)
        fig = self.plot.reset(); ax = fig.add_subplot(121)
        x = df["thickness_mm"]
        col = ACCENT
        ax.semilogy(x, df["B_montecarlo"], "o-", color=col, lw=2, ms=5,
                    label="Monte Carlo")
        ax.semilogy(x, df["B_narrow"], "--", color="#999", lw=1.4,
                    label="narrow beam - physical floor")
        if self.show_archer.get() and df["B_archer"].notna().all():
            ax.semilogy(x, df["B_archer"], "s:", color=WARN, lw=1.6, ms=4,
                        label="Archer empirical")
        ax.grid(True, which="both", alpha=0.25); ax.set_axisbelow(True)
        ax.set_xlabel("thickness (mm)", fontsize=9)
        ax.set_ylabel("transmission B", fontsize=9)
        ax.set_title("%s, %s" % (MATERIAL_LABEL.get(mat, mat),
                                 sp.model if sp.model.startswith("isotope:")
                                 else "%s kVp" % fmt(sp.kvp)),
                     fontsize=10, color=ACCENT, fontweight="bold")
        if sel:
            ax.axvline(sel[0], color="#111", lw=0.9, ls=":", alpha=0.7)
        ax.legend(fontsize=8); ax.tick_params(labelsize=8)

        ax2 = fig.add_subplot(122)
        inc = sp.N / max(float(np.max(sp.N)), 1e-30)
        ax2.plot(sp.E, inc, color="#999", lw=1.3, label="incident")
        if sel is not None and sel[2].sum() > 0:
            Ht, He, Hw = sel
            Hw = Hw / max(float(np.max(Hw)), 1e-30)
            ax2.fill_between(He, Hw, color=col, alpha=0.22)
            ax2.plot(He, Hw, color=col, lw=1.8,
                     label="transmitted, %s mm" % fmt(Ht))
            mean_out = float(np.sum(sel[1] * sel[2]))
            ax2.axvline(mean_out, color="black", ls="--", lw=1.0)
            ax2.annotate("mean %.1f keV" % mean_out,
                         (mean_out, 0.9), fontsize=7, color="black",
                         xytext=(4, 0), textcoords="offset points")
        ax2.set_xlim(0, sp.kvp); ax2.set_ylim(0, 1.45)
        ax2.grid(alpha=0.25); ax2.set_axisbelow(True)
        ax2.set_xlabel("photon energy (keV)", fontsize=9)
        ax2.set_ylabel("relative fluence per keV", fontsize=9)
        ax2.set_title("Spectrum reaching the detector, %s mm"
                      % (fmt(sel[0]) if sel else "-"), fontsize=10,
                      color=ACCENT, fontweight="bold")
        ax2.tick_params(labelsize=8)

        # Attenuation coefficients on their own axis. This is a genuine
        # use for a twin axis - fluence and cm^-1 are different
        # quantities - and it shows why the transmitted peak sits where
        # it does: the photoelectric minimum just below the K edge.
        if self.show_mu.get():
            ax3 = ax2.twinx()
            Eg = np.linspace(max(sp.E.min(), 5.0), sp.kvp, 600)
            # Total is the SUM, so it is always the topmost curve. Drawn
            # first as a thick pale envelope: a thin dotted total sits
            # almost exactly on photoelectric and is easily mistaken for
            # the flat Compton line further down.
            ax3.semilogy(Eg, XS.mu(mat, Eg, "total", dens),
                         color="#444444", lw=4.0, alpha=0.25,
                         label="total (sum)", solid_capstyle="round")
            for kind, colr, style, lbl in (
                    ("photo", "#c1666b", ":", "photoelectric"),
                    ("coh", "#1a7f4b", ":", "coherent"),
                    ("incoh", "#7a5195", ":", "Compton")):
                ax3.semilogy(Eg, XS.mu(mat, Eg, kind, dens),
                             color=colr, ls=style, lw=1.3, alpha=0.95,
                             label=lbl)
            ax3.set_ylabel("attenuation coefficient (cm$^{-1}$)",
                           fontsize=9)
            ax3.tick_params(labelsize=8)
            h1, l1 = ax2.get_legend_handles_labels()
            h2, l2 = ax3.get_legend_handles_labels()
            ax2.legend(h1 + h2, l1 + l2, fontsize=6.5, loc="upper center",
                       framealpha=0.95, ncol=3)
        else:
            ax2.legend(fontsize=7, loc="upper center", framealpha=0.95, ncol=2)
        self.plot.done()

    def run(self):
        r = self.results; r.clear()
        try:
            is_iso = self.source.get() == "Isotope"
            beyond = []
            if is_iso:
                iso_name = self.isotope.get()
                sp = make_isotope_spectrum(iso_name)
                beyond = isotope_lines_beyond_xsdata(iso_name)
            else:
                sp = self.spectrum_tab.spectrum
                if sp is None:
                    self.spectrum_tab.run()
                    sp = self.spectrum_tab.spectrum
                if sp is None:
                    raise ValueError("Generate a spectrum first.")
            mat = mat_key(self.material.get())
            t0, t1 = self.tmin.get(), self.tmax.get()
            ns = max(int(self.nsteps.get(7)), 1)
            nh = max(int(self.nhist.get(20000)), 100)
            dens = self.density.get(0.0) or None
            acc = min(max(self.accept.get(90.0), 0.5), 90.0)
            ethr = max(self.ethr.get(0.0), 0.0)
            if t0 <= 0 or t1 < t0:
                raise ValueError("Thickness range is not valid.")
            ts = np.linspace(t0, t1, ns) if ns > 1 else np.array([t0])

            self.btn.configure(state="disabled")
            sc = ("kerma" if self.score.get().startswith("Air")
                  else "fluence")
            rows = []
            specs = []
            for i, t in enumerate(ts):
                self.prog["value"] = 100.0 * i / len(ts)
                self.update_idletasks()
                # Coherent scatter is sampled from the real atomic form
                # factor (xsdata.npz's formfactor_x/*__ffcum), not bracketed
                # between forward-only and Thomson - see the Notes tab.
                res = run_monte_carlo(sp, mat, float(t), n_hist=nh,
                                      density=dens, score=sc,
                                      acceptance_deg=acc, e_threshold=ethr,
                                      fluorescence=True)
                ar = None
                if self.show_archer.get():
                    if is_iso:
                        ar = archer_reference_nuclide(float(t), iso_name, mat)
                    else:
                        ar = archer_reference(float(t), sp.kvp, mat)
                specs.append((float(t), res["hist_E"], res["hist_norm"]))
                rows.append({"thickness_mm": round(float(t), 4),
                             "B_narrow": res["B_narrow"],
                             "B_montecarlo": res["B_mc"],
                             "buildup": res["buildup"],
                             "noise_pc": 100 * res["rel_std_err"],
                             "B_archer": ar,
                             "MC_over_Archer": (res["B_mc"] / ar
                                                if ar else None)})
            self.prog["value"] = 100
            self.btn.configure(state="normal")
            df = pd.DataFrame(rows)

            r.head("Monte Carlo transmission through %s"
                   % MATERIAL_LABEL.get(mat, mat))
            if is_iso:
                r.row("Source", "isotope: %s" % iso_name)
                r.row("Lines", ", ".join(
                    "%s keV (%.1f%%)" % (fmt(e), 100 * i)
                    for e, i in ISOTOPE_LINES[iso_name]))
            else:
                r.row("Spectrum", "%s kVp, %s" % (fmt(sp.kvp), sp.model))
            r.row("Mean energy", fmt(sp.mean_energy), "keV")
            r.row("Histories per point", "%d" % nh)
            r.row("Coherent scatter", "atomic form factors")
            r.row("Detector acceptance", fmt(acc), "deg half-angle")
            r.row("Detector threshold", fmt(ethr), "keV")
            r.row("K fluorescence", "tracked")
            r.row("Score", self.score.get())
            r.row("Density", fmt(dens if dens else XS.density(mat)), "g/cm3")
            if is_iso and self.show_archer.get():
                basis, unit = nuclide_info(iso_name, mat)
                if basis:
                    r.row("Archer overlay", "radionuclide fit (%s, source "
                          "units: %s)" % (basis, unit))
                    if basis == "inferred":
                        r.warn("This nuclide's Archer fit is flagged "
                               "'inferred' at the source - both unit "
                               "readings were arguable and the one "
                               "consistent with %s's photon energies was "
                               "taken. Worth re-checking before relying on "
                               "it." % iso_name)
            if beyond:
                r.warn("%s line(s) at %s keV exceed xsdata.npz's 800 keV "
                       "cross-section range (the Elam tables it comes from "
                       "are flagged unreliable above there too). Those "
                       "lines are transported using the frozen 800 keV "
                       "attenuation coefficient, which is measurably wrong "
                       "- real mu keeps falling with energy in this "
                       "Compton-dominated regime, so transmission for "
                       "those lines is UNDERESTIMATED here. Do not use "
                       "this result for a Co-60 shielding decision; a "
                       "validated code (OpenMC, EGSnrc, PENELOPE) or a "
                       "published Co-60 buildup factor table is needed."
                       % (len(beyond), ", ".join(fmt(e) for e in beyond)))
            r.head("Results")
            disp = df.copy()
            for c in ("B_narrow", "B_montecarlo", "B_archer"):
                disp[c] = disp[c].map(lambda v: fmt(v, 4)
                                      if v is not None else "-")
            for c in ("buildup", "MC_over_Archer"):
                disp[c] = disp[c].map(lambda v: fmt(v, 3)
                                      if v is not None else "-")
            disp["noise_pc"] = disp["noise_pc"].map(lambda v: "%.1f" % v)
            r.frame(disp)
            r._df = df

            r.head("Interpretation")
            r.row("Buildup at thinnest point", fmt(df["buildup"].iloc[0], 3))
            r.row("Buildup at thickest point", fmt(df["buildup"].iloc[-1], 3))
            if df["MC_over_Archer"].notna().all():
                r.row("MC / Archer, mean",
                      fmt(df["MC_over_Archer"].mean(), 3))
                r.row("MC / Archer, range",
                      "%s to %s" % (fmt(df["MC_over_Archer"].min(), 3),
                                    fmt(df["MC_over_Archer"].max(), 3)))
            below = df[df["B_archer"] < df["B_narrow"]] \
                if df["B_archer"].notna().all() else df.iloc[0:0]
            if len(below):
                r.warn("Beyond %s mm Archer drops below the no-scatter floor."
                       % fmt(below["thickness_mm"].iloc[0]))
                al = mu_top = None
                if is_iso:
                    p = nuclide_params(iso_name, mat)
                    if p:
                        al = p[0]
                        mu_top = float(XS.mu(mat, float(sp.mean_energy),
                                             "total")) / 10.0
                        r.row("Archer alpha", fmt(al, 4), "/mm")
                        r.row("mu at mean energy (%s keV)"
                              % fmt(sp.mean_energy), fmt(mu_top, 4), "/mm")
                else:
                    tbl = ARCHER_REF.get(mat, {})
                    if tbl:
                        ks = sorted(tbl)
                        kv = min(max(sp.kvp, ks[0]), ks[-1])
                        near = min(ks, key=lambda k: abs(k - kv))
                        al = tbl[near][0]
                        mu_top = float(XS.mu(mat, float(sp.kvp),
                                             "total")) / 10.0
                        r.row("Archer alpha", fmt(al, 4), "/mm")
                        r.row("mu at %s keV" % fmt(sp.kvp), fmt(mu_top, 4),
                              "/mm")
                if al is not None and mu_top is not None and al > mu_top:
                    r.note("Alpha is the fit's asymptotic decay rate. It "
                           "exceeds mu at the beam's representative "
                           "energy, so the fit decays faster than the "
                           "most penetrating photons in the beam "
                           "physically can. Crossing the floor is then "
                           "unavoidable at some depth.")
                    r.note("That is not an error in Archer's measurements. "
                           "It is a three-parameter form being asked to "
                           "describe attenuation far beyond the data it "
                           "was fitted to. Ignore the comparison there.")
            worst = df["noise_pc"].max()
            r.row("Worst random noise", "%.1f %%" % worst)
            if worst > 15:
                r.warn("Noise above 15 %. Increase histories or drop the "
                       "thickest points.")
            r.top()

            self._state = dict(df=df, sp=sp, mat=mat, dens=dens,
                               specs=specs)
            n = max(len(specs) - 1, 0)
            self.slider.configure(from_=0, to=n)
            self.slider_var.set(n)
            self._draw(n)
        except ValueError as e:
            self.btn.configure(state="normal")
            r.head("Input error"); r.warn(str(e))
        except Exception as e:
            self.btn.configure(state="normal")
            r.head("Error"); r.warn("%s: %s" % (type(e).__name__, e))


# =============================================================================
#  SECTION 6 -- IMAGE QUALITY: MTF, NNPS, DQE (IEC 62220-1 methodology)
#
#  MTF: slanted-edge method. The edge is deliberately tilted 1-3 degrees so
#  that successive rows (or columns) sample the edge profile at slightly
#  different sub-pixel offsets - collecting all of them together builds an
#  oversampled edge spread function from an ordinary pixel grid, without
#  needing a genuinely sub-pixel detector. Differentiate -> line spread
#  function, window, FFT, normalise to MTF(0)=1.
#
#  NNPS: flat-field image split into many overlapping ROIs, each detrended
#  (removes heel effect / large-scale non-uniformity), windowed, FFT'd and
#  averaged into a 2D noise power spectrum, then normalised by the square
#  of the mean signal.
#
#  DQE(f) = MTF(f)^2 / (q * K * NNPS(f))
#  where q = photon fluence per unit air kerma (mm^-2.uGy^-1) - this reuses
#  RQA_SNR_IN2 from earlier, since that is exactly this quantity - and K is
#  the air kerma (uGy) at the detector for the flat-field exposure. Verified
#  by dimensional analysis: q*K is the incident photon fluence (mm^-2) for
#  that exposure, NNPS carries mm^2, so q*K*NNPS is dimensionless, matching
#  MTF^2.
#
#  SIMPLIFICATIONS (see the Notes tab):
#  * one edge per direction, not the two-edge average IEC allows
#  * no rebinning onto IEC's exact reported frequency spacing
#  * detector response assumed linear through zero (single flat-field
#    exposure, no multi-point signal-transfer curve)
# =============================================================================

def load_dicom_pixels(path):
    """
    Read a DICOM image and return (pixel_array float32 with rescale
    slope/intercept applied, pixel_spacing_mm, pydicom Dataset).

    float32 not float64: this is 16-bit integer detector data (max value
    65535), and float32 has 24 bits of mantissa - room for integers up to
    ~16.7 million exactly. Zero precision loss, half the memory, and
    measurably faster to convert (verified: ~3x faster for a 4316x4316
    image in testing, both from the smaller dtype and from skipping the
    rescale multiply entirely when slope=1/intercept=0, a very common
    case that was previously costing a full-array multiply-add for
    nothing).
    """
    if not HAVE_PYDICOM:
        raise ValueError("pydicom is not installed - pip install pydicom")
    ds = pydicom.dcmread(path)
    raw = ds.pixel_array
    slope = float(getattr(ds, "RescaleSlope", 1.0))
    intercept = float(getattr(ds, "RescaleIntercept", 0.0))
    if slope == 1.0 and intercept == 0.0:
        arr = raw.astype(np.float32)
    else:
        arr = raw.astype(np.float32) * np.float32(slope) + np.float32(intercept)
    px_mm = None
    for attr in ("PixelSpacing", "ImagerPixelSpacing"):
        if hasattr(ds, attr):
            px_mm = float(getattr(ds, attr)[0])
            break
    if px_mm is None:
        raise ValueError("No PixelSpacing or ImagerPixelSpacing tag found - "
                         "cannot convert pixels to mm.")
    return arr, px_mm, ds


def slanted_edge_mtf(roi, px_mm, oversample=8, edge_axis="vertical"):
    """
    IEC-style slanted-edge MTF from a single ROI containing one edge.

    edge_axis: 'vertical' - the edge separates left/right (scan each row
    for the crossing); 'horizontal' - separates top/bottom (scan each
    column). Polarity-agnostic: works off the ROI's own 2nd/98th
    percentile levels, not an assumption about which side is "bright".

    Returns dict with freq (cy/mm), mtf, angle_deg, esf, lsf, bin_mm.
    """
    work = roi.T if edge_axis == "horizontal" else roi
    n_rows, n_cols = work.shape
    lo, hi = np.percentile(work, [2, 98])
    if hi <= lo:
        raise ValueError("No contrast found in this ROI - check the edge "
                         "is actually inside the selected region.")
    mid = 0.5 * (lo + hi)

    crossings = []
    for i in range(n_rows):
        line = work[i]
        above = line > mid
        diffs = np.diff(above.astype(int))
        idx = np.where(diffs != 0)[0]
        if len(idx) == 0:
            continue
        j = idx[len(idx) // 2]
        y0, y1 = line[j], line[j + 1]
        if y1 == y0:
            continue
        crossings.append((i, j + (mid - y0) / (y1 - y0)))
    crossings = np.array(crossings)
    frac_ok = len(crossings) / n_rows if n_rows else 0.0
    if frac_ok < 0.5:
        wrong_axis = "vertical" if edge_axis == "horizontal" else "horizontal"
        raise ValueError(
            "Only %.0f%% of %s found a clean crossing (need at least "
            "50%%). Contrast range in this ROI: %.0f to %.0f. Most likely "
            "causes: (1) the ROI doesn't actually straddle an edge - it "
            "may sit entirely on one side; (2) it includes a corner of "
            "the aperture, where the edge direction changes; (3) this ROI "
            "actually has a %s edge in it, which belongs in the OTHER "
            "MTF direction, not this one." %
            (100 * frac_ok, "rows" if edge_axis == "vertical" else "columns",
             lo, hi, wrong_axis))

    a, b = np.polyfit(crossings[:, 0], crossings[:, 1], 1)
    angle_deg = float(np.degrees(np.arctan(a)))

    rows_idx, cols_idx = np.meshgrid(np.arange(n_rows), np.arange(n_cols),
                                     indexing="ij")
    dist_px = (cols_idx - (a * rows_idx + b)).ravel()
    v = work.ravel()

    bin_w = 1.0 / oversample
    lo_s, hi_s = np.percentile(dist_px, [1, 99])
    bins = np.arange(lo_s, hi_s, bin_w)
    bin_idx = np.digitize(dist_px, bins)
    esf = np.full(len(bins) - 1, np.nan)
    for k in range(1, len(bins)):
        sel = bin_idx == k
        if sel.any():
            esf[k - 1] = v[sel].mean()
    valid = ~np.isnan(esf)
    if valid.sum() < 16:
        raise ValueError("Too few valid ESF bins - ROI may be too small.")
    esf = esf[valid]

    lsf = np.gradient(esf)
    window = np.hanning(len(lsf))
    lsf_w = lsf * window

    mtf = np.abs(np.fft.rfft(lsf_w))
    mtf = mtf / mtf[0] if mtf[0] != 0 else mtf
    freq = np.fft.rfftfreq(len(lsf_w), d=bin_w * px_mm)

    return {"freq": freq, "mtf": mtf, "angle_deg": angle_deg,
           "esf": esf, "lsf": lsf, "bin_mm": bin_w * px_mm}


def _detrend_2d(roi, order=2):
    """Subtract a low-order 2D polynomial fit (removes heel effect etc.)."""
    ny, nx = roi.shape
    y, x = np.mgrid[0:ny, 0:nx]
    x = x.ravel(); y = y.ravel(); z = roi.ravel()
    terms = [(x ** i) * (y ** j) for i in range(order + 1)
            for j in range(order + 1 - i)]
    A = np.vstack(terms).T
    coef, *_ = np.linalg.lstsq(A, z, rcond=None)
    return roi - (A @ coef).reshape(ny, nx)


def flat_field_nnps(image, px_mm, roi_size=256, step=128, mask_below=1.0,
                    progress=None):
    """
    2D NNPS from a flat-field image: overlapping detrended, windowed,
    FFT'd ROIs, averaged, normalised by the square of the mean signal.

    mask_below: ROIs containing any pixel at or below this value are
    skipped (avoids the masked-out / unexposed border seen at the edges
    of the true collimated field).

    Returns (nnps_2d, freq_axis_cy_per_mm, mean_signal, n_rois_used).
    """
    ny, nx = image.shape
    win = np.outer(np.hanning(roi_size), np.hanning(roi_size))
    win_mean_sq = np.mean(win ** 2)
    accum = None
    count = 0
    positions = [(y0, x0) for y0 in range(0, ny - roi_size, step)
                for x0 in range(0, nx - roi_size, step)]
    valid_means = []
    for i, (y0, x0) in enumerate(positions):
        block = image[y0:y0 + roi_size, x0:x0 + roi_size]
        if (block <= mask_below).any():
            continue
        valid_means.append(block.mean())
        d = _detrend_2d(block, order=2) * win
        F = np.fft.fftshift(np.fft.fft2(d))
        ps = (np.abs(F) ** 2) * (px_mm ** 2) / (roi_size * roi_size) \
            / win_mean_sq
        accum = ps if accum is None else accum + ps
        count += 1
        if progress and i % 20 == 0:
            progress(i / len(positions))
    if count == 0:
        raise ValueError("No usable ROIs - flat field may be smaller than "
                         "expected, or mostly masked-out background.")
    nnps2d = (accum / count) / (np.mean(valid_means) ** 2)
    freq = np.fft.fftshift(np.fft.fftfreq(roi_size, d=px_mm))
    return nnps2d, freq, float(np.mean(valid_means)), count


def nnps_profile(nnps2d, freq, axis="h", band=3):
    """1D NNPS profile: a narrow band of rows/columns straddling the
    relevant axis, excluding the true zero-frequency bin itself."""
    n = nnps2d.shape[0]
    c = n // 2
    prof = (nnps2d[c - band:c + band, :].mean(axis=0) if axis == "h"
           else nnps2d[:, c - band:c + band].mean(axis=1))
    return freq, prof


def compute_dqe(freq_mtf, mtf, freq_nnps, nnps, q, K_uGy):
    """
    DQE(f) = MTF(f)^2 / (q * K * NNPS(f)), interpolating MTF and NNPS onto
    a shared frequency grid (the positive-frequency half of freq_nnps).
    """
    pos = freq_nnps >= 0
    fN = freq_nnps[pos]
    nN = nnps[pos]
    fcommon = fN[fN <= freq_mtf.max()]
    mtf_i = np.interp(fcommon, freq_mtf, mtf)
    nnps_i = np.interp(fcommon, fN, nN)
    with np.errstate(divide="ignore", invalid="ignore"):
        dqe = np.where(nnps_i > 0, mtf_i ** 2 / (q * K_uGy * nnps_i), np.nan)
    return fcommon, dqe


def simple_open_dicom_dialog(parent, initialdir=None):
    """
    A file picker built entirely from Tkinter widgets - deliberately NOT
    using filedialog.askopenfilename(). On some Windows machines that
    native dialog takes 30+ seconds to first appear (Windows initialising
    its own shell/COM subsystem - icon handlers, thumbnail providers,
    cloud-storage shell extensions - the first time any such dialog is
    shown in a process). That cost lives entirely inside Windows' own
    dialog implementation, invoked as one blocking call with nothing of
    this app's code running during it - not something reachable or
    speedable from here. This picker only ever calls os.listdir(), so
    there is no shell extension in the loop to be slow.

    Returns the selected file path, or None if cancelled.
    """
    top = tk.Toplevel(parent)
    top.title("Select DICOM file")
    top.geometry("640x460")
    top.transient(parent)
    top.grab_set()
    result = {"path": None}

    row1 = ttk.Frame(top, padding=6); row1.pack(fill="x")
    ttk.Label(row1, text="Folder:").pack(side="left")
    path_var = tk.StringVar(value=initialdir or os.path.expanduser("~"))
    path_entry = ttk.Entry(row1, textvariable=path_var)
    path_entry.pack(side="left", fill="x", expand=True, padx=(4, 4))

    if os.name == "nt":
        import string
        drives = [d + ":\\" for d in string.ascii_uppercase
                 if os.path.exists(d + ":\\")]
        drive_var = tk.StringVar(value="")
        drive_combo = ttk.Combobox(row1, textvariable=drive_var, values=drives,
                                   width=5, state="readonly")
        drive_combo.pack(side="left", padx=(0, 4))

        def _go_drive(_=None):
            if drive_var.get():
                refresh(drive_var.get())
        drive_combo.bind("<<ComboboxSelected>>", _go_drive)

    list_frame = ttk.Frame(top, padding=(6, 0)); list_frame.pack(fill="both", expand=True)
    scrollbar = ttk.Scrollbar(list_frame, orient="vertical")
    lb = tk.Listbox(list_frame, yscrollcommand=scrollbar.set, font=("Consolas", 10))
    scrollbar.configure(command=lb.yview)
    lb.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")

    entries = []   # parallel to listbox rows: (full_path, is_dir)

    def refresh(path):
        try:
            names = sorted(os.listdir(path), key=lambda s: s.lower())
        except Exception as e:
            messagebox.showerror("Cannot open folder", str(e), parent=top)
            return
        entries.clear()
        lb.delete(0, "end")
        parent_dir = os.path.dirname(path.rstrip("\\/"))
        if parent_dir and parent_dir != path.rstrip("\\/"):
            entries.append((parent_dir, True))
            lb.insert("end", "[..]")
        dirs, files = [], []
        for name in names:
            full = os.path.join(path, name)
            try:
                if os.path.isdir(full):
                    dirs.append((name, full))
                elif name.lower().endswith(".dcm"):
                    files.append((name, full))
            except OSError:
                continue
        for name, full in dirs:
            entries.append((full, True)); lb.insert("end", "[%s]" % name)
        for name, full in files:
            entries.append((full, False)); lb.insert("end", name)
        path_var.set(path)

    def open_selected():
        sel = lb.curselection()
        if not sel:
            return
        full, is_dir = entries[sel[0]]
        if is_dir:
            refresh(full)
        else:
            result["path"] = full
            top.destroy()

    def go_to_typed_path():
        p = path_var.get()
        if os.path.isdir(p):
            refresh(p)
        elif os.path.isfile(p):
            result["path"] = p
            top.destroy()
        else:
            messagebox.showerror("Not found", "No such folder or file:\n%s" % p,
                                 parent=top)

    lb.bind("<Double-Button-1>", lambda e: open_selected())
    path_entry.bind("<Return>", lambda e: go_to_typed_path())
    ttk.Button(row1, text="Go", command=go_to_typed_path).pack(side="left")

    row2 = ttk.Frame(top, padding=6); row2.pack(fill="x")
    ttk.Label(row2, text="Double-click a folder to open it, or a .dcm "
                        "file to select it.", foreground="#777",
              font=("Segoe UI", 8)).pack(side="left")
    ttk.Button(row2, text="Open", command=open_selected).pack(side="right")
    ttk.Button(row2, text="Cancel", command=top.destroy).pack(side="right",
                                                               padx=(0, 6))

    refresh(path_var.get())
    top.wait_window()
    return result["path"]


class ImageQualityTab(ttk.Frame):
    title = "  MTF / NNPS / DQE  "

    figsize = (7.6, 6.4)

    def __init__(self, parent):
        super().__init__(parent, padding=10)
        if not HAVE_PYDICOM:
            ttk.Label(self, text="pydicom is not installed - this tab needs "
                                 "it to read DICOM images.\npip install "
                                 "pydicom",
                      foreground=WARN, justify="left").pack(pady=20, padx=10,
                                                            anchor="w")
            return

        ttk.Label(self, text=(
            "MTF (slanted-edge), NNPS (flat field) and DQE, per IEC 62220-1 "
            "methodology."),
            foreground="#555", wraplength=1200, justify="left").pack(
                fill="x", pady=(0, 8))
        body = ttk.Frame(self); body.pack(fill="both", expand=True)

        # The control panel (4 sections, the last with 3 text boxes) can be
        # taller than the available window height on a smaller screen - it
        # was silently getting clipped with no way to reach whatever fell
        # below the bottom edge, e.g. the DQE factory spec box. Wrapping it
        # in a scrollable canvas means everything stays reachable regardless
        # of window size, rather than depending on the screen being big
        # enough.
        left_container = ttk.Frame(body)
        left_container.pack(side="left", fill="y", padx=(0, 12))
        left_canvas = tk.Canvas(left_container, borderwidth=0,
                                highlightthickness=0, width=290)
        left_vsb = ttk.Scrollbar(left_container, orient="vertical",
                                 command=left_canvas.yview)
        left_canvas.configure(yscrollcommand=left_vsb.set)
        left_canvas.pack(side="left", fill="both", expand=True)
        left_vsb.pack(side="right", fill="y")
        left = ttk.Frame(left_canvas)
        left_window = left_canvas.create_window((0, 0), window=left,
                                                 anchor="nw")

        def _on_left_configure(event):
            left_canvas.configure(scrollregion=left_canvas.bbox("all"))
        left.bind("<Configure>", _on_left_configure)

        def _on_canvas_configure(event):
            left_canvas.itemconfig(left_window, width=event.width)
        left_canvas.bind("<Configure>", _on_canvas_configure)

        def _on_mousewheel(event):
            left_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        def _bind_wheel(event):
            left_canvas.bind_all("<MouseWheel>", _on_mousewheel)
        def _unbind_wheel(event):
            left_canvas.unbind_all("<MouseWheel>")
        left_canvas.bind("<Enter>", _bind_wheel)
        left_canvas.bind("<Leave>", _unbind_wheel)

        s1 = ttk.LabelFrame(left, text="1. Edge image (MTF)", padding=8)
        s1.pack(fill="x", pady=(0, 8))
        ttk.Button(s1, text="Load MTF tool image",
                  command=self.load_edge_image).pack(fill="x")
        self.edge_file_lbl = ttk.Label(s1, text="(none loaded)",
                                       foreground="#777", font=("Segoe UI", 8))
        self.edge_file_lbl.pack(anchor="w", pady=(2, 6))
        self.btn_show_edge = ttk.Button(s1, text="View edge image / select ROIs",
                                        command=self._show_edge_image)
        self.btn_show_edge.pack(fill="x", pady=(0, 6))
        self.roi_mode = LabelledCombo(
            s1, "Now selecting",
            ["Horizontal MTF - drag a HORIZONTAL edge",
             "Vertical MTF - drag a VERTICAL edge"],
            "Horizontal MTF - drag a HORIZONTAL edge", width=24)
        self.roi_mode.pack(fill="x", pady=2)
        self.roi_h_lbl = ttk.Label(s1, text="Horizontal-MTF ROI: not set",
                                   foreground="#777", font=("Segoe UI", 8))
        self.roi_h_lbl.pack(anchor="w")
        self.roi_v_lbl = ttk.Label(s1, text="Vertical-MTF ROI: not set",
                                   foreground="#777", font=("Segoe UI", 8))
        self.roi_v_lbl.pack(anchor="w", pady=(0, 4))
        row_h = ttk.Frame(s1); row_h.pack(fill="x", pady=(4, 2))
        self.btn_mtf_h = ttk.Button(row_h, text="Calculate horizontal MTF",
                                    command=self.calc_mtf_h)
        self.btn_mtf_h.pack(side="left", fill="x", expand=True)
        self.btn_copy_mtf_h = ttk.Button(row_h, text="Copy", width=6,
                                         command=lambda: self.copy_curve(
                                             "res_h", 103, "Horizontal MTF"))
        self.btn_copy_mtf_h.pack(side="left", padx=(4, 0))
        row_v = ttk.Frame(s1); row_v.pack(fill="x")
        self.btn_mtf_v = ttk.Button(row_v, text="Calculate vertical MTF",
                                    command=self.calc_mtf_v)
        self.btn_mtf_v.pack(side="left", fill="x", expand=True)
        self.btn_copy_mtf_v = ttk.Button(row_v, text="Copy", width=6,
                                         command=lambda: self.copy_curve(
                                             "res_v", 103, "Vertical MTF"))
        self.btn_copy_mtf_v.pack(side="left", padx=(4, 0))
        self.mtf_status_lbl = ttk.Label(s1, text="MTF: neither done yet",
                                        foreground="#777", font=("Segoe UI", 8))
        self.mtf_status_lbl.pack(anchor="w", pady=(4, 0))

        s2 = ttk.LabelFrame(left, text="2. Flat field (NNPS)", padding=8)
        s2.pack(fill="x", pady=(0, 8))
        ttk.Button(s2, text="Load NNPS image",
                  command=self.load_flat_image).pack(fill="x")
        self.flat_file_lbl = ttk.Label(s2, text="(none loaded)",
                                       foreground="#777", font=("Segoe UI", 8))
        self.flat_file_lbl.pack(anchor="w", pady=(2, 6))
        self.btn_nnps = ttk.Button(s2, text="Calculate NNPS",
                                   command=self.calc_nnps)
        self.btn_nnps.pack(fill="x")
        row_nnps_copy = ttk.Frame(s2); row_nnps_copy.pack(fill="x", pady=(4, 0))
        self.btn_copy_nnps_h = ttk.Button(row_nnps_copy, text="Copy NNPS (h)",
                                          command=lambda: self.copy_nnps("h"))
        self.btn_copy_nnps_h.pack(side="left", fill="x", expand=True)
        self.btn_copy_nnps_v = ttk.Button(row_nnps_copy, text="Copy NNPS (v)",
                                          command=lambda: self.copy_nnps("v"))
        self.btn_copy_nnps_v.pack(side="left", fill="x", expand=True, padx=(4, 0))
        self.nnps_status_lbl = ttk.Label(s2, text="NNPS: not done yet",
                                         foreground="#777", font=("Segoe UI", 8))
        self.nnps_status_lbl.pack(anchor="w", pady=(4, 0))

        s3 = ttk.LabelFrame(left, text="3. DQE", padding=8)
        s3.pack(fill="x", pady=(0, 8))
        self.kerma = LabelledEntry(s3, "Air kerma at detector", 0.0, "uGy")
        self.kerma.pack(fill="x", pady=2)
        self.quality = LabelledCombo(
            s3, "RQA quality",
            sorted(RQA_SNR_IN2, key=lambda q: RQA_ED2[q][0]),
            "RQA5", width=12)
        self.quality.pack(fill="x", pady=2)
        ttk.Label(s3, text="q (SNR_in^2) comes straight from the IEC "
                          "62220-1 Table 3 data already in the app -\n"
                          "only RQA3/5/7/9 carry a mandatory value.",
                  foreground="#777", font=("Segoe UI", 8),
                  justify="left").pack(anchor="w", pady=(2, 4))
        self.btn_dqe = ttk.Button(s3, text="Calculate DQE",
                                  style="Accent.TButton",
                                  command=self.calc_dqe)
        self.btn_dqe.pack(fill="x", ipady=4)
        self.btn_copy_dqe = ttk.Button(s3, text="Copy DQE",
                                       command=self.copy_dqe)
        self.btn_copy_dqe.pack(fill="x", pady=(4, 0))

        s4 = ttk.LabelFrame(left, text="4. Factory spec (optional)", padding=8)
        s4.pack(fill="x", pady=(0, 8))
        ttk.Label(s4, text="One 'lp/mm, value' pair per line, for each:",
                  font=("Segoe UI", 8), foreground="#777",
                  wraplength=200, justify="left").pack(anchor="w")
        ttk.Label(s4, text="MTF (%):", font=("Segoe UI", 8)).pack(anchor="w", pady=(4, 0))
        self.spec_text_mtf = tk.Text(s4, height=4, width=24, font=("Consolas", 9))
        self.spec_text_mtf.pack(fill="x", pady=(0, 2))
        ttk.Label(s4, text="NNPS:", font=("Segoe UI", 8)).pack(anchor="w", pady=(4, 0))
        self.spec_text_nnps = tk.Text(s4, height=4, width=24, font=("Consolas", 9))
        self.spec_text_nnps.pack(fill="x", pady=(0, 2))
        ttk.Label(s4, text="DQE (%):", font=("Segoe UI", 8)).pack(anchor="w", pady=(4, 0))
        self.spec_text_dqe = tk.Text(s4, height=4, width=24, font=("Consolas", 9))
        self.spec_text_dqe.pack(fill="x", pady=(0, 2))

        self.prog = ttk.Progressbar(left, mode="determinate", maximum=100)
        self.prog.pack(fill="x", pady=(6, 0))

        pane = ttk.PanedWindow(body, orient="vertical")
        pane.pack(side="left", fill="both", expand=True)
        self.results = ResultPane(pane)
        self.plot = PlotPane(pane, self.figsize)
        pane.add(self.results, weight=2); pane.add(self.plot, weight=3)

        self.edge_arr = self.edge_px_mm = self.edge_ds = None
        self.flat_arr = self.flat_px_mm = self.flat_ds = None
        self.roi_h = self.roi_v = None
        self.res_h = self.res_v = None
        self.nnps_data = None      # dict: nnps2d, freq, mean_sig, n_used, fH, profH, fV, profV
        self.dqe_data = None       # dict: fDh, dqe_h, fDv, dqe_v
        self._selector = None
        self._edge_ax = None
        self._roi_patches = []
        self._task_queue = queue.Queue()
        self._last_dir = None
        self._last_dialog_time = None
        self.after(50, self._poll_tasks)

    def _poll_tasks(self):
        """
        Runs on the main thread only, via self.after - the single place
        that drains results background worker threads have put on the
        queue. Worker threads never touch Tkinter/matplotlib directly,
        and never call self.after themselves: calling after() from a
        non-main thread is not reliably safe in Tkinter (it can raise
        "main thread is not in main loop" depending on Python/Tk version),
        so everything crosses back over via this queue instead.
        """
        try:
            while True:
                tag, args = self._task_queue.get_nowait()
                if tag == "edge_loaded":
                    self._edge_loaded(*args)
                elif tag == "flat_loaded":
                    self._flat_loaded(*args)
                elif tag == "progress":
                    self.prog.configure(value=args[0])
                elif tag == "mtf_h_done":
                    self._mtf_h_done(*args)
                elif tag == "mtf_v_done":
                    self._mtf_v_done(*args)
                elif tag == "nnps_done":
                    self._nnps_done(*args)
                elif tag == "dqe_done":
                    self._dqe_done(*args)
        except queue.Empty:
            pass
        self.after(50, self._poll_tasks)

    def _busy(self, on):
        try:
            self.winfo_toplevel().configure(cursor="watch" if on else "")
        except Exception:
            pass
        self.update_idletasks()

    def load_edge_image(self):
        path = simple_open_dicom_dialog(self, initialdir=self._last_dir)
        if not path:
            return
        self._start_edge_load(path)

    def _start_edge_load(self, path):
        self._last_dir = os.path.dirname(path)
        self._busy(True)
        self.edge_file_lbl.configure(text="Loading...")

        def worker():
            try:
                arr, px_mm, ds = load_dicom_pixels(path)
                self._task_queue.put(("edge_loaded", (path, arr, px_mm, ds, None)))
            except Exception as e:
                self._task_queue.put(("edge_loaded", (path, None, None, None, e)))
        threading.Thread(target=worker, daemon=True).start()

    def _edge_loaded(self, path, arr, px_mm, ds, err):
        self._busy(False)
        if err is not None:
            self.edge_file_lbl.configure(text="(none loaded)")
            messagebox.showerror("Load failed", str(err))
            return
        self.edge_arr, self.edge_px_mm, self.edge_ds = arr, px_mm, ds
        self.edge_file_lbl.configure(text=os.path.basename(path))
        # Only the ROI selections reset here - they're pixel coordinates
        # tied to whichever image is currently loaded, so they can't
        # carry over to a new one. Any MTF already calculated (res_h /
        # res_v) is a finished result, not tied to the loaded image, so
        # it's kept: load a horizontal-edge image, measure it, then load
        # a separate vertical-edge image and measure that too, without
        # losing the first result.
        self.roi_h = None; self.roi_v = None
        self.roi_h_lbl.configure(text="Horizontal ROI: not set")
        self.roi_v_lbl.configure(text="Vertical ROI: not set")
        self.mtf_status_lbl.configure(text=self._mtf_status_text())
        self._show_edge_image()

    def load_flat_image(self):
        path = simple_open_dicom_dialog(self, initialdir=self._last_dir)
        if not path:
            return
        self._start_flat_load(path)

    def _start_flat_load(self, path):
        self._last_dir = os.path.dirname(path)
        self._busy(True)
        self.flat_file_lbl.configure(text="Loading...")

        def worker():
            t0 = time.time()
            try:
                arr, px_mm, ds = load_dicom_pixels(path)
                elapsed = time.time() - t0
                print("[timing] load_dicom_pixels (flat): %.2fs" % elapsed,
                     flush=True)
                self._task_queue.put(("flat_loaded", (path, arr, px_mm, ds, None)))
            except Exception as e:
                self._task_queue.put(("flat_loaded", (path, None, None, None, e)))
        threading.Thread(target=worker, daemon=True).start()

    def _flat_loaded(self, path, arr, px_mm, ds, err):
        self._busy(False)
        if err is not None:
            self.flat_file_lbl.configure(text="(none loaded)")
            messagebox.showerror("Load failed", str(err))
            return
        self.flat_arr, self.flat_px_mm, self.flat_ds = arr, px_mm, ds
        self.flat_file_lbl.configure(text=os.path.basename(path))
        self.nnps_data = None
        self.dqe_data = None
        self.nnps_status_lbl.configure(text="NNPS: not done yet")

    def _show_edge_image(self):
        if self.edge_arr is None:
            self.results.clear()
            self.results.head("Input error")
            self.results.warn("Load the edge image first.")
            self.results.top()
            return
        fig = self.plot.reset()
        ax = fig.add_subplot(111)
        sign = getattr(self.edge_ds, "PixelIntensityRelationshipSign", 1)
        cmap = "gray_r" if sign == -1 else "gray"
        ny, nx = self.edge_arr.shape
        ds_factor = max(1, ny // 1000)
        disp = self.edge_arr[::ds_factor, ::ds_factor]
        ax.imshow(disp, cmap=cmap, extent=[0, nx, ny, 0], aspect="equal")
        ax.set_title("Drag a rectangle around a straight section of edge "
                     "(avoid corners)", fontsize=9)
        ax.tick_params(labelsize=7)
        self._edge_ax = ax
        self._roi_patches = []

        def onselect(eclick, erelease):
            x0, x1 = sorted([eclick.xdata, erelease.xdata])
            y0, y1 = sorted([eclick.ydata, erelease.ydata])
            roi = (int(y0), int(y1), int(x0), int(x1))
            if "Horizontal" in self.roi_mode.get():
                self.roi_h = roi
                self.roi_h_lbl.configure(
                    text="Horizontal ROI: rows %d-%d, cols %d-%d" % roi)
            else:
                self.roi_v = roi
                self.roi_v_lbl.configure(
                    text="Vertical ROI: rows %d-%d, cols %d-%d" % roi)
            self._redraw_roi_overlays()

        self._selector = RectangleSelector(
            ax, onselect, useblit=True, button=[1], interactive=True,
            minspanx=10, minspany=10, spancoords="data")
        self._redraw_roi_overlays()
        self.plot.done()

    def _redraw_roi_overlays(self):
        """
        Both ROIs stay visible at once, in different colours, once set -
        RectangleSelector on its own only ever shows the most recent drag,
        so without this a previously-selected ROI becomes invisible (even
        though it's still stored correctly) as soon as you draw the other
        one, making it easy to not notice a bad placement until Calculate
        fails.
        """
        ax = getattr(self, "_edge_ax", None)
        if ax is None:
            return
        for p in list(getattr(self, "_roi_patches", [])):
            p.remove()
        self._roi_patches = []
        for roi, color, label in ((self.roi_h, ACCENT, "horizontal ROI"),
                                  (self.roi_v, "#c1666b", "vertical ROI")):
            if roi is None:
                continue
            r0, r1, c0, c1 = roi
            patch = mpatches.Rectangle((c0, r0), c1 - c0, r1 - r0,
                                       fill=False, edgecolor=color,
                                       linewidth=1.5, linestyle="-",
                                       label=label)
            ax.add_patch(patch)
            self._roi_patches.append(patch)
        if self.roi_h is not None or self.roi_v is not None:
            ax.legend(fontsize=7, loc="upper right")
        self.plot.canvas.draw_idle()

    def _get_spec_pts(self, text_widget):
        spec_pts = []
        for line in text_widget.get("1.0", "end").splitlines():
            line = line.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.replace("\t", ",").split(",")
                    if p.strip()]
            if len(parts) >= 2:
                try:
                    spec_pts.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    pass
        return spec_pts

    def _mtf_status_text(self):
        h = "done" if self.res_h is not None else "not done"
        v = "done" if self.res_v is not None else "not done"
        return "MTF: horizontal %s, vertical %s" % (h, v)

    def calc_mtf_h(self):
        r = self.results; r.clear()
        if self.edge_arr is None:
            r.head("Input error"); r.warn("Load the edge image first.")
            return
        if self.roi_h is None:
            r.head("Input error")
            r.warn("Select a horizontal-MTF ROI (drag a horizontal edge) "
                   "on the edge image first.")
            return
        self._busy(True)
        roi, arr, px_mm = self.roi_h, self.edge_arr, self.edge_px_mm

        def worker():
            try:
                r0, r1, c0, c1 = roi
                res = slanted_edge_mtf(arr[r0:r1, c0:c1], px_mm,
                                       edge_axis="horizontal")
                self._task_queue.put(("mtf_h_done", (res, None)))
            except Exception as e:
                self._task_queue.put(("mtf_h_done", (None, e)))
        threading.Thread(target=worker, daemon=True).start()

    def _mtf_h_done(self, res, err):
        self._busy(False)
        r = self.results
        if err is not None:
            r.head("Error"); r.warn("%s: %s" % (type(err).__name__, err))
            return
        self.res_h = res
        self.mtf_status_lbl.configure(text=self._mtf_status_text())
        r.head("Horizontal MTF")
        r.row("Edge angle", fmt(res["angle_deg"]), "deg")
        if not (0.3 <= abs(res["angle_deg"]) <= 5):
            r.warn("An edge angle outside about 0.3-5 degrees is outside "
                   "the range this tilt-sampling technique is designed "
                   "for - check the ROI straddles a single straight edge "
                   "segment, away from a corner.")
        for f0 in (0.5, 1.0, 2.0, 3.0, 4.0):
            m = np.interp(f0, res["freq"], res["mtf"])
            r.row("%.1f lp/mm" % f0, "%.3f" % m)
        r.top()
        self._draw_mtf_single(res, "horizontal")

    def calc_mtf_v(self):
        r = self.results; r.clear()
        if self.edge_arr is None:
            r.head("Input error"); r.warn("Load the edge image first.")
            return
        if self.roi_v is None:
            r.head("Input error")
            r.warn("Select a vertical-MTF ROI (drag a vertical edge) "
                   "on the edge image first.")
            return
        self._busy(True)
        roi, arr, px_mm = self.roi_v, self.edge_arr, self.edge_px_mm

        def worker():
            try:
                r0, r1, c0, c1 = roi
                res = slanted_edge_mtf(arr[r0:r1, c0:c1], px_mm,
                                       edge_axis="vertical")
                self._task_queue.put(("mtf_v_done", (res, None)))
            except Exception as e:
                self._task_queue.put(("mtf_v_done", (None, e)))
        threading.Thread(target=worker, daemon=True).start()

    def _mtf_v_done(self, res, err):
        self._busy(False)
        r = self.results
        if err is not None:
            r.head("Error"); r.warn("%s: %s" % (type(err).__name__, err))
            return
        self.res_v = res
        self.mtf_status_lbl.configure(text=self._mtf_status_text())
        r.head("Vertical MTF")
        r.row("Edge angle", fmt(res["angle_deg"]), "deg")
        if not (0.3 <= abs(res["angle_deg"]) <= 5):
            r.warn("An edge angle outside about 0.3-5 degrees is outside "
                   "the range this tilt-sampling technique is designed "
                   "for - check the ROI straddles a single straight edge "
                   "segment, away from a corner.")
        for f0 in (0.5, 1.0, 2.0, 3.0, 4.0):
            m = np.interp(f0, res["freq"], res["mtf"])
            r.row("%.1f lp/mm" % f0, "%.3f" % m)
        r.top()
        self._draw_mtf_single(res, "vertical")

    def calc_nnps(self):
        r = self.results; r.clear()
        if self.flat_arr is None:
            r.head("Input error"); r.warn("Load the flat-field image first.")
            return
        self._busy(True)
        self.prog["value"] = 0
        arr, px_mm = self.flat_arr, self.flat_px_mm

        def progress(f):
            self._task_queue.put(("progress", (100 * f,)))

        def worker():
            try:
                nnps2d, freqN, mean_sig, n_used = flat_field_nnps(
                    arr, px_mm, progress=progress)
                fH, profH = nnps_profile(nnps2d, freqN, axis="h")
                fV, profV = nnps_profile(nnps2d, freqN, axis="v")
                self._task_queue.put(("nnps_done", (
                    nnps2d, freqN, mean_sig, n_used, fH, profH, fV, profV,
                    None)))
            except Exception as e:
                self._task_queue.put(("nnps_done", (
                    None, None, None, None, None, None, None, None, e)))
        threading.Thread(target=worker, daemon=True).start()

    def _nnps_done(self, nnps2d, freqN, mean_sig, n_used, fH, profH, fV,
                   profV, err):
        self._busy(False)
        r = self.results
        if err is not None:
            r.head("Error"); r.warn("%s: %s" % (type(err).__name__, err))
            return
        self.prog["value"] = 100
        self.nnps_data = {"nnps2d": nnps2d, "freq": freqN,
                          "mean_sig": mean_sig, "n_used": n_used,
                          "fH": fH, "profH": profH, "fV": fV, "profV": profV}
        self.nnps_status_lbl.configure(text="NNPS: done (%d ROIs)" % n_used)
        r.head("NNPS")
        r.row("ROIs used", "%d" % n_used)
        r.row("Mean signal", fmt(mean_sig), "DN")
        for f0 in (0.5, 1.0, 2.0, 3.0):
            nh = np.interp(f0, fH[fH >= 0], profH[fH >= 0])
            nv = np.interp(f0, fV[fV >= 0], profV[fV >= 0])
            r.row("%.1f lp/mm" % f0, "h=%.3e  v=%.3e" % (nh, nv))
        r.top()
        self._draw_nnps(fH, profH, fV, profV)

    def calc_dqe(self):
        r = self.results; r.clear()
        try:
            if self.res_h is None or self.res_v is None:
                raise ValueError("Calculate both horizontal and vertical "
                                 "MTF first.")
            if self.nnps_data is None:
                raise ValueError("Calculate NNPS first.")
            K = self.kerma.get(0.0)
            if K <= 0:
                raise ValueError("Enter the measured air kerma at the "
                                 "detector for the flat-field exposure.")
            q_name = self.quality.get()
            q = RQA_SNR_IN2.get(q_name)
            if q is None:
                raise ValueError("%s has no mandatory SNR_in^2 value in "
                                 "Table 3 (only RQA3/5/7/9 do)." % q_name)
        except ValueError as e:
            r.head("Input error"); r.warn(str(e))
            return

        res_h, res_v = self.res_h, self.res_v
        nd = self.nnps_data
        spec_pts = self._get_spec_pts(self.spec_text_dqe)
        self._busy(True)

        def worker():
            try:
                fDh, dqe_h = compute_dqe(res_h["freq"], res_h["mtf"],
                                         nd["fH"], nd["profH"], q, K)
                fDv, dqe_v = compute_dqe(res_v["freq"], res_v["mtf"],
                                         nd["fV"], nd["profV"], q, K)
                self._task_queue.put(("dqe_done", (
                    fDh, dqe_h, fDv, dqe_v, q_name, q, K, spec_pts, None)))
            except Exception as e:
                self._task_queue.put(("dqe_done", (
                    None, None, None, None, None, None, None, None, e)))
        threading.Thread(target=worker, daemon=True).start()

    def _dqe_done(self, fDh, dqe_h, fDv, dqe_v, q_name, q, K, spec_pts, err):
        self._busy(False)
        r = self.results
        if err is not None:
            r.head("Error"); r.warn("%s: %s" % (type(err).__name__, err))
            return
        self.dqe_data = {"fDh": fDh, "dqe_h": dqe_h, "fDv": fDv, "dqe_v": dqe_v}
        r.head("DQE")
        r.row("Air kerma", fmt(K), "uGy")
        r.row("q (%s, Table 3)" % q_name, fmt(q, 1), "mm^-2.uGy^-1")
        for f0 in (0.5, 1.0, 2.0, 3.0):
            dh = np.interp(f0, fDh, dqe_h)
            dv = np.interp(f0, fDv, dqe_v)
            r.row("%.1f lp/mm" % f0, "h=%.3f  v=%.3f" % (dh, dv))
        r.note("Single edge per direction (not IEC's optional two-edge "
               "average), and detector response assumed linear through "
               "zero from this one flat-field exposure - see the Notes "
               "tab.")
        r.top()
        self._draw_dqe(fDh, dqe_h, fDv, dqe_v)

    def _draw_mtf_single(self, res, direction):
        fig = self.plot.reset()
        ax = fig.add_subplot(111)
        ax.plot(res["freq"], res["mtf"], color=ACCENT, label=direction)
        spec_pts = self._get_spec_pts(self.spec_text_mtf)
        if spec_pts:
            sx = [p[0] for p in spec_pts]
            sy = [p[1] / 100.0 for p in spec_pts]
            ax.plot(sx, sy, "ko--", ms=4, lw=1, label="factory spec")
        ax.set_xlim(0, min(res["freq"].max(), 6)); ax.set_ylim(0, 1.05)
        ax.set_xlabel("lp/mm"); ax.set_ylabel("MTF")
        ax.legend(fontsize=8); ax.grid(alpha=0.25)
        ax.set_title("%s MTF (edge angle %.2f deg)"
                     % (direction.capitalize(), res["angle_deg"]),
                     color=ACCENT, fontweight="bold")
        self.plot.done()

    def _draw_nnps(self, fH, profH, fV, profV):
        fig = self.plot.reset()
        ax = fig.add_subplot(111)
        ax.semilogy(fH[fH >= 0], profH[fH >= 0], color=ACCENT, label="horizontal")
        ax.semilogy(fV[fV >= 0], profV[fV >= 0], color="#c1666b", label="vertical")
        spec_pts = self._get_spec_pts(self.spec_text_nnps)
        if spec_pts:
            sx = [p[0] for p in spec_pts]
            sy = [p[1] for p in spec_pts]
            ax.semilogy(sx, sy, "ko--", ms=4, lw=1, label="factory spec")
        ax.set_xlim(0, min(fH.max(), 6))
        ax.set_xlabel("lp/mm"); ax.set_ylabel("NNPS (mm^2)")
        ax.legend(fontsize=8); ax.grid(alpha=0.25, which="both")
        ax.set_title("NNPS", color=ACCENT, fontweight="bold")
        self.plot.done()

    def _draw_dqe(self, fDh, dqe_h, fDv, dqe_v):
        fig = self.plot.reset()
        ax = fig.add_subplot(111)
        ax.plot(fDh, dqe_h, color=ACCENT, label="horizontal")
        ax.plot(fDv, dqe_v, color="#c1666b", label="vertical")
        spec_pts = self._get_spec_pts(self.spec_text_dqe)
        if spec_pts:
            sx = [p[0] for p in spec_pts]
            sy = [p[1] / 100.0 for p in spec_pts]
            ax.plot(sx, sy, "ko--", ms=4, lw=1, label="factory spec")
        ax.set_xlim(0, min(fDh.max(), 6))
        ax.set_ylim(0, max(0.05, 1.05 * np.nanmax(dqe_h)))
        ax.set_xlabel("lp/mm"); ax.set_ylabel("DQE")
        ax.legend(fontsize=8); ax.grid(alpha=0.25)
        ax.set_title("DQE", color=ACCENT, fontweight="bold")
        self.plot.done()

    def _copy_to_clipboard(self, grid, values, label, n_rows, fmax_cap=None):
        r = self.results
        valid = np.isfinite(values)
        if valid.sum() < 2:
            r.clear(); r.head("Copy"); r.warn("%s has no usable data yet." % label)
            return
        fmax = float(np.max(grid[valid]))
        if fmax_cap is not None:
            fmax = min(fmax, fmax_cap)
        out_grid = np.linspace(0.0, fmax, n_rows)
        out_vals = np.interp(out_grid, grid[valid], values[valid])
        text = "\n".join("%.6g\t%.6g" % (f, v)
                         for f, v in zip(out_grid, out_vals))
        self.clipboard_clear()
        self.clipboard_append(text)
        self.update()
        r.clear()
        r.head("Copy")
        r.row("Copied", "%s - %d rows, 0 to %.3f lp/mm" % (label, n_rows, fmax))
        r.top()

    def copy_curve(self, attr_name, n_rows, label):
        res = getattr(self, attr_name, None)
        if res is None:
            self.results.clear(); self.results.head("Copy")
            self.results.warn("%s hasn't been calculated yet." % label)
            return
        # Cap at the detector's real Nyquist (1/(2*pixel pitch)) - the MTF
        # array itself extends to the OVERSAMPLED grid's Nyquist (~8x
        # higher with the default oversample=8), which is a real array of
        # numbers but meaningless past the detector's actual resolution
        # limit: an artifact of the oversampling method, not real detector
        # response, and not something that belongs in an exported table.
        nyquist = 0.5 / self.edge_px_mm if self.edge_px_mm else None
        self._copy_to_clipboard(res["freq"], res["mtf"], label, n_rows,
                                fmax_cap=nyquist)

    def copy_nnps(self, direction):
        if self.nnps_data is None:
            self.results.clear(); self.results.head("Copy")
            self.results.warn("NNPS hasn't been calculated yet.")
            return
        nd = self.nnps_data
        if direction == "h":
            grid, vals = nd["fH"], nd["profH"]
        else:
            grid, vals = nd["fV"], nd["profV"]
        pos = grid >= 0
        self._copy_to_clipboard(grid[pos], vals[pos],
                                "NNPS (%s)" % direction, 63)

    def copy_dqe(self):
        """
        A single DQE(f) curve, combining horizontal and vertical - per
        IEC 62220-1, DQE(u) and DQE(v) may be averaged into one value
        ONLY where they agree closely (the standard's own criterion:
        DQE(u)/DQE(v) between 0.9 and 1.1). Outside that band the two
        directions are telling you something genuinely different
        (anisotropy), so averaging there would hide real information
        rather than simplify it - flagged rather than silently done.
        """
        if self.dqe_data is None:
            self.results.clear(); self.results.head("Copy")
            self.results.warn("DQE hasn't been calculated yet.")
            return
        dd = self.dqe_data
        fDh, dqe_h, fDv, dqe_v = dd["fDh"], dd["dqe_h"], dd["fDv"], dd["dqe_v"]

        fmax = min(np.nanmax(fDh[np.isfinite(dqe_h)]),
                  np.nanmax(fDv[np.isfinite(dqe_v)]))
        n_rows = 63
        grid = np.linspace(0.0, fmax, n_rows)
        h_i = np.interp(grid, fDh, dqe_h)
        v_i = np.interp(grid, fDv, dqe_v)
        avg = 0.5 * (h_i + v_i)

        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(v_i > 0, h_i / v_i, np.nan)
        out_of_range = np.isfinite(ratio) & ((ratio < 0.9) | (ratio > 1.1))

        text = "\n".join("%.6g\t%.6g" % (f, v) for f, v in zip(grid, avg))
        self.clipboard_clear()
        self.clipboard_append(text)
        self.update()

        r = self.results; r.clear()
        r.head("Copy DQE (h/v average)")
        r.row("Copied", "%d rows, 0 to %.3f lp/mm" % (n_rows, fmax))
        if out_of_range.any():
            bad_freqs = grid[out_of_range]
            r.warn("DQE(h)/DQE(v) falls outside IEC 62220-1's 0.9-1.1 "
                   "averaging criterion at %d of %d frequencies (from "
                   "%.2f lp/mm up) - the two directions diverge there, "
                   "so the averaged value copied at those frequencies is "
                   "less meaningful than at the lower frequencies where "
                   "they agree. Copied anyway since you asked for one "
                   "curve; worth checking the DQE plot's h/v split "
                   "before relying on the tail of this table."
                   % (out_of_range.sum(), n_rows, float(bad_freqs.min())))
        else:
            r.note("DQE(h) and DQE(v) agree within IEC's 0.9-1.1 "
                   "averaging criterion across the whole copied range.")
        r.top()


class NotesTab(ttk.Frame):
    title = "  Notes  "
    TEXT = """
ASSUMPTIONS AND LIMITATIONS
===========================
Read this before quoting any number from this program.

WHAT IS AND IS NOT MONTE CARLO
  Spectrum generation is NOT Monte Carlo. SpekCalc and SpekPy are
  semi-analytical. Poludniowski used Monte Carlo once, upstream, to derive
  the electron energy and depth distribution in a tungsten target, then
  fitted an analytical bremsstrahlung expression to it. The tool evaluates
  that expression in milliseconds.

  Filtering a spectrum is NOT Monte Carlo either. It is Beer-Lambert applied
  bin by bin, and it gives NARROW-BEAM attenuation.

  Only the transport step is Monte Carlo, and it exists to recover the
  scattered photons that narrow-beam attenuation discards.

VALIDATION PERFORMED
  With scattering disabled, the transport estimator reproduces exp(-mu.d)
  to machine precision at 60, 100 and 150 keV and to 3 mm of lead. The
  Klein-Nishina sampler satisfies Compton kinematics to 1e-6 and respects
  the backscatter energy limit. Buildup exceeds 1 everywhere and increases
  monotonically with thickness, as it must.

COHERENT SCATTER USES THE REAL ATOMIC FORM FACTOR
  Rayleigh scattering angle is sampled from the tabulated cumulative of
  x.F(x)^2, using the real atomic form factor from xsdata.npz, not a
  bracketing approximation. Validated against exact integration to better
  than 0.0021 in mean cos(theta) across four materials and 40-120 keV.
  The remaining approximation is the free-electron Klein-Nishina treatment
  of incoherent scattering below - no binding correction - which is poor
  below about 30 keV. If that matters for your purpose, this program is
  not accurate enough and you want a validated code such as OpenMC,
  EGSnrc or PENELOPE.

OTHER SIMPLIFICATIONS
  * Incoherent scattering uses the free-electron Klein-Nishina cross
    section with no binding correction. This is reasonable above about
    30 keV and poor below it.
  * K-shell fluorescence is tracked (toggle: fluorescence=True/False) for
    materials with real K-edge data - originally lead, lead glass and
    tungsten, now every added element with a K-edge inside 1-800 keV too.
    For the other curated compounds (concrete, gypsum, steel, glass,
    borosilicate, PMMA, water, tissue, air) photoelectric absorption is
    still terminal, since no K data was generated for those mixtures.
  * Geometry is an infinite slab at normal incidence. No obliquity, no
    finite field, no room return.
  * HVL is air-kerma defined throughout, which is the conventional
    definition and what a survey meter reports. An energy-fluence HVL is
    available in the API but differs by tens of per cent and must not be
    compared against a measurement.
  * Transmission can be scored as a true AIR-KERMA ratio, using mass
    energy-absorption coefficients for dry air taken verbatim from the NIST
    Hubbell and Seltzer tables. An energy-fluence score is also offered for
    comparison. Air kerma is the quantity shielding calculations use.

CROSS SECTION DATA
  xsdata.npz holds photoelectric, coherent and incoherent mass attenuation
  coefficients for 15 curated shielding/filtration materials plus every
  element Z=1-92, all from 1 to 800 keV, derived from the Elam tables via
  xraydb. Verified to 0.011 per cent against xraydb and against published
  NIST values - water gives 0.1707 cm2/g at 100 keV, lead 1.77 cm^-1 at
  511 keV. The Elam tables are flagged unreliable above 800 keV, which is
  why the table stops there.

  The 15 curated materials have absorption edges explicitly resolved by
  extra points on the energy grid; the added elements sit on that same
  grid without extra points at their own edges, so their K-edge is a
  little more smoothed by interpolation than lead's or tungsten's is.
  K-fluorescence yield for the added elements is a jump-ratio
  approximation (P_K x omega_K) rather than whatever compilation produced
  the lead/tungsten values, reproducing those two to about 3 per cent, not
  exactly - a second-order effect on an already-optional channel.

REGULATORY STATUS
  None. Shielding guidance built on NCRP 147 and BIR 2012 requires the
  Archer empirical model and its published attenuation data. This program
  is a cross-check and a teaching tool. The shielding calculator remains
  the authority for design work.
"""

    def __init__(self, parent):
        super().__init__(parent, padding=10)
        t = tk.Text(self, wrap="word", font=MONO, background=CARD,
                    relief="solid", borderwidth=1, padx=14, pady=10)
        sb = ttk.Scrollbar(self, orient="vertical", command=t.yview)
        t.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y"); t.pack(fill="both", expand=True)
        t.insert("1.0", self.TEXT)
        t.configure(state="disabled")


def size_window(root, want_w=1480, want_h=1000, min_w=1100, min_h=740):
    """
    Open at a comfortable size, but never larger than the screen.

    A fixed geometry string is fine until it meets a 1366x768 laptop, where
    it silently puts the bottom of the window under the taskbar. Clamping to
    the actual screen and centring avoids that.
    """
    root.update_idletasks()
    sw = root.winfo_screenwidth()
    sh = root.winfo_screenheight()
    w = min(want_w, max(sw - 80, 640))
    h = min(want_h, max(sh - 120, 480))
    x = max((sw - w) // 2, 0)
    y = max((sh - h) // 3, 0)
    root.geometry("%dx%d+%d+%d" % (w, h, x, y))
    root.minsize(min(min_w, w), min(min_h, h))
    return w, h


def _set_windows_app_id(app_id="spectrum.montecarlo.transmission"):
    if sys.platform != "win32":
        return
    try:
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(app_id)
    except Exception:
        pass


def _apply_icon(root):
    for name in ("icon.ico", "shielding.ico"):
        p = os.path.join(_HERE, name)
        if os.path.exists(p):
            try:
                root.iconbitmap(default=p)
                return
            except tk.TclError:
                pass


def main():
    global XS
    _set_windows_app_id()
    root = tk.Tk()
    root.title("%s  v%s" % (APP_TITLE, APP_VERSION))
    size_window(root, want_w=1480, want_h=1000, min_w=1100, min_h=740)
    root.configure(background=BG)
    _apply_icon(root)
    # Windows creates the taskbar button lazily - setting the icon once
    # before the window is mapped often reaches the title bar fine but
    # is ignored for the taskbar button specifically. Reapplying shortly
    # after the window is actually on screen fixes that in practice.
    root.after(250, lambda: _apply_icon(root))

    try:
        XS = CrossSections()
    except FileNotFoundError as e:
        messagebox.showerror("Missing data file", str(e))
        sys.exit(1)

    # An older xsdata.npz will not have every material this build offers.
    # Trim the lists to what is actually present and say so, rather than
    # letting a KeyError surface when someone picks the missing one.
    missing = ([m for m in FILTER_MATERIALS if m not in XS.materials]
               + [m for m in BARRIER_MATERIALS if m not in XS.materials])
    FILTER_MATERIALS[:] = [m for m in FILTER_MATERIALS if m in XS.materials]
    BARRIER_MATERIALS[:] = [m for m in BARRIER_MATERIALS if m in XS.materials]
    if missing:
        messagebox.showwarning(
            "Data file is out of date",
            "xsdata.npz does not contain:\n\n    "
            + "\n    ".join(sorted(MATERIAL_LABEL.get(m, m)
                                    for m in set(missing)))
            + "\n\nThose materials have been hidden. Replace xsdata.npz "
              "with the version that came with this build to enable them.")

    style = ttk.Style(root)
    for th in ("vista", "winnative", "clam", "default"):
        if th in style.theme_names():
            style.theme_use(th); break
    for w in ("TFrame", "TLabelframe", "TLabel", "TCheckbutton", "TRadiobutton"):
        try:
            style.configure(w, background=BG)
        except tk.TclError:
            pass
    style.configure("TLabelframe.Label", background=BG,
                    font=("Segoe UI", 9, "bold"), foreground=ACCENT)
    style.configure("Accent.TButton", font=("Segoe UI", 10, "bold"))

    bar = ttk.Frame(root, padding=(10, 4)); bar.pack(fill="x")
    ttk.Label(bar, text="spectrum model: %s     transport: %s     "
                        "cross sections: Elam tables, 1-800 keV"
                        % ("SpekPy" if HAVE_SPEKPY else "built-in approximation",
                           "numba" if HAVE_NUMBA else "pure python (slow)"),
              foreground="#666").pack(side="left")

    nb = ttk.Notebook(root); nb.pack(fill="both", expand=True, padx=8, pady=(4, 4))
    stab = SpectrumTab(nb)
    nb.add(stab, text=SpectrumTab.title)
    nb.add(MonteCarloTab(nb, stab), text=MonteCarloTab.title)
    nb.add(ImageQualityTab(nb), text=ImageQualityTab.title)
    nb.add(NotesTab(nb), text=NotesTab.title)

    ttk.Label(root, anchor="w", foreground="#666",
              text="Research and cross-check tool. Not for design "
                   "submissions - see the Notes tab.").pack(fill="x", padx=12,
                                                            pady=(0, 6))

    def crash(exc, val, tb):
        traceback.print_exception(exc, val, tb)
        messagebox.showerror("Unexpected error", "%s: %s" % (exc.__name__, val))
    root.report_callback_exception = crash
    root.mainloop()


if __name__ == "__main__":
    main()