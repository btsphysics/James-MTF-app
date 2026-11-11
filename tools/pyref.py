"""
Load the computational parts of app.py (sections 1-3 and 6) without its
Tk/matplotlib GUI, so the Rust port can be checked against the original.
"""
import math, os, sys, types
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = open(os.path.join(HERE, "..", "app.py"), encoding="utf-8").read().splitlines(True)

def _slice(start_marker, end_marker):
    a = next(i for i, l in enumerate(SRC) if start_marker in l)
    b = next(i for i, l in enumerate(SRC) if end_marker in l and i > a)
    return "".join(SRC[a:b])

code = (_slice("# numpy 2.0 renamed trapz", "#  SECTION 4 -- GUI")
        + _slice("def load_dicom_pixels(path):", "def simple_open_dicom_dialog"))
g = {"__name__": "pyref", "__file__": os.path.join(HERE, "..", "app.py"),
     "np": np, "math": math, "os": os, "sys": sys}
def njit(*a, **k):
    if a and callable(a[0]):
        return a[0]
    return lambda f: f
g["njit"] = njit
try:
    import spekpy as sp
    g["sp"] = sp; g["HAVE_SPEKPY"] = True
except Exception:
    g["HAVE_SPEKPY"] = False
try:
    import pydicom
    g["pydicom"] = pydicom; g["HAVE_PYDICOM"] = True
except Exception:
    g["HAVE_PYDICOM"] = False
exec(compile(code, "app.py(extract)", "exec"), g)
g["XS"] = g["CrossSections"](os.path.join(HERE, "..", "xsdata.npz"))
ref = types.SimpleNamespace(**g)
