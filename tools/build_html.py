#!/usr/bin/env python3
"""
Build the single-file browser app (app.html) from the Rust crate and web/.

    rustup target add wasm32-unknown-unknown
    python3 tools/build_html.py            # builds the wasm, writes app.html

The page embeds the WebAssembly module and xsdata.npz (base64), so it opens
straight from disk with no server and no network access. Pass --no-npz to
leave xsdata.npz out; the page then asks for it with a file picker.
"""
import base64, os, subprocess, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
WEB = os.path.join(ROOT, "web")
JS_ORDER = ["plot.js", "util.js", "engine.js", "notes.js", "app_spectrum.js", "app_mc.js", "app_iq.js", "main.js"]

if "--skip-cargo" not in sys.argv:
    subprocess.run(["cargo", "build", "--release", "--target", "wasm32-unknown-unknown", "--lib"],
                   cwd=os.path.join(ROOT, "rust"), check=True)
wasm = open(os.path.join(ROOT, "rust", "target", "wasm32-unknown-unknown", "release", "xray_app.wasm"), "rb").read()
npz = b"" if "--no-npz" in sys.argv else open(os.path.join(ROOT, "xsdata.npz"), "rb").read()

js = "\n".join(open(os.path.join(WEB, "src", f), encoding="utf-8").read() for f in JS_ORDER)
js = ('const WASM_B64 = "%s";\nconst NPZ_B64 = "%s";\n' % (base64.b64encode(wasm).decode(), base64.b64encode(npz).decode())) + js
css = open(os.path.join(WEB, "style.css"), encoding="utf-8").read()
html = open(os.path.join(WEB, "template.html"), encoding="utf-8").read()
html = html.replace("/*__CSS__*/", css).replace("/*__JS__*/", js.replace("</script", "<\\/script"))
out = os.path.join(ROOT, "app.html")
with open(out, "w", encoding="utf-8") as f:
    f.write(html)
print("wrote %s (%.1f MB: wasm %.0f kB, xsdata.npz %.0f kB)" % (out, os.path.getsize(out) / 1e6, len(wasm) / 1e3, len(npz) / 1e3))
