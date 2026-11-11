/* Application shell - port of main() */
"use strict";

const APP_TITLE = "X-ray Spectrum and Monte Carlo Transmission", APP_VERSION = "1.0";

class App {
  async start() {
    const status = document.getElementById("status");
    status.textContent = "starting the Rust/WebAssembly engine...";
    this.engine = new Engine(b64ToBytes(WASM_B64));
    let nw;
    try { nw = await this.engine.start(); }
    catch (e) { status.textContent = "Could not start the engine: " + e.message; throw e; }
    this.buildTabs();
    const tables = await this.engine.call({ cmd: "tables" });
    for (const t of [this.stab, this.mtab, this.itab]) t.setTables(tables);
    this.spekRange = tables.spek_range;
    this.mtab.workersHint.textContent = "Transport runs in Rust/WebAssembly on " + nw + " parallel worker" + (nw > 1 ? "s" : "") + " - one thickness point per worker.";
    const npz = (typeof NPZ_B64 === "string" && NPZ_B64.length) ? b64ToBytes(NPZ_B64) : null;
    if (npz) await this.loadXS(npz, "built-in xsdata.npz");
    else this.banner("xsdata.npz is not embedded in this page. Use \"Load xsdata.npz\" (top right) to load it.", true);
    this.statusText();
  }
  statusText() {
    const sr = this.spekRange, nw = this.engine.workers.length;
    document.getElementById("status").textContent =
      "spectrum model: SpekPy (pre-computed, " + sr[0] + "-" + sr[1] + " kVp)     transport: Rust/WebAssembly" +
      (this.engine.inline ? " (single thread)" : ", " + nw + " workers") + "     cross sections: " + (this.xsName || "not loaded") +
      (this.xs ? ", Elam tables, " + Number(this.xs.emin.toPrecision(6)) + "-" + Number(this.xs.emax.toPrecision(6)) + " keV" : "");
  }
  banner(msg, warn) {
    const b = document.getElementById("banner"); b.innerHTML = "";
    if (!msg) { b.style.display = "none"; return; }
    b.append(el("span", null, msg), button("Dismiss", () => { b.style.display = "none"; }, "small"));
    b.className = warn ? "warn" : ""; b.style.display = "";
  }
  async loadXS(bytes, name) {
    try {
      const info = await this.engine.loadXS(bytes);
      this.xs = info; this.xsName = name;
      const have = new Set(info.materials.map((m) => m.key));
      const all = _BASE_MATERIALS.concat(_ELEMENTS_BY_LABEL);
      const missing = all.filter((m) => !have.has(m));
      FILTER_MATERIALS = all.filter((m) => have.has(m)); BARRIER_MATERIALS = all.filter((m) => have.has(m));
      this.stab.setMaterials(); this.mtab.setMaterials();
      if (missing.length) this.banner("xsdata.npz does not contain: " + missing.map(LBL).sort().join(", ") + ". Those materials have been hidden. Replace xsdata.npz with the version that came with this build to enable them.", true);
      else this.banner(null);
      this.statusText();
      this.mtab.updateArcher();
      await this.stab.run();
    } catch (e) { this.banner("Could not load " + name + ": " + e.message, true); }
  }
  buildTabs() {
    const tabs = [["spectrum", "Spectrum"], ["mc", "Monte Carlo"], ["iq", "MTF / NNPS / DQE"], ["notes", "Notes"]];
    const bar = document.getElementById("tabs"), host = document.getElementById("pages");
    const pages = {};
    for (const [k, t] of tabs) {
      const b = el("button", { class: "tab", role: "tab", "data-k": k }, t);
      b.onclick = () => this.show(k);
      bar.append(b);
      pages[k] = el("section", { class: "page", role: "tabpanel", "data-k": k }); host.append(pages[k]);
    }
    this.pages = pages;
    this.stab = new SpectrumTab(pages.spectrum, this);
    this.mtab = new MonteCarloTab(pages.mc, this, this.stab);
    this.itab = new ImageQualityTab(pages.iq, this);
    pages.notes.append(el("pre", { class: "notes" }, NOTES_TEXT));
    let last = "spectrum"; try { last = localStorage.getItem("xray.tab") || "spectrum"; } catch (e) { /* ignore */ }
    this.show(pages[last] ? last : "spectrum");
    const inp = el("input", { type: "file", accept: ".npz", style: "display:none" });
    inp.onchange = async () => { const f = inp.files[0]; if (f) await this.loadXS(new Uint8Array(await f.arrayBuffer()), f.name); inp.value = ""; };
    const lb = button("Load xsdata.npz", () => inp.click(), "small");
    lb.title = "Replace the built-in cross-section data with another xsdata.npz";
    document.getElementById("hdr-actions").append(inp, lb);
  }
  show(k) {
    for (const b of document.querySelectorAll(".tab")) b.setAttribute("aria-selected", b.dataset.k === k ? "true" : "false");
    for (const [kk, p] of Object.entries(this.pages)) p.style.display = kk === k ? "" : "none";
    try { localStorage.setItem("xray.tab", k); } catch (e) { /* ignore */ }
  }
}

window.addEventListener("error", (e) => { const s = document.getElementById("status"); if (s) s.textContent = "Unexpected error: " + e.message; });
window.addEventListener("DOMContentLoaded", () => { document.title = APP_TITLE; window.__app = new App(); window.__app.start(); });
