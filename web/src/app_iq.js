/* MTF / NNPS / DQE tab - port of ImageQualityTab */
"use strict";

class ImageQualityTab {
  constructor(root, app) {
    this.app = app;
    root.append(el("p", { class: "intro" }, "MTF (slanted-edge), NNPS (flat field) and DQE, per IEC 62220-1 methodology."));
    const body = el("div", { class: "body" }); root.append(body);
    const left = el("div", { class: "left" }); const right = el("div", { class: "right" }); body.append(left, right);

    this.edgeInput = el("input", { type: "file", accept: ".dcm,.dicom,application/dicom,*/*", style: "display:none" });
    this.edgeInput.onchange = () => { if (this.edgeInput.files[0]) this.loadEdge(this.edgeInput.files[0]); this.edgeInput.value = ""; };
    this.edgeLbl = el("div", { class: "hint" }, "(none loaded)");
    this.roiMode = new Combo("Now selecting", [["h", "Horizontal MTF - drag a HORIZONTAL edge"], ["v", "Vertical MTF - drag a VERTICAL edge"]], "h", { stack: true });
    this.roiHLbl = el("div", { class: "hint" }, "Horizontal-MTF ROI: not set");
    this.roiVLbl = el("div", { class: "hint" }, "Vertical-MTF ROI: not set");
    this.mtfStatus = el("div", { class: "hint" }, "MTF: neither done yet");
    left.append(section("1. Edge image (MTF)", this.edgeInput, button("Load MTF tool image", () => this.edgeInput.click(), "wide"), this.edgeLbl,
      button("View edge image / select ROIs", () => this.showEdge(), "wide"), this.roiMode.root,
      hint("Left-drag to draw an ROI. Wheel to zoom; right-drag or Shift+drag to pan."), this.roiHLbl, this.roiVLbl,
      el("div", { class: "btnrow" }, button("Calculate horizontal MTF", () => this.calcMTF("h")), button("Copy", () => this.copyCurve("h"), "small")),
      el("div", { class: "btnrow" }, button("Calculate vertical MTF", () => this.calcMTF("v")), button("Copy", () => this.copyCurve("v"), "small")),
      this.mtfStatus));

    this.flatInput = el("input", { type: "file", accept: ".dcm,.dicom,application/dicom,*/*", style: "display:none" });
    this.flatInput.onchange = () => { if (this.flatInput.files[0]) this.loadFlat(this.flatInput.files[0]); this.flatInput.value = ""; };
    this.flatLbl = el("div", { class: "hint" }, "(none loaded)");
    this.nnpsStatus = el("div", { class: "hint" }, "NNPS: not done yet");
    left.append(section("2. Flat field (NNPS)", this.flatInput, button("Load NNPS image", () => this.flatInput.click(), "wide"), this.flatLbl,
      button("Calculate NNPS", () => this.calcNNPS(), "wide"),
      el("div", { class: "btnrow" }, button("Copy NNPS (h)", () => this.copyNNPS("h")), button("Copy NNPS (v)", () => this.copyNNPS("v"))),
      this.nnpsStatus));

    this.kerma = new Entry("Air kerma at detector", 0.0, "uGy");
    this.quality = new Combo("RQA quality", [], "RQA5");
    left.append(section("3. DQE", this.kerma.root, this.quality.root,
      hint("q (SNR_in^2) comes straight from the IEC 62220-1 Table 3 data already in the app - only RQA3/5/7/9 carry a mandatory value."),
      button("Calculate DQE", () => this.calcDQE(), "accent big"), button("Copy DQE", () => this.copyDQE(), "wide")));

    const ta = () => el("textarea", { rows: "4", spellcheck: "false" });
    this.specMTF = ta(); this.specNNPS = ta(); this.specDQE = ta();
    left.append(section("4. Factory spec (optional)", hint("One 'lp/mm, value' pair per line, for each:"),
      el("div", { class: "sublbl" }, "MTF (%):"), this.specMTF, el("div", { class: "sublbl" }, "NNPS:"), this.specNNPS,
      el("div", { class: "sublbl" }, "DQE (%):"), this.specDQE));
    this.prog = el("progress", { max: "100", value: "0" }); left.append(this.prog);

    const pane = makePane(right);
    this.results = new ResultPane(pane.top, "Results", "image_quality.csv");
    this.plotHost = pane.bot;
    this.figHost = el("div", { class: "fill" }); this.ivHost = el("div", { class: "fill", style: "display:none" });
    this.plotHost.append(this.figHost, this.ivHost);
    this.fig = new Figure(this.figHost);
    this.viewer = new ImageViewer(this.ivHost, app.engine);
    this.viewer.onROI = (roi) => this.onROI(roi);

    this.edge = null; this.flat = null; this.roiH = null; this.roiV = null; this.resH = null; this.resV = null;
    this.nnps = null; this.dqe = null;
    app.engine.onRestart0 = () => {
      this.edge = null; this.flat = null; this.edgeLbl.textContent = "(none loaded)"; this.flatLbl.textContent = "(none loaded)";
      this.results.clear(); this.results.head("Error"); this.results.warn("The analysis engine hit an internal error and was restarted. Please load the images again.");
    };
  }
  setTables(t) { this.quality.setValues(t.snr_in2.map((x) => x[0]).sort((a, b) => t.rqa.find((r) => r[0] === a)[1] - t.rqa.find((r) => r[0] === b)[1]), "RQA5"); this.tables = t; }
  showFig() { this.ivHost.style.display = "none"; this.figHost.style.display = ""; this.fig.draw(); }
  showViewer() { this.figHost.style.display = "none"; this.ivHost.style.display = ""; }
  busy(on) { document.body.classList.toggle("busy", on); }

  async loadFile(file, slot, lbl) {
    lbl.textContent = "Loading...";
    this.busy(true);
    try {
      const buf = await file.arrayBuffer();
      const meta = await this.app.engine.call({ cmd: "dicom_load", slot }, buf);
      meta.name = file.name; lbl.textContent = file.name + "  (" + meta.w + " x " + meta.h + ", " + meta.px_mm + " mm)";
      return meta;
    } catch (e) {
      lbl.textContent = "(none loaded)";
      this.results.clear(); this.results.head("Load failed"); this.results.warn(e.message); return null;
    } finally { this.busy(false); }
  }
  dicomInfo(meta, title) {
    const r = this.results; r.clear(); r.head(title);
    r.row("File", meta.name);
    for (const [k, v] of meta.info) r.row(k, v);
    if (meta.sign === -1) r.note("Pixel Intensity Relationship Sign is -1, so the image is displayed inverted (as app.py does).");
  }
  async loadEdge(file) {
    const meta = await this.loadFile(file, "edge", this.edgeLbl);
    if (!meta) return;
    this.edge = meta;
    this.roiH = null; this.roiV = null;
    this.roiHLbl.textContent = "Horizontal ROI: not set"; this.roiVLbl.textContent = "Vertical ROI: not set";
    this.mtfStatus.textContent = this.mtfStatusText();
    this.dicomInfo(meta, "Edge image loaded");
    this.showEdge();
  }
  async loadFlat(file) {
    const meta = await this.loadFile(file, "flat", this.flatLbl);
    if (!meta) return;
    this.flat = meta; this.nnps = null; this.dqe = null; this.nnpsStatus.textContent = "NNPS: not done yet";
    this.dicomInfo(meta, "Flat-field image loaded");
  }
  showEdge() {
    if (!this.edge) { this.results.clear(); this.results.head("Input error"); this.results.warn("Load the edge image first."); this.results.top(); return; }
    this.showViewer();
    this.updateOverlays();
    this.viewer.show("edge", this.edge, "Drag a rectangle around a straight section of edge (avoid corners)");
  }
  updateOverlays() {
    this.viewer.rois = [{ roi: this.roiH, color: ACCENT, label: "horizontal ROI" }, { roi: this.roiV, color: "#c1666b", label: "vertical ROI" }];
  }
  onROI(roi) {
    const t = "rows " + roi[0] + "-" + roi[1] + ", cols " + roi[2] + "-" + roi[3];
    if (this.roiMode.get() === "h") { this.roiH = roi; this.roiHLbl.textContent = "Horizontal ROI: " + t; }
    else { this.roiV = roi; this.roiVLbl.textContent = "Vertical ROI: " + t; }
    this.updateOverlays(); this.viewer.draw();
  }
  mtfStatusText() { return "MTF: horizontal " + (this.resH ? "done" : "not done") + ", vertical " + (this.resV ? "done" : "not done"); }
  specPts(ta) {
    const out = [];
    for (let line of ta.value.split(/\r?\n/)) {
      line = line.trim(); if (!line) continue;
      const parts = line.replace(/\t/g, ",").split(",").map((p) => p.trim()).filter((p) => p);
      if (parts.length >= 2) { const a = Number(parts[0]), b = Number(parts[1]); if (isFinite(a) && isFinite(b) && parts[0] !== "" && parts[1] !== "") out.push([a, b]); }
    }
    return out;
  }
  async calcMTF(dir) {
    const r = this.results; r.clear();
    const H = dir === "h", roi = H ? this.roiH : this.roiV;
    if (!this.edge) { r.head("Input error"); r.warn("Load the edge image first."); return; }
    if (!roi) { r.head("Input error"); r.warn(H ? "Select a horizontal-MTF ROI (drag a horizontal edge) on the edge image first." : "Select a vertical-MTF ROI (drag a vertical edge) on the edge image first."); return; }
    this.busy(true);
    try {
      const res = await this.app.engine.call({ cmd: "mtf", slot: "edge", roi, horizontal: H, oversample: 8 });
      res.freq = nums(res.freq); res.mtf = nums(res.mtf);
      if (H) this.resH = res; else this.resV = res;
      this.mtfStatus.textContent = this.mtfStatusText();
      r.head(H ? "Horizontal MTF" : "Vertical MTF");
      r.row("Edge angle", fmt(res.angle_deg), "deg");
      if (!(Math.abs(res.angle_deg) >= 0.3 && Math.abs(res.angle_deg) <= 5)) r.warn("An edge angle outside about 0.3-5 degrees is outside the range this tilt-sampling technique is designed for - check the ROI straddles a single straight edge segment, away from a corner.");
      for (const f0 of [0.5, 1.0, 2.0, 3.0, 4.0]) r.row(f0.toFixed(1) + " lp/mm", interp1(f0, res.freq, res.mtf).toFixed(3));
      r.setDF(res.freq.map((f, i) => ({ freq_lp_per_mm: f, mtf: res.mtf[i] })), ["freq_lp_per_mm", "mtf"]);
      r.top();
      this.drawMTF(res, H ? "horizontal" : "vertical");
    } catch (e) { r.head("Error"); r.warn(e.message); } finally { this.busy(false); }
  }
  async calcNNPS() {
    const r = this.results; r.clear();
    if (!this.flat) { r.head("Input error"); r.warn("Load the flat-field image first."); return; }
    this.busy(true); this.prog.value = 0;
    try {
      const d = await this.app.engine.callP({ cmd: "nnps", slot: "flat" }, (f) => { this.prog.value = 100 * f; });
      this.prog.value = 100;
      const f = nums(d.freq);
      this.nnps = { fH: f, profH: nums(d.profH), fV: f, profV: nums(d.profV), mean_sig: d.mean_sig, n_used: d.n_used };
      this.nnpsStatus.textContent = "NNPS: done (" + d.n_used + " ROIs)";
      r.head("NNPS"); r.row("ROIs used", String(d.n_used)); r.row("Mean signal", fmt(d.mean_sig), "DN");
      const pos = (fa, pa) => { const x = [], y = []; fa.forEach((v, i) => { if (v >= 0) { x.push(v); y.push(pa[i]); } }); return [x, y]; };
      const [hx, hy] = pos(this.nnps.fH, this.nnps.profH), [vx, vy] = pos(this.nnps.fV, this.nnps.profV);
      for (const f0 of [0.5, 1.0, 2.0, 3.0]) r.row(f0.toFixed(1) + " lp/mm", "h=" + pyExp(interp1(f0, hx, hy), 3) + "  v=" + pyExp(interp1(f0, vx, vy), 3));
      r.setDF(hx.map((v, i) => ({ freq_lp_per_mm: v, nnps_h: hy[i], nnps_v: vy[i] })), ["freq_lp_per_mm", "nnps_h", "nnps_v"]);
      r.top();
      this.drawNNPS(hx, hy, vx, vy);
    } catch (e) { r.head("Error"); r.warn(e.message); } finally { this.busy(false); }
  }
  async calcDQE() {
    const r = this.results; r.clear();
    let K, q, qn;
    try {
      if (!this.resH || !this.resV) throw new UserError("Calculate both horizontal and vertical MTF first.");
      if (!this.nnps) throw new UserError("Calculate NNPS first.");
      K = this.kerma.get(0.0);
      if (K <= 0) throw new UserError("Enter the measured air kerma at the detector for the flat-field exposure.");
      qn = this.quality.get(); const e = this.tables.snr_in2.find((x) => x[0] === qn); q = e ? e[1] : null;
      if (q === null) throw new UserError(qn + " has no mandatory SNR_in^2 value in Table 3 (only RQA3/5/7/9 do).");
    } catch (e) { r.head("Input error"); r.warn(e.message); return; }
    this.busy(true);
    try {
      const n = this.nnps;
      const d = await this.app.engine.call({ cmd: "dqe", q, K, fmh: this.resH.freq, mtfh: this.resH.mtf, fH: n.fH, profH: n.profH,
        fmv: this.resV.freq, mtfv: this.resV.mtf, fV: n.fV, profV: n.profV });
      this.dqe = { fDh: nums(d.fDh), dqe_h: nums(d.dqe_h), fDv: nums(d.fDv), dqe_v: nums(d.dqe_v) };
      const dd = this.dqe;
      r.head("DQE"); r.row("Air kerma", fmt(K), "uGy"); r.row("q (" + qn + ", Table 3)", fmt(q, 1), "mm^-2.uGy^-1");
      for (const f0 of [0.5, 1.0, 2.0, 3.0]) r.row(f0.toFixed(1) + " lp/mm", "h=" + interp1(f0, dd.fDh, dd.dqe_h).toFixed(3) + "  v=" + interp1(f0, dd.fDv, dd.dqe_v).toFixed(3));
      r.note("Single edge per direction (not IEC's optional two-edge average), and detector response assumed linear through zero from this one flat-field exposure - see the Notes tab.");
      r.setDF(dd.fDh.map((f, i) => ({ freq_lp_per_mm: f, dqe_h: dd.dqe_h[i], dqe_v: dd.dqe_v[i] })), ["freq_lp_per_mm", "dqe_h", "dqe_v"]);
      r.top();
      this.drawDQE();
    } catch (e) { r.head("Error"); r.warn(e.message); } finally { this.busy(false); }
  }
  drawMTF(res, direction) {
    this.showFig();
    const [ax] = this.fig.reset(1);
    ax.plot(res.freq, res.mtf, { color: ACCENT, label: direction });
    const sp = this.specPts(this.specMTF);
    if (sp.length) ax.plot(sp.map((p) => p[0]), sp.map((p) => p[1] / 100), { color: "black", lw: 1, dash: [5, 4], marker: "o", ms: 4, label: "factory spec" });
    ax.xlim = [0, Math.min(amax(res.freq), 6)]; ax.ylim = [0, 1.05];
    ax.xlabel = "lp/mm"; ax.ylabel = "MTF"; ax.legend({ size: 10 });
    ax.title = direction[0].toUpperCase() + direction.slice(1) + " MTF (edge angle " + res.angle_deg.toFixed(2) + " deg)";
    this.fig.draw();
  }
  drawNNPS(hx, hy, vx, vy) {
    this.showFig();
    const [ax] = this.fig.reset(1); ax.ylog = true;
    ax.plot(hx, hy, { color: ACCENT, label: "horizontal" }); ax.plot(vx, vy, { color: "#c1666b", label: "vertical" });
    const sp = this.specPts(this.specNNPS);
    if (sp.length) ax.plot(sp.map((p) => p[0]), sp.map((p) => p[1]), { color: "black", lw: 1, dash: [5, 4], marker: "o", ms: 4, label: "factory spec" });
    ax.xlim = [0, Math.min(amax(this.nnps.fH), 6)];
    ax.xlabel = "lp/mm"; ax.ylabel = "NNPS (mm²)"; ax.title = "NNPS"; ax.legend({ size: 10 });
    this.fig.draw();
  }
  drawDQE() {
    this.showFig();
    const d = this.dqe; const [ax] = this.fig.reset(1);
    ax.plot(d.fDh, d.dqe_h, { color: ACCENT, label: "horizontal" }); ax.plot(d.fDv, d.dqe_v, { color: "#c1666b", label: "vertical" });
    const sp = this.specPts(this.specDQE);
    if (sp.length) ax.plot(sp.map((p) => p[0]), sp.map((p) => p[1] / 100), { color: "black", lw: 1, dash: [5, 4], marker: "o", ms: 4, label: "factory spec" });
    ax.xlim = [0, Math.min(amax(d.fDh), 6)]; ax.ylim = [0, Math.max(0.05, 1.05 * nanmax(d.dqe_h))];
    ax.xlabel = "lp/mm"; ax.ylabel = "DQE"; ax.title = "DQE"; ax.legend({ size: 10 });
    this.fig.draw();
  }
  copyToClipboard(grid, values, label, nRows, fmaxCap) {
    const r = this.results;
    const gx = [], gy = []; grid.forEach((g, i) => { if (isFinite(values[i])) { gx.push(g); gy.push(values[i]); } });
    if (gx.length < 2) { r.clear(); r.head("Copy"); r.warn(label + " has no usable data yet."); return; }
    let fmax = amax(gx); if (fmaxCap != null) fmax = Math.min(fmax, fmaxCap);
    const og = linspace(0, fmax, nRows);
    copyText(og.map((f) => pyG(f) + "\t" + pyG(interp1(f, gx, gy))).join("\n"));
    r.clear(); r.head("Copy"); r.row("Copied", label + " - " + nRows + " rows, 0 to " + fmax.toFixed(3) + " lp/mm"); r.top();
  }
  copyCurve(dir) {
    const res = dir === "h" ? this.resH : this.resV, label = dir === "h" ? "Horizontal MTF" : "Vertical MTF";
    if (!res) { this.results.clear(); this.results.head("Copy"); this.results.warn(label + " hasn't been calculated yet."); return; }
    const nyq = res.px_mm ? 0.5 / res.px_mm : null;
    this.copyToClipboard(res.freq, res.mtf, label, 103, nyq);
  }
  copyNNPS(dir) {
    if (!this.nnps) { this.results.clear(); this.results.head("Copy"); this.results.warn("NNPS hasn't been calculated yet."); return; }
    const g = dir === "h" ? this.nnps.fH : this.nnps.fV, v = dir === "h" ? this.nnps.profH : this.nnps.profV;
    const x = [], y = []; g.forEach((f, i) => { if (f >= 0) { x.push(f); y.push(v[i]); } });
    this.copyToClipboard(x, y, "NNPS (" + dir + ")", 63);
  }
  copyDQE() {
    if (!this.dqe) { this.results.clear(); this.results.head("Copy"); this.results.warn("DQE hasn't been calculated yet."); return; }
    const d = this.dqe;
    const fin = (f, v) => f.filter((_, i) => isFinite(v[i]));
    const fmax = Math.min(amax(fin(d.fDh, d.dqe_h)), amax(fin(d.fDv, d.dqe_v)));
    const n = 63, grid = linspace(0, fmax, n);
    const hi = grid.map((f) => interp1(f, d.fDh, d.dqe_h)), vi = grid.map((f) => interp1(f, d.fDv, d.dqe_v));
    const avg = hi.map((h, i) => 0.5 * (h + vi[i]));
    const ratio = hi.map((h, i) => (vi[i] > 0 ? h / vi[i] : NaN));
    const bad = grid.filter((_, i) => isFinite(ratio[i]) && (ratio[i] < 0.9 || ratio[i] > 1.1));
    copyText(grid.map((f, i) => pyG(f) + "\t" + pyG(avg[i])).join("\n"));
    const r = this.results; r.clear(); r.head("Copy DQE (h/v average)");
    r.row("Copied", n + " rows, 0 to " + fmax.toFixed(3) + " lp/mm");
    if (bad.length) r.warn("DQE(h)/DQE(v) falls outside IEC 62220-1's 0.9-1.1 averaging criterion at " + bad.length + " of " + n + " frequencies (from " + Math.min(...bad).toFixed(2) + " lp/mm up) - the two directions diverge there, so the averaged value copied at those frequencies is less meaningful than at the lower frequencies where they agree. Copied anyway since you asked for one curve; worth checking the DQE plot's h/v split before relying on the tail of this table.");
    else r.note("DQE(h) and DQE(v) agree within IEC's 0.9-1.1 averaging criterion across the whole copied range.");
    r.top();
  }
}
