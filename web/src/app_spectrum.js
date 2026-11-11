/* Spectrum tab - port of SpectrumTab */
"use strict";

const MATERIAL_LABEL = {
  lead: "Lead", concrete: "Concrete", gypsum: "Gypsum", steel: "Steel", glass: "Glass", borosilicate: "Borosilicate",
  leadglass: "Lead glass", aluminium: "Aluminium", copper: "Copper", water: "Water", pmma: "PMMA", tungsten: "Tungsten",
  beryllium: "Beryllium", air: "Air", tissue: "Tissue",
};
const ELEMENT_LABEL = {
  hydrogen: "Hydrogen", helium: "Helium", lithium: "Lithium", boron: "Boron", carbon: "Carbon", nitrogen: "Nitrogen",
  oxygen: "Oxygen", fluorine: "Fluorine", neon: "Neon", sodium: "Sodium", magnesium: "Magnesium", silicon: "Silicon",
  phosphorus: "Phosphorus", sulfur: "Sulfur", chlorine: "Chlorine", argon: "Argon", potassium: "Potassium", calcium: "Calcium",
  scandium: "Scandium", titanium: "Titanium", vanadium: "Vanadium", chromium: "Chromium", manganese: "Manganese", iron: "Iron",
  cobalt: "Cobalt", nickel: "Nickel", zinc: "Zinc", gallium: "Gallium", germanium: "Germanium", arsenic: "Arsenic",
  selenium: "Selenium", bromine: "Bromine", krypton: "Krypton", rubidium: "Rubidium", strontium: "Strontium", yttrium: "Yttrium",
  zirconium: "Zirconium", niobium: "Niobium", molybdenum: "Molybdenum", technetium: "Technetium", ruthenium: "Ruthenium",
  rhodium: "Rhodium", palladium: "Palladium", silver: "Silver", cadmium: "Cadmium", indium: "Indium", tin: "Tin",
  antimony: "Antimony", tellurium: "Tellurium", iodine: "Iodine", xenon: "Xenon", cesium: "Caesium", barium: "Barium",
  lanthanum: "Lanthanum", cerium: "Cerium", praseodymium: "Praseodymium", neodymium: "Neodymium", promethium: "Promethium",
  samarium: "Samarium", europium: "Europium", gadolinium: "Gadolinium", terbium: "Terbium", dysprosium: "Dysprosium",
  holmium: "Holmium", erbium: "Erbium", thulium: "Thulium", ytterbium: "Ytterbium", lutetium: "Lutetium", hafnium: "Hafnium",
  tantalum: "Tantalum", rhenium: "Rhenium", osmium: "Osmium", iridium: "Iridium", platinum: "Platinum", gold: "Gold",
  mercury: "Mercury", thallium: "Thallium", bismuth: "Bismuth", polonium: "Polonium", astatine: "Astatine", radon: "Radon",
  francium: "Francium", radium: "Radium", actinium: "Actinium", thorium: "Thorium", protactinium: "Protactinium", uranium: "Uranium",
};
Object.assign(MATERIAL_LABEL, ELEMENT_LABEL);
const LBL = (k) => MATERIAL_LABEL[k] || k;
const _ELEMENTS_BY_LABEL = Object.keys(ELEMENT_LABEL).sort((a, b) => (ELEMENT_LABEL[a] < ELEMENT_LABEL[b] ? -1 : ELEMENT_LABEL[a] > ELEMENT_LABEL[b] ? 1 : 0));
const _BASE_MATERIALS = ["aluminium", "copper", "pmma", "beryllium", "borosilicate", "glass", "water", "tissue", "air", "lead",
  "leadglass", "concrete", "gypsum", "steel", "tungsten"];
let FILTER_MATERIALS = _BASE_MATERIALS.concat(_ELEMENTS_BY_LABEL);
let BARRIER_MATERIALS = _BASE_MATERIALS.concat(_ELEMENTS_BY_LABEL);
const matOpts = (keys) => keys.map((k) => [k, LBL(k)]);
const filtStr = (fl) => fl.map((f) => LBL(f[0]) + " " + fmt(f[1]) + " mm").join(", ");

const CALIBRATION = { scale: null, kvp: null, filters: null, measured: null, distance: null, engine: null, spek_ratio: null };
function clearCalibration() { for (const k in CALIBRATION) CALIBRATION[k] = null; }
function outputBasis() {
  if (CALIBRATION.scale !== null) return "calibrated: " + pyRepr(CALIBRATION.measured) + " uGy/mAs at " + pyRepr(CALIBRATION.distance) + " m, " + pyRepr(CALIBRATION.kvp) + " kVp";
  return "SpekPy absolute scale";
}
function pyRepr(v) { return Number.isInteger(v) ? v.toFixed(1) : String(v); }

class SpectrumTab {
  constructor(root, app) {
    this.app = app; this.shapeSet = false; this.spectrum = null; this.params = null;
    root.append(el("p", { class: "intro" }, "Generate a tube spectrum and filter it. Filtering is Beer-Lambert per energy bin, which gives narrow-beam attenuation."));
    const body = el("div", { class: "body" }); root.append(body);
    const left = el("div", { class: "left" }); const right = el("div", { class: "right" }); body.append(left, right);

    this.kvp = new Entry("Tube potential", 125, "kVp");
    this.angle = new Entry("Anode angle", 12, "deg");
    this.engine = new Combo("Spectrum model", [["spekpy", "SpekPy"], ["builtin", "Built-in (approx)"]], "spekpy");
    left.append(section("Tube", this.kvp.root, this.angle.root, this.engine.root));

    this.filters = [];
    const s2 = section("Total filtration (inherent + added)", hint("2.5 mm Al is the regulatory minimum total above 70 kV; typical general radiography is 2.5-3.5 mm. Tube labels usually quote total, not inherent - inherent alone is 0.5-1.5 mm."));
    [["aluminium", 2.5], ["copper", 0.0], ["water", 0.0]].forEach(([m, mm], i) => {
      const c = new Combo("Filter " + (i + 1), matOpts(FILTER_MATERIALS), m, { compact: true });
      const e = new Entry("", mm, "mm", { compact: true });
      s2.append(el("div", { class: "row2" }, c.root, e.root));
      this.filters.push([c, e]);
    });
    left.append(s2);

    this.hvlMat = new Combo("Measure HVL in", matOpts(FILTER_MATERIALS), "aluminium");
    this.hvlMeas = new Entry("Measured HVL", "", "mm");
    this.measKvp = new Entry("measured at", "", "kVp", { title: "Blank uses the tube potential above." });
    this.meas = new Entry("Measured output", "", "uGy/mAs");
    this.measD = new Entry("at distance", 1.0, "m");
    this.calLbl = el("div", { class: "hint" }, "not calibrated");
    left.append(section("Calibration", hint("Measured on the tube. Enter either or both, then press Calibrate once."),
      this.hvlMat.root, this.hvlMeas.root, this.measKvp.root, this.meas.root, this.measD.root,
      el("div", { class: "btnrow" }, button("Calibrate", () => this.calibrate()), button("Clear", () => this.uncalibrate())),
      this.calLbl,
      hint("HVL sets the shape by solving filtration. Output sets the scale. Both are properties of the tube, so once solved at the kVp measured at, the calibration applies at every other kVp.")));

    left.append(button("Generate", () => this.run(), "accent big"));
    left.append(button("Sweep kVp vs Archer 1994", () => this.sweep(), "wide"));

    this.quality = new Combo("Quality", [], "RQA5");
    left.append(section("RQA beam qualities", this.quality.root,
      el("div", { class: "btnrow" }, button("Identify beam", () => this.rqaIdentify()), button("Solve for target", () => this.rqaSolve())),
      hint("IEC 61267 Ed.2 (2005). Qualities are defined by HVL, not by the nominal added filtration.")));

    const pane = makePane(right);
    this.results = new ResultPane(pane.top, "Results", "spectrum.csv");
    this.fig = new Figure(pane.bot);
  }
  setTables(t) {
    this.tables = t;
    this.quality.setValues(t.rqa.slice().sort((a, b) => a[1] - b[1]).map((r) => r[0]), "RQA5");
  }
  setMaterials() {
    for (const [c] of this.filters) { const v = c.get(); c.setValues(matOpts(FILTER_MATERIALS), FILTER_MATERIALS.includes(v) ? v : FILTER_MATERIALS[0]); }
    const v = this.hvlMat.get(); this.hvlMat.setValues(matOpts(FILTER_MATERIALS), FILTER_MATERIALS.includes(v) ? v : "aluminium");
  }
  currentFilters() {
    const out = [];
    for (const [c, e] of this.filters) { let mm = 0; try { mm = e.get(0.0); } catch (x) { mm = 0; } if (mm > 0) out.push([c.get(), mm]); }
    return out;
  }
  calLabel() {
    if (CALIBRATION.scale === null && !this.shapeSet) { this.calLbl.textContent = "not calibrated"; this.calLbl.style.color = ""; return; }
    const bits = []; if (this.shapeSet) bits.push("shape"); if (CALIBRATION.scale !== null) bits.push("scale");
    let k = 0; try { k = this.kvp.get(0); } catch (e) { k = 0; }
    this.calLbl.textContent = bits.join(" and ") + " set at " + fmt(k) + " kVp"; this.calLbl.style.color = OK;
  }
  useSpek() { return this.engine.get() === "spekpy"; }
  async busy(p) { document.body.classList.add("busy"); try { return await p; } finally { document.body.classList.remove("busy"); } }

  async calibrate() {
    const r = this.results;
    try {
      const kvpShow = this.kvp.get();
      const kvp = this.measKvp.get(0.0) || kvpShow;
      const ang = this.angle.get(12.0), hvlT = this.hvlMeas.get(0.0), outT = this.meas.get(0.0), dist = this.measD.get(1.0);
      const engine = this.engine.get();
      if (hvlT <= 0 && outT <= 0) throw new UserError("Enter a measured HVL, a measured output, or both.");
      const report = [];
      let shapeChanged = false;
      if (hvlT > 0) {
        const s = await this.busy(this.app.engine.call({ cmd: "solve_hvl", engine, kvp, hvl: hvlT, angle: ang }));
        const mm = s.mm === null ? null : num(s.mm), ach = num(s.achieved);
        if (mm === null) {
          report.push(["warn", "Cannot reach an HVL of " + fmt(hvlT) + " mm. With no filtration at all the model already gives " + fmt(ach) + " mm, so it is making a harder beam than your tube."]);
          report.push(["note", "On the built-in model, select SpekPy. On SpekPy, check the anode angle."]);
          this.shapeSet = false;
        } else {
          this.filters[0][0].set("aluminium"); this.filters[0][1].set(Math.round(mm * 1000) / 1000);
          for (const [, e] of this.filters.slice(1)) e.set(0.0);
          this.shapeSet = true; shapeChanged = true;
          report.push(["row", ["Measured HVL", fmt(hvlT), "mm"]]);
          report.push(["big", ["Total filtration", fmt(mm), "mm Al equivalent"]]);
          report.push(["row", ["Reproduced HVL", fmt(ach), "mm"]]);
          report.push(["row", ["Deviation", signed(100 * (ach - hvlT) / hvlT, 3) + " %", ""]]);
        }
      } else report.push(["note", "No HVL given - filtration left as entered above, so the spectral shape is only as good as that figure."]);
      if (shapeChanged && outT <= 0 && CALIBRATION.scale !== null) {
        clearCalibration();
        report.push(["warn", "The filtration changed, so the previous scale factor referred to a different beam and has been cleared. Re-enter your measured output to restore it."]);
      }
      const measFilters = this.currentFilters();
      await this.run();
      if (outT > 0) {
        const c = await this.busy(this.app.engine.call({ cmd: "set_calibration", engine, kvp, angle: ang, filters: measFilters, measured: outT, distance: dist }));
        Object.assign(CALIBRATION, { scale: num(c.scale), kvp: c.kvp, filters: measFilters, measured: outT, distance: dist, engine: c.model,
          spek_ratio: c.spek_ratio === null ? null : num(c.spek_ratio) });
        await this.run();
        report.push(["row", ["Measured output", fmt(outT) + " uGy/mAs at " + fmt(dist) + " m", ""]]);
        report.push(["row", ["Scale factor", fmt(CALIBRATION.scale, 6), ""]]);
        if (CALIBRATION.spek_ratio) {
          report.push(["row", ["Correction to SpekPy", "x " + fmt(CALIBRATION.spek_ratio, 4), ""]]);
          report.push(["note", "A factor near 1 means SpekPy already agreed with your measurement."]);
        }
      } else report.push(["note", "No output given - beam quality is calibrated but dose cannot be quoted in uGy."]);
      this.calLabel();
      r.head("Calibration from measurements at " + fmt(kvp) + " kVp");
      if (Math.abs(kvp - kvpShow) > 1e-9) r.note("Tube potential on display is " + fmt(kvpShow) + " kVp. Filtration and scale were derived at " + fmt(kvp) + " kVp and carried across, which is what makes them tube properties rather than settings.");
      for (const [kind, p] of report) {
        if (kind === "row") r.row(p[0], p[1], p[2]); else if (kind === "big") r.big(p[0], p[1], p[2]);
        else if (kind === "note") r.note(p); else r.warn(p);
      }
      r.head("What was set");
      r.row("Shape", this.shapeSet ? "from measured HVL" : "from the filtration entered above");
      r.row("Scale", CALIBRATION.scale !== null ? "from measured output" : "not set");
      r.head("How to check it");
      r.note("Change the kVp and press Generate. Both the HVL and the output should MOVE. If the output stays put, the scale is not being applied.");
      r.note("Add filtration and press Generate. HVL should rise and output should fall.");
      r.note("Press Sweep to see the whole kVp range on your scale.");
      r.top();
    } catch (e) { r.clear(); r.head("Calibration error"); r.warn(e.message); }
  }
  async uncalibrate() { clearCalibration(); this.shapeSet = false; this.calLabel(); await this.run(); }

  async currentSpectrum() {
    if (!this.spectrum) await this.run();
    if (!this.spectrum) throw new UserError("Generate a spectrum first.");
    return this.params;
  }
  calScale() { return CALIBRATION.scale; }

  async rqaIdentify() {
    const r = this.results;
    try {
      const p = await this.currentSpectrum();
      r.clear();
      const d = await this.busy(this.app.engine.call(Object.assign({ cmd: "rqa_identify", cal_scale: this.calScale() }, p)));
      const ed = this.tables.rqa_edition, rows = d.rows, hvl = num(d.hvl);
      r.head("Your beam against " + ed);
      r.row("Tube potential", fmt(d.kvp), "kVp");
      r.row("Filtration", filtStr(d.filters) || "none");
      r.big("Measured HVL", fmt(hvl), "mm Al");
      r.row("Mean energy", fmt(d.mean), "keV");
      const best = rows[0];
      r.head("Closest match");
      if (Math.abs(best.deviation_pc) < 5) {
        r.big("Nearest quality", best.quality); r.row("Target HVL", fmt(best.target_HVL_mm), "mm Al"); r.row("Deviation", signed(best.deviation_pc, 2) + " %");
      } else {
        r.row("Nearest quality", best.quality + " (" + best.deviation_pc.toFixed(1) + " % off)");
        r.warn("Your beam is not close to any RQA quality. RQA beams are heavily hardened - RQA5 needs about 21 mm of added aluminium on top of a normal tube.");
      }
      if (best.confidence !== "verified") r.warn("The " + best.quality + " entry could not be confirmed from an accessible published source. Check it against the standard before relying on it.");
      if (best.kV !== d.kvp) r.note("Note the kV differs: " + best.quality + " is specified at " + best.kV + " kV and your spectrum is at " + fmt(d.kvp) + " kVp. HVL alone does not make a beam an RQA quality - the tube potential must match too.");
      const cols = ["quality", "kV", "nominal_added_Al_mm", "target_HVL_mm", "your_HVL_mm", "deviation_pc", "confidence"];
      r.head("All qualities, sorted by closeness");
      r.frame(rows, cols);
      r.head("Photon fluence per air kerma (IEC 62220-1 Table 3)");
      if (d.q === null) {
        r.note("Not available - this spectrum has no absolute scale. Press Calibrate with a measured output, or use SpekPy's own normalisation, to get a real photons/mm2/uGy figure. HVL alone does not fix this: two beams can share an RQA's target HVL while differing in Q/K_air, since the standard constrains the beam by HVL and kV, not by fluence.");
      } else {
        r.big("Your Q/K_air (SNR_in^2)", fmt(num(d.q), 1), "mm^-2.uGy^-1");
        const fb = d.frows[0];
        r.row("Closest mandatory quality", fb.quality + " (" + signed(fb.deviation_pc, 2) + " %)");
        r.note("Only RQA3/5/7/9 carry a mandatory value in Table 3 - this is not compared against RQA2/4/6/8/10.");
        const df = r.df; r.frame(d.frows, ["quality", "target_SNRin2", "your_SNRin2", "deviation_pc"]); r.df = df;
      }
      r.top();
      const srt = rows.slice().sort((a, b) => a.kV - b.kV);
      const [ax] = this.fig.reset(1);
      ax.plot(srt.map((x) => x.kV), srt.map((x) => x.target_HVL_mm), { color: WARN, lw: 1.6, marker: "s", ms: 6, label: ed + " targets" });
      for (const x of srt) ax.text(x.kV, x.target_HVL_mm, x.quality, { color: WARN, size: 10, dx: 3, dy: -10 });
      ax.plot([d.kvp], [hvl], { color: ACCENT, lw: 0, marker: "o", ms: 11, label: "your beam" });
      ax.xlabel = "tube potential (kV)"; ax.ylabel = "HVL (mm Al)"; ax.title = "Beam quality against the RQA series"; ax.legend();
      this.fig.draw();
    } catch (e) { r.head("Input error"); r.warn(e.message); }
  }

  async rqaSolve() {
    const r = this.results;
    try {
      const p = await this.currentSpectrum();
      r.clear();
      const q = this.quality.get();
      const d = await this.busy(this.app.engine.call(Object.assign({ cmd: "rqa_solve", quality: q }, p)));
      const ed = this.tables.rqa_edition, tgt = d.target, mm = d.mm === null ? null : num(d.mm), mmc = d.mm_cu === null ? null : num(d.mm_cu);
      r.head("Reaching " + q + " from your beam");
      r.row("Standard", ed); r.row("Specified at", d.kv + " kV"); r.row("Target HVL", fmt(tgt), "mm Al");
      r.row("Nominal added Al", fmt(d.nominal), "mm"); r.row("Your starting HVL", fmt(num(d.start_hvl)), "mm Al");
      if (mm === null) r.warn("Pure Al: cannot reach the target. Your beam is either already harder than the target, or 80 mm of aluminium is not enough. Check the tube potential matches.");
      else {
        r.big("Add (pure Al)", fmt(mm), "mm Al"); r.row("Achieved HVL", fmt(num(d.achieved)), "mm Al");
        r.row("Deviation from target", signed(100 * (num(d.achieved) - tgt) / tgt, 3) + " %");
        r.row("Versus nominal", signed(mm - d.nominal, 2) + " mm");
        r.note("The nominal figure assumes the standard's reference tube. Yours differs, which is exactly why the standard specifies HVL rather than filtration.");
      }
      r.head("Alternative: 0.5 mm Cu + Al");
      if (mmc === null) r.warn("0.5 mm Cu alone already exceeds the target HVL - no amount of added Al can bring it back down. This Cu thickness cannot reach " + q + " from your beam.");
      else {
        r.big("Add (on top of 0.5 mm Cu)", fmt(mmc), "mm Al"); r.row("Achieved HVL", fmt(num(d.achieved_cu)), "mm Al");
        r.row("Deviation from target", signed(100 * (num(d.achieved_cu) - tgt) / tgt, 3) + " %");
        r.note("Same target HVL, different filter composition. A beam that matches by HVL alone is not guaranteed to match in photon fluence per unit air kerma - see the fluence comparison in Identify beam. This is offered as a genuinely different route to the same HVL, not a more accurate version of the pure-Al answer above.");
      }
      if (d.sp.kvp !== d.kv) r.warn("Your spectrum is at " + fmt(d.sp.kvp) + " kVp but " + q + " is specified at " + d.kv + " kV. Matching HVL alone does not reproduce the quality.");
      if (d.conf !== "verified") r.warn("This entry is unverified - confirm against the standard before use.");
      r.top();
      const [ax] = this.fig.reset(1);
      const nrm = (a) => { const m = Math.max(amax(a), 1e-30); return a.map((v) => v / m); };
      ax.plot(d.sp.E, nrm(nums(d.sp.N)), { color: "#999", lw: 1.4, label: "your beam" });
      if (d.hard) {
        ax.plot(d.hard.E, nrm(nums(d.hard.N)), { color: ACCENT, lw: 2, label: q + "  (+" + fmt(mm) + " mm Al)" });
        ax.vline(d.hard.mean, { color: OK, lw: 1, dash: [5, 4], label: "mean " + d.hard.mean.toFixed(1) + " keV" });
      }
      if (d.hard_cu) ax.plot(d.hard_cu.E, nrm(nums(d.hard_cu.N)), { color: "#c1666b", lw: 2, dash: [6, 4], label: q + "  (0.5 mm Cu + " + fmt(mmc) + " mm Al)" });
      ax.xlim = [0, d.sp.kvp * 1.02]; ax.ylim = [0, 1.15];
      ax.xlabel = "photon energy (keV)"; ax.ylabel = "relative fluence per keV"; ax.title = "Hardening to " + q; ax.legend();
      this.fig.draw();
    } catch (e) { r.head("Input error"); r.warn(e.message); }
  }

  drawSpectrum(sp, kvp, mean) {
    const [ax] = this.fig.reset(1);
    const N = nums(sp.N), peak = N.length ? amax(N) : 1;
    ax.fill(sp.E, N, { color: ACCENT, alpha: 0.22 });
    ax.plot(sp.E, N, { color: ACCENT, lw: 1.9 });
    ax.vline(mean, { color: OK, lw: 1.1, dash: [5, 4] });
    ax.text(mean, peak * 0.94, "mean " + mean.toFixed(1) + " keV", { color: OK, size: 11, dx: 5 });
    ax.xlim = [0, kvp]; ax.ylim = [0, peak > 0 ? peak * 1.15 : 1];
    ax.xlabel = "photon energy (keV)"; ax.ylabel = "fluence per keV";
    ax.title = fmt(kvp) + " kVp, " + (filtStr(sp.filters) || "no filtration");
    this.fig.draw();
  }

  async sweep() {
    const r = this.results; r.clear();
    try {
      const ang = this.angle.get(12.0), fl = this.currentFilters(), engine = this.engine.get();
      const d = await this.busy(this.app.engine.call({ cmd: "sweep", engine, angle: ang, filters: fl, cal_scale: this.calScale() }));
      const rows = d.rows.map((x) => ({ kVp: x.kVp, K_app_mGy_per_mAmin: x.K_app === null ? null : num(x.K_app), K_Archer_mGy_per_mAmin: x.K_archer,
        HVL_Al_mm: num(x.HVL), mean_keV: x.mean }));
      const haveAbs = rows.every((x) => x.K_app_mGy_per_mAmin !== null);
      r.head("Output against kVp, filtration held fixed");
      r.row("Filtration", filtStr(fl) || "none"); r.row("Anode angle", fmt(ang), "deg");
      const cols = ["kVp", "K_app_mGy_per_mAmin", "K_Archer_mGy_per_mAmin", "HVL_Al_mm", "mean_keV"];
      if (!haveAbs) r.warn("No absolute scale, so only HVL and mean energy are meaningful here. Select SpekPy, or press Calibrate with a measured output, to compare output.");
      else {
        for (const x of rows) x.ratio = x.K_app_mGy_per_mAmin / x.K_Archer_mGy_per_mAmin;
        cols.push("ratio");
        const k100 = rows.find((x) => x.kVp === 100).K_app_mGy_per_mAmin, k125 = rows.find((x) => x.kVp === 125).K_app_mGy_per_mAmin;
        const ak = (kv) => 1.222170186 - 0.056644073 * kv + 0.001227127 * kv * kv - 3.13579e-6 * kv * kv * kv;
        r.head("The test that actually discriminates");
        r.row("app K(125)/K(100)", fmt(k125 / k100, 4)); r.row("Archer tabulated 7.170/4.722", fmt(7.170 / 4.722, 4));
        r.row("Archer fitted ratio", fmt(ak(125) / ak(100), 4));
        r.note("This ratio is independent of absolute normalisation, so it is insensitive to anode angle and generator waveform. Matching it is evidence about FILTRATION. Matching a single output value is not - that is one free parameter fitted to one number.");
        const rat = rows.map((x) => x.ratio);
        r.row("Spread of app/Archer ratio", fmt(Math.min(...rat), 4) + " to " + fmt(Math.max(...rat), 4));
        r.note("A flat ratio near 1 across the whole range means the filtration is right. A ratio that trends with kVp means it is compensating for something else.");
        r.note("The reference is a fit to the measured data in Archer 1994. NCRP 147 carries the same coefficients with one sign wrong, so the paper is the source used here.");
      }
      r.head("Sweep");
      const disp = rows.map((x) => { const o = { kVp: x.kVp }; for (const c of cols.slice(1)) o[c] = x[c] === null ? "-" : fmt(x[c], 4); return o; });
      r.frame(rows, cols, disp);
      r.head("Archer 1994 tabulated anchor values");
      r.row("100 kVp", "4.722", "mGy per mA.min at 1 m"); r.row("125 kVp", "7.170", "mGy per mA.min at 1 m");
      r.note("These sit alongside the lead parameters a=2.500, b=15.28, g=0.7557 - the same pair your shielding calculator uses at 100 kVp. The curve above is the fit to this data, so the two are one source, not two.");
      r.top();
      const [a1, a2] = this.fig.reset(2);
      const kv = rows.map((x) => x.kVp);
      if (haveAbs) a1.plot(kv, rows.map((x) => x.K_app_mGy_per_mAmin), { color: ACCENT, lw: 1.8, marker: "o", ms: 4, label: "this app" });
      a1.plot(kv, rows.map((x) => x.K_Archer_mGy_per_mAmin), { color: WARN, lw: 1.5, dash: [6, 4], marker: "s", ms: 4, label: "Archer 1994 (fitted)" });
      a1.xlabel = "kVp"; a1.ylabel = "mGy per mA.min at 1 m"; a1.title = "Output"; a1.legend({ size: 10 });
      if (haveAbs) { a2.plot(kv, rows.map((x) => x.ratio), { color: ACCENT, lw: 1.8, marker: "o", ms: 4 }); a2.hline(1.0, { color: "#c33", dash: [5, 4] }); a2.ylabel = "app / Archer"; a2.title = "Ratio - flat is good"; }
      else { a2.plot(kv, rows.map((x) => x.HVL_Al_mm), { color: OK, lw: 1.8, marker: "o", ms: 4 }); a2.ylabel = "HVL in Al (mm)"; a2.title = "Beam quality"; }
      a2.xlabel = "kVp";
      this.fig.draw();
    } catch (e) { r.head("Input error"); r.warn(e.message); }
  }

  async run() {
    const r = this.results; r.clear();
    try {
      const kvpIn = this.kvp.get(), ang = this.angle.get(12.0), fl = this.currentFilters(), engine = this.engine.get();
      const useSpek = engine === "spekpy";
      const hm = this.hvlMat.get();
      const p = { kvp: kvpIn, angle: ang, engine, filters: fl, hvl_mat: hm };
      const d = await this.busy(this.app.engine.call(Object.assign({ cmd: "spectrum", cal_scale: this.calScale() }, p)));
      const sp = d.sp, unf = d.unf, kvp = sp.kvp;
      this.spectrum = sp; this.params = p; this.result = d;
      r.head("Spectrum at " + fmt(kvp) + " kVp");
      r.row("Model", useSpek ? "SpekPy (pre-computed table)" : sp.model);
      if (useSpek && Math.abs(kvp - kvpIn) > 1e-9) r.note("The SpekPy table holds SpekPy's exact output at whole kV, so " + fmt(kvpIn) + " kVp was run as " + fmt(kvp) + " kVp.");
      if (!useSpek) r.note("SpekPy is available - select it for better accuracy.");
      r.row("Filtration", filtStr(fl) || "none");
      r.head("Beam quality");
      r.big("Mean energy", fmt(d.mean), "keV");
      r.row("Effective energy", fmt(d.eff), "keV");
      r.row("HVL 1 in " + LBL(hm), fmt(num(d.hvl1)), "mm");
      r.row("HVL 2 (0.5 -> 0.25)", fmt(num(d.hvl2)), "mm");
      r.row("Homogeneity coefficient", fmt(d.homog));
      r.note("HVL is air-kerma defined, matching what an ion chamber or a RaySafe measures.");
      r.row("Photon fluence (rel.)", fmt(d.fluence)); r.row("Energy fluence (rel.)", fmt(d.efluence));
      if (d.unf_fluence > 0 && fl.length) {
        r.head("What the filtration removed");
        r.row("Photons", (100 * (1 - d.fluence / d.unf_fluence)).toFixed(1) + " %");
        if (d.k0 > 0) r.row("Air kerma", (100 * (1 - d.k1 / d.k0)).toFixed(1) + " %");
        r.row("Mean energy raised", fmt(d.unf_mean) + " -> " + fmt(d.mean) + " keV");
        r.note("The plot shows the beam as it is. An unfiltered spectrum peaks far to the left and far higher, so drawing both on one scale hides the real beam - and no tube emits an unfiltered spectrum anyway.");
      }
      r.head("Output at 1 m");
      const k = d.output === null ? null : num(d.output);
      if (k === null) {
        r.warn("No absolute scale. The built-in spectrum model is in arbitrary units, so dose cannot be quoted. Either select SpekPy, or measure your tube output once and press Calibrate - one measurement fixes the scale at every kVp and filtration.");
        r.note("Beam quality figures above are unaffected - HVL, mean and effective energy are all ratios and need no absolute scale.");
      } else {
        const kmin = k * 60.0 / 1000.0;
        r.row("Basis", outputBasis());
        if (CALIBRATION.scale !== null) {
          if (CALIBRATION.kvp !== kvp) r.note("Calibrated at " + fmt(CALIBRATION.kvp) + " kVp, now showing " + fmt(kvp) + " kVp - the model is carrying the scale across, which is the point of calibrating.");
          if (CALIBRATION.spek_ratio) { r.row("Correction to SpekPy", "x " + fmt(CALIBRATION.spek_ratio, 4)); r.note("A factor near 1 means SpekPy already agreed with your measurement."); }
        }
        r.big("Air kerma", fmt(k), "uGy per mAs at 1 m");
        r.big("Same as", fmt(kmin), "mGy per mA.min at 1 m");
        const ncrp = d.archer_K;
        r.row("Archer 1994", fmt(ncrp), "mGy per mA.min");
        if (ncrp > 0) r.row("Difference", signed(100 * (kmin - ncrp) / ncrp, 1) + " %");
        r.note("Archer measured a particular beam quality, so agreement depends on matching its filtration. A difference of tens of per cent is normal and does not imply either is wrong.");
        r.note("K = SUM Phi(E) . E . (mu_en/rho)_air(E), with air coefficients from the NIST Hubbell and Seltzer tables. No Monte Carlo is involved - this is an integral over the spectrum.");
      }
      const rows = sp.E.map((e, i) => ({ E_keV: e, unfiltered: num(unf.N[i]), filtered: num(sp.N[i]) }));
      r.head("Spectrum data");
      const step = Math.max(1, Math.floor(rows.length / 25));
      r.frame(rows, ["E_keV", "unfiltered", "filtered"], rows.filter((_, i) => i % step === 0));
      r.note("Full spectrum available via Export.");
      r.top();
      this.drawSpectrum(sp, kvp, d.mean);
    } catch (e) { r.head("Input error"); r.warn(e.message); }
  }
}
