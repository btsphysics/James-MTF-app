/* Monte Carlo tab - port of MonteCarloTab */
"use strict";

class MonteCarloTab {
  constructor(root, app, spectrumTab) {
    this.app = app; this.spectrumTab = spectrumTab; this.state = null; this.running = false;
    root.append(el("p", { class: "intro" }, "Monte Carlo photon transport through a slab."));
    const body = el("div", { class: "body" }); root.append(body);
    const left = el("div", { class: "left" }); const right = el("div", { class: "right" }); body.append(left, right);

    this.source = new Combo("Source", ["X-ray spectrum", "Isotope"], "X-ray spectrum");
    this.isotope = new Combo("Isotope", [], "Cs-137");
    this.isoRow = this.isotope.root;
    left.append(section("Source", this.source.root, this.isoRow));
    this.source.select.onchange = () => this.sourceChanged();
    this.isotope.select.onchange = () => this.updateArcher();

    this.material = new Combo("Material", matOpts(BARRIER_MATERIALS), "lead");
    this.material.select.onchange = () => this.updateArcher();
    this.tmin = new Entry("Thickness from", 0.25, "mm"); this.tmax = new Entry("Thickness to", 3.0, "mm");
    this.nsteps = new Entry("Number of steps", 7, ""); this.density = new Entry("Density override", "", "g/cm3");
    left.append(section("Barrier", this.material.root, this.tmin.root, this.tmax.root, this.nsteps.root, this.density.root));

    this.nhist = new Entry("Histories per point", 1000000, "");
    this.score = new Combo("Score", ["Air kerma", "Energy fluence"], "Air kerma");
    this.accept = new Entry("Detector acceptance angle", 90, "deg");
    this.ethr = new Entry("Detector energy threshold", 0, "keV");
    this.showMu = el("input", { type: "checkbox", checked: "checked" });
    this.showArcher = el("input", { type: "checkbox", checked: "checked" });
    this.archerLbl = el("label", { class: "chk" }, this.showArcher, " overlay Archer reference");
    this.showMu.onchange = () => this.redraw(); this.showArcher.onchange = () => this.redraw();
    left.append(section("Simulation", this.nhist.root, hint("Coherent scatter uses real atomic form factors, so it is computed, not bracketed."),
      this.score.root, this.accept.root, this.ethr.root,
      el("label", { class: "chk" }, this.showMu, " overlay attenuation coefficients"), this.archerLbl));
    this.workersHint = hint(""); left.append(this.workersHint);
    this.btn = button("Run", () => this.run(), "accent big"); left.append(this.btn);
    this.prog = el("progress", { max: "100", value: "0" }); left.append(this.prog);
    this.progLbl = el("div", { class: "hint" }, ""); left.append(this.progLbl);

    this.slider = el("input", { type: "range", min: "0", max: "0", step: "1", value: "0", "aria-label": "thickness shown" });
    this.sliderLbl = el("div", { class: "hint" }, "run first");
    this.slider.oninput = () => this.sliderMoved();
    left.append(section("Spectrum shown", this.slider, this.sliderLbl));

    const pane = makePane(right);
    this.results = new ResultPane(pane.top, "Results", "montecarlo.csv");
    this.fig = new Figure(pane.bot);
    this.sourceChanged();
  }
  setTables(t) {
    this.tables = t;
    const names = t.isotopes.map((x) => x[0]).sort();
    this.isotope.setValues(names, "Cs-137");
    this.isoLines = Object.fromEntries(t.isotopes);
  }
  setMaterials() { const v = this.material.get(); this.material.setValues(matOpts(BARRIER_MATERIALS), BARRIER_MATERIALS.includes(v) ? v : BARRIER_MATERIALS[0]); }
  sourceChanged() { this.isoRow.style.display = this.source.get() === "Isotope" ? "" : "none"; this.updateArcher(); }
  async updateArcher() {
    if (this.source.get() !== "Isotope") { this.showArcher.disabled = false; this.archerLbl.classList.remove("disabled"); return; }
    try {
      const d = await this.app.engine.call({ cmd: "archer", mat: this.material.get(), x: [], iso: this.isotope.get(), mean_energy: 0 });
      this.showArcher.disabled = !d.has_fit; this.showArcher.checked = d.has_fit;
      this.archerLbl.classList.toggle("disabled", !d.has_fit);
    } catch (e) { /* engine not ready yet */ }
  }
  sliderMoved() { if (this.state && this.state.specs.length) this.draw(Math.max(0, Math.min(Number(this.slider.value), this.state.specs.length - 1))); }
  redraw() { if (this.state) this.draw(Number(this.slider.value)); }

  draw(idx) {
    const st = this.state, { rows, sp, mat, dens, specs, mu } = st;
    idx = specs.length ? Math.max(0, Math.min(idx, specs.length - 1)) : 0;
    const sel = specs.length ? specs[idx] : null;
    if (sel) { this.sliderLbl.textContent = fmt(sel.t) + " mm of " + LBL(mat); this.sliderLbl.style.color = OK; }
    const [ax, ax2] = this.fig.reset(2);
    const x = rows.map((r) => r.thickness_mm);
    ax.ylog = true;
    ax.plot(x, rows.map((r) => r.B_montecarlo), { color: ACCENT, lw: 2, marker: "o", ms: 5, label: "Monte Carlo" });
    ax.plot(x, rows.map((r) => r.B_narrow), { color: "#999", lw: 1.4, dash: [6, 4], label: "narrow beam - physical floor" });
    if (this.showArcher.checked && rows.every((r) => r.B_archer !== null)) ax.plot(x, rows.map((r) => r.B_archer), { color: WARN, lw: 1.6, dash: [2, 3], marker: "s", ms: 4, label: "Archer empirical" });
    ax.xlabel = "thickness (mm)"; ax.ylabel = "transmission B";
    ax.title = LBL(mat) + ", " + (sp.model.startsWith("isotope:") ? sp.model : fmt(sp.kvp) + " kVp");
    if (sel) ax.vline(sel.t, { color: "#111", lw: 0.9, dash: [2, 3], alpha: 0.7 });
    ax.legend({ size: 10 });

    const N = nums(sp.N), m = Math.max(amax(N), 1e-30);
    ax2.plot(sp.E, N.map((v) => v / m), { color: "#999", lw: 1.3, label: "incident" });
    if (sel && sel.w.reduce((a, b) => a + b, 0) > 0) {
      const mw = Math.max(amax(sel.w), 1e-30), Hw = sel.w.map((v) => v / mw);
      ax2.fill(sel.E, Hw, { color: ACCENT, alpha: 0.22 });
      ax2.plot(sel.E, Hw, { color: ACCENT, lw: 1.8, label: "transmitted, " + fmt(sel.t) + " mm" });
      let meanOut = 0; for (let i = 0; i < sel.E.length; i++) meanOut += sel.E[i] * sel.w[i];
      ax2.vline(meanOut, { color: "black", lw: 1, dash: [5, 4] });
      ax2.text(meanOut, 0.9, "mean " + meanOut.toFixed(1) + " keV", { color: "black", size: 10, dx: 4 });
    }
    ax2.xlim = [0, sp.kvp]; ax2.ylim = [0, 1.45];
    ax2.xlabel = "photon energy (keV)"; ax2.ylabel = "relative fluence per keV";
    ax2.title = "Spectrum reaching the detector, " + (sel ? fmt(sel.t) : "-") + " mm";
    if (this.showMu.checked && mu) {
      const a3 = ax2.twinx(); a3.ylog = true;
      a3.plot(mu.E, nums(mu.total), { color: "#444444", lw: 4, alpha: 0.25, label: "total (sum)", cap: "round" });
      for (const [k, c, l] of [["photo", "#c1666b", "photoelectric"], ["coh", "#1a7f4b", "coherent"], ["incoh", "#7a5195", "Compton"]])
        a3.plot(mu.E, nums(mu[k]), { color: c, lw: 1.3, dash: [2, 2], alpha: 0.95, label: l });
      a3.ylabel = "attenuation coefficient (cm⁻¹)";
      ax2.legend({ loc: "upper center", ncol: 3, size: 9 });
    } else ax2.legend({ loc: "upper center", ncol: 2, size: 10 });
    this.fig.draw();
  }

  async run() {
    if (this.running) return;
    const r = this.results; r.clear();
    try {
      const isIso = this.source.get() === "Isotope";
      let sp, beyond = [], isoName = null;
      if (isIso) {
        isoName = this.isotope.get();
        sp = await this.app.engine.call({ cmd: "isotope", name: isoName });
        beyond = nums(sp.beyond);
      } else {
        if (!this.spectrumTab.spectrum) await this.spectrumTab.run();
        sp = this.spectrumTab.spectrum;
        if (!sp) throw new UserError("Generate a spectrum first.");
      }
      const mat = this.material.get();
      const t0 = this.tmin.get(), t1 = this.tmax.get();
      const ns = Math.max(Math.trunc(this.nsteps.get(7)), 1), nh = Math.max(Math.trunc(this.nhist.get(20000)), 100);
      const dens = this.density.get(0.0) || null;
      const acc = Math.min(Math.max(this.accept.get(90.0), 0.5), 90.0), ethr = Math.max(this.ethr.get(0.0), 0.0);
      if (t0 <= 0 || t1 < t0) throw new UserError("Thickness range is not valid.");
      const ts = ns > 1 ? linspace(t0, t1, ns) : [t0];
      const sc = this.score.get().startsWith("Air") ? "kerma" : "fluence";
      const meanE = isIso ? sp.mean : (this.spectrumTab.result ? this.spectrumTab.result.mean : 0);

      this.running = true; this.btn.disabled = true; this.prog.value = 0;
      const t_start = performance.now();
      const prog = ts.map(() => 0);
      const upd = () => { const f = prog.reduce((a, b) => a + b, 0) / ts.length; this.prog.value = 100 * f;
        this.progLbl.textContent = (100 * f).toFixed(0) + " % of " + (nh * ts.length).toLocaleString() + " histories"; };
      const spj = { E: sp.E, N: sp.N, kvp: sp.kvp, model: sp.model };
      const jobs = ts.map((t, i) => this.app.engine.mc({ sp: spj, material: mat, thickness: t, n_hist: nh, density: dens, seed: 1234, batch: 20000,
        score: sc, acceptance: acc, e_threshold: ethr, fluorescence: true }, (f) => { prog[i] = f; upd(); }));
      const archer = this.showArcher.checked
        ? this.app.engine.call(isIso ? { cmd: "archer", mat, x: ts, iso: isoName, mean_energy: meanE } : { cmd: "archer", mat, x: ts, kvp: sp.kvp })
        : Promise.resolve(null);
      const muP = this.app.engine.call({ cmd: "mu_curves", mat, density: dens, emin: Math.max(Math.min(...sp.E), 5.0), emax: sp.kvp, n: 600 });
      const dP = this.app.engine.call({ cmd: "density", mat });
      const res = await Promise.all(jobs);
      const ar = await archer, mu = await muP, rho = (await dP).density;
      const secs = (performance.now() - t_start) / 1000;
      this.prog.value = 100; this.progLbl.textContent = "done in " + secs.toFixed(1) + " s";

      const rows = [], specs = [];
      res.forEach((m, i) => {
        const a = ar ? (ar.B[i] === null ? null : num(ar.B[i])) : null;
        const Bmc = num(m.B_mc);
        specs.push({ t: ts[i], E: nums(m.hist_E), w: nums(m.hist_norm) });
        rows.push({ thickness_mm: Math.round(ts[i] * 1e4) / 1e4, B_narrow: num(m.B_narrow), B_montecarlo: Bmc, buildup: num(m.buildup),
          noise_pc: 100 * num(m.rel_std_err), B_archer: a, MC_over_Archer: a ? Bmc / a : null });
      });
      const cols = ["thickness_mm", "B_narrow", "B_montecarlo", "buildup", "noise_pc", "B_archer", "MC_over_Archer"];

      r.head("Monte Carlo transmission through " + LBL(mat));
      if (isIso) { r.row("Source", "isotope: " + isoName); r.row("Lines", this.isoLines[isoName].map(([e, i]) => fmt(e) + " keV (" + (100 * i).toFixed(1) + "%)").join(", ")); }
      else r.row("Spectrum", fmt(sp.kvp) + " kVp, " + sp.model);
      r.row("Mean energy", fmt(meanE), "keV");
      r.row("Histories per point", String(nh));
      r.row("Coherent scatter", "atomic form factors");
      r.row("Detector acceptance", fmt(acc), "deg half-angle");
      r.row("Detector threshold", fmt(ethr), "keV");
      r.row("K fluorescence", "tracked");
      r.row("Score", this.score.get());
      r.row("Density", fmt(dens ? dens : rho), "g/cm3");
      r.row("Run time", secs.toFixed(1) + " s on " + this.app.engine.workers.length + " worker(s)");
      if (isIso && this.showArcher.checked && ar && ar.basis) {
        r.row("Archer overlay", "radionuclide fit (" + ar.basis + ", source units: " + ar.unit + ")");
        if (ar.basis === "inferred") r.warn("This nuclide's Archer fit is flagged 'inferred' at the source - both unit readings were arguable and the one consistent with " + isoName + "'s photon energies was taken. Worth re-checking before relying on it.");
      }
      if (beyond.length) r.warn(beyond.length + " line(s) at " + beyond.map((e) => fmt(e)).join(", ") + " keV exceed xsdata.npz's 800 keV cross-section range (the Elam tables it comes from are flagged unreliable above there too). Those lines are transported using the frozen 800 keV attenuation coefficient, which is measurably wrong - real mu keeps falling with energy in this Compton-dominated regime, so transmission for those lines is UNDERESTIMATED here. Do not use this result for a Co-60 shielding decision; a validated code (OpenMC, EGSnrc, PENELOPE) or a published Co-60 buildup factor table is needed.");
      r.head("Results");
      const disp = rows.map((x) => ({ thickness_mm: x.thickness_mm,
        B_narrow: fmt(x.B_narrow, 4), B_montecarlo: fmt(x.B_montecarlo, 4), buildup: fmt(x.buildup, 3), noise_pc: x.noise_pc.toFixed(1),
        B_archer: x.B_archer !== null ? fmt(x.B_archer, 4) : "-", MC_over_Archer: x.MC_over_Archer !== null ? fmt(x.MC_over_Archer, 3) : "-" }));
      r.frame(rows, cols, disp);
      r.head("Interpretation");
      r.row("Buildup at thinnest point", fmt(rows[0].buildup, 3));
      r.row("Buildup at thickest point", fmt(rows[rows.length - 1].buildup, 3));
      const allA = rows.every((x) => x.MC_over_Archer !== null);
      if (allA) {
        const v = rows.map((x) => x.MC_over_Archer);
        r.row("MC / Archer, mean", fmt(v.reduce((a, b) => a + b, 0) / v.length, 3));
        r.row("MC / Archer, range", fmt(Math.min(...v), 3) + " to " + fmt(Math.max(...v), 3));
      }
      const below = rows.every((x) => x.B_archer !== null) ? rows.filter((x) => x.B_archer < x.B_narrow) : [];
      if (below.length) {
        r.warn("Beyond " + fmt(below[0].thickness_mm) + " mm Archer drops below the no-scatter floor.");
        let al = null, muTop = null;
        if (ar && ar.alpha !== null && ar.mu_top !== null) {
          al = num(ar.alpha); muTop = num(ar.mu_top);
          r.row("Archer alpha", fmt(al, 4), "/mm");
          if (isIso) r.row("mu at mean energy (" + fmt(meanE) + " keV)", fmt(muTop, 4), "/mm");
          else r.row("mu at " + fmt(sp.kvp) + " keV", fmt(muTop, 4), "/mm");
        }
        if (al !== null && muTop !== null && al > muTop) {
          r.note("Alpha is the fit's asymptotic decay rate. It exceeds mu at the beam's representative energy, so the fit decays faster than the most penetrating photons in the beam physically can. Crossing the floor is then unavoidable at some depth.");
          r.note("That is not an error in Archer's measurements. It is a three-parameter form being asked to describe attenuation far beyond the data it was fitted to. Ignore the comparison there.");
        }
      }
      const worst = Math.max(...rows.map((x) => x.noise_pc));
      r.row("Worst random noise", worst.toFixed(1) + " %");
      if (worst > 15) r.warn("Noise above 15 %. Increase histories or drop the thickest points.");
      r.top();
      this.state = { rows, sp: Object.assign({}, sp), mat, dens, specs, mu };
      const n = Math.max(specs.length - 1, 0);
      this.slider.max = String(n); this.slider.value = String(n);
      this.draw(n);
    } catch (e) {
      r.head(e instanceof UserError ? "Input error" : "Error"); r.warn(e.message);
    } finally { this.running = false; this.btn.disabled = false; }
  }
}
