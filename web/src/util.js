/* Python-compatible formatting, form widgets, results pane, table formatter. */
"use strict";

const BG = "#f4f6f8", CARD = "#ffffff", ACCENT = "#1f5f8b", OK = "#1a7f4b", WARN = "#b3541e";

function num(v) { if (v === null || v === undefined) return NaN; if (v === "inf") return Infinity; if (v === "-inf") return -Infinity; return Number(v); }
function nums(a) { return a.map(num); }

/* Python "%.*e" (two-digit exponent) */
function pyExp(v, d) {
  let s = v.toExponential(d);
  return s.replace(/e([+-])(\d)$/, "e$10$2");
}
/* Python "%.6g" */
function pyG(v, p = 6) {
  if (!isFinite(v)) return isNaN(v) ? "nan" : (v > 0 ? "inf" : "-inf");
  if (v === 0) return "0";
  const e = Math.floor(Math.log10(Math.abs(Number(v.toPrecision(p)))));
  if (e < -4 || e >= p) {
    let s = pyExp(v, p - 1);
    s = s.replace(/\.?0+e/, "e");
    return s;
  }
  let s = v.toFixed(Math.max(0, p - 1 - e));
  if (s.indexOf(".") >= 0) s = s.replace(/\.?0+$/, "");
  return s;
}
/* app.py's fmt() */
function fmt(v, sig = 4) {
  if (v === null || v === undefined) return "n/a";
  if (typeof v === "string") return v;
  v = Number(v);
  if (isNaN(v)) return "n/a";
  if (!isFinite(v)) return "inf";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e5 || a < 1e-3) return pyExp(v, sig - 1);
  if (a >= 100) return v.toFixed(1);
  if (a >= 10) return v.toFixed(2);
  if (a >= 1) return v.toFixed(3);
  return v.toFixed(sig + 1);
}
function signed(v, d) { return (v >= 0 ? "+" : "") + v.toFixed(d); }
function linspace(a, b, n) { if (n === 1) return [a]; const o = []; for (let i = 0; i < n; i++) o.push(i === n - 1 ? b : a + (b - a) * i / (n - 1)); return o; }
function interp1(x, xp, fp) {
  const n = xp.length; if (!n) return NaN;
  if (x <= xp[0]) return fp[0]; if (x >= xp[n - 1]) return fp[n - 1];
  let lo = 0, hi = n - 1; while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xp[m] <= x) lo = m; else hi = m; }
  return fp[lo] + (fp[hi] - fp[lo]) * (x - xp[lo]) / (xp[hi] - xp[lo]);
}
function amax(a) { let m = -Infinity; for (const v of a) if (v > m) m = v; return m; }
function nanmax(a) { let m = -Infinity; for (const v of a) if (isFinite(v) && v > m) m = v; return m === -Infinity ? NaN : m; }
function el(tag, attrs, ...kids) {
  const e = document.createElement(tag);
  if (attrs) for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") e.className = v; else if (k === "style") e.style.cssText = v;
    else if (k.startsWith("on")) e.addEventListener(k.slice(2), v); else if (v !== null && v !== undefined && v !== false) e.setAttribute(k, v);
  }
  for (const k of kids) if (k !== null && k !== undefined) e.append(k.nodeType ? k : document.createTextNode(String(k)));
  return e;
}

/* --- form widgets (LabelledEntry / LabelledCombo) --- */
class Entry {
  constructor(label, def, unit, opts = {}) {
    this.input = el("input", { type: "text", value: def === "" ? "" : String(def), spellcheck: "false", inputmode: "decimal", "aria-label": label || unit || "value" });
    if (opts.width) this.input.style.width = opts.width;
    this.root = el("div", { class: "row" + (opts.compact ? " compact" : "") },
      label !== "" ? el("label", { class: "lbl" }, label) : null, this.input, unit ? el("span", { class: "unit" }, unit) : null);
    if (opts.title) this.root.title = opts.title;
  }
  get(def) {
    const t = this.input.value.trim();
    if (!t) { if (def === undefined) throw new UserError("A required value is blank."); return def; }
    const v = Number(t);
    if (!isFinite(v) || !/^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/.test(t)) throw new UserError("'" + t + "' is not a number.");
    return v;
  }
  set(v) { this.input.value = String(v); }
}
class Combo {
  constructor(label, values, def, opts = {}) {
    this.select = el("select", { "aria-label": label || "choice" });
    this.setValues(values, def);
    if (opts.width) this.select.style.width = opts.width;
    this.root = el("div", { class: "row" + (opts.compact ? " compact" : "") + (opts.stack ? " stack" : "") }, label !== "" ? el("label", { class: "lbl" }, label) : null, this.select);
  }
  setValues(values, def) {
    this.select.innerHTML = "";
    for (const v of values) { const [val, txt] = Array.isArray(v) ? v : [v, v]; this.select.append(el("option", { value: val }, txt)); }
    if (def !== undefined) this.select.value = def;
  }
  get() { return this.select.value; }
  set(v) { this.select.value = v; }
}
class UserError extends Error {}
function section(title, ...kids) { return el("fieldset", { class: "card" }, el("legend", null, title), ...kids); }
function hint(text) { return el("div", { class: "hint" }, text); }
function button(text, onclick, cls) { return el("button", { class: "btn" + (cls ? " " + cls : ""), onclick }, text); }

/* --- results pane (ResultPane) --- */
class ResultPane {
  constructor(host, title = "Results", exportName = "results.csv") {
    this.exportName = exportName; this.df = null;
    this.exportBtn = button("Export", () => this.export()); this.exportBtn.disabled = true;
    this.body = el("div", { class: "res-body", tabindex: "0" });
    host.append(el("div", { class: "res-bar" }, el("span", { class: "res-title" }, title),
      el("span", { class: "spacer" }), this.exportBtn, button("Copy", () => this.copy())), this.body);
  }
  clear() { this.body.innerHTML = ""; this.df = null; this.exportBtn.disabled = true; }
  line(cls, ...parts) { const d = el("div", { class: cls }, ...parts); this.body.append(d); return d; }
  head(s) { this.line("r-h", s); }
  row(k, v, u = "", w = 30) { this.line("r-row", "  " + String(k).padEnd(w), el("b", null, String(v)), u ? " " + u : ""); }
  big(k, v, u = "") { this.line("r-row", "  " + k.padEnd(30), el("span", { class: "r-big" }, String(v) + (u ? " " + u : ""))); }
  note(s) { this.line("r-note", s); }
  warn(s) { this.line("r-warn", "! " + s); }
  frame(rows, cols, displayRows) {
    this.df = { rows, cols }; this.exportBtn.disabled = false;
    const txt = frameToString(displayRows || rows, cols);
    this.line("r-frame", txt);
  }
  setDF(rows, cols) { this.df = { rows, cols }; this.exportBtn.disabled = false; }
  top() { this.body.scrollTop = 0; }
  text() {
    const out = [];
    for (const d of this.body.children) {
      if (d.classList.contains("r-h")) out.push("", d.textContent, "-".repeat(Math.max(d.textContent.length, 48)));
      else if (d.classList.contains("r-note")) out.push("  " + d.textContent);
      else if (d.classList.contains("r-warn")) out.push("  " + d.textContent);
      else if (d.classList.contains("r-frame")) out.push(...d.textContent.split("\n").map((l) => "  " + l));
      else out.push(d.textContent);
    }
    return out.join("\n").replace(/^\n/, "");
  }
  copy() { copyText(this.text()); }
  export() {
    if (!this.df) return;
    const { rows, cols } = this.df;
    const esc = (v) => { const s = v === null || v === undefined ? "" : (typeof v === "number" ? (isFinite(v) ? String(v) : "") : String(v)); return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; };
    const csv = [cols.join(",")].concat(rows.map((r) => cols.map((c) => esc(r[c])).join(","))).join("\n");
    downloadText(this.exportName, csv, "text/csv");
  }
}

/* pandas DataFrame.to_string(index=False), approximately */
function frameToString(rows, cols) {
  const colStr = cols.map((c) => {
    const vals = rows.map((r) => r[c]);
    const allNum = vals.every((v) => typeof v === "number" || v === null || v === undefined);
    let strs;
    if (allNum && vals.some((v) => typeof v === "number")) {
      const isInt = vals.every((v) => v === null || v === undefined || (Number.isInteger(v) && Math.abs(v) < 1e15));
      if (isInt) strs = vals.map((v) => v === null || v === undefined ? "None" : String(v));
      else {
        const fin = vals.filter((v) => typeof v === "number" && isFinite(v));
        const big = fin.some((v) => Math.abs(v) >= 1e6 || (v !== 0 && Math.abs(v) < 1e-4));
        if (big) strs = vals.map((v) => typeof v === "number" ? (isFinite(v) ? pyExp(v, 6) : (isNaN(v) ? "NaN" : "inf")) : "None");
        else {
          let dec = 1;
          for (const v of fin) { const s = pyG(v, 7); const i = s.indexOf("."); if (i >= 0) dec = Math.max(dec, Math.min(6, s.length - i - 1)); }
          strs = vals.map((v) => typeof v === "number" ? (isFinite(v) ? v.toFixed(dec) : (isNaN(v) ? "NaN" : "inf")) : "None");
        }
      }
    } else strs = vals.map((v) => v === null || v === undefined ? "None" : String(v));
    const w = Math.max(c.length, ...strs.map((s) => s.length));
    return { w, head: c.padStart(w), strs: strs.map((s) => s.padStart(w)) };
  });
  const lines = [colStr.map((c) => c.head).join(" ")];
  for (let i = 0; i < rows.length; i++) lines.push(colStr.map((c) => c.strs[i]).join(" "));
  return lines.join("\n");
}

function copyText(text) {
  const fallback = () => {
    const ta = el("textarea", { style: "position:fixed;left:-9999px;top:0" }); ta.value = text; document.body.append(ta);
    ta.select(); try { document.execCommand("copy"); } catch (e) { /* ignore */ } ta.remove();
  };
  if (navigator.clipboard && window.isSecureContext) navigator.clipboard.writeText(text).catch(fallback); else fallback();
  toast("Copied to clipboard");
}
function downloadText(name, text, type) {
  const a = el("a", { download: name, href: URL.createObjectURL(new Blob([text], { type })) });
  document.body.append(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
}
let toastTimer = null;
function toast(msg) {
  let t = document.getElementById("toast");
  if (!t) { t = el("div", { id: "toast" }); document.body.append(t); }
  t.textContent = msg; t.classList.add("show"); clearTimeout(toastTimer); toastTimer = setTimeout(() => t.classList.remove("show"), 1600);
}

/* vertical split between results (top) and plot (bottom), like ttk.PanedWindow */
function makePane(host, topWeight = 2, bottomWeight = 3) {
  const top = el("div", { class: "pane-top" }), grip = el("div", { class: "pane-grip", title: "Drag to resize" }), bot = el("div", { class: "pane-bot" });
  host.classList.add("pane"); host.append(top, grip, bot);
  top.style.flex = topWeight + " 1 0"; bot.style.flex = bottomWeight + " 1 0";
  grip.addEventListener("pointerdown", (e) => {
    e.preventDefault(); grip.setPointerCapture(e.pointerId);
    const H = host.clientHeight, y0 = e.clientY, h0 = top.getBoundingClientRect().height;
    const mv = (ev) => { const h = Math.min(Math.max(h0 + ev.clientY - y0, 60), H - 120); top.style.flex = "0 0 " + h + "px"; bot.style.flex = "1 1 0"; };
    const up = () => { grip.removeEventListener("pointermove", mv); grip.removeEventListener("pointerup", up); };
    grip.addEventListener("pointermove", mv); grip.addEventListener("pointerup", up);
  });
  return { top, bot };
}

/* --- DICOM image viewer with zoom / pan and ROI rectangle selection --- */
class ImageViewer {
  constructor(host, engine) {
    this.host = host; this.engine = engine; host.classList.add("iv-host");
    this.canvas = el("canvas"); this.info = el("div", { class: "iv-info" });
    this.tools = el("div", { class: "iv-tools" },
      button("Fit", () => { this.fit(); this.request(); }),
      el("label", null, "Contrast ", this.contrast = el("select", { "aria-label": "contrast" }, el("option", { value: "full" }, "full range"), el("option", { value: "p" }, "1-99 %"))));
    this.contrast.onchange = () => this.request();
    host.append(this.canvas, this.tools, this.info);
    this.slot = null; this.meta = null; this.view = null; this.imgData = null; this.rois = []; this.onROI = null; this.title = "";
    this.drag = null; this.reqSeq = 0;
    new ResizeObserver(() => { if (this.meta) { this.fitIfNeeded(); this.request(); } }).observe(host);
    const c = this.canvas;
    c.addEventListener("wheel", (e) => this.onWheel(e), { passive: false });
    c.addEventListener("pointerdown", (e) => this.onDown(e));
    c.addEventListener("pointermove", (e) => this.onMove(e));
    c.addEventListener("pointerup", (e) => this.onUp(e));
    c.addEventListener("contextmenu", (e) => e.preventDefault());
  }
  show(slot, meta, title) { this.slot = slot; this.meta = meta; this.title = title || ""; this.view = null; this.fitIfNeeded(); this.request(); }
  geom() { const w = this.host.clientWidth, h = this.host.clientHeight - 30; return { w: Math.max(w, 10), h: Math.max(h, 10), top: 30 }; }
  fit() {
    const { w, h } = this.geom(), m = this.meta; const s = Math.max(m.w / w, m.h / h);
    const cx = m.w / 2, cy = m.h / 2; this.view = { x0: cx - w * s / 2, y0: cy - h * s / 2, s };
  }
  fitIfNeeded() { if (!this.view) this.fit(); }
  toImg(px, py) { const v = this.view; return [v.x0 + px * v.s, v.y0 + py * v.s]; }
  toScr(x, y) { const v = this.view; return [(x - v.x0) / v.s, (y - v.y0) / v.s]; }
  async request() {
    if (!this.meta || !this.view) return;
    const { w, h } = this.geom(), v = this.view, dpr = Math.min(window.devicePixelRatio || 1, 2);
    const ow = Math.round(w * dpr), oh = Math.round(h * dpr), seq = ++this.reqSeq;
    const m = this.meta, pc = this.contrast.value === "p";
    try {
      const r = await this.engine.call({ cmd: "image_view", slot: this.slot, x0: v.x0, y0: v.y0, x1: v.x0 + w * v.s, y1: v.y0 + h * v.s,
        out_w: ow, out_h: oh, vmin: pc ? m.p1 : m.vmin, vmax: pc ? m.p99 : m.vmax, invert: m.sign === -1 });
      if (seq !== this.reqSeq) return;
      const rgba = new Uint8ClampedArray(ow * oh * 4);
      for (let i = 0, j = 0; i < r.bin.length; i++, j += 4) { const g = r.bin[i]; rgba[j] = rgba[j + 1] = rgba[j + 2] = g; rgba[j + 3] = 255; }
      this.imgData = new ImageData(rgba, ow, oh); this.dpr = dpr; this.draw();
    } catch (e) { this.info.textContent = String(e.message || e); }
  }
  draw() {
    const { w, h, top } = this.geom(), c = this.canvas, dpr = this.dpr || 1;
    c.width = Math.round(w * dpr); c.height = Math.round((h + top) * dpr); c.style.width = w + "px"; c.style.height = (h + top) + "px";
    const g = c.getContext("2d");
    g.setTransform(1, 0, 0, 1, 0, 0); g.fillStyle = "#fff"; g.fillRect(0, 0, c.width, c.height);
    if (this.imgData) g.putImageData(this.imgData, 0, Math.round(top * dpr));
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.font = "12px " + PLOT_FONT; g.fillStyle = "#222"; g.textBaseline = "middle"; g.textAlign = "left";
    g.fillText(this.title, 8, 14);
    g.save(); g.translate(0, top); g.beginPath(); g.rect(0, 0, w, h); g.clip();
    const legend = [];
    for (const r of this.rois) {
      if (!r.roi) continue;
      const [r0, r1, c0, c1] = r.roi; const [x0, y0] = this.toScr(c0, r0), [x1, y1] = this.toScr(c1, r1);
      g.strokeStyle = r.color; g.lineWidth = 1.5; g.strokeRect(x0, y0, x1 - x0, y1 - y0); legend.push(r);
    }
    if (this.drag && this.drag.mode === "roi" && this.drag.cur) {
      const d = this.drag; const [x0, y0] = this.toScr(d.ix, d.iy), [x1, y1] = this.toScr(d.cur[0], d.cur[1]);
      g.setLineDash([5, 4]); g.strokeStyle = "#ffcc00"; g.lineWidth = 1.5; g.strokeRect(x0, y0, x1 - x0, y1 - y0); g.setLineDash([]);
    }
    if (legend.length) {
      g.font = "11px " + PLOT_FONT; const bw = 120, bh = legend.length * 16 + 8;
      g.fillStyle = "rgba(255,255,255,0.9)"; g.fillRect(w - bw - 8, 8, bw, bh);
      legend.forEach((r, i) => { g.strokeStyle = r.color; g.lineWidth = 2; g.beginPath(); g.moveTo(w - bw, 20 + i * 16); g.lineTo(w - bw + 20, 20 + i * 16); g.stroke();
        g.fillStyle = "#222"; g.fillText(r.label, w - bw + 26, 20 + i * 16); });
    }
    g.restore();
  }
  pos(e) { const r = this.canvas.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top - this.geom().top]; }
  onWheel(e) {
    if (!this.view) return; e.preventDefault();
    const [px, py] = this.pos(e); const [ix, iy] = this.toImg(px, py);
    const f = Math.exp(e.deltaY * 0.0015); const v = this.view; const s = Math.min(Math.max(v.s * f, 0.02), 100);
    this.view = { x0: ix - px * s, y0: iy - py * s, s }; this.draw(); clearTimeout(this.wt); this.wt = setTimeout(() => this.request(), 60);
  }
  onDown(e) {
    if (!this.view) return; this.canvas.setPointerCapture(e.pointerId);
    const [px, py] = this.pos(e);
    if (e.button === 1 || e.button === 2 || e.shiftKey) this.drag = { mode: "pan", px, py, v: Object.assign({}, this.view) };
    else { const [ix, iy] = this.toImg(px, py); this.drag = { mode: "roi", ix, iy, cur: null }; }
  }
  onMove(e) {
    if (!this.view) return;
    const [px, py] = this.pos(e); const [ix, iy] = this.toImg(px, py);
    this.info.textContent = (ix >= 0 && iy >= 0 && ix < this.meta.w && iy < this.meta.h) ? "col " + Math.floor(ix) + ", row " + Math.floor(iy) : "";
    if (!this.drag) return;
    if (this.drag.mode === "pan") { const d = this.drag; this.view = { x0: d.v.x0 - (px - d.px) * d.v.s, y0: d.v.y0 - (py - d.py) * d.v.s, s: d.v.s }; this.draw(); clearTimeout(this.wt); this.wt = setTimeout(() => this.request(), 60); }
    else { this.drag.cur = [ix, iy]; this.draw(); }
  }
  onUp() {
    const d = this.drag; this.drag = null; if (!d) return;
    if (d.mode === "pan") { this.request(); return; }
    if (!d.cur) { this.draw(); return; }
    const x0 = Math.min(d.ix, d.cur[0]), x1 = Math.max(d.ix, d.cur[0]), y0 = Math.min(d.iy, d.cur[1]), y1 = Math.max(d.iy, d.cur[1]);
    if (x1 - x0 < 10 || y1 - y0 < 10) { this.draw(); return; }   // minspanx / minspany = 10
    const clamp = (v, n) => Math.max(0, Math.min(n, v));
    const roi = [Math.trunc(clamp(y0, this.meta.h)), Math.trunc(clamp(y1, this.meta.h)), Math.trunc(clamp(x0, this.meta.w)), Math.trunc(clamp(x1, this.meta.w))];
    if (this.onROI) this.onROI(roi);
    this.draw();
  }
}
