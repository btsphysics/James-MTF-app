/* Small dependency-free canvas plotting library (lines, fills, markers,
   log axes, twin y axis, annotations, legends, PNG export). */
"use strict";

const PLOT_FONT = '"Segoe UI", system-ui, -apple-system, Roboto, Arial, sans-serif';
const SUP = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
  "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const sup = (s) => String(s).split("").map((c) => SUP[c] || c).join("");

function niceStep(span, target) {
  const raw = span / Math.max(target, 1);
  const p = Math.pow(10, Math.floor(Math.log10(raw)));
  const f = raw / p;
  return (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p;
}
function linTicks(lo, hi, target) {
  if (!(hi > lo)) return [lo];
  const st = niceStep(hi - lo, target);
  const out = [];
  for (let v = Math.ceil(lo / st - 1e-9) * st; v <= hi + st * 1e-9; v += st) out.push(Math.abs(v) < st * 1e-9 ? 0 : v);
  return out;
}
function tickLabel(v, step) {
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e5 || a < 1e-3) return v.toExponential(1).replace("e+", "e");
  const dec = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
  return v.toFixed(Math.min(dec, 6));
}
function finite(v) { return typeof v === "number" && isFinite(v); }

class Axes {
  constructor(fig) {
    this.fig = fig; this.items = []; this.title = ""; this.xlabel = ""; this.ylabel = "";
    this.xlim = null; this.ylim = null; this.ylog = false; this.grid = true; this.legendOpt = null;
    this.twin = null; this.isTwin = false; this.titleColor = "#1f5f8b";
  }
  plot(x, y, o = {}) { this.items.push(Object.assign({ type: "line", x, y, lw: 1.5, color: "#1f5f8b" }, o)); return this; }
  fill(x, y, o = {}) { this.items.push(Object.assign({ type: "fill", x, y, color: "#1f5f8b", alpha: 0.22 }, o)); return this; }
  vline(x, o = {}) { this.items.push(Object.assign({ type: "vline", x, color: "#333", lw: 1, dash: [4, 3] }, o)); return this; }
  hline(y, o = {}) { this.items.push(Object.assign({ type: "hline", y, color: "#333", lw: 1, dash: [4, 3] }, o)); return this; }
  text(x, y, s, o = {}) { this.items.push(Object.assign({ type: "text", x, y, s, color: "#333", size: 11, dx: 0, dy: 0 }, o)); return this; }
  legend(o = {}) { this.legendOpt = Object.assign({ loc: "best", ncol: 1, size: 11 }, o); return this; }
  twinx() { this.twin = new Axes(this.fig); this.twin.isTwin = true; this.twin.grid = false; return this.twin; }

  dataRange(axis) {
    let lo = Infinity, hi = -Infinity;
    for (const it of this.items) {
      if (it.type === "line" || it.type === "fill") {
        const arr = axis === "x" ? it.x : it.y;
        for (let i = 0; i < arr.length; i++) {
          const v = arr[i], o = axis === "x" ? it.y[i] : it.x[i];
          if (!finite(v) || !finite(o)) continue;
          if (axis === "y" && this.ylog && v <= 0) continue;
          if (v < lo) lo = v; if (v > hi) hi = v;
        }
        if (it.type === "fill" && axis === "y" && !this.ylog) { lo = Math.min(lo, 0); }
      } else if (it.type === "vline" && axis === "x" && finite(it.x)) { lo = Math.min(lo, it.x); hi = Math.max(hi, it.x); }
      else if (it.type === "hline" && axis === "y" && finite(it.y)) { lo = Math.min(lo, it.y); hi = Math.max(hi, it.y); }
    }
    return [lo, hi];
  }
  limits() {
    let [x0, x1] = this.xlim || this.dataRange("x");
    if (!this.xlim) { if (!finite(x0)) { x0 = 0; x1 = 1; } const m = (x1 - x0) * 0.05 || Math.abs(x0) * 0.05 || 1; x0 -= m; x1 += m; }
    let y0, y1;
    if (this.ylim) [y0, y1] = this.ylim;
    else {
      [y0, y1] = this.dataRange("y");
      if (!finite(y0)) { y0 = this.ylog ? 0.1 : 0; y1 = this.ylog ? 10 : 1; }
      if (this.ylog) { const l0 = Math.log10(y0), l1 = Math.log10(y1); const m = Math.max((l1 - l0) * 0.05, 0.02); y0 = Math.pow(10, l0 - m); y1 = Math.pow(10, l1 + m); }
      else { const m = (y1 - y0) * 0.05 || Math.abs(y0) * 0.05 || 1; y0 -= m; y1 += m; }
    }
    if (x1 === x0) x1 = x0 + 1;
    if (y1 === y0) y1 = y0 + (this.ylog ? y0 : 1);
    return { x0, x1, y0, y1 };
  }
}

class Figure {
  constructor(host) {
    this.host = host; host.classList.add("fig-host");
    this.canvas = document.createElement("canvas"); host.appendChild(this.canvas);
    this.readout = document.createElement("div"); this.readout.className = "fig-readout"; host.appendChild(this.readout);
    const btn = document.createElement("button"); btn.className = "fig-save"; btn.textContent = "Save PNG"; btn.title = "Save the plot as a PNG image";
    btn.onclick = () => this.savePNG(); host.appendChild(btn); this.saveBtn = btn;
    this.axes = []; this.layout = [];
    new ResizeObserver(() => this.draw()).observe(host);
    this.canvas.addEventListener("mousemove", (e) => this.onMove(e));
    this.canvas.addEventListener("mouseleave", () => { this.readout.textContent = ""; });
  }
  reset(ncols = 1) { this.axes = []; for (let i = 0; i < ncols; i++) this.axes.push(new Axes(this)); return this.axes; }
  savePNG() {
    const a = document.createElement("a"); a.download = "plot.png"; a.href = this.canvas.toDataURL("image/png"); a.click();
  }
  onMove(e) {
    const r = this.canvas.getBoundingClientRect(); const px = e.clientX - r.left, py = e.clientY - r.top;
    for (const L of this.layout) {
      if (px >= L.l && px <= L.r && py >= L.t && py <= L.b) {
        const x = L.xinv(px), y = L.yinv(py);
        let s = "x=" + fmtShort(x) + "  y=" + fmtShort(y);
        if (L.twin) s += "  y2=" + fmtShort(L.twin.yinv(py));
        this.readout.textContent = s; return;
      }
    }
    this.readout.textContent = "";
  }
  draw() {
    const w = this.host.clientWidth, h = this.host.clientHeight;
    if (w < 20 || h < 20) return;
    const dpr = window.devicePixelRatio || 1;
    const c = this.canvas; c.width = Math.round(w * dpr); c.height = Math.round(h * dpr);
    c.style.width = w + "px"; c.style.height = h + "px";
    const g = c.getContext("2d"); g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.fillStyle = "#ffffff"; g.fillRect(0, 0, w, h);
    this.layout = [];
    const n = this.axes.length; if (!n) return;
    const gap = 14, cellW = (w - gap * (n - 1)) / n;
    this.axes.forEach((ax, i) => this.drawAxes(g, ax, i * (cellW + gap), 0, cellW, h));
  }
  drawAxes(g, ax, X, Y, W, H) {
    const hasTwin = !!ax.twin;
    const ml = 62, mr = hasTwin ? 62 : 14, mt = ax.title ? 30 : 12, mb = 44;
    const l = X + ml, r = X + W - mr, t = Y + mt, b = Y + H - mb;
    if (r - l < 30 || b - t < 30) return;
    const lim = ax.limits();
    const xs = (v) => l + (v - lim.x0) / (lim.x1 - lim.x0) * (r - l);
    const mk = (A, L) => {
      if (A.ylog) {
        const a = Math.log10(L.y0), bb = Math.log10(L.y1);
        return { f: (v) => b - (Math.log10(v) - a) / (bb - a) * (b - t), inv: (p) => Math.pow(10, a + (b - p) / (b - t) * (bb - a)) };
      }
      return { f: (v) => b - (v - L.y0) / (L.y1 - L.y0) * (b - t), inv: (p) => L.y0 + (b - p) / (b - t) * (L.y1 - L.y0) };
    };
    const ym = mk(ax, lim);
    const entry = { l, r, t, b, xinv: (p) => lim.x0 + (p - l) / (r - l) * (lim.x1 - lim.x0), yinv: ym.inv };
    this.layout.push(entry);

    g.save(); g.font = "11px " + PLOT_FONT; g.fillStyle = "#333"; g.strokeStyle = "#333";
    // ticks & grid
    const xt = linTicks(lim.x0, lim.x1, Math.max(3, Math.floor((r - l) / 75)));
    const xstep = xt.length > 1 ? xt[1] - xt[0] : 1;
    g.textAlign = "center"; g.textBaseline = "top";
    for (const v of xt) {
      const p = xs(v); if (p < l - 0.5 || p > r + 0.5) continue;
      if (ax.grid) { g.strokeStyle = "rgba(0,0,0,0.09)"; g.beginPath(); g.moveTo(p, t); g.lineTo(p, b); g.stroke(); }
      g.strokeStyle = "#333"; g.beginPath(); g.moveTo(p, b); g.lineTo(p, b + 4); g.stroke();
      g.fillText(tickLabel(v, xstep), p, b + 6);
    }
    this.yTicks(g, ax, lim, ym.f, l, r, t, b, "left");
    // clip & draw items
    g.save(); g.beginPath(); g.rect(l, t, r - l, b - t); g.clip(); g.__clipR = r;
    this.drawItems(g, ax, xs, ym.f, t, b);
    let twinEntry = null, tl = null, tym = null;
    if (hasTwin) {
      const tw = ax.twin; tw.xlim = [lim.x0, lim.x1]; tl = tw.limits(); tym = mk(tw, tl);
      this.drawItems(g, tw, xs, tym.f, t, b);
      twinEntry = { yinv: tym.inv }; entry.twin = twinEntry;
    }
    g.restore(); g.__clipR = null;
    if (hasTwin) this.yTicks(g, ax.twin, tl, tym.f, l, r, t, b, "right");
    // frame
    g.strokeStyle = "#555"; g.lineWidth = 1; g.strokeRect(l + 0.5, t + 0.5, r - l - 1, b - t - 1);
    // labels
    g.fillStyle = "#222"; g.font = "12px " + PLOT_FONT; g.textAlign = "center"; g.textBaseline = "bottom";
    if (ax.xlabel) g.fillText(ax.xlabel, (l + r) / 2, Y + H - 4);
    if (ax.ylabel) { g.save(); g.translate(X + 13, (t + b) / 2); g.rotate(-Math.PI / 2); g.textBaseline = "middle"; g.fillText(ax.ylabel, 0, 0); g.restore(); }
    if (hasTwin && ax.twin.ylabel) { g.save(); g.translate(X + W - 10, (t + b) / 2); g.rotate(-Math.PI / 2); g.textBaseline = "middle"; g.fillText(ax.twin.ylabel, 0, 0); g.restore(); }
    if (ax.title) {
      g.font = "bold 13px " + PLOT_FONT; g.fillStyle = ax.titleColor; g.textBaseline = "middle";
      let tt = ax.title; while (g.measureText(tt).width > W - 8 && tt.length > 4) tt = tt.slice(0, -2);
      if (tt !== ax.title) tt = tt.slice(0, -1) + "…";
      g.fillText(tt, X + W / 2, Y + 15);
    }
    // legend (both axes' entries)
    if (ax.legendOpt) {
      const ents = ax.items.filter((i) => i.label && i.type !== "text");
      if (hasTwin) for (const i of ax.twin.items) if (i.label) ents.push(i);
      if (ents.length) this.drawLegend(g, ax, ents, l, r, t, b, xs, ym.f);
    }
    g.restore();
  }
  yTicks(g, ax, lim, yf, l, r, t, b, side) {
    g.font = "11px " + PLOT_FONT; g.textBaseline = "middle"; g.textAlign = side === "left" ? "right" : "left";
    const xa = side === "left" ? l : r, dir = side === "left" ? -1 : 1;
    let ticks = [], labels = [], minor = [];
    if (ax.ylog) {
      const a = Math.floor(Math.log10(lim.y0)), bb = Math.ceil(Math.log10(lim.y1));
      const span = bb - a, every = Math.max(1, Math.ceil(span / Math.max(2, Math.floor((b - t) / 28))));
      for (let e = a; e <= bb; e++) {
        if (e % every === 0) { ticks.push(Math.pow(10, e)); labels.push(e === 0 ? "1" : e === 1 ? "10" : "10" + sup(e)); }
        if (span <= 6) for (let m = 2; m <= 9; m++) minor.push(m * Math.pow(10, e));
      }
    } else {
      ticks = linTicks(lim.y0, lim.y1, Math.max(3, Math.floor((b - t) / 40)));
      const st = ticks.length > 1 ? ticks[1] - ticks[0] : 1;
      labels = ticks.map((v) => tickLabel(v, st));
    }
    for (const v of minor) {
      const p = yf(v); if (p < t || p > b) continue;
      if (ax.grid) { g.strokeStyle = "rgba(0,0,0,0.04)"; g.beginPath(); g.moveTo(l, p); g.lineTo(r, p); g.stroke(); }
      g.strokeStyle = "#555"; g.beginPath(); g.moveTo(xa, p); g.lineTo(xa + dir * 2, p); g.stroke();
    }
    ticks.forEach((v, i) => {
      const p = yf(v); if (p < t - 0.5 || p > b + 0.5) return;
      if (ax.grid) { g.strokeStyle = "rgba(0,0,0,0.09)"; g.beginPath(); g.moveTo(l, p); g.lineTo(r, p); g.stroke(); }
      g.strokeStyle = "#333"; g.beginPath(); g.moveTo(xa, p); g.lineTo(xa + dir * 4, p); g.stroke();
      g.fillStyle = "#333"; g.fillText(labels[i], xa + dir * 6, p);
    });
  }
  drawItems(g, ax, xs, yf, t, b) {
    const ok = (v) => finite(v) && (!ax.ylog || v > 0);
    for (const it of ax.items) {
      g.save(); g.globalAlpha = it.alpha != null ? it.alpha : 1;
      if (it.type === "fill") {
        g.fillStyle = it.color; g.beginPath(); let started = false, lastX = null;
        const base = ax.ylog ? b : yf(0);
        for (let i = 0; i < it.x.length; i++) {
          if (!ok(it.y[i]) || !finite(it.x[i])) continue;
          const px = xs(it.x[i]), py = yf(it.y[i]);
          if (!started) { g.moveTo(px, base); started = true; }
          g.lineTo(px, py); lastX = px;
        }
        if (started) { g.lineTo(lastX, base); g.closePath(); g.fill(); }
      } else if (it.type === "line") {
        g.strokeStyle = it.color; g.lineWidth = it.lw; g.setLineDash(it.dash || []); g.lineJoin = "round"; g.lineCap = it.cap || "butt";
        if (it.lw > 0) {
          g.beginPath(); let pen = false;
          for (let i = 0; i < it.x.length; i++) {
            if (!ok(it.y[i]) || !finite(it.x[i])) { pen = false; continue; }
            const px = xs(it.x[i]), py = yf(it.y[i]);
            if (!pen) { g.moveTo(px, py); pen = true; } else g.lineTo(px, py);
          }
          g.stroke();
        }
        if (it.marker) {
          g.setLineDash([]); g.fillStyle = it.mfc || it.color; const s = it.ms || 5;
          for (let i = 0; i < it.x.length; i++) {
            if (!ok(it.y[i]) || !finite(it.x[i])) continue;
            drawMarker(g, it.marker, xs(it.x[i]), yf(it.y[i]), s);
          }
        }
      } else if (it.type === "vline") {
        g.strokeStyle = it.color; g.lineWidth = it.lw; g.setLineDash(it.dash || []);
        const p = xs(it.x); g.beginPath(); g.moveTo(p, t); g.lineTo(p, b); g.stroke();
      } else if (it.type === "hline") {
        if (ok(it.y)) { g.strokeStyle = it.color; g.lineWidth = it.lw; g.setLineDash(it.dash || []);
          const p = yf(it.y); g.beginPath(); g.moveTo(-1e4, p); g.lineTo(1e5, p); g.stroke(); }
      } else if (it.type === "text") {
        if (ok(it.y) && finite(it.x)) {
          g.fillStyle = it.color; g.font = it.size + "px " + PLOT_FONT; g.textBaseline = "middle";
          const px = xs(it.x), tw = g.measureText(it.s).width, right = g.canvas.width / (window.devicePixelRatio || 1);
          const clip = g.__clipR || right;
          if (px + it.dx + tw > clip - 2) { g.textAlign = "right"; g.fillText(it.s, px - Math.abs(it.dx), yf(it.y) - it.dy); }
          else { g.textAlign = "left"; g.fillText(it.s, px + it.dx, yf(it.y) - it.dy); }
        }
      }
      g.restore();
    }
  }
  drawLegend(g, ax, ents, l, r, t, b, xs, yf) {
    const o = ax.legendOpt, fs = o.size || 11;
    g.font = fs + "px " + PLOT_FONT;
    const ncol = Math.max(1, Math.min(o.ncol || 1, ents.length)), nrow = Math.ceil(ents.length / ncol);
    const sw = 26, pad = 6, lh = fs + 6;
    const colW = []; for (let c = 0; c < ncol; c++) { let m = 0; for (let k = c; k < ents.length; k += ncol) m = Math.max(m, g.measureText(ents[k].label).width); colW.push(sw + 6 + m + 10); }
    const W = colW.reduce((a, v) => a + v, 0) + pad, H = nrow * lh + pad * 2 - 2;
    const cands = { "upper right": [r - W - 6, t + 6], "upper left": [l + 6, t + 6], "lower right": [r - W - 6, b - H - 6],
      "lower left": [l + 6, b - H - 6], "upper center": [(l + r - W) / 2, t + 6], "center right": [r - W - 6, (t + b - H) / 2] };
    let pos = cands[o.loc];
    if (!pos) {
      let best = null, bestN = Infinity;
      for (const key of ["upper right", "upper left", "lower left", "lower right", "upper center", "center right"]) {
        const [x0, y0] = cands[key]; let cnt = 0;
        for (const it of ax.items) if (it.type === "line" || it.type === "fill") for (let i = 0; i < it.x.length; i++) {
          const v = it.y[i]; if (!finite(v) || (ax.ylog && v <= 0)) continue;
          const px = xs(it.x[i]), py = yf(v); if (px >= x0 && px <= x0 + W && py >= y0 && py <= y0 + H) cnt++;
        }
        if (cnt < bestN) { bestN = cnt; best = cands[key]; }
      }
      pos = best;
    }
    const [x0, y0] = pos;
    g.save(); g.globalAlpha = 0.93; g.fillStyle = "#fff"; g.fillRect(x0, y0, W, H); g.globalAlpha = 1;
    g.strokeStyle = "#ccc"; g.lineWidth = 1; g.strokeRect(x0 + 0.5, y0 + 0.5, W - 1, H - 1);
    g.textBaseline = "middle"; g.textAlign = "left";
    ents.forEach((it, k) => {
      const c = k % ncol, rr = Math.floor(k / ncol);
      const cx = x0 + pad + colW.slice(0, c).reduce((a, v) => a + v, 0), cy = y0 + pad + rr * lh + lh / 2 - 1;
      g.save(); g.globalAlpha = it.alpha != null ? Math.max(it.alpha, 0.35) : 1;
      if (it.type === "fill") { g.fillStyle = it.color; g.fillRect(cx, cy - 5, sw, 10); }
      else {
        g.strokeStyle = it.color; g.lineWidth = Math.min(it.lw || 1.5, 5); g.setLineDash(it.dash || []);
        if ((it.lw || 0) > 0) { g.beginPath(); g.moveTo(cx, cy); g.lineTo(cx + sw, cy); g.stroke(); }
        if (it.marker) { g.setLineDash([]); g.fillStyle = it.mfc || it.color; drawMarker(g, it.marker, cx + sw / 2, cy, it.ms || 5); }
      }
      g.restore();
      g.fillStyle = "#222"; g.fillText(it.label, cx + sw + 6, cy);
    });
    g.restore();
  }
}

function drawMarker(g, m, x, y, s) {
  g.beginPath();
  if (m === "s") g.rect(x - s / 2, y - s / 2, s, s);
  else g.arc(x, y, s / 2, 0, Math.PI * 2);
  g.fill();
}
function fmtShort(v) {
  if (!finite(v)) return "-";
  const a = Math.abs(v);
  if (a !== 0 && (a >= 1e5 || a < 1e-3)) return v.toExponential(3);
  return Number(v.toPrecision(5)).toString();
}
