/* Worker pool around the Rust/WebAssembly module. Every worker holds its own
   instance with xsdata.npz loaded; worker 0 also holds the DICOM images.
   Falls back to running on the page itself if workers are unavailable. */
"use strict";

function workerMain(self) {
  let X = null, curId = null;
  const enc = new TextEncoder(), dec = new TextDecoder();
  function put(bytes) { const p = X.alloc(bytes.length); new Uint8Array(X.memory.buffer, p, bytes.length).set(bytes); return p; }
  function call(obj, blob) {
    const r = enc.encode(JSON.stringify(obj));
    const rp = put(r);
    let bp = 0, bl = 0;
    if (blob) { bl = blob.length; bp = put(blob); }
    let n;
    try { n = X.call(rp, r.length, bp, bl); }
    catch (e) {
      let msg = "";
      try { msg = dec.decode(new Uint8Array(X.memory.buffer, X.panic_ptr(), X.panic_len())); } catch (e2) { /* ignore */ }
      const err = new Error(msg || ("WebAssembly error: " + e.message)); err.fatal = true; throw err;
    }
    const out = JSON.parse(dec.decode(new Uint8Array(X.memory.buffer, X.out_ptr(), n)));
    if (out.err !== undefined) throw new Error(out.err);
    return out.ok;
  }
  self.onmessage = async (ev) => {
    const m = ev.data;
    if (m.type === "init") {
      try {
        const res = await WebAssembly.instantiate(m.wasm, { env: { js_progress: (f) => self.postMessage({ id: curId, progress: f }) } });
        X = res.instance.exports;
        self.postMessage({ id: m.id, ok: true });
      } catch (e) { self.postMessage({ id: m.id, err: String(e.message || e) }); }
      return;
    }
    curId = m.id;
    try {
      if (m.type === "mc") {
        call(Object.assign({ cmd: "mc_start" }, m.req));
        let st, last = 0;
        do {
          st = call({ cmd: "mc_step", batches: 1 });
          const now = Date.now();
          if (now - last > 120 || st.finished) { self.postMessage({ id: m.id, progress: st.done / st.n_hist }); last = now; }
        } while (!st.finished);
        self.postMessage({ id: m.id, ok: call({ cmd: "mc_result" }) });
      } else {
        const r = call(m.req, m.blob ? new Uint8Array(m.blob) : null);
        let bin = null;
        if (r && r.bin) bin = new Uint8Array(X.memory.buffer, X.bin_ptr(), X.bin_len()).slice();
        self.postMessage({ id: m.id, ok: r, bin }, bin ? [bin.buffer] : []);
      }
    } catch (e) { self.postMessage({ id: m.id, err: String(e.message || e), fatal: !!e.fatal }); }
  };
}

class Engine {
  constructor(wasmBytes) {
    this.wasm = wasmBytes; this.nextId = 1; this.pending = new Map(); this.workers = []; this.inline = false;
    this.xsBytes = null; this.busy = [];
  }
  async start() {
    const hc = navigator.hardwareConcurrency || 2;
    const n = Math.max(1, Math.min(hc > 2 ? hc - 1 : hc, 12));
    let url = null;
    try {
      url = URL.createObjectURL(new Blob(["(" + workerMain.toString() + ")(self);"], { type: "text/javascript" }));
      const w0 = this.makeWorker(url);
      this.workers.push(w0);
      await this.send(0, { type: "init", wasm: this.wasm.slice(0) });
      for (let i = 1; i < n; i++) this.workers.push(this.makeWorker(url));
      await Promise.all(this.workers.slice(1).map((_, k) => this.send(k + 1, { type: "init", wasm: this.wasm.slice(0) })));
    } catch (e) {
      console.warn("Web Workers unavailable, running on the page:", e);
      this.workers = []; this.inline = true;
      this.workers.push(this.makeInline());
      await this.send(0, { type: "init", wasm: this.wasm.slice(0) });
    }
    this.url = url;
    this.busy = this.workers.map(() => 0);
    return this.workers.length;
  }
  makeWorker(url) {
    const w = new Worker(url);
    w.onmessage = (ev) => this.onMsg(ev.data);
    w.onerror = (e) => { console.error(e); };
    return w;
  }
  makeInline() {
    const self = { postMessage: (msg) => setTimeout(() => this.onMsg(msg), 0) };
    workerMain(self);
    return { postMessage: (m) => setTimeout(() => self.onmessage({ data: m }), 0), terminate() {} };
  }
  onMsg(m) {
    const p = this.pending.get(m.id);
    if (!p) return;
    if (m.progress !== undefined && m.ok === undefined && m.err === undefined) { if (p.onProgress) p.onProgress(m.progress); return; }
    this.pending.delete(m.id);
    this.busy[p.w] = Math.max(0, (this.busy[p.w] || 1) - 1);
    if (m.err !== undefined) { const e = new Error(m.err); e.fatal = m.fatal; if (m.fatal) this.restart(p.w); p.reject(e); }
    else p.resolve(m.bin ? { res: m.ok, bin: m.bin } : m.ok);
  }
  send(wi, msg, onProgress, transfer) {
    const id = this.nextId++; msg.id = id;
    return new Promise((resolve, reject) => {
      this.pending.set(id, { resolve, reject, onProgress, w: wi });
      this.busy[wi] = (this.busy[wi] || 0) + 1;
      this.workers[wi].postMessage(msg, transfer || []);
    });
  }
  async restart(wi) {
    try {
      this.workers[wi].terminate();
      this.workers[wi] = this.inline ? this.makeInline() : this.makeWorker(this.url);
      await this.send(wi, { type: "init", wasm: this.wasm.slice(0) });
      if (this.xsBytes) await this.send(wi, { type: "call", req: { cmd: "load_xs" }, blob: this.xsBytes.slice(0) });
      if (wi === 0 && this.onRestart0) this.onRestart0();
    } catch (e) { console.error("restart failed", e); }
  }
  /* run a command on worker 0 (the one holding images) */
  call(req, blob) { return this.send(0, { type: "call", req, blob }, null, blob ? [blob] : []); }
  callP(req, onProgress) { return this.send(0, { type: "call", req }, onProgress); }
  async loadXS(bytes) {
    this.xsBytes = bytes.slice(0);
    const res = await Promise.all(this.workers.map((_, i) => this.send(i, { type: "call", req: { cmd: "load_xs" }, blob: bytes.slice(0) })));
    return res[0];
  }
  /* run a Monte Carlo point on the least busy worker */
  mc(req, onProgress) {
    let wi = 0;
    for (let i = 0; i < this.workers.length; i++) if (this.busy[i] < this.busy[wi] || (this.busy[i] === this.busy[wi] && i > wi && this.busy[wi] > 0)) wi = i;
    if (this.workers.length > 1 && this.busy[wi] > 0) { // prefer a non-zero idle worker
      for (let i = this.workers.length - 1; i >= 0; i--) if (this.busy[i] === 0) { wi = i; break; }
    }
    return this.send(wi, { type: "mc", req }, onProgress);
  }
}

function b64ToBytes(b64) {
  const bin = atob(b64); const out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}
