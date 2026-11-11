//! Native driver for validation: reads JSON requests (one per line) on stdin
//! and writes JSON replies. A request may carry "file": path, whose bytes are
//! passed as the blob (for load_xs / dicom_load).
use std::io::{BufRead, Write};
use xray_app::api::{dispatch, State};

fn main() {
    let mut st = State::new(Box::new(|_| {}));
    let stdin = std::io::stdin();
    let mut out = std::io::stdout();
    for line in stdin.lock().lines() {
        let line = line.unwrap();
        if line.trim().is_empty() {
            continue;
        }
        let v: serde_json::Value = serde_json::from_str(&line).expect("bad json");
        let blob = v.get("file").and_then(|p| p.as_str()).map(|p| std::fs::read(p).expect("read file"));
        let t = std::time::Instant::now();
        let r = match dispatch(&mut st, &v, blob) {
            Ok(r) => serde_json::json!({"ok": r, "secs": t.elapsed().as_secs_f64()}),
            Err(e) => serde_json::json!({"err": e}),
        };
        writeln!(out, "{}", r).unwrap();
    }
}
