//! X-ray Spectrum and Monte Carlo Transmission - Rust port of app.py.
//!
//! The physics and image analysis live here; `web/` holds the browser front
//! end that drives this crate compiled to WebAssembly.

pub mod api;
pub mod archer;
pub mod dicom;
pub mod iq;
pub mod mc;
pub mod npz;
pub mod num;
pub mod spectrum;
pub mod xs;

#[cfg(target_arch = "wasm32")]
mod wasm {
    use crate::api::{dispatch, State};
    use serde_json::{json, Value};

    #[link(wasm_import_module = "env")]
    extern "C" {
        fn js_progress(f: f64);
    }

    static mut STATE: Option<State> = None;
    static mut OUT: Vec<u8> = Vec::new();
    static mut PANIC_MSG: Vec<u8> = Vec::new();

    #[allow(static_mut_refs)]
    fn state() -> &'static mut State {
        unsafe {
            if STATE.is_none() {
                std::panic::set_hook(Box::new(|info| {
                    PANIC_MSG = format!("internal error: {}", info).into_bytes();
                }));
                STATE = Some(State::new(Box::new(|f| js_progress(f))));
            }
            STATE.as_mut().unwrap()
        }
    }

    #[no_mangle]
    pub extern "C" fn alloc(n: usize) -> *mut u8 {
        let mut v: Vec<u8> = Vec::with_capacity(n.max(1));
        let p = v.as_mut_ptr();
        std::mem::forget(v);
        p
    }

    #[no_mangle]
    pub unsafe extern "C" fn dealloc(p: *mut u8, n: usize) {
        drop(Vec::from_raw_parts(p, 0, n.max(1)));
    }

    /// Run a JSON command. `req` is consumed; `blob` (may be null) is an
    /// optional byte buffer allocated with `alloc` and consumed too.
    /// Returns the length of the JSON reply, readable at `out_ptr()`.
    #[no_mangle]
    #[allow(static_mut_refs)]
    pub unsafe extern "C" fn call(req: *mut u8, req_len: usize, blob: *mut u8, blob_len: usize) -> usize {
        let st = state();
        let rq = Vec::from_raw_parts(req, req_len, req_len.max(1));
        let bl = if blob.is_null() { None } else { Some(Vec::from_raw_parts(blob, blob_len, blob_len.max(1))) };
        let reply = match serde_json::from_slice::<Value>(&rq) {
            Ok(v) => match dispatch(st, &v, bl) {
                Ok(r) => json!({"ok": r}),
                Err(e) => json!({"err": e}),
            },
            Err(e) => json!({"err": format!("bad request: {}", e)}),
        };
        OUT = serde_json::to_vec(&reply).unwrap_or_default();
        OUT.len()
    }

    #[no_mangle]
    #[allow(static_mut_refs)]
    pub extern "C" fn out_ptr() -> *const u8 {
        unsafe { OUT.as_ptr() }
    }

    #[no_mangle]
    pub extern "C" fn bin_ptr() -> *const u8 {
        state().bin.as_ptr()
    }

    #[no_mangle]
    pub extern "C" fn bin_len() -> usize {
        state().bin.len()
    }

    #[no_mangle]
    #[allow(static_mut_refs)]
    pub extern "C" fn panic_ptr() -> *const u8 {
        unsafe { PANIC_MSG.as_ptr() }
    }

    #[no_mangle]
    #[allow(static_mut_refs)]
    pub extern "C" fn panic_len() -> usize {
        unsafe { PANIC_MSG.len() }
    }
}
