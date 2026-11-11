//! Minimal reader for NumPy `.npz` archives (a zip of `.npy` files).
//! Supports stored and deflated members and the numeric dtypes
//! `<f8 <f4 <i8 <i4 <u4 <i2 <u2 |u1 |i1`, which is everything xsdata.npz uses.

use std::collections::BTreeMap;

fn u16le(b: &[u8], o: usize) -> Result<u16, String> {
    b.get(o..o + 2).map(|s| u16::from_le_bytes([s[0], s[1]])).ok_or_else(|| "npz: truncated".to_string())
}
fn u32le(b: &[u8], o: usize) -> Result<u32, String> {
    b.get(o..o + 4).map(|s| u32::from_le_bytes([s[0], s[1], s[2], s[3]])).ok_or_else(|| "npz: truncated".to_string())
}

/// Read every member of a zip archive into memory: name -> bytes.
pub fn read_zip(buf: &[u8]) -> Result<BTreeMap<String, Vec<u8>>, String> {
    // locate end-of-central-directory record
    if buf.len() < 22 {
        return Err("Not a zip/npz file (too short).".into());
    }
    let mut eocd = None;
    let start = buf.len().saturating_sub(22 + 65535);
    let mut i = buf.len() - 22;
    loop {
        if &buf[i..i + 4] == b"PK\x05\x06" {
            eocd = Some(i);
            break;
        }
        if i == start {
            break;
        }
        i -= 1;
    }
    let eocd = eocd.ok_or("Not a zip/npz file (no end-of-central-directory).")?;
    let mut n_entries = u16le(buf, eocd + 10)? as u64;
    let mut cd_off = u32le(buf, eocd + 16)? as u64;
    // zip64?
    if n_entries == 0xFFFF || cd_off == 0xFFFF_FFFF {
        if eocd >= 20 && &buf[eocd - 20..eocd - 16] == b"PK\x06\x07" {
            let z64 = u64::from_le_bytes(buf[eocd - 12..eocd - 4].try_into().unwrap()) as usize;
            if buf.get(z64..z64 + 4) == Some(b"PK\x06\x06") {
                n_entries = u64::from_le_bytes(buf[z64 + 32..z64 + 40].try_into().unwrap());
                cd_off = u64::from_le_bytes(buf[z64 + 48..z64 + 56].try_into().unwrap());
            }
        }
    }
    let mut out = BTreeMap::new();
    let mut p = cd_off as usize;
    for _ in 0..n_entries {
        if buf.get(p..p + 4) != Some(b"PK\x01\x02") {
            return Err("npz: bad central directory".into());
        }
        let method = u16le(buf, p + 10)?;
        let mut csize = u32le(buf, p + 20)? as u64;
        let mut usize_ = u32le(buf, p + 24)? as u64;
        let nlen = u16le(buf, p + 28)? as usize;
        let xlen = u16le(buf, p + 30)? as usize;
        let clen = u16le(buf, p + 32)? as usize;
        let mut loff = u32le(buf, p + 42)? as u64;
        let name = String::from_utf8_lossy(&buf[p + 46..p + 46 + nlen]).to_string();
        // zip64 extra field
        let mut q = p + 46 + nlen;
        let xend = q + xlen;
        while q + 4 <= xend {
            let id = u16le(buf, q)?;
            let sz = u16le(buf, q + 2)? as usize;
            if id == 1 {
                let mut r = q + 4;
                if usize_ == 0xFFFF_FFFF {
                    usize_ = u64::from_le_bytes(buf[r..r + 8].try_into().unwrap());
                    r += 8;
                }
                if csize == 0xFFFF_FFFF {
                    csize = u64::from_le_bytes(buf[r..r + 8].try_into().unwrap());
                    r += 8;
                }
                if loff == 0xFFFF_FFFF {
                    loff = u64::from_le_bytes(buf[r..r + 8].try_into().unwrap());
                }
            }
            q += 4 + sz;
        }
        p = xend + clen;

        let lo = loff as usize;
        if buf.get(lo..lo + 4) != Some(b"PK\x03\x04") {
            return Err("npz: bad local header".into());
        }
        let ln = u16le(buf, lo + 26)? as usize;
        let lx = u16le(buf, lo + 28)? as usize;
        let d0 = lo + 30 + ln + lx;
        let data = buf.get(d0..d0 + csize as usize).ok_or("npz: truncated member")?;
        let bytes = match method {
            0 => data.to_vec(),
            8 => miniz_oxide::inflate::decompress_to_vec(data).map_err(|e| format!("npz: inflate failed for {}: {:?}", name, e))?,
            m => return Err(format!("npz: unsupported compression method {} in {}", m, name)),
        };
        let _ = usize_;
        out.insert(name, bytes);
    }
    Ok(out)
}

/// Parse a `.npy` blob into a flat f64 vector (C order) and its shape.
pub fn parse_npy(b: &[u8]) -> Result<(Vec<f64>, Vec<usize>), String> {
    if b.len() < 10 || &b[0..6] != b"\x93NUMPY" {
        return Err("not an .npy array".into());
    }
    let major = b[6];
    let (hlen, hstart) = if major == 1 {
        (u16le(b, 8)? as usize, 10)
    } else {
        (u32le(b, 8)? as usize, 12)
    };
    let header = std::str::from_utf8(&b[hstart..hstart + hlen]).map_err(|_| "npy: bad header")?;
    let descr = extract(header, "descr").ok_or("npy: no descr")?;
    let descr = descr.trim().trim_matches(|c| c == '\'' || c == '"').to_string();
    let fortran = extract(header, "fortran_order").map(|s| s.trim().starts_with("True")).unwrap_or(false);
    let shape_s = extract(header, "shape").ok_or("npy: no shape")?;
    let shape: Vec<usize> = shape_s
        .trim()
        .trim_start_matches('(')
        .trim_end_matches(')')
        .split(',')
        .filter_map(|t| t.trim().parse().ok())
        .collect();
    let n: usize = shape.iter().product();
    let data = &b[hstart + hlen..];
    let (sz, conv): (usize, fn(&[u8]) -> f64) = match descr.as_str() {
        "<f8" => (8, |s| f64::from_le_bytes(s.try_into().unwrap())),
        "<f4" => (4, |s| f32::from_le_bytes(s.try_into().unwrap()) as f64),
        "<i8" => (8, |s| i64::from_le_bytes(s.try_into().unwrap()) as f64),
        "<i4" => (4, |s| i32::from_le_bytes(s.try_into().unwrap()) as f64),
        "<u4" => (4, |s| u32::from_le_bytes(s.try_into().unwrap()) as f64),
        "<i2" => (2, |s| i16::from_le_bytes(s.try_into().unwrap()) as f64),
        "<u2" => (2, |s| u16::from_le_bytes(s.try_into().unwrap()) as f64),
        "|u1" => (1, |s| s[0] as f64),
        "|i1" => (1, |s| s[0] as i8 as f64),
        ">f8" => (8, |s| f64::from_be_bytes(s.try_into().unwrap())),
        ">f4" => (4, |s| f32::from_be_bytes(s.try_into().unwrap()) as f64),
        d => return Err(format!("npy: unsupported dtype {}", d)),
    };
    if data.len() < n * sz {
        return Err("npy: truncated data".into());
    }
    let mut v: Vec<f64> = (0..n).map(|i| conv(&data[i * sz..i * sz + sz])).collect();
    if fortran && shape.len() == 2 {
        let (r, c) = (shape[0], shape[1]);
        let mut t = vec![0.0; n];
        for i in 0..r {
            for j in 0..c {
                t[i * c + j] = v[j * r + i];
            }
        }
        v = t;
    }
    Ok((v, shape))
}

fn extract<'a>(header: &'a str, key: &str) -> Option<&'a str> {
    let k1 = format!("'{}':", key);
    let k2 = format!("\"{}\":", key);
    let pos = header.find(&k1).map(|p| p + k1.len()).or_else(|| header.find(&k2).map(|p| p + k2.len()))?;
    let rest = &header[pos..];
    let rest_trim = rest.trim_start();
    if rest_trim.starts_with('(') {
        let end = rest_trim.find(')')?;
        Some(&rest_trim[..=end])
    } else {
        let end = rest.find(',').unwrap_or(rest.len());
        Some(&rest[..end])
    }
}

/// Load an npz into name (without `.npy`) -> flat f64 data.
pub fn load_npz(buf: &[u8]) -> Result<BTreeMap<String, Vec<f64>>, String> {
    let z = read_zip(buf)?;
    let mut out = BTreeMap::new();
    for (name, bytes) in z {
        let key = name.strip_suffix(".npy").unwrap_or(&name).to_string();
        let (v, _shape) = parse_npy(&bytes).map_err(|e| format!("{}: {}", key, e))?;
        out.insert(key, v);
    }
    Ok(out)
}
