//! A self-contained DICOM reader for the image-quality tab (replaces pydicom).
//!
//! Handles Part-10 files (and preamble-less raw datasets), implicit/explicit
//! VR little endian, explicit VR big endian, deflated explicit VR little
//! endian, RLE Lossless, and encapsulated JPEG (baseline and lossless
//! process 14) via the pure-Rust jpeg-decoder crate. JPEG-LS and JPEG 2000
//! are reported as unsupported rather than mis-read.
//!
//! As with pydicom's `pixel_array`, values are returned with Rescale
//! Slope/Intercept applied, as float32. Only top-level attributes are used
//! (attributes nested in sequences are skipped, as `getattr(ds, ...)` does).

use crate::iq::Image;
use std::collections::HashMap;

pub struct Dicom {
    pub image: Image,
    pub px_mm: f64,
    pub sign: i32,
    pub photometric: String,
    pub transfer_syntax: String,
    pub frames: usize,
    pub info: Vec<(String, String)>,
}

#[derive(Clone, Copy, PartialEq)]
enum Syntax {
    ImplicitLE,
    ExplicitLE,
    ExplicitBE,
}

struct Reader<'a> {
    b: &'a [u8],
    p: usize,
    syn: Syntax,
}

const UNDEF: u32 = 0xFFFF_FFFF;

impl<'a> Reader<'a> {
    fn u16(&mut self) -> Result<u16, String> {
        let s = self.b.get(self.p..self.p + 2).ok_or("DICOM: unexpected end of file")?;
        self.p += 2;
        Ok(if self.syn == Syntax::ExplicitBE { u16::from_be_bytes([s[0], s[1]]) } else { u16::from_le_bytes([s[0], s[1]]) })
    }
    fn u32(&mut self) -> Result<u32, String> {
        let s = self.b.get(self.p..self.p + 4).ok_or("DICOM: unexpected end of file")?;
        self.p += 4;
        let a = [s[0], s[1], s[2], s[3]];
        Ok(if self.syn == Syntax::ExplicitBE { u32::from_be_bytes(a) } else { u32::from_le_bytes(a) })
    }
    fn eof(&self) -> bool {
        self.p + 8 > self.b.len()
    }

    /// Read one element header: (group, element, vr, length).
    fn header(&mut self) -> Result<(u16, u16, [u8; 2], u32), String> {
        // group 0002 is always explicit LE
        let save = self.syn;
        let g_peek = u16::from_le_bytes([self.b[self.p], self.b[self.p + 1]]);
        if g_peek == 0x0002 {
            self.syn = Syntax::ExplicitLE;
        }
        let g = self.u16()?;
        let e = self.u16()?;
        let r;
        if g == 0xFFFE {
            // item / delimiters: never have a VR
            let len = self.u32()?;
            self.syn = save;
            return Ok((g, e, *b"  ", len));
        }
        if self.syn == Syntax::ImplicitLE {
            let len = self.u32()?;
            r = (g, e, *b"  ", len);
        } else {
            let vr = [self.b[self.p], self.b[self.p + 1]];
            self.p += 2;
            let long = matches!(&vr, b"OB" | b"OW" | b"OF" | b"SQ" | b"UT" | b"UN" | b"OD" | b"OL" | b"OV" | b"UC" | b"UR" | b"SV" | b"UV");
            let len = if long {
                self.p += 2;
                self.u32()?
            } else {
                self.u16()? as u32
            };
            r = (g, e, vr, len);
        }
        self.syn = save;
        Ok(r)
    }

    /// Skip a sequence / item body with undefined length.
    fn skip_undefined(&mut self) -> Result<(), String> {
        loop {
            if self.eof() {
                return Ok(());
            }
            let (g, e, vr, len) = self.header()?;
            if g == 0xFFFE && e == 0xE0DD {
                return Ok(()); // sequence delimiter
            }
            if g == 0xFFFE && e == 0xE00D {
                return Ok(()); // item delimiter
            }
            if g == 0xFFFE && e == 0xE000 {
                if len == UNDEF {
                    self.skip_undefined()?;
                } else {
                    self.p += len as usize;
                }
                continue;
            }
            if len == UNDEF {
                let _ = vr;
                self.skip_undefined()?;
            } else {
                self.p += len as usize;
            }
        }
    }
}

fn ts_name(uid: &str) -> &'static str {
    match uid {
        "1.2.840.10008.1.2" => "Implicit VR Little Endian",
        "1.2.840.10008.1.2.1" => "Explicit VR Little Endian",
        "1.2.840.10008.1.2.1.99" => "Deflated Explicit VR Little Endian",
        "1.2.840.10008.1.2.2" => "Explicit VR Big Endian",
        "1.2.840.10008.1.2.4.50" => "JPEG Baseline",
        "1.2.840.10008.1.2.4.51" => "JPEG Extended",
        "1.2.840.10008.1.2.4.57" => "JPEG Lossless",
        "1.2.840.10008.1.2.4.70" => "JPEG Lossless SV1",
        "1.2.840.10008.1.2.4.80" => "JPEG-LS Lossless",
        "1.2.840.10008.1.2.4.81" => "JPEG-LS Near-lossless",
        "1.2.840.10008.1.2.4.90" => "JPEG 2000 Lossless",
        "1.2.840.10008.1.2.4.91" => "JPEG 2000",
        "1.2.840.10008.1.2.5" => "RLE Lossless",
        _ => "unknown",
    }
}

fn clean(s: &[u8]) -> String {
    String::from_utf8_lossy(s).trim_matches(|c: char| c == '\0' || c.is_whitespace()).to_string()
}

pub fn read_dicom(buf: &[u8]) -> Result<Dicom, String> {
    let mut start = 0usize;
    let has_preamble = buf.len() > 132 && &buf[128..132] == b"DICM";
    if has_preamble {
        start = 132;
    }
    // --- meta group (explicit LE) ---
    let mut ts = String::from("1.2.840.10008.1.2");
    let mut r = Reader { b: buf, p: start, syn: Syntax::ExplicitLE };
    if has_preamble {
        while r.p + 8 <= buf.len() {
            let g = u16::from_le_bytes([buf[r.p], buf[r.p + 1]]);
            if g != 0x0002 {
                break;
            }
            let (_g, e, _vr, len) = r.header()?;
            let v = buf.get(r.p..r.p + len as usize).ok_or("DICOM: truncated meta header")?;
            if e == 0x0010 {
                ts = clean(v);
            }
            r.p += len as usize;
        }
    } else {
        // no preamble: guess explicit vs implicit from the first element
        if buf.len() > 6 && buf[4].is_ascii_uppercase() && buf[5].is_ascii_uppercase() {
            ts = "1.2.840.10008.1.2.1".into();
        }
    }
    let ds_start = r.p;
    let inflated;
    let (body, syn): (&[u8], Syntax) = match ts.as_str() {
        "1.2.840.10008.1.2" => (&buf[ds_start..], Syntax::ImplicitLE),
        "1.2.840.10008.1.2.2" => (&buf[ds_start..], Syntax::ExplicitBE),
        "1.2.840.10008.1.2.1.99" => {
            inflated = miniz_oxide::inflate::decompress_to_vec(&buf[ds_start..]).map_err(|e| format!("DICOM: deflate failed: {:?}", e))?;
            (&inflated[..], Syntax::ExplicitLE)
        }
        _ => (&buf[ds_start..], Syntax::ExplicitLE),
    };
    let encapsulated = !matches!(ts.as_str(), "1.2.840.10008.1.2" | "1.2.840.10008.1.2.1" | "1.2.840.10008.1.2.2" | "1.2.840.10008.1.2.1.99");

    // --- dataset: collect top-level attributes we need ---
    let mut rd = Reader { b: body, p: 0, syn };
    let mut vals: HashMap<(u16, u16), Vec<u8>> = HashMap::new();
    let mut pixel: Option<(usize, u32)> = None; // offset, length (UNDEF = encapsulated)
    while !rd.eof() {
        let (g, e, vr, len) = rd.header()?;
        if g == 0x7FE0 && e == 0x0010 {
            pixel = Some((rd.p, len));
            break;
        }
        if len == UNDEF {
            rd.skip_undefined()?; // sequences (and anything else undefined)
            continue;
        }
        if &vr == b"SQ" {
            rd.p += len as usize;
            continue;
        }
        let end = (rd.p + len as usize).min(body.len());
        if matches!(g, 0x0008 | 0x0010 | 0x0018 | 0x0020 | 0x0028) {
            vals.insert((g, e), body[rd.p..end].to_vec());
        }
        rd.p = end;
    }
    let be = syn == Syntax::ExplicitBE;
    let us = |g: u16, e: u16| -> Option<u32> {
        let v = vals.get(&(g, e))?;
        if v.len() >= 2 {
            Some(if be { u16::from_be_bytes([v[0], v[1]]) } else { u16::from_le_bytes([v[0], v[1]]) } as u32)
        } else {
            None
        }
    };
    let ss = |g: u16, e: u16| -> Option<i32> { us(g, e).map(|v| v as u16 as i16 as i32) };
    let text = |g: u16, e: u16| -> Option<String> { vals.get(&(g, e)).map(|v| clean(v)) };
    let ds = |g: u16, e: u16| -> Option<f64> { text(g, e).and_then(|s| s.split('\\').next().and_then(|t| t.trim().parse().ok())) };

    let rows = us(0x0028, 0x0010).ok_or("DICOM: no Rows attribute - is this an image?")? as usize;
    let cols = us(0x0028, 0x0011).ok_or("DICOM: no Columns attribute - is this an image?")? as usize;
    let bits_alloc = us(0x0028, 0x0100).unwrap_or(16) as usize;
    let bits_stored = us(0x0028, 0x0101).unwrap_or(bits_alloc as u32) as usize;
    let pix_rep = us(0x0028, 0x0103).unwrap_or(0);
    let spp = us(0x0028, 0x0002).unwrap_or(1);
    let photometric = text(0x0028, 0x0004).unwrap_or_else(|| "MONOCHROME2".into());
    let frames: usize = text(0x0028, 0x0008).and_then(|s| s.parse().ok()).unwrap_or(1).max(1);
    let slope = ds(0x0028, 0x1053).unwrap_or(1.0);
    let intercept = ds(0x0028, 0x1052).unwrap_or(0.0);
    let sign = ss(0x0028, 0x1041).unwrap_or(1);
    if spp != 1 {
        return Err(format!("DICOM: {} samples per pixel (colour) - a greyscale detector image is needed.", spp));
    }
    let px_mm = ds(0x0028, 0x0030)
        .or_else(|| ds(0x0018, 0x1164))
        .ok_or("No PixelSpacing or ImagerPixelSpacing tag found - cannot convert pixels to mm.")?;
    let (poff, plen) = pixel.ok_or("DICOM: no Pixel Data element found.")?;
    let npx = rows * cols;

    // --- decode the first frame into raw integer samples ---
    let raw: Vec<f64> = if !encapsulated && plen != UNDEF {
        let data = &body[poff..(poff + plen as usize).min(body.len())];
        decode_native(data, npx, bits_alloc, bits_stored, pix_rep, be)?
    } else {
        let frags = fragments(body, poff)?;
        let fr0 = first_frame(&frags, frames);
        match ts.as_str() {
            "1.2.840.10008.1.2.5" => decode_rle(&fr0, npx, bits_alloc, bits_stored, pix_rep)?,
            "1.2.840.10008.1.2.4.50" | "1.2.840.10008.1.2.4.51" | "1.2.840.10008.1.2.4.57" | "1.2.840.10008.1.2.4.70" => {
                decode_jpeg(&fr0, npx, bits_stored, pix_rep)?
            }
            other => {
                return Err(format!(
                    "This DICOM uses the {} transfer syntax ({}), which the browser build cannot decode. \
                     Re-export the image uncompressed (Explicit VR Little Endian) from the modality or PACS.",
                    ts_name(other),
                    other
                ))
            }
        }
    };
    let data: Vec<f32> = if slope == 1.0 && intercept == 0.0 {
        raw.iter().map(|&v| v as f32).collect()
    } else {
        raw.iter().map(|&v| (v as f32) * (slope as f32) + (intercept as f32)).collect()
    };
    let mut info = vec![
        ("Transfer syntax".to_string(), format!("{} ({})", ts_name(&ts), ts)),
        ("Size".to_string(), format!("{} x {}", cols, rows)),
        ("Bits stored".to_string(), format!("{} of {}", bits_stored, bits_alloc)),
        ("Pixel spacing".to_string(), format!("{} mm", px_mm)),
        ("Photometric".to_string(), photometric.clone()),
    ];
    if frames > 1 {
        info.push(("Frames".to_string(), format!("{} (first frame used)", frames)));
    }
    if let Some(m) = text(0x0008, 0x0070) {
        info.push(("Manufacturer".to_string(), m));
    }
    if let Some(m) = text(0x0008, 0x1090) {
        info.push(("Model".to_string(), m));
    }
    Ok(Dicom { image: Image { w: cols, h: rows, data }, px_mm, sign, photometric, transfer_syntax: ts, frames, info })
}

fn sign_mask(v: u32, bits_stored: usize, pix_rep: u32) -> f64 {
    let bs = bits_stored.min(32);
    if bs == 0 || bs >= 32 {
        return if pix_rep == 1 { v as i32 as f64 } else { v as f64 };
    }
    let m = (1u32 << bs) - 1;
    let x = v & m;
    if pix_rep == 1 && (x >> (bs - 1)) & 1 == 1 {
        (x as i64 - (1i64 << bs)) as f64
    } else {
        x as f64
    }
}

fn decode_native(d: &[u8], npx: usize, ba: usize, bs: usize, pr: u32, be: bool) -> Result<Vec<f64>, String> {
    let bytes = ba / 8;
    if ba == 1 || bytes == 0 {
        return Err("DICOM: 1-bit images are not supported.".into());
    }
    if d.len() < npx * bytes {
        return Err("DICOM: pixel data is shorter than Rows x Columns.".into());
    }
    let mut out = Vec::with_capacity(npx);
    for i in 0..npx {
        let s = &d[i * bytes..i * bytes + bytes];
        let v: u32 = match (bytes, be) {
            (1, _) => s[0] as u32,
            (2, false) => u16::from_le_bytes([s[0], s[1]]) as u32,
            (2, true) => u16::from_be_bytes([s[0], s[1]]) as u32,
            (4, false) => u32::from_le_bytes([s[0], s[1], s[2], s[3]]),
            (4, true) => u32::from_be_bytes([s[0], s[1], s[2], s[3]]),
            _ => return Err(format!("DICOM: {} bits allocated is not supported.", ba)),
        };
        out.push(sign_mask(v, bs, pr));
    }
    Ok(out)
}

/// Encapsulated pixel data fragments (after the basic offset table).
fn fragments(b: &[u8], mut p: usize) -> Result<Vec<Vec<u8>>, String> {
    let mut out = vec![];
    let mut first = true;
    while p + 8 <= b.len() {
        let g = u16::from_le_bytes([b[p], b[p + 1]]);
        let e = u16::from_le_bytes([b[p + 2], b[p + 3]]);
        let len = u32::from_le_bytes([b[p + 4], b[p + 5], b[p + 6], b[p + 7]]) as usize;
        p += 8;
        if g != 0xFFFE {
            break;
        }
        if e == 0xE0DD {
            break;
        }
        let end = (p + len).min(b.len());
        if !first {
            out.push(b[p..end].to_vec());
        }
        first = false;
        p = end;
    }
    if out.is_empty() {
        return Err("DICOM: encapsulated pixel data has no fragments.".into());
    }
    Ok(out)
}

fn first_frame(frags: &[Vec<u8>], frames: usize) -> Vec<u8> {
    if frames <= 1 {
        frags.concat()
    } else {
        frags[0].clone()
    }
}

fn decode_rle(d: &[u8], npx: usize, ba: usize, bs: usize, pr: u32) -> Result<Vec<f64>, String> {
    if d.len() < 64 {
        return Err("DICOM: bad RLE header.".into());
    }
    let nseg = u32::from_le_bytes([d[0], d[1], d[2], d[3]]) as usize;
    let bytes = ba / 8;
    if nseg != bytes {
        return Err("DICOM: RLE segment count does not match bits allocated (colour RLE is not supported).".into());
    }
    let offs: Vec<usize> = (0..nseg).map(|i| u32::from_le_bytes(d[4 + 4 * i..8 + 4 * i].try_into().unwrap()) as usize).collect();
    let mut segs = vec![];
    for i in 0..nseg {
        let a = offs[i];
        let b = if i + 1 < nseg { offs[i + 1] } else { d.len() };
        let src = &d[a..b.min(d.len())];
        let mut out = Vec::with_capacity(npx);
        let mut k = 0;
        while k < src.len() && out.len() < npx {
            let n = src[k] as i8;
            k += 1;
            if n >= 0 {
                let c = n as usize + 1;
                out.extend_from_slice(&src[k..(k + c).min(src.len())]);
                k += c;
            } else if n != -128 {
                let c = (-(n as i32)) as usize + 1;
                if k < src.len() {
                    out.extend(std::iter::repeat(src[k]).take(c));
                }
                k += 1;
            }
        }
        out.resize(npx, 0);
        segs.push(out);
    }
    // segments are most-significant byte first
    Ok((0..npx)
        .map(|i| {
            let mut v = 0u32;
            for s in &segs {
                v = (v << 8) | s[i] as u32;
            }
            sign_mask(v, bs, pr)
        })
        .collect())
}

fn decode_jpeg(d: &[u8], npx: usize, bs: usize, pr: u32) -> Result<Vec<f64>, String> {
    let mut dec = jpeg_decoder::Decoder::new(d);
    let px = dec.decode().map_err(|e| format!("DICOM: JPEG decode failed: {}", e))?;
    let info = dec.info().ok_or("DICOM: JPEG has no header")?;
    let out: Vec<f64> = match info.pixel_format {
        jpeg_decoder::PixelFormat::L8 => px.iter().map(|&v| sign_mask(v as u32, bs.min(8), pr)).collect(),
        jpeg_decoder::PixelFormat::L16 => px.chunks_exact(2).map(|c| sign_mask(u16::from_ne_bytes([c[0], c[1]]) as u32, bs, pr)).collect(),
        _ => return Err("DICOM: colour JPEG is not supported.".into()),
    };
    if out.len() < npx {
        return Err("DICOM: decoded JPEG is smaller than Rows x Columns.".into());
    }
    Ok(out[..npx].to_vec())
}
