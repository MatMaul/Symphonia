use super::entdec::EcDec;
use super::intrin::ec_ilog;

pub type CeltVal = f32;

include!("generated/cwrs_tables.rs");

fn pvq_u_row(n: i32, k: i32) -> u32 {
    let n = n as usize;
    let k = k as usize;
    let base = CELT_PVQ_U_ROW_OFFSETS[n] + k;
    debug_assert!(base < CELT_PVQ_U_DATA.len());
    CELT_PVQ_U_DATA[base]
}

fn pvq_u(n: i32, k: i32) -> u32 {
    if n < k {
        pvq_u_row(n, k)
    } else {
        pvq_u_row(k, n)
    }
}

fn pvq_v(n: i32, k: i32) -> u32 {
    pvq_u(n, k) + pvq_u(n, k + 1)
}

pub fn log2_frac(mut val: u32, mut frac: i32) -> i32 {
    let mut l = ec_ilog(val);
    if val & (val - 1) != 0 {
        if l > 16 {
            val = ((val - 1) >> (l - 16)) + 1;
        } else {
            val <<= 16 - l;
        }
        l = (l - 1) << frac;
        loop {
            let b = (val >> 16) as i32;
            l += b << frac;
            val = (val + b as u32) >> (b as u32);
            val = (val * val + 0x7fff) >> 15;
            if frac == 0 {
                break;
            }
            frac -= 1;
        }
        if val > 0x8000 { l + 1 } else { l }
    } else {
        (l - 1) << frac
    }
}

pub fn get_required_bits(bits: &mut [i16], n: i32, maxk: i32, frac: i32) {
    debug_assert!(maxk > 0);
    bits[0] = 0;
    for k in 1..=maxk {
        let v = pvq_v(n, k);
        bits[k as usize] = log2_frac(v, frac) as i16;
    }
}

fn icwrs(n: i32, y: &[i32]) -> u32 {
    debug_assert!(n >= 2);
    let mut j = (n - 1) as usize;
    let mut i = if y[j] < 0 { 1 } else { 0 };
    let mut k = y[j].abs();
    while j > 0 {
        j -= 1;
        i += pvq_u((n as usize - j) as i32, k);
        k += y[j].abs();
        if y[j] < 0 {
            i += pvq_u((n as usize - j) as i32, k + 1);
        }
    }
    i
}

fn cwrsi(n: i32, mut k: i32, mut i: u32, y: &mut [i32]) -> CeltVal {
    let mut n = n;
    let mut yy: CeltVal = 0.0;
    debug_assert!(k > 0);
    debug_assert!(n > 1);
    let mut idx = 0usize;
    while n > 2 {
        if k >= n {
            let p = pvq_u_row(n, k + 1);
            let mut s = if i >= p { -1 } else { 0 };
            if s < 0 { i -= p; }
            let k0 = k;
            let q = pvq_u_row(n, n);
            if q > i {
                k = n;
                loop {
                    k -= 1;
                    if pvq_u_row(k, n) <= i {
                        break;
                    }
                }
            } else {
                while pvq_u_row(n, k) > i {
                    k -= 1;
                }
            }
            i -= pvq_u_row(n, k);
            let val = (k0 - k + s) ^ s;
            y[idx] = val;
            yy += (val * val) as CeltVal;
        } else {
            let p = pvq_u_row(k, n);
            let q = pvq_u_row(k + 1, n);
            if p <= i && i < q {
                i -= p;
                y[idx] = 0;
            } else {
                let mut s = if i >= q { -1 } else { 0 };
                if s < 0 { i -= q; }
                let k0 = k;
                loop {
                    k -= 1;
                    if pvq_u_row(k, n) <= i {
                        break;
                    }
                }
                i -= pvq_u_row(k, n);
                let val = (k0 - k + s) ^ s;
                y[idx] = val;
                yy += (val * val) as CeltVal;
            }
        }
        n -= 1;
        idx += 1;
    }
    let p = 2 * k + 1;
    let mut s = if i >= p as u32 { -1 } else { 0 };
    if s < 0 { i -= p as u32; }
    let k0 = k;
    k = ((i + 1) >> 1) as i32;
    if k > 0 {
        i -= (2 * k - 1) as u32;
    }
    let val = (k0 - k + s) ^ s;
    y[idx] = val;
    yy += (val * val) as CeltVal;
    idx += 1;
    s = -(i as i32);
    let val = (k + s) ^ s;
    y[idx] = val;
    yy += (val * val) as CeltVal;
    yy
}

pub fn decode_pulses(y: &mut [i32], n: i32, k: i32, dec: &mut EcDec<'_>) -> CeltVal {
    debug_assert!(k > 0);
    let index = dec.dec_uint(pvq_v(n, k));
    cwrsi(n, k, index, y)
}
