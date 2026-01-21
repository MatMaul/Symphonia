use super::entcode::BITRES;
use super::intrin::{celt_sudiv, celt_udiv, ec_ilog};
use super::math::celt_exp2_db;
use super::math::{celt_rsqrt, celt_sqrt};
use super::{renormalise_vector, CeltGlog, CeltMode, CeltNorm, CeltSig, E_MEANS};
use alloc::vec;

const ORDERY_TABLE: [i32; 30] = [
    1, 0, 3, 0, 2, 1, 7, 0, 4, 3, 6, 1, 5, 2, 15, 0, 8, 7, 12, 3, 11, 4, 14, 1, 9, 6,
    13, 2, 10, 5,
];

#[inline]
fn frac_mul16(a: i32, b: i32) -> i32 {
    (16384 + a * b) >> 15
}

pub fn hysteresis_decision(val: f32, thresholds: &[f32], hysteresis: &[f32], prev: i32) -> i32 {
    let mut i = 0i32;
    while (i as usize) < thresholds.len() {
        if val < thresholds[i as usize] {
            break;
        }
        i += 1;
    }
    if i > prev && val < thresholds[prev as usize] + hysteresis[prev as usize] {
        i = prev;
    }
    if i < prev && prev > 0 && val > thresholds[(prev - 1) as usize] - hysteresis[(prev - 1) as usize] {
        i = prev;
    }
    i
}

pub fn celt_lcg_rand(seed: u32) -> u32 {
    seed.wrapping_mul(1664525).wrapping_add(1013904223)
}

pub fn bitexact_cos(x: i16) -> i16 {
    let tmp = (4096 + (x as i32 * x as i32)) >> 13;
    debug_assert!(tmp <= 32767);
    let mut x2 = tmp;
    x2 = (32767 - x2)
        + frac_mul16(x2, (-7651 + frac_mul16(x2, (8277 + frac_mul16(-626, x2)))));
    debug_assert!(x2 <= 32766);
    (1 + x2) as i16
}

pub fn bitexact_log2tan(isin: i32, icos: i32) -> i32 {
    let lc = ec_ilog(icos as u32);
    let ls = ec_ilog(isin as u32);
    let icos = icos << (15 - lc);
    let isin = isin << (15 - ls);
    (ls - lc) * (1 << 11)
        + frac_mul16(isin, frac_mul16(isin, -2597) + 7932)
        - frac_mul16(icos, frac_mul16(icos, -2597) + 7932)
}

pub fn denormalise_bands(
    mode: &CeltMode,
    x: &[CeltNorm],
    freq: &mut [CeltSig],
    band_log_e: &[CeltGlog],
    start: i32,
    end: i32,
    m: i32,
    downsample: i32,
    silence: bool,
) {
    let n = (m as usize) * (mode.short_mdct_size as usize);
    debug_assert!(x.len() >= n);
    debug_assert!(freq.len() >= n);

    if silence {
        freq[..n].fill(0.0);
        return;
    }

    let ebands = mode.ebands;
    let mut bound = (m * ebands[end as usize] as i32) as usize;
    if downsample != 1 {
        bound = bound.min(n / downsample as usize);
    }

    let start_idx = (m * ebands[start as usize] as i32) as usize;
    if start_idx > 0 {
        freq[..start_idx].fill(0.0);
    }

    let mut x_idx = start_idx;
    for band in start..end {
        let band_start = (m * ebands[band as usize] as i32) as usize;
        let band_end = (m * ebands[band as usize + 1] as i32) as usize;
        debug_assert_eq!(x_idx, band_start);
        let lg = band_log_e[band as usize] + E_MEANS[band as usize];
        let g = celt_exp2_db(lg.min(32.0));
        for j in band_start..band_end {
            freq[j] = x[x_idx] * g;
            x_idx += 1;
        }
    }

    if bound < n {
        freq[bound..n].fill(0.0);
    }
}

pub fn anti_collapse(
    mode: &CeltMode,
    x: &mut [CeltNorm],
    collapse_masks: &[u8],
    lm: i32,
    channels: i32,
    size: i32,
    start: i32,
    end: i32,
    log_e: &[CeltGlog],
    prev1_log_e: &[CeltGlog],
    prev2_log_e: &[CeltGlog],
    pulses: &[i32],
    mut seed: u32,
    encode: bool,
) {
    let ebands = mode.ebands;
    let nb_ebands = mode.nb_ebands as usize;
    for band in start..end {
        let n0 = ebands[band as usize + 1] - ebands[band as usize];
        let depth = (celt_udiv((1 + pulses[band as usize]) as u32, n0 as u32) as i32) >> lm;
        let thresh = 0.5 * celt_exp2_db(-0.125 * depth as f32);
        let sqrt_1 = celt_rsqrt(((n0 as i32) << lm) as f32);

        for c in 0..channels {
            let mut prev1 = prev1_log_e[c as usize * nb_ebands + band as usize];
            let mut prev2 = prev2_log_e[c as usize * nb_ebands + band as usize];
            if !encode && channels == 1 {
                prev1 = prev1.max(prev1_log_e[nb_ebands + band as usize]);
                prev2 = prev2.max(prev2_log_e[nb_ebands + band as usize]);
            }
            let ediff = (log_e[c as usize * nb_ebands + band as usize] - prev1.min(prev2)).max(0.0);
            let mut r = 2.0 * celt_exp2_db(-ediff);
            if lm == 3 {
                r *= 1.41421356;
            }
            if r > thresh {
                r = thresh;
            }
            r *= sqrt_1;

            let base =
                c as usize * size as usize + ((ebands[band as usize] as i32) << lm) as usize;
            let lm_usize = lm as usize;
            let blocks = 1usize << lm_usize;
            let mut renormalize = false;
            for k in 0..blocks {
                let mask = 1u8 << k;
                if collapse_masks[band as usize * channels as usize + c as usize] & mask == 0 {
                    for j in 0..(n0 as usize) {
                        seed = celt_lcg_rand(seed);
                        let val = if (seed & 0x8000) != 0 { r } else { -r };
                        x[base + (j << lm_usize) + k] = val;
                    }
                    renormalize = true;
                }
            }
            if renormalize {
                renormalise_vector(
                    &mut x[base..base + ((n0 as i32) << lm) as usize],
                    (n0 as i32) << lm,
                    1.0,
                );
            }
        }
    }
}

pub fn deinterleave_hadamard(x: &mut [CeltNorm], n0: i32, stride: i32, hadamard: bool) {
    let n = (n0 * stride) as usize;
    let mut tmp = vec![0.0f32; n];
    if hadamard {
        let offset = (stride as usize).saturating_sub(2);
        let ordery = &ORDERY_TABLE[offset..offset + stride as usize];
        for i in 0..(stride as usize) {
            for j in 0..(n0 as usize) {
                tmp[(ordery[i] as usize) * n0 as usize + j] = x[j * stride as usize + i];
            }
        }
    } else {
        for i in 0..(stride as usize) {
            for j in 0..(n0 as usize) {
                tmp[i * n0 as usize + j] = x[j * stride as usize + i];
            }
        }
    }
    x[..n].copy_from_slice(&tmp[..n]);
}

pub fn interleave_hadamard(x: &mut [CeltNorm], n0: i32, stride: i32, hadamard: bool) {
    let n = (n0 * stride) as usize;
    let mut tmp = vec![0.0f32; n];
    if hadamard {
        let offset = (stride as usize).saturating_sub(2);
        let ordery = &ORDERY_TABLE[offset..offset + stride as usize];
        for i in 0..(stride as usize) {
            for j in 0..(n0 as usize) {
                tmp[j * stride as usize + i] = x[(ordery[i] as usize) * n0 as usize + j];
            }
        }
    } else {
        for i in 0..(stride as usize) {
            for j in 0..(n0 as usize) {
                tmp[j * stride as usize + i] = x[i * n0 as usize + j];
            }
        }
    }
    x[..n].copy_from_slice(&tmp[..n]);
}

pub fn haar1(x: &mut [CeltNorm], n0: i32, stride: i32) {
    let n0 = (n0 >> 1) as usize;
    let stride = stride as usize;
    for i in 0..stride {
        for j in 0..n0 {
            let idx0 = stride * 2 * j + i;
            let idx1 = stride * (2 * j + 1) + i;
            let tmp1 = 0.70710678 * x[idx0];
            let tmp2 = 0.70710678 * x[idx1];
            x[idx0] = tmp1 + tmp2;
            x[idx1] = tmp1 - tmp2;
        }
    }
}

pub fn compute_qn(n: i32, b: i32, offset: i32, pulse_cap: i32, stereo: bool) -> i32 {
    const EXP2_TABLE8: [i16; 8] = [16384, 17866, 19483, 21247, 23170, 25267, 27554, 30048];
    let mut n2 = 2 * n - 1;
    if stereo && n == 2 {
        n2 -= 1;
    }
    let bitres = BITRES as i32;
    let mut qb = celt_sudiv(b + n2 * offset, n2);
    qb = qb.min(b - pulse_cap - (4 << bitres));
    qb = qb.min(8 << bitres);

    if qb < (1 << bitres >> 1) {
        1
    } else {
        let idx = (qb & 0x7) as usize;
        let qn = (EXP2_TABLE8[idx] as i32) >> (14 - (qb >> bitres));
        ((qn + 1) >> 1) << 1
    }
}

pub fn intensity_stereo(
    mode: &CeltMode,
    x: &mut [CeltNorm],
    y: &[CeltNorm],
    band_e: &[CeltNorm],
    band_id: i32,
    n: i32,
) {
    let i = band_id as usize;
    let left = band_e[i];
    let right = band_e[i + mode.nb_ebands as usize];
    let norm = 1.0e-15 + celt_sqrt(1.0e-15 + left * left + right * right);
    let a1 = left / norm;
    let a2 = right / norm;
    for j in 0..(n as usize) {
        x[j] = a1 * x[j] + a2 * y[j];
    }
}

pub fn stereo_split(x: &mut [CeltNorm], y: &mut [CeltNorm], n: i32) {
    for j in 0..(n as usize) {
        let l = 0.70710678 * x[j];
        let r = 0.70710678 * y[j];
        x[j] = l + r;
        y[j] = r - l;
    }
}

pub fn stereo_merge(x: &mut [CeltNorm], y: &mut [CeltNorm], mid: f32, n: i32) {
    let mut xp = 0.0f32;
    let mut side = 0.0f32;
    for j in 0..(n as usize) {
        xp += y[j] * x[j];
        side += y[j] * y[j];
    }
    xp *= mid;
    let el = mid * mid + side - 2.0 * xp;
    let er = mid * mid + side + 2.0 * xp;
    if el < 6.0e-4 || er < 6.0e-4 {
        y[..(n as usize)].copy_from_slice(&x[..(n as usize)]);
        return;
    }
    let lgain = 1.0 / celt_sqrt(el);
    let rgain = 1.0 / celt_sqrt(er);
    for j in 0..(n as usize) {
        let l = mid * x[j];
        let r = y[j];
        x[j] = lgain * (l - r);
        y[j] = rgain * (l + r);
    }
}
