use super::cwrs::decode_pulses;
use super::entdec::EcDec;
use super::intrin::celt_udiv;
use super::math::{celt_cos_norm, celt_div, celt_sqrt};
use super::types::{CeltNorm, Q15_ONE};
use alloc::vec;

pub const SPREAD_NONE: i32 = 0;
pub const SPREAD_LIGHT: i32 = 1;
pub const SPREAD_NORMAL: i32 = 2;
pub const SPREAD_AGGRESSIVE: i32 = 3;

fn exp_rotation1(x: &mut [CeltNorm], len: usize, stride: usize, c: f32, s: f32) {
    let ms = -s;
    let mut idx = 0usize;
    for _ in 0..(len.saturating_sub(stride)) {
        let x1 = x[idx];
        let x2 = x[idx + stride];
        x[idx + stride] = c * x2 + s * x1;
        x[idx] = c * x1 + ms * x2;
        idx += 1;
    }

    if len < 2 * stride + 1 {
        return;
    }
    let mut idx = len - 2 * stride - 1;
    loop {
        let x1 = x[idx];
        let x2 = x[idx + stride];
        x[idx + stride] = c * x2 + s * x1;
        x[idx] = c * x1 + ms * x2;
        if idx == 0 {
            break;
        }
        idx -= 1;
    }
}

pub fn exp_rotation(x: &mut [CeltNorm], len: i32, dir: i32, stride: i32, k: i32, spread: i32) {
    const SPREAD_FACTOR: [i32; 3] = [15, 10, 5];
    if 2 * k >= len || spread == SPREAD_NONE {
        return;
    }
    let factor = SPREAD_FACTOR[(spread - 1) as usize];
    let gain = celt_div(Q15_ONE * len as f32, (len + factor * k) as f32);
    let theta = 0.5 * gain * gain;
    let c = celt_cos_norm(theta);
    let s = celt_cos_norm(Q15_ONE - theta);

    let mut stride2 = 0;
    if len >= 8 * stride {
        stride2 = 1;
        while (stride2 * stride2 + stride2) * stride + (stride >> 2) < len {
            stride2 += 1;
        }
    }

    let len = celt_udiv(len as u32, stride as u32) as usize;
    for i in 0..(stride as usize) {
        let start = i * len;
        let band = &mut x[start..start + len];
        if dir < 0 {
            if stride2 != 0 {
                exp_rotation1(band, len, stride2 as usize, s, c);
            }
            exp_rotation1(band, len, 1, c, s);
        } else {
            exp_rotation1(band, len, 1, c, -s);
            if stride2 != 0 {
                exp_rotation1(band, len, stride2 as usize, s, -c);
            }
        }
    }
}

fn normalise_residual(iy: &[i32], x: &mut [CeltNorm], n: usize, ryy: f32, gain: f32) {
    let scale = if ryy > 0.0 { gain / celt_sqrt(ryy) } else { 0.0 };
    for i in 0..n {
        x[i] = iy[i] as f32 * scale;
    }
}

fn extract_collapse_mask(iy: &[i32], n: i32, b: i32) -> u32 {
    if b <= 1 {
        return 1;
    }
    let n0 = celt_udiv(n as u32, b as u32) as usize;
    let mut collapse_mask = 0u32;
    for i in 0..(b as usize) {
        let mut tmp = 0;
        for j in 0..n0 {
            tmp |= iy[i * n0 + j];
        }
        collapse_mask |= ((tmp != 0) as u32) << i;
    }
    collapse_mask
}

pub fn alg_unquant(
    x: &mut [CeltNorm],
    n: i32,
    k: i32,
    spread: i32,
    b: i32,
    dec: &mut EcDec<'_>,
    gain: f32,
) -> u32 {
    let mut iy = vec![0i32; n as usize];
    let ryy = decode_pulses(&mut iy, n, k, dec);
    normalise_residual(&iy, x, n as usize, ryy, gain);
    exp_rotation(x, n, -1, b, k, spread);
    extract_collapse_mask(&iy, n, b)
}

pub fn renormalise_vector(x: &mut [CeltNorm], n: i32, gain: f32) {
    let mut energy = 1.0e-15f32;
    for i in 0..(n as usize) {
        energy += x[i] * x[i];
    }
    let scale = gain / celt_sqrt(energy);
    for i in 0..(n as usize) {
        x[i] *= scale;
    }
}
