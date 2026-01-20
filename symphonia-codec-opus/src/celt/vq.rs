// Symphonia
// Copyright (c) 2024 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! Vector Quantization (VQ) decoder for CELT.
//!
//! This module implements the pyramid vector quantization (PVQ) decoding
//! and spreading operations used in CELT.

use crate::celt::cwrs::decode_pulses;
use crate::celt::tables::{SPREAD_AGGRESSIVE, SPREAD_LIGHT, SPREAD_NONE, SPREAD_NORMAL};
use crate::entropy::RangeDecoder;
use crate::util::math::{celt_ilog2, celt_rsqrt_norm, celt_udiv};

/// Spread factors for different spreading modes.
const SPREAD_FACTOR: [i32; 3] = [15, 10, 5];

/// Q15 constant for 1.0.
const Q15_ONE: i32 = 32767;

/// Helper: Multiply-accumulate in Q15.
#[inline]
fn mult16_16_q15(a: i32, b: i32) -> i32 {
    (a * b) >> 15
}

/// Helper: Multiply-accumulate with rounding in Q15.
#[inline]
fn mult16_16_p15(a: i32, b: i32) -> i32 {
    (16384 + a * b) >> 15
}

/// Helper: Shift right with rounding.
#[inline]
fn pshr32(a: i32, shift: i32) -> i32 {
    (a + (1 << (shift - 1))) >> shift
}

/// Helper: Variable shift right.
#[inline]
fn vshr32(x: i32, shift: i32) -> i32 {
    if shift > 0 {
        x >> shift
    } else {
        x << -shift
    }
}

/// Helper: Negate 16-bit value.
#[inline]
fn neg16(a: i32) -> i32 {
    -a
}

/// Helper: Extract 16-bit value (saturate to 16 bits if needed).
#[inline]
fn extract16(x: i32) -> i32 {
    x.clamp(i16::MIN as i32, i16::MAX as i32)
}

/// Helper: Half of a 16-bit value.
#[inline]
fn half16(x: i32) -> i32 {
    x >> 1
}

/// Compute approximate cosine using Q15 fixed-point.
///
/// Input is in Q15 format where 1.0 = 32768.
/// Uses a polynomial approximation.
#[inline]
fn celt_cos_norm(x: i32) -> i32 {
    // Approximation: cos(pi/2 * x) ≈ 1 - x^2/2 for small x
    // But we use a more accurate polynomial:
    // cos(x) ≈ 1 - x²/2 + x⁴/24 in normalized form
    let x2 = mult16_16_q15(x, x);
    let x4 = mult16_16_q15(x2, x2);

    // Coefficients tuned for Q15 input
    // Using: 32767 - x²/2 + x⁴/24
    Q15_ONE - (x2 >> 1) + (x4 >> 4)
}

/// Perform a single exponential rotation step.
///
/// This applies a rotation to spread energy across the vector.
fn exp_rotation1(x: &mut [i32], len: usize, stride: usize, c: i32, s: i32) {
    let ms = neg16(s);

    // Forward pass
    for i in 0..(len - stride) {
        let x1 = x[i];
        let x2 = x[i + stride];
        x[i + stride] = extract16(pshr32(mult16_16_q15(c, x2) + mult16_16_q15(s, x1), 0));
        x[i] = extract16(pshr32(mult16_16_q15(c, x1) + mult16_16_q15(ms, x2), 0));
    }

    // Backward pass
    for i in (0..=(len - 2 * stride - 1)).rev() {
        let x1 = x[i];
        let x2 = x[i + stride];
        x[i + stride] = extract16(pshr32(mult16_16_q15(c, x2) + mult16_16_q15(s, x1), 0));
        x[i] = extract16(pshr32(mult16_16_q15(c, x1) + mult16_16_q15(ms, x2), 0));
    }
}

/// Apply exponential rotation (spreading) to a vector.
///
/// This operation spreads the energy of pulses across the vector
/// to prevent them from clustering together.
///
/// # Arguments
/// * `x` - Vector to rotate (modified in place)
/// * `len` - Length of vector
/// * `dir` - Direction: 1 for forward (encode), -1 for inverse (decode)
/// * `stride` - Stride between blocks
/// * `k` - Number of pulses
/// * `spread` - Spreading mode (SPREAD_NONE, SPREAD_LIGHT, SPREAD_NORMAL, SPREAD_AGGRESSIVE)
pub fn exp_rotation(x: &mut [i32], len: usize, dir: i32, stride: usize, k: usize, spread: i32) {
    if 2 * k >= len || spread == SPREAD_NONE {
        return;
    }

    let factor = SPREAD_FACTOR[(spread - 1) as usize];

    // Compute gain = len / (len + factor * k)
    let gain = celt_udiv(Q15_ONE * len as i32, len as i32 + factor * k as i32);

    // theta = gain^2 / 2
    let theta = half16(mult16_16_q15(gain, gain));

    // Compute cos and sin
    let c = celt_cos_norm(theta);
    let s = celt_cos_norm(Q15_ONE - theta);

    // Compute secondary stride if len is large enough
    let mut stride2 = 0usize;
    if len >= 8 * stride {
        stride2 = 1;
        while (stride2 * stride2 + stride2) * stride + (stride >> 2) < len {
            stride2 += 1;
        }
    }

    let block_len = celt_udiv(len as i32, stride as i32) as usize;

    for i in 0..stride {
        let x_block = &mut x[i * block_len..(i + 1) * block_len];

        if dir < 0 {
            // Inverse rotation (decoding)
            if stride2 != 0 {
                exp_rotation1(x_block, block_len, stride2, s, c);
            }
            exp_rotation1(x_block, block_len, 1, c, s);
        } else {
            // Forward rotation (encoding)
            exp_rotation1(x_block, block_len, 1, c, neg16(s));
            if stride2 != 0 {
                exp_rotation1(x_block, block_len, stride2, s, neg16(c));
            }
        }
    }
}

/// Normalize a residual vector after PVQ decoding.
///
/// # Arguments
/// * `iy` - Input pulse vector (integer)
/// * `x` - Output normalized vector (fixed-point)
/// * `n` - Vector length
/// * `ryy` - Sum of squared pulses
/// * `gain` - Target gain
pub fn normalise_residual(iy: &[i32], x: &mut [i32], n: usize, ryy: i32, gain: i32) {
    let k = celt_ilog2(ryy as u32) >> 1;
    let t = vshr32(ryy, 2 * (k as i32 - 7));
    let g = mult16_16_p15(celt_rsqrt_norm(t), gain);

    for i in 0..n {
        x[i] = extract16(pshr32(mult16_16_q15(g, iy[i]), k as i32 + 1));
    }
}

/// Extract the collapse mask from a decoded pulse vector.
///
/// The collapse mask indicates which sub-bands have non-zero energy.
///
/// # Arguments
/// * `iy` - Pulse vector
/// * `n` - Total vector length
/// * `b` - Number of sub-bands
///
/// # Returns
/// A bitmask where bit i is set if sub-band i has non-zero pulses
pub fn extract_collapse_mask(iy: &[i32], n: usize, b: usize) -> u32 {
    if b <= 1 {
        return 1;
    }

    let n0 = celt_udiv(n as i32, b as i32) as usize;
    let mut collapse_mask = 0u32;

    for i in 0..b {
        let mut tmp = 0i32;
        for j in 0..n0 {
            tmp |= iy[i * n0 + j];
        }
        if tmp != 0 {
            collapse_mask |= 1 << i;
        }
    }

    collapse_mask
}

/// Unquantize (decode) a PVQ vector.
///
/// This is the main decoding function that:
/// 1. Decodes the pulse codeword from the bitstream
/// 2. Normalizes the result
/// 3. Applies inverse spreading rotation
///
/// # Arguments
/// * `x` - Output vector (modified in place)
/// * `n` - Vector length
/// * `k` - Number of pulses
/// * `spread` - Spreading mode
/// * `b` - Number of blocks for collapse mask
/// * `dec` - Range decoder
/// * `gain` - Target gain
///
/// # Returns
/// Collapse mask indicating which sub-bands are non-zero
pub fn alg_unquant(
    x: &mut [i32],
    n: usize,
    k: usize,
    spread: i32,
    b: usize,
    dec: &mut RangeDecoder<'_>,
    gain: i32,
) -> u32 {
    debug_assert!(k > 0, "alg_unquant() needs at least one pulse");
    debug_assert!(n > 1, "alg_unquant() needs at least two dimensions");

    // Allocate temporary buffer for pulses
    let mut iy = vec![0i32; n];

    // Decode pulses using CWRS
    let ryy = decode_pulses(&mut iy, n, k, dec);

    // Normalize the residual
    normalise_residual(&iy, x, n, ryy, gain);

    // Apply inverse rotation
    exp_rotation(x, n, -1, b, k, spread);

    // Extract collapse mask
    extract_collapse_mask(&iy, n, b)
}

/// Renormalize a vector to unit energy with specified gain.
///
/// # Arguments
/// * `x` - Vector to normalize (modified in place)
/// * `n` - Vector length
/// * `gain` - Target gain in Q15
pub fn renormalise_vector(x: &mut [i32], n: usize, gain: i32) {
    const EPSILON: i32 = 1;

    // Compute energy E = sum(x[i]^2)
    let mut e = EPSILON;
    for i in 0..n {
        e += mult16_16_q15(x[i], x[i]);
    }

    // Normalize
    let k = celt_ilog2(e as u32) >> 1;
    let t = vshr32(e, 2 * (k as i32 - 7));
    let g = mult16_16_p15(celt_rsqrt_norm(t), gain);

    for i in 0..n {
        x[i] = extract16(pshr32(mult16_16_q15(g, x[i]), k as i32 + 1));
    }
}

/// Compute stereo theta angle from mid/side energies.
///
/// # Arguments
/// * `x` - First channel
/// * `y` - Second channel
/// * `stereo` - Whether to compute as M/S stereo
/// * `n` - Vector length
///
/// # Returns
/// Theta angle in Q15 format
pub fn stereo_itheta(x: &[i32], y: &[i32], stereo: bool, n: usize) -> i32 {
    const EPSILON: i32 = 1;

    let mut emid = EPSILON;
    let mut eside = EPSILON;

    if stereo {
        // Compute M/S energies
        for i in 0..n {
            let m = (x[i] >> 1) + (y[i] >> 1);
            let s = (x[i] >> 1) - (y[i] >> 1);
            emid += mult16_16_q15(m, m);
            eside += mult16_16_q15(s, s);
        }
    } else {
        // Compute L/R energies
        for i in 0..n {
            emid += mult16_16_q15(x[i], x[i]);
            eside += mult16_16_q15(y[i], y[i]);
        }
    }

    // Compute sqrt of energies
    let mid = celt_sqrt(emid);
    let side = celt_sqrt(eside);

    // Compute atan2(side, mid) scaled to [0, 16384]
    // Using approximation: 2/pi * atan2(y, x) ≈ y / (|x| + |y|) for rough estimate
    // More accurate: use CORDIC or polynomial
    let itheta = mult16_16_q15(20861, celt_atan2p(side, mid)); // 20861 ≈ 2/pi in Q15

    itheta
}

/// Approximate integer square root.
#[inline]
fn celt_sqrt(x: i32) -> i32 {
    if x <= 0 {
        return 0;
    }

    // Newton-Raphson iteration
    let mut guess = 1i32 << ((32 - x.leading_zeros()) / 2);
    guess = (guess + x / guess) >> 1;
    guess = (guess + x / guess) >> 1;
    guess = (guess + x / guess) >> 1;

    guess
}

/// Approximate atan2 for positive values.
#[inline]
fn celt_atan2p(y: i32, x: i32) -> i32 {
    if x == 0 {
        return if y > 0 { 16384 } else { 0 }; // pi/2 or 0 in Q14
    }

    // Simple approximation using ratio
    // atan(y/x) ≈ y/x for small angles, with correction
    let ratio = if y < x {
        (y << 14) / x
    } else {
        16384 - ((x << 14) / y)
    };

    ratio.clamp(0, 16384)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mult16_16_q15() {
        // 0.5 * 0.5 = 0.25
        let half = 16384; // 0.5 in Q15
        assert_eq!(mult16_16_q15(half, half), 8192); // 0.25 in Q15
    }

    #[test]
    fn test_extract_collapse_mask() {
        // All zeros -> mask is 0 (except minimum is 1 for b<=1)
        let iy = [0, 0, 0, 0];
        assert_eq!(extract_collapse_mask(&iy, 4, 2), 0);

        // First half has pulse
        let iy = [1, 0, 0, 0];
        assert_eq!(extract_collapse_mask(&iy, 4, 2), 1);

        // Second half has pulse
        let iy = [0, 0, 1, 0];
        assert_eq!(extract_collapse_mask(&iy, 4, 2), 2);

        // Both halves have pulses
        let iy = [1, 0, 1, 0];
        assert_eq!(extract_collapse_mask(&iy, 4, 2), 3);
    }

    #[test]
    fn test_celt_sqrt() {
        assert_eq!(celt_sqrt(0), 0);
        assert_eq!(celt_sqrt(1), 1);
        assert_eq!(celt_sqrt(4), 2);
        assert_eq!(celt_sqrt(9), 3);
        assert_eq!(celt_sqrt(100), 10);
    }

    #[test]
    fn test_exp_rotation_noop() {
        // With SPREAD_NONE, rotation should be a no-op
        let mut x = [100, 200, 300, 400];
        let original = x;
        exp_rotation(&mut x, 4, 1, 1, 2, SPREAD_NONE);
        assert_eq!(x, original);
    }
}
