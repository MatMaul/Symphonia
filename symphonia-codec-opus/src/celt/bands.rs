// Symphonia
// Copyright (c) 2024 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! Band processing for CELT decoder.
//!
//! This module handles frequency band operations including:
//! - Band energy denormalization
//! - Anti-collapse processing
//! - The main band quantization/dequantization loop

#![allow(dead_code)]

#[cfg(not(feature = "std"))]
use alloc::{vec, vec::Vec};
#[cfg(feature = "std")]
use std::{vec, vec::Vec};

use crate::celt::constants::{BITRES, DB_SHIFT, Q15ONE, SPREAD_AGGRESSIVE};
use crate::celt::mode::CeltMode;
use crate::celt::rate::{bits2pulses, get_pulses, pulses2bits};
use crate::celt::tables::{E_MEANS, ORDERY_TABLE};
use crate::celt::vq::{alg_unquant, renormalise_vector};
use crate::entropy::RangeDecoder;
use crate::util::math::{celt_exp2, celt_ilog2, celt_rsqrt_norm, celt_sqrt, celt_udiv, mult16_16_q15, shl32, vshr32};

/// Linear congruential generator for pseudo-random numbers.
#[inline]
pub fn celt_lcg_rand(seed: u32) -> u32 {
    seed.wrapping_mul(1664525).wrapping_add(1013904223)
}

/// Denormalize band energies from quantized log-domain to linear scale.
///
/// Converts quantized log-energy values back to linear scale MDCT coefficients.
pub fn denormalise_bands(
    mode: &CeltMode,
    x: &[i32],
    freq: &mut [i32],
    band_log_e: &[i32],
    start: usize,
    end: usize,
    m: usize,
    downsample: usize,
    silence: bool,
) {
    let n = m * mode.short_mdct_size;
    let mut bound = m * mode.ebands[end] as usize;

    if downsample != 1 {
        bound = bound.min(n / downsample);
    }

    let (start, end, bound) = if silence {
        (0, 0, 0)
    } else {
        (start, end, bound)
    };

    // Zero out frequencies before start band
    let start_bin = m * mode.ebands[start] as usize;
    for f in freq.iter_mut().take(start_bin) {
        *f = 0;
    }

    let mut x_idx = start_bin;
    let mut f_idx = start_bin;

    for i in start..end {
        let band_end = m * mode.ebands[i + 1] as usize;

        // Compute linear gain from log energy
        let lg = band_log_e[i] + shl32(E_MEANS[i] as i32, 6);
        let shift = 16 - (lg >> DB_SHIFT);

        let g = if shift > 31 {
            0
        } else {
            celt_exp2_frac(lg & ((1 << DB_SHIFT) - 1))
        };

        // Apply gain to each coefficient in the band
        if shift < 0 {
            let shift = if shift < -2 {
                -2
            } else {
                shift
            };
            while f_idx < band_end {
                freq[f_idx] = shr32(mult16_16(x[x_idx], g.min(32767)), -shift);
                x_idx += 1;
                f_idx += 1;
            }
        } else {
            while f_idx < band_end {
                freq[f_idx] = shr32(mult16_16(x[x_idx], g), shift);
                x_idx += 1;
                f_idx += 1;
            }
        }
    }

    // Zero out remaining frequencies
    for f in freq.iter_mut().take(n).skip(bound) {
        *f = 0;
    }
}

/// Compute the fractional part of 2^x (Q15 output).
#[inline]
fn celt_exp2_frac(x: i32) -> i32 {
    let frac = x << 4;
    16383 + mult16_16_q15(frac, 22804 + mult16_16_q15(frac, 14819 + mult16_16_q15(10204, frac)))
}

/// Shift right helper.
#[inline]
fn shr32(a: i32, shift: i32) -> i32 {
    if shift >= 0 {
        a >> shift
    } else {
        a << -shift
    }
}

/// Multiply 16x16 helper.
#[inline]
fn mult16_16(a: i32, b: i32) -> i32 {
    a * b
}

/// Anti-collapse processing to prevent bands from going silent.
///
/// When bands have very low energy due to quantization, this adds
/// noise to prevent artifacts.
pub fn anti_collapse(
    mode: &CeltMode,
    x: &mut [Vec<i32>],
    collapse_masks: &[i16],
    lm: usize,
    channels: usize,
    _size: usize,
    start: usize,
    end: usize,
    log_e: &[i32],
    prev1_log_e: &[i32],
    prev2_log_e: &[i32],
    pulses: &[i32],
    seed: &mut u32,
) {
    for i in start..end {
        let n0 = (mode.ebands[i + 1] - mode.ebands[i]) as usize;

        // Compute depth (bits per sample)
        let depth = celt_udiv(1 + pulses[i], (mode.ebands[i + 1] - mode.ebands[i]) as i32) >> lm;

        // Compute threshold
        let thresh32 = shr32(celt_exp2(-shl32(depth as i32, 10 - BITRES)), 1);
        let thresh = mult16_32_q15(16384, thresh32.min(32767)); // 0.5 in Q15

        // Compute sqrt_1 for normalization
        let t = (n0 << lm) as i32;
        let shift = celt_ilog2(t as u32) >> 1;
        let t_scaled = shl32(t, (7 - shift as i32) * 2);
        let sqrt_1 = celt_rsqrt_norm(t_scaled);

        for c in 0..channels {
            let prev1 = if channels == 1 {
                prev1_log_e[i].max(prev1_log_e[mode.nb_ebands + i])
            } else {
                prev1_log_e[c * mode.nb_ebands + i]
            };

            let prev2 = if channels == 1 {
                prev2_log_e[i].max(prev2_log_e[mode.nb_ebands + i])
            } else {
                prev2_log_e[c * mode.nb_ebands + i]
            };

            let e_diff = (log_e[c * mode.nb_ebands + i] - prev1.min(prev2)).max(0);

            let r = if e_diff < 16384 {
                let r32 = shr32(celt_exp2(-e_diff), 1);
                (2 * r32.min(16383)) as i32
            } else {
                0
            };

            let r = if lm == 3 {
                mult16_16_q14(23170, r.min(23169)) // Scale for longest frames
            } else {
                r
            };

            let r = shr32(thresh.min(r), 1);
            let r = shr32(mult16_16_q15(sqrt_1, r), shift as i32);

            let x_start = (mode.ebands[i] as usize) << lm;
            let mut renormalize = false;

            for k in 0..(1 << lm) {
                // Check if this sub-band collapsed
                if (collapse_masks[i * channels + c] as u32) & (1 << k) == 0 {
                    // Fill with noise
                    let xk = x_start + k;
                    for j in 0..n0 {
                        *seed = celt_lcg_rand(*seed);
                        let val = if (*seed & 0x8000) != 0 { r } else { -r };
                        x[c][xk + (j << lm)] = val;
                    }
                    renormalize = true;
                }
            }

            if renormalize {
                let x_slice = &mut x[c][x_start..x_start + (n0 << lm)];
                renormalise_vector(x_slice, n0 << lm, Q15ONE);
            }
        }
    }
}

/// Q14 multiply helper.
#[inline]
fn mult16_16_q14(a: i32, b: i32) -> i32 {
    (a * b) >> 14
}

/// Q15 multiply with 32-bit result.
#[inline]
fn mult16_32_q15(a: i32, b: i32) -> i32 {
    ((a as i64 * b as i64) >> 15) as i32
}

/// Hadamard de-interleaving for band processing.
pub fn deinterleave_hadamard(x: &mut [i32], n0: usize, stride: usize, hadamard: bool) {
    let n = n0 * stride;
    let mut tmp = vec![0i32; n];

    if hadamard {
        let ordery = stride - 2;
        for i in 0..stride {
            for j in 0..n0 {
                tmp[ORDERY_TABLE[ordery + i] as usize * n0 + j] = x[j * stride + i];
            }
        }
    } else {
        for i in 0..stride {
            for j in 0..n0 {
                tmp[i * n0 + j] = x[j * stride + i];
            }
        }
    }

    x[..n].copy_from_slice(&tmp);
}

/// Hadamard interleaving for band processing.
pub fn interleave_hadamard(x: &mut [i32], n0: usize, stride: usize, hadamard: bool) {
    let n = n0 * stride;
    let mut tmp = vec![0i32; n];

    if hadamard {
        let ordery = stride - 2;
        for i in 0..stride {
            for j in 0..n0 {
                tmp[j * stride + i] = x[ORDERY_TABLE[ordery + i] as usize * n0 + j];
            }
        }
    } else {
        for i in 0..stride {
            for j in 0..n0 {
                tmp[j * stride + i] = x[i * n0 + j];
            }
        }
    }

    x[..n].copy_from_slice(&tmp);
}

/// Haar wavelet transform (single step).
pub fn haar1(x: &mut [i32], n0: usize, stride: usize) {
    let n0 = n0 >> 1;
    const SQRT_HALF: i32 = 23170; // sqrt(0.5) in Q15

    for i in 0..stride {
        for j in 0..n0 {
            let idx = i + stride * 2 * j;
            let tmp1 = mult16_16(SQRT_HALF, x[idx]);
            let tmp2 = mult16_16(SQRT_HALF, x[idx + stride]);
            x[idx] = pshr32(tmp1 + tmp2, 15);
            x[idx + stride] = pshr32(tmp1 - tmp2, 15);
        }
    }
}

/// Shift right with rounding.
#[inline]
fn pshr32(a: i32, shift: i32) -> i32 {
    (a + (1 << (shift - 1))) >> shift
}

/// Compute theta quantization step size.
fn compute_qn(n: usize, b: i32, offset: i32, pulse_cap: i32, stereo: bool) -> i32 {
    static EXP2_TABLE8: [i16; 8] = [16384, 17866, 19483, 21247, 23170, 25267, 27554, 30048];

    let n2 = 2 * n as i32 - 1 - if stereo && n == 2 { 1 } else { 0 };
    let mut qb = celt_udiv(b + n2 * offset, n2);
    qb = qb.min(b - pulse_cap - (4 << BITRES));
    qb = qb.min(8 << BITRES);

    if qb >= (1 << BITRES) >> 1 {
        let qn = (EXP2_TABLE8[(qb & 0x7) as usize] as i32) >> (14 - (qb >> BITRES));
        ((qn + 1) >> 1) << 1
    } else {
        1
    }
}

/// Band context for quant_band operations.
pub struct BandCtx<'a, 'b> {
    pub encode: bool,
    pub mode: &'a CeltMode,
    pub band_idx: usize,
    pub intensity: usize,
    pub spread: i32,
    pub tf_change: i32,
    pub dec: &'b mut RangeDecoder<'a>,
    pub remaining_bits: i32,
    pub seed: u32,
}

/// Split context for band splitting.
pub struct SplitCtx {
    pub inv: bool,
    pub imid: i32,
    pub iside: i32,
    pub delta: i32,
    pub itheta: i32,
    pub qalloc: i32,
}

/// Decode a single band with N=1 sample.
fn quant_band_n1(
    ctx: &mut BandCtx<'_, '_>,
    x: &mut [i32],
    y: Option<&mut [i32]>,
    b: i32,
    lowband_out: Option<&mut [i32]>,
) -> u32 {
    let stereo = y.is_some();

    let mut process = |out: &mut [i32]| {
        let sign = if ctx.remaining_bits >= 1 << BITRES {
            ctx.remaining_bits -= 1 << BITRES;
            ctx.dec.dec_bit_logp(1)
        } else {
            0
        };

        out[0] = if sign != 0 { -16384 } else { 16384 }; // NORM_SCALING
    };

    process(x);
    if let Some(y_slice) = y {
        process(y_slice);
    }

    if let Some(lb_out) = lowband_out {
        lb_out[0] = x[0] >> 4;
    }

    1
}

/// Quantize a single partition recursively.
fn quant_partition(
    ctx: &mut BandCtx<'_, '_>,
    x: &mut [i32],
    n: usize,
    b: i32,
    b_blocks: usize,
    lowband: Option<&[i32]>,
    lm: i32,
    gain: i32,
    fill: u32,
) -> u32 {
    let mode = ctx.mode;
    let i = ctx.band_idx;
    let cache = &mode.cache.bits;
    let cache_idx = mode.cache.index[((lm + 1) as usize) * mode.nb_ebands + i] as usize;

    // Check if we need to split
    if lm != -1
        && b > cache[cache_idx + cache[cache_idx] as usize] as i32 + 12
        && n > 2
    {
        // Split the band in two
        let n_half = n >> 1;
        let lm = lm - 1;
        let b_half = (b_blocks + 1) >> 1;
        let fill_half = if b_blocks == 1 {
            (fill & 1) | (fill << 1)
        } else {
            fill
        };

        // For MVP, use simpler bit allocation without theta coding
        let mbits = (b / 2).max(0);
        let sbits = b - mbits;

        // Recursively quantize both halves
        let (x_lo, x_hi) = x.split_at_mut(n_half);
        let lb_lo = lowband.map(|lb| &lb[..n_half]);
        let lb_hi = lowband.map(|lb| &lb[n_half..]);

        let cm_lo = quant_partition(ctx, x_lo, n_half, mbits, b_half, lb_lo, lm, gain, fill_half & ((1 << b_half) - 1));
        let cm_hi = quant_partition(ctx, x_hi, n_half, sbits, b_half, lb_hi, lm, gain, fill_half >> b_half);

        cm_lo | (cm_hi << b_half)
    } else {
        // Base case: decode pulses
        let q = bits2pulses(mode, i, lm, b);
        let curr_bits = pulses2bits(mode, i, lm, q);
        ctx.remaining_bits -= curr_bits;

        // Ensure we don't bust the budget
        let mut q = q;
        while ctx.remaining_bits < 0 && q > 0 {
            ctx.remaining_bits += curr_bits;
            q -= 1;
            let curr_bits = pulses2bits(mode, i, lm, q);
            ctx.remaining_bits -= curr_bits;
        }

        if q != 0 {
            let k = get_pulses(q);
            alg_unquant(x, n, k as usize, ctx.spread, b_blocks, ctx.dec, gain)
        } else {
            // No pulses, fill with noise or fold from lowband
            let cm_mask = (1u32 << b_blocks) - 1;
            let fill = fill & cm_mask;

            if fill == 0 {
                for val in x.iter_mut().take(n) {
                    *val = 0;
                }
            } else if let Some(lb) = lowband {
                // Folded spectrum with small noise
                for j in 0..n {
                    ctx.seed = celt_lcg_rand(ctx.seed);
                    let tmp = 128; // Small noise ~48dB below normal
                    let noise = if (ctx.seed & 0x8000) != 0 { tmp } else { -tmp };
                    x[j] = lb[j] + noise;
                }
                renormalise_vector(x, n, gain);
                cm_mask
            } else {
                // Pure noise
                for j in 0..n {
                    ctx.seed = celt_lcg_rand(ctx.seed);
                    x[j] = (ctx.seed as i32) >> 20;
                }
                renormalise_vector(x, n, gain);
                cm_mask
            }
        }
    }
}

/// Decode a complete band (mono).
pub fn quant_band(
    ctx: &mut BandCtx<'_, '_>,
    x: &mut [i32],
    n: usize,
    b: i32,
    b_blocks: usize,
    lowband: Option<&[i32]>,
    lm: i32,
    lowband_out: Option<&mut [i32]>,
    gain: i32,
    lowband_scratch: Option<&mut [i32]>,
    fill: u32,
) -> u32 {
    let n0 = n;
    let b0 = b_blocks;
    let mut n_b = n / b_blocks;
    let mut time_divide = 0;
    let mut recombine = 0;
    let long_blocks = b_blocks == 1;
    let tf_change = ctx.tf_change;

    // Handle N=1 special case
    if n == 1 {
        return quant_band_n1(ctx, x, None, b, lowband_out);
    }

    // Band recombining for TF resolution changes
    if tf_change > 0 {
        recombine = tf_change as usize;
    }

    // Copy lowband to scratch if needed
    let lowband = if let (Some(scratch), Some(lb)) = (lowband_scratch, lowband) {
        if recombine != 0 || (n_b & 1) == 0 && tf_change < 0 || b0 > 1 {
            scratch[..n].copy_from_slice(&lb[..n]);
            Some(&scratch[..n] as &[i32])
        } else {
            Some(lb)
        }
    } else {
        lowband
    };

    // Apply recombine transforms
    let mut b_blocks = b_blocks;
    for _ in 0..recombine {
        if let Some(lb) = lowband {
            // haar1 on lowband - would need mutable access
        }
        b_blocks >>= 1;
        n_b <<= 1;
    }

    // Time resolution changes
    while (n_b & 1) == 0 && tf_change < 0 {
        b_blocks <<= 1;
        n_b >>= 1;
        time_divide += 1;
    }

    let b0_final = b_blocks;
    let n_b0 = n_b;

    // Reorganize samples
    if b0_final > 1 {
        deinterleave_hadamard(x, n_b >> recombine, b0_final << recombine, long_blocks);
    }

    // Do the actual quantization
    let mut cm = quant_partition(ctx, x, n, b, b_blocks, lowband, lm, gain, fill);

    // Undo sample reorganization
    if b0_final > 1 {
        interleave_hadamard(x, n_b0 >> recombine, b0_final << recombine, long_blocks);
    }

    // Undo time-freq changes
    let mut b_blocks = b0_final;
    let mut n_b = n_b0;
    for _ in 0..time_divide {
        b_blocks >>= 1;
        n_b <<= 1;
        cm |= cm >> b_blocks;
        haar1(x, n_b, b_blocks);
    }

    // Undo recombine
    for _ in 0..recombine {
        haar1(x, n0 >> recombine, 1 << recombine);
    }

    // Scale output for folding
    if let Some(lb_out) = lowband_out {
        let norm = celt_sqrt(shl32(n0 as i32, 22));
        for j in 0..n0 {
            lb_out[j] = mult16_16_q15(norm, x[j]);
        }
    }

    cm & ((1 << b_blocks) - 1)
}

/// Decode all bands in a frame.
pub fn quant_all_bands(
    mode: &CeltMode,
    start: usize,
    end: usize,
    x: &mut [i32],
    y: Option<&mut [i32]>,
    collapse_masks: &mut [i16],
    pulses: &[i32],
    short_blocks: bool,
    spread: i32,
    dual_stereo: bool,
    intensity: usize,
    tf_res: &[i32],
    total_bits: i32,
    mut balance: i32,
    dec: &mut RangeDecoder<'_>,
    lm: usize,
    coded_bands: usize,
    seed: &mut u32,
) {
    let m = 1usize << lm;
    let b_blocks = if short_blocks { m } else { 1 };
    let channels = if y.is_some() { 2 } else { 1 };

    let norm_offset = m * mode.ebands[start] as usize;
    let norm_size = m * mode.ebands[mode.nb_ebands - 1] as usize - norm_offset;
    let mut norm = vec![0i32; channels * norm_size];

    let lowband_scratch_offset = m * mode.ebands[mode.nb_ebands - 1] as usize;
    let mut lowband_offset = 0usize;
    let mut update_lowband = true;

    // For dual_stereo, split bits between channels
    let mut dual_stereo = dual_stereo;

    for i in start..end {
        let last = i == end - 1;
        let x_start = m * mode.ebands[i] as usize;
        let n = m * (mode.ebands[i + 1] - mode.ebands[i]) as usize;

        let tell = dec.tell_frac();
        if i != start {
            balance -= tell;
        }
        let remaining_bits = total_bits - tell - 1;

        // Compute bits for this band
        let b = if i <= coded_bands - 1 {
            let curr_balance = celt_udiv(balance, 3.min((coded_bands - i) as i32));
            (remaining_bits + 1).max(0).min(16383).min(pulses[i] + curr_balance)
        } else {
            0
        };

        // Determine effective lowband
        let effective_lowband = if lowband_offset != 0
            && (spread != SPREAD_AGGRESSIVE || b_blocks > 1 || tf_res[i] < 0)
        {
            (m * mode.ebands[lowband_offset] as usize).saturating_sub(norm_offset + n).max(0)
        } else {
            0
        };

        // Build collapse mask from previous bands
        let (mut x_cm, mut y_cm) = if effective_lowband > 0 {
            let mut xcm = 0i64;
            let mut ycm = 0i64;
            // Simplified: just use full mask
            xcm = (1 << b_blocks) - 1;
            ycm = (1 << b_blocks) - 1;
            (xcm as u32, ycm as u32)
        } else {
            ((1u32 << b_blocks) - 1, (1u32 << b_blocks) - 1)
        };

        // Handle intensity stereo transition
        if dual_stereo && i == intensity {
            dual_stereo = false;
        }

        let tf_change = tf_res[i];

        // Create band context
        let mut ctx = BandCtx {
            encode: false,
            mode,
            band_idx: i,
            intensity,
            spread,
            tf_change,
            dec,
            remaining_bits,
            seed: *seed,
        };

        // Get lowband for folding
        let lowband = if effective_lowband > 0 {
            Some(&norm[..effective_lowband + n])
        } else {
            None
        };

        // Decode the band
        let x_slice = &mut x[x_start..x_start + n];

        if channels == 1 || !dual_stereo {
            // Mono or joint stereo
            x_cm = quant_band(
                &mut ctx,
                x_slice,
                n,
                b,
                b_blocks,
                lowband.map(|lb| &lb[effective_lowband..]),
                lm as i32,
                if !last { Some(&mut norm[m * mode.ebands[i] as usize - norm_offset..]) } else { None },
                Q15ONE,
                None,
                x_cm | y_cm,
            );
            y_cm = x_cm;
        } else {
            // Dual stereo - decode each channel separately
            x_cm = quant_band(
                &mut ctx,
                x_slice,
                n,
                b / 2,
                b_blocks,
                lowband.map(|lb| &lb[effective_lowband..]),
                lm as i32,
                if !last { Some(&mut norm[m * mode.ebands[i] as usize - norm_offset..]) } else { None },
                Q15ONE,
                None,
                x_cm,
            );

            if let Some(y_data) = y.as_deref() {
                // Would need mutable y here for real stereo
            }
        }

        *seed = ctx.seed;

        // Update collapse masks
        collapse_masks[i * channels] = x_cm as i16;
        if channels == 2 {
            collapse_masks[i * channels + 1] = y_cm as i16;
        }

        balance += pulses[i] + tell;

        // Update lowband offset for folding
        if m * mode.ebands[i] as usize >= norm_offset + n
            && (update_lowband || lowband_offset == 0)
        {
            lowband_offset = i;
        }

        update_lowband = b > (n as i32) << BITRES;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_celt_lcg_rand() {
        let mut seed = 12345u32;
        seed = celt_lcg_rand(seed);
        assert_ne!(seed, 12345);

        // Check sequence is deterministic
        let seed1 = celt_lcg_rand(0);
        let seed2 = celt_lcg_rand(0);
        assert_eq!(seed1, seed2);
    }

    #[test]
    fn test_haar1() {
        let mut x = vec![100, 100, 0, 0];
        haar1(&mut x, 4, 1);
        // After Haar transform, energy should be preserved but redistributed
        assert!(x[0] != 100 || x[1] != 100);
    }
}
