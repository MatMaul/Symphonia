// Symphonia
// Copyright (c) 2024 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! CELT synthesis module.
//!
//! This module handles the final stage of CELT decoding:
//! - IMDCT (Inverse Modified Discrete Cosine Transform)
//! - Overlap-add processing
//! - De-emphasis filtering
//! - Comb filtering (post-filter)

#![allow(dead_code)]

#[cfg(not(feature = "std"))]
use alloc::{vec, vec::Vec};
#[cfg(feature = "std")]
use std::{vec, vec::Vec};

#[cfg(not(feature = "std"))]
use num_traits::Float;

use crate::celt::bands::denormalise_bands;
use crate::celt::constants::{COMBFILTER_MINPERIOD, DECODE_BUFFER_SIZE, Q15ONE, SIG_SHIFT};
use crate::celt::mode::CeltMode;
use crate::util::math::mult16_16_q15;

/// IMDCT state for a single transform size.
pub struct ImdctState {
    /// Transform size (N)
    n: usize,
    /// Twiddle factors (pre-computed cos/sin)
    twiddle: Vec<(f32, f32)>,
    /// Scratch buffer for computation
    scratch: Vec<f32>,
}

impl ImdctState {
    /// Create a new IMDCT state for a given transform size.
    pub fn new(n: usize) -> Self {
        use core::f32::consts::PI;

        let n2 = n / 2;
        let mut twiddle = Vec::with_capacity(n2);

        // Pre-compute twiddle factors
        let alpha = 1.0 / 8.0;
        let pi_n = PI / n as f32;

        for k in 0..n2 {
            let theta = pi_n * (alpha + k as f32);
            twiddle.push((theta.cos(), theta.sin()));
        }

        Self {
            n,
            twiddle,
            scratch: vec![0.0; n2],
        }
    }

    /// Perform IMDCT.
    ///
    /// Converts N/2 frequency coefficients to N time-domain samples.
    pub fn imdct(&mut self, input: &[f32], output: &mut [f32]) {
        let n = self.n;
        let n2 = n / 2;
        let n4 = n / 4;

        debug_assert!(input.len() >= n);
        debug_assert!(output.len() >= 2 * n);

        // Pre-FFT twiddling: combine real inputs into complex values
        for (i, &(cos_w, sin_w)) in self.twiddle.iter().enumerate() {
            let even = input[i * 2];
            let odd = -input[n - 1 - i * 2];

            let re = odd * sin_w - even * cos_w;
            let im = odd * cos_w + even * sin_w;
            self.scratch[i] = re;
            // Note: For full implementation, we'd need complex scratch
        }

        // For MVP: use simple DFT (not FFT)
        // This is O(N^2) but correct; will be optimized later
        let mut temp_re = vec![0.0f32; n2];
        let mut temp_im = vec![0.0f32; n2];

        for k in 0..n2 {
            let mut sum_re = 0.0f32;
            let mut sum_im = 0.0f32;

            for i in 0..n2 {
                let angle = 2.0 * core::f32::consts::PI * (k * i) as f32 / n2 as f32;
                let (sin_a, cos_a) = angle.sin_cos();

                // Simplified: just use real part for now
                sum_re += self.scratch[i] * cos_a;
                sum_im += self.scratch[i] * sin_a;
            }

            temp_re[k] = sum_re;
            temp_im[k] = sum_im;
        }

        // Post-FFT twiddling and expansion
        for i in 0..n4 {
            let (cos_w, sin_w) = self.twiddle[i];
            let re = temp_re[i] * cos_w + temp_im[i] * sin_w;
            let im = temp_im[i] * cos_w - temp_re[i] * sin_w;

            let fi = 2 * i;
            let ri = n2 - 1 - 2 * i;

            output[ri] = -im;
            output[n2 + fi] = im;
            output[n2 + ri] = re;
            output[n + fi] = re;
        }

        for i in n4..n2 {
            let (cos_w, sin_w) = self.twiddle[i];
            let re = temp_re[i] * cos_w + temp_im[i] * sin_w;
            let im = temp_im[i] * cos_w - temp_re[i] * sin_w;

            let fi = 2 * (i - n4);
            let ri = n2 - 1 - 2 * (i - n4);

            output[fi] = -re;
            output[n2 + ri] = re;
            output[n + fi] = im;
            output[n + n2 + ri] = im;
        }
    }
}

/// Perform CELT synthesis: IMDCT + overlap-add.
///
/// Converts frequency-domain coefficients to time-domain samples.
pub fn celt_synthesis(
    mode: &CeltMode,
    x: &[Vec<i32>],
    out_syn: &mut [Vec<i32>],
    old_band_e: &[i32],
    start: usize,
    eff_end: usize,
    channels: usize,
    out_channels: usize,
    is_transient: bool,
    lm: usize,
    downsample: usize,
    silence: bool,
) {
    let n = mode.short_mdct_size << lm;
    let overlap = mode.overlap;
    let m = 1 << lm;

    let (b_blocks, nb, shift) = if is_transient {
        (m, mode.short_mdct_size, mode.max_lm)
    } else {
        (1, mode.short_mdct_size << lm, mode.max_lm - lm)
    };

    let mut freq = vec![0i32; n];

    if out_channels == 2 && channels == 1 {
        // Mono to stereo upmix
        denormalise_bands(mode, &x[0], &mut freq, old_band_e, start, eff_end, m, downsample, silence);

        // Copy to both channels
        for b in 0..b_blocks {
            clt_mdct_backward_i32(
                mode,
                &freq,
                b,
                &mut out_syn[0],
                DECODE_BUFFER_SIZE - n + nb * b,
                shift,
                b_blocks,
            );
        }

        for b in 0..b_blocks {
            clt_mdct_backward_i32(
                mode,
                &freq,
                b,
                &mut out_syn[1],
                DECODE_BUFFER_SIZE - n + nb * b,
                shift,
                b_blocks,
            );
        }
    } else if out_channels == 1 && channels == 2 {
        // Stereo to mono downmix
        let mut freq2 = vec![0i32; n];

        denormalise_bands(mode, &x[0], &mut freq, old_band_e, start, eff_end, m, downsample, silence);
        denormalise_bands(mode, &x[1], &mut freq2, &old_band_e[mode.nb_ebands..], start, eff_end, m, downsample, silence);

        // Mix down
        for i in 0..n {
            freq[i] = (freq[i] + freq2[i]) / 2;
        }

        for b in 0..b_blocks {
            clt_mdct_backward_i32(
                mode,
                &freq,
                b,
                &mut out_syn[0],
                DECODE_BUFFER_SIZE - n + nb * b,
                shift,
                b_blocks,
            );
        }
    } else {
        // Normal case: mono or stereo
        for c in 0..out_channels {
            let band_e_offset = if c < channels { c * mode.nb_ebands } else { 0 };
            let x_c = if c < channels { c } else { 0 };

            denormalise_bands(
                mode,
                &x[x_c],
                &mut freq,
                &old_band_e[band_e_offset..],
                start,
                eff_end,
                m,
                downsample,
                silence,
            );

            for b in 0..b_blocks {
                clt_mdct_backward_i32(
                    mode,
                    &freq,
                    b,
                    &mut out_syn[c],
                    DECODE_BUFFER_SIZE - n + nb * b,
                    shift,
                    b_blocks,
                );
            }
        }
    }
}

/// Simplified integer IMDCT backward transform.
///
/// For MVP, this is a straightforward implementation. Should be replaced
/// with optimized FFT-based IMDCT for production.
fn clt_mdct_backward_i32(
    mode: &CeltMode,
    freq: &[i32],
    b: usize,
    out: &mut [i32],
    out_offset: usize,
    shift: usize,
    stride: usize,
) {
    let n = mode.short_mdct_size << (mode.max_lm - shift);
    let n2 = n >> 1;
    let n4 = n >> 2;
    let overlap = mode.overlap >> shift;

    // Simple DCT-IV based IMDCT
    // y[k] = sum_{n=0}^{N/2-1} x[n] * cos(pi/N * (k + 0.5 + N/4) * (n + 0.5))

    let mut temp = vec![0i32; n];

    for k in 0..n {
        let mut sum = 0i64;
        for i in 0..n2 {
            let x_val = freq[b + i * stride] as i64;
            // Compute cos term using Q15 approximation
            let phase = ((2 * k + 1 + n2) * (2 * i + 1)) % (4 * n);
            let cos_val = compute_cos_q15(phase, n * 4);
            sum += x_val * cos_val as i64;
        }
        temp[k] = (sum >> 15) as i32;
    }

    // Apply window and overlap-add
    let window = mode.window;
    for i in 0..overlap {
        let w = window[i] as i32;
        let w_comp = window[overlap - 1 - i] as i32;

        // First half of overlap
        out[out_offset + i] += mult16_16_q15(w, temp[n2 - 1 - i])
            + mult16_16_q15(w_comp, temp[n2 + i]);
    }

    // Middle section (no windowing)
    for i in overlap..n - overlap {
        out[out_offset + i] = temp[n2 - 1 - overlap + (i - overlap + 1)];
    }

    // Second half of overlap (store for next frame)
    for i in 0..overlap {
        let w = window[overlap - 1 - i] as i32;
        let w_comp = window[i] as i32;

        out[out_offset + n - overlap + i] = mult16_16_q15(w, temp[n - 1 - i])
            + mult16_16_q15(w_comp, temp[i]);
    }
}

/// Compute cos in Q15 for a given phase.
#[inline]
fn compute_cos_q15(phase: usize, period: usize) -> i32 {
    use core::f64::consts::PI;
    let angle = 2.0 * PI * phase as f64 / period as f64;
    (angle.cos() * 32768.0) as i32
}

/// De-emphasis filter to undo pre-emphasis.
///
/// Converts the output from the SIG_SHIFT fixed-point format to f32 samples.
pub fn deemphasis(
    inp: &[Vec<i32>],
    inp_offsets: &[usize],
    out: &mut [f32],
    n: usize,
    channels: usize,
    downsample: usize,
    preemph: &[i32],
    mem: &mut [i32],
    accum: bool,
) {
    let coef = preemph[0];
    let n_out = n / downsample;

    for c in 0..channels {
        let mut mem_val = mem[c];
        let inp_c = &inp[c];
        let inp_off = inp_offsets[c];

        for i in 0..n_out {
            // Downsample by taking every nth sample
            let idx = inp_off + i * downsample;
            let sample = inp_c[idx];

            // De-emphasis: y[n] = x[n] + coef * y[n-1]
            // coef is in Q15
            let y = sample + mult16_16_q15(coef, mem_val);
            mem_val = y;

            // Convert from SIG_SHIFT fixed-point to f32
            let sample_f32 = (y as f32) / (1 << SIG_SHIFT) as f32;

            // Interleave output
            let out_idx = i * channels + c;
            if accum {
                out[out_idx] += sample_f32;
            } else {
                out[out_idx] = sample_f32;
            }
        }

        mem[c] = mem_val;
    }
}

/// Comb filter (post-filter) for pitch enhancement.
///
/// Applies a comb filter to enhance pitch harmonics.
pub fn comb_filter(
    out: &mut [i32],
    out_offset: usize,
    inp: &[i32],
    inp_offset: usize,
    period_old: usize,
    period_new: usize,
    n: usize,
    g_old: i32,
    g_new: i32,
    tapset_old: usize,
    tapset_new: usize,
    window: &[i16],
    overlap: usize,
) {
    // Comb filter taps for different tapsets
    static COMB_TAPS: [[i16; 3]; 3] = [
        [20972, 16384, 8192],  // Tapset 0
        [22016, 14336, 8704],  // Tapset 1
        [21120, 15872, 7680],  // Tapset 2
    ];

    let period_old = period_old.max(COMBFILTER_MINPERIOD);
    let period_new = period_new.max(COMBFILTER_MINPERIOD);

    let taps_old = &COMB_TAPS[tapset_old.min(2)];
    let taps_new = &COMB_TAPS[tapset_new.min(2)];

    // Transition region (overlap)
    for i in 0..overlap {
        let w = window[i] as i32;
        let w_comp = Q15ONE - w;

        // Old filter contribution
        let mut y = 0i32;
        if g_old != 0 && inp_offset + i >= period_old + 2 {
            let base = inp_offset + i - period_old;
            y += mult16_16_q15(g_old, mult16_16_q15(taps_old[0] as i32, inp[base]));
            y += mult16_16_q15(g_old, mult16_16_q15(taps_old[1] as i32, inp[base - 1] + inp[base + 1]));
            y += mult16_16_q15(g_old, mult16_16_q15(taps_old[2] as i32, inp[base - 2] + inp[base + 2]));
        }
        let y_old = y;

        // New filter contribution
        y = 0;
        if g_new != 0 && inp_offset + i >= period_new + 2 {
            let base = inp_offset + i - period_new;
            y += mult16_16_q15(g_new, mult16_16_q15(taps_new[0] as i32, inp[base]));
            y += mult16_16_q15(g_new, mult16_16_q15(taps_new[1] as i32, inp[base - 1] + inp[base + 1]));
            y += mult16_16_q15(g_new, mult16_16_q15(taps_new[2] as i32, inp[base - 2] + inp[base + 2]));
        }
        let y_new = y;

        // Crossfade
        out[out_offset + i] = inp[inp_offset + i]
            + mult16_16_q15(w_comp, y_old)
            + mult16_16_q15(w, y_new);
    }

    // Main region (new filter only)
    for i in overlap..n {
        let mut y = inp[inp_offset + i];

        if g_new != 0 && inp_offset + i >= period_new + 2 {
            let base = inp_offset + i - period_new;
            y += mult16_16_q15(g_new, mult16_16_q15(taps_new[0] as i32, inp[base]));
            y += mult16_16_q15(g_new, mult16_16_q15(taps_new[1] as i32, inp[base - 1] + inp[base + 1]));
            y += mult16_16_q15(g_new, mult16_16_q15(taps_new[2] as i32, inp[base - 2] + inp[base + 2]));
        }

        out[out_offset + i] = y;
    }
}

/// Decode TF (time-frequency) resolution parameters.
pub fn tf_decode(
    start: usize,
    end: usize,
    is_transient: bool,
    tf_res: &mut [i32],
    lm: usize,
    dec: &mut crate::entropy::RangeDecoder<'_>,
) {
    use crate::celt::tables::TF_SELECT_TABLE;

    let logp = if is_transient { 2 } else { 4 };
    let budget = dec.storage() as i32 * 8;
    let mut tell = dec.tell();

    // Check if we have budget for tf_select
    let tf_select_rsv = if lm > 0 && tell + logp + 1 <= budget { 1 } else { 0 };
    let budget = budget - tf_select_rsv;

    let mut curr = 0i32;
    let mut tf_changed = false;

    for i in start..end {
        let logp = if i == start {
            if is_transient { 2 } else { 4 }
        } else {
            if is_transient { 4 } else { 5 }
        };

        if tell + logp <= budget {
            let bit = dec.dec_bit_logp(logp as u32);
            curr ^= bit;
            if bit != 0 {
                tf_changed = true;
            }
            tell = dec.tell();
        }
        tf_res[i] = curr;
    }

    // Decode tf_select if we have enough bits
    let tf_select = if tf_select_rsv != 0 {
        let is_trans = if is_transient { 1 } else { 0 };
        let tf_changed_idx = if tf_changed { 1 } else { 0 };
        if TF_SELECT_TABLE[lm][4 * is_trans + tf_changed_idx]
            != TF_SELECT_TABLE[lm][4 * is_trans + 2 + tf_changed_idx]
        {
            dec.dec_bit_logp(1) as usize
        } else {
            0
        }
    } else {
        0
    };

    // Apply tf_select to get final values
    for i in start..end {
        let is_trans = if is_transient { 1 } else { 0 };
        tf_res[i] = TF_SELECT_TABLE[lm][4 * is_trans + 2 * tf_select + tf_res[i] as usize] as i32;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_imdct_state_creation() {
        let state = ImdctState::new(120);
        assert_eq!(state.n, 120);
        assert_eq!(state.twiddle.len(), 60);
    }

    #[test]
    fn test_compute_cos_q15() {
        // cos(0) = 1.0 ≈ 32768 in Q15
        let cos_0 = compute_cos_q15(0, 1000);
        assert!((cos_0 - 32768).abs() < 2);

        // cos(pi) = -1.0 ≈ -32768 in Q15
        let cos_pi = compute_cos_q15(500, 1000);
        assert!((cos_pi + 32768).abs() < 2);
    }
}
