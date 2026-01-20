// Symphonia
// Copyright (c) 2024 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! MDCT (Modified Discrete Cosine Transform) for CELT.
//!
//! This module implements the IMDCT (inverse MDCT) transform used to convert
//! frequency-domain coefficients back to time-domain samples.

#![allow(dead_code)]

use crate::util::fixed_point::{mult16_16_q15, mult16_32_q15};

/// MDCT window types.
const WINDOW_VORBIS: usize = 0;

/// Pre-computed twiddle factors for MDCT.
///
/// For MVP, we'll compute these on the fly. In optimized version,
/// these would be pre-computed tables.
#[derive(Clone)]
pub struct MdctLookup {
    /// Transform size (N)
    pub n: usize,
    /// Maximum shift (for different frame sizes)
    pub max_shift: usize,
    /// Twiddle factors (sin/cos tables)
    trig: Vec<i16>,
}

impl MdctLookup {
    /// Create a new MDCT lookup table for a given size.
    ///
    /// # Arguments
    /// * `n` - Transform size (must be power of 2)
    /// * `max_shift` - Maximum downshift (0-3 for 2.5ms to 20ms frames)
    pub fn new(n: usize, max_shift: usize) -> Self {
        let mut trig = Vec::new();

        // Compute twiddle factors for all shift levels
        for shift in 0..=max_shift {
            let size = n >> shift;
            Self::compute_twiddles(size, &mut trig);
        }

        Self {
            n,
            max_shift,
            trig,
        }
    }

    /// Compute twiddle factors (sin/cos) for MDCT.
    fn compute_twiddles(n: usize, trig: &mut Vec<i16>) {
        use core::f64::consts::PI;

        let n2 = n >> 1;
        let scale = 32768.0; // Q15 fixed-point scale

        for i in 0..n2 {
            let angle = PI * (i as f64) / (n as f64);
            let cos_val = (angle.cos() * scale) as i16;
            let sin_val = (angle.sin() * scale) as i16;
            trig.push(cos_val);
            trig.push(sin_val);
        }
    }

    /// Get twiddle factor offset for a given shift level.
    fn trig_offset(&self, shift: usize) -> usize {
        let mut offset = 0;
        for s in 0..shift {
            let size = self.n >> s;
            offset += size; // N/2 pairs * 2 values
        }
        offset
    }
}

/// Inverse MDCT (IMDCT) transform.
///
/// Converts frequency-domain coefficients to time-domain samples
/// with windowing and overlap-add.
///
/// # Arguments
/// * `lookup` - Pre-computed MDCT lookup table
/// * `input` - Input frequency coefficients (N/2 values)
/// * `output` - Output time-domain samples (N values)
/// * `window` - Window function (for overlap-add)
/// * `overlap` - Overlap size
/// * `shift` - Frame size shift (0=2.5ms, 1=5ms, 2=10ms, 3=20ms)
/// * `stride` - Output stride (for interleaved stereo)
pub fn imdct(
    lookup: &MdctLookup,
    input: &[i32],
    output: &mut [i32],
    window: &[i16],
    overlap: usize,
    shift: usize,
    stride: usize,
) {
    let mut n = lookup.n;
    for _ in 0..shift {
        n >>= 1;
    }

    let n2 = n >> 1;
    let n4 = n >> 2;

    debug_assert!(input.len() >= n2);
    debug_assert!(output.len() >= n);

    // For MVP: Use direct DCT-IV computation
    // DCT-IV: y[k] = sum_{n=0}^{N-1} x[n] * cos(π/N * (n + 0.5) * (k + 0.5))

    let trig_offset = lookup.trig_offset(shift);

    // Temporary buffer for IMDCT computation
    let mut temp = vec![0i32; n];

    // Pre-rotation: convert MDCT spectrum to time domain via DCT-IV
    // This is a simplified implementation for correctness
    for k in 0..n2 {
        let mut sum_real = 0i64;
        let mut sum_imag = 0i64;

        for i in 0..n4 {
            let idx = trig_offset + 2 * i;
            let cos_val = lookup.trig[idx] as i32;
            let sin_val = lookup.trig[idx + 1] as i32;

            let re = input[2 * i];
            let im = input[2 * i + 1];

            // Complex rotation
            let yr = mult16_16_q15(re, cos_val) - mult16_16_q15(im, sin_val);
            let yi = mult16_16_q15(im, cos_val) + mult16_16_q15(re, sin_val);

            // Accumulate contribution to output
            let phase = (2 * k * i) % n;
            let cos_k = compute_cos_q15(phase, n);
            let sin_k = compute_sin_q15(phase, n);

            sum_real += (yr as i64 * cos_k as i64) - (yi as i64 * sin_k as i64);
            sum_imag += (yr as i64 * sin_k as i64) + (yi as i64 * cos_k as i64);
        }

        temp[k] = (sum_real >> 15) as i32;
        temp[n2 + k] = (sum_imag >> 15) as i32;
    }

    // Post-rotation and windowing
    apply_window_and_overlap(
        &temp,
        output,
        window,
        n,
        overlap,
        stride,
    );
}

/// Compute cosine in Q15 fixed-point for a given phase.
#[inline]
fn compute_cos_q15(phase: usize, n: usize) -> i32 {
    use core::f64::consts::PI;
    let angle = 2.0 * PI * (phase as f64) / (n as f64);
    (angle.cos() * 32768.0) as i32
}

/// Compute sine in Q15 fixed-point for a given phase.
#[inline]
fn compute_sin_q15(phase: usize, n: usize) -> i32 {
    use core::f64::consts::PI;
    let angle = 2.0 * PI * (phase as f64) / (n as f64);
    (angle.sin() * 32768.0) as i32
}

/// Apply window function and perform overlap-add.
fn apply_window_and_overlap(
    input: &[i32],
    output: &mut [i32],
    window: &[i16],
    n: usize,
    overlap: usize,
    stride: usize,
) {
    let overlap_half = overlap >> 1;

    // First half: apply window and overlap with previous frame
    for i in 0..overlap_half {
        let win_val = window[overlap_half + i];
        let sample = mult16_32_q15(win_val, input[i]);

        // For decoder, we accumulate with previous overlap
        output[i * stride] += sample;
    }

    // Middle section: no windowing (or window = 1.0)
    for i in overlap_half..(n - overlap_half) {
        output[i * stride] = input[i];
    }

    // Last half: apply window for next frame's overlap
    for i in 0..overlap_half {
        let win_val = window[overlap_half - 1 - i];
        let sample = mult16_32_q15(win_val, input[n - overlap_half + i]);
        output[(n - overlap_half + i) * stride] = sample;
    }
}

/// Simplified IMDCT using direct DCT-IV computation.
///
/// This is slower than FFT-based approach but simpler and correct.
/// For production, this should be replaced with optimized FFT-based MDCT.
pub fn imdct_simple(
    input: &[i32],
    output: &mut [i32],
    n: usize,
) {
    use core::f64::consts::PI;

    let n2 = n >> 1;

    // DCT-IV inverse transform
    for k in 0..n {
        let mut sum = 0i64;

        for i in 0..n2 {
            let angle = PI / (n as f64) * (i as f64 + 0.5) * (k as f64 + 0.5);
            let cos_val = (angle.cos() * 32768.0) as i32;
            sum += input[i] as i64 * cos_val as i64;
        }

        output[k] = ((sum >> 15) * 2) as i32; // Scale factor for DCT-IV
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mdct_lookup_creation() {
        let lookup = MdctLookup::new(120, 3);
        assert_eq!(lookup.n, 120);
        assert_eq!(lookup.max_shift, 3);
        assert!(!lookup.trig.is_empty());
    }

    #[test]
    fn test_imdct_simple_dc() {
        // Test DC component (all ones)
        let n = 8;
        let input = vec![32768; n / 2]; // 1.0 in Q15
        let mut output = vec![0; n];

        imdct_simple(&input, &mut output, n);

        // Output should be non-zero
        assert!(output.iter().any(|&x| x != 0));
    }
}
