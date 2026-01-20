// Symphonia
// Copyright (c) 2024 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! Band energy quantization for CELT.
//!
//! This module handles the decoding of coarse and fine energy values
//! for each frequency band.

use crate::celt::constants::DB_SHIFT;
use crate::celt::laplace::ec_laplace_decode;
use crate::celt::mode::{CeltMode, MAX_FINE_BITS};
use crate::celt::tables::{BETA_COEF, BETA_INTRA, E_PROB_MODEL, PRED_COEF, SMALL_ENERGY_ICDF};
use crate::entropy::RangeDecoder;
use crate::util::math::{pshr, shl32};

/// Maximum negative energy value (in Q10 dB).
const MINUS_9DB: i32 = -(((9.0 * (1 << DB_SHIFT) as f64) + 0.5) as i32);
const MINUS_28DB: i32 = -(((28.0 * (1 << DB_SHIFT) as f64) + 0.5) as i32);

/// Decode coarse energy for all bands.
///
/// Coarse energy is coded using Laplace distribution for most bands,
/// with fallback to simpler coding for very low bit rates.
///
/// # Arguments
/// * `m` - CELT mode configuration
/// * `start` - First band to decode
/// * `end` - Last band + 1 to decode
/// * `old_ebands` - Output array for decoded energies (also used for prediction)
/// * `intra` - Whether to use intra (non-predictive) mode
/// * `dec` - Range decoder
/// * `channels` - Number of channels (1 or 2)
/// * `lm` - Log2 of frame size multiplier (0=2.5ms, 1=5ms, 2=10ms, 3=20ms)
pub fn unquant_coarse_energy(
    m: &CeltMode,
    start: usize,
    end: usize,
    old_ebands: &mut [i32],
    intra: bool,
    dec: &mut RangeDecoder<'_>,
    channels: usize,
    lm: usize,
) {
    let prob_model = &E_PROB_MODEL[lm][if intra { 1 } else { 0 }];
    let mut prev = [0i32; 2];

    let (coef, beta) = if intra { (0, BETA_INTRA) } else { (PRED_COEF[lm], BETA_COEF[lm]) };

    let budget = dec.storage() as i32 * 8;

    for i in start..end {
        for c in 0..channels {
            let tell = dec.tell();
            let qi: i32;

            if budget - tell >= 15 {
                // Use Laplace coding
                let pi = 2 * i.min(20);
                qi = ec_laplace_decode(
                    dec,
                    (prob_model[pi] as u32) << 7,
                    (prob_model[pi + 1] as i32) << 6,
                );
            }
            else if budget - tell >= 2 {
                // Use small energy coding
                let idx = dec.dec_icdf(SMALL_ENERGY_ICDF, 2);
                qi = (idx >> 1) ^ -(idx & 1);
            }
            else if budget - tell >= 1 {
                // Single bit
                qi = -dec.dec_bit_logp(1);
            }
            else {
                // No bits left, assume -1
                qi = -1;
            }

            let q = shl32(qi, DB_SHIFT);

            // Apply prediction and update state
            let idx = i + c * m.nb_ebands;
            old_ebands[idx] = old_ebands[idx].max(MINUS_9DB);

            let tmp = pshr(mult16_16(coef, old_ebands[idx]), 8) + prev[c] + shl32(q, 7);
            let tmp = tmp.max(MINUS_28DB << 7);

            old_ebands[idx] = pshr(tmp, 7);
            prev[c] = prev[c] + shl32(q, 7) - mult16_16(beta, pshr(q, 8));
        }
    }
}

/// Decode fine energy for all bands.
///
/// Fine energy provides additional precision beyond coarse energy.
///
/// # Arguments
/// * `m` - CELT mode configuration
/// * `start` - First band to decode
/// * `end` - Last band + 1 to decode
/// * `old_ebands` - Energy array to update
/// * `fine_quant` - Number of fine bits per band
/// * `dec` - Range decoder
/// * `channels` - Number of channels
pub fn unquant_fine_energy(
    m: &CeltMode,
    start: usize,
    end: usize,
    old_ebands: &mut [i32],
    fine_quant: &[i32],
    dec: &mut RangeDecoder<'_>,
    channels: usize,
) {
    for i in start..end {
        if fine_quant[i] <= 0 {
            continue;
        }

        for c in 0..channels {
            let q2 = dec.dec_bits(fine_quant[i]);
            let offset = shr32(shl32(q2, DB_SHIFT) + half_db(), fine_quant[i]) - half_db();

            let idx = i + c * m.nb_ebands;
            old_ebands[idx] += offset;
        }
    }
}

/// Decode final fine energy bits (leftover bits).
///
/// After all bands have been coded, any remaining bits are used to
/// improve the precision of fine energy values.
///
/// # Arguments
/// * `m` - CELT mode configuration
/// * `start` - First band
/// * `end` - Last band + 1
/// * `old_ebands` - Energy array to update
/// * `fine_quant` - Number of fine bits per band
/// * `fine_priority` - Priority for extra bits (0 or 1)
/// * `bits_left` - Number of remaining bits
/// * `dec` - Range decoder
/// * `channels` - Number of channels
pub fn unquant_energy_finalise(
    m: &CeltMode,
    start: usize,
    end: usize,
    old_ebands: &mut [i32],
    fine_quant: &[i32],
    fine_priority: &[i32],
    mut bits_left: i32,
    dec: &mut RangeDecoder<'_>,
    channels: usize,
) {
    let c = channels as i32;

    for prio in 0..2 {
        for i in start..end {
            if bits_left < c {
                break;
            }
            if fine_quant[i] >= MAX_FINE_BITS || fine_priority[i] != prio {
                continue;
            }

            for ch in 0..channels {
                let q2 = dec.dec_bits(1);
                let offset = (shl32(q2, DB_SHIFT) - half_db()) >> (fine_quant[i] + 1);

                let idx = i + ch * m.nb_ebands;
                old_ebands[idx] += offset;
                bits_left -= 1;
            }
        }
    }
}

/// Helper: 0.5 in DB_SHIFT fixed point.
#[inline]
fn half_db() -> i32 {
    1 << (DB_SHIFT - 1)
}

/// Helper: Multiply two 16-bit values.
#[inline]
fn mult16_16(a: i32, b: i32) -> i32 {
    a * b
}

/// Helper: Shift right with rounding.
#[inline]
fn shr32(a: i32, shift: i32) -> i32 {
    if shift >= 0 { a >> shift } else { a << -shift }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constants() {
        // Verify the Q10 fixed-point constants
        assert!(MINUS_9DB < 0);
        assert!(MINUS_28DB < MINUS_9DB);
    }

    #[test]
    fn test_half_db() {
        assert_eq!(half_db(), 1 << (DB_SHIFT - 1));
    }
}
