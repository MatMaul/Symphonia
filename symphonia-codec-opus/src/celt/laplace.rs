// Symphonia
// Copyright (c) 2024 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! Laplace distribution decoder for CELT coarse energy.
//!
//! The Laplace distribution is used for efficient entropy coding of
//! coarse energy values which tend to cluster around zero.

use crate::entropy::RangeDecoder;

/// Constants for Laplace coding
const LAPLACE_LOG_MINP: i32 = 0;
const LAPLACE_MINP: i32 = 1 << LAPLACE_LOG_MINP;
const LAPLACE_NMIN: i32 = 16;

/// Compute the frequency of the first symbol after zero.
///
/// Given the probability of zero (fs0) and decay factor, compute the
/// probability of ±1.
#[inline]
fn ec_laplace_get_freq1(fs0: u32, decay: i32) -> u32 {
    let ft = cap_to_u32(32768 - LAPLACE_MINP * (2 * LAPLACE_NMIN) - fs0 as i32);
    cap_to_u32((ft as i32 * (16384 - decay)) >> 15)
}

/// Decode a value from a Laplace distribution.
///
/// The Laplace distribution has probability concentrated around zero,
/// with exponentially decaying probability for larger magnitudes.
///
/// # Arguments
/// * `dec` - Range decoder
/// * `fs` - Initial probability of zero (in Q15)
/// * `decay` - Decay factor for probability of larger values
///
/// # Returns
/// The decoded value (can be positive, negative, or zero)
pub fn ec_laplace_decode(dec: &mut RangeDecoder<'_>, fs: u32, decay: i32) -> i32 {
    let mut val = 0i32;
    let fm = dec.decode_bin(15);
    let mut fl: u32 = 0;
    let mut fs = fs;

    if fm >= fs as u32 {
        val += 1;
        fl = fs;
        fs = ec_laplace_get_freq1(fs, decay) + LAPLACE_MINP as u32;

        // Search the decaying part of the PDF
        while fs > LAPLACE_MINP as u32 && fm >= fl + 2 * fs {
            fs *= 2;
            fl = cap_to_u32(fl as i32 + fs as i32);
            fs = cap_to_u32(((fs as i32 - 2 * LAPLACE_MINP) * decay) >> 15) + LAPLACE_MINP as u32;
            val += 1;
        }

        // Everything beyond that has probability LAPLACE_MINP
        if fs <= LAPLACE_MINP as u32 {
            let di = ((fm - fl) >> (LAPLACE_LOG_MINP + 1)) as i32;
            val += di;
            fl = cap_to_u32(fl as i32 + cap_to_u32(2 * di * LAPLACE_MINP) as i32);
        }

        // Determine sign
        if fm < fl + fs {
            val = -val;
        } else {
            fl = cap_to_u32(fl as i32 + fs as i32);
        }
    }

    debug_assert!(fl < 32768);
    debug_assert!(fs > 0);
    debug_assert!(fl as u32 <= fm);

    dec.dec_update(fl as u32, (fl + fs).min(32768), 32768);
    val
}

/// Clamp a value to u32, treating negative values as 0.
#[inline]
fn cap_to_u32(val: i32) -> u32 {
    val.max(0) as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ec_laplace_get_freq1() {
        // Test with typical parameters
        let fs0 = 72u32 << 7; // probability of zero
        let decay = 127 << 6; // decay factor
        let freq1 = ec_laplace_get_freq1(fs0, decay);
        // Frequency should be positive and less than remaining probability
        assert!(freq1 > 0);
        assert!(freq1 < 32768 - fs0);
    }

    #[test]
    fn test_cap_to_u32() {
        assert_eq!(cap_to_u32(100), 100);
        assert_eq!(cap_to_u32(0), 0);
        assert_eq!(cap_to_u32(-1), 0);
        assert_eq!(cap_to_u32(-1000), 0);
    }
}
