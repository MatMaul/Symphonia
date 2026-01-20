// Symphonia
// Copyright (c) 2024 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! Integer math utilities for Opus decoder.
//!
//! This module provides integer-only mathematical operations used throughout
//! the CELT decoder, including logarithms and square roots.

/// Computes the integer log base 2 of a value.
/// Returns floor(log2(x)). Only defined for strictly positive values.
///
/// # Panics
/// Panics if x <= 0 (in debug builds)
#[inline]
pub fn ilog2(x: i32) -> i32 {
    debug_assert!(x > 0, "ilog2() only defined for strictly positive numbers");
    ec_ilog(x as u32) as i32 - 1
}

/// Computes the integer log base 2, returning 0 for non-positive values.
#[inline]
pub fn zlog2(x: i32) -> i32 {
    if x <= 0 {
        0
    }
    else {
        ilog2(x)
    }
}

/// Computes the base-2 integer logarithm of a 32-bit value.
/// This is the number of bits required to represent the value.
///
/// Returns 1 for x == 0 (by convention for entropy coding).
#[inline]
fn ec_ilog(x: u32) -> u32 {
    if x == 0 {
        return 1;
    }

    // Count leading zeros and subtract from 32
    32 - x.leading_zeros()
}

/// Computes fixed-point log2 in Q10 format (10 fractional bits).
///
/// Returns approximately `1024 * log2(x)`.
pub fn celt_log2(x: i32) -> i32 {
    if x == 0 {
        return -32767;
    }

    let i = ilog2(x);
    let n = vshr32(x, i - 15) - 32768 - 16384;

    // Polynomial approximation for the fractional part
    const LOG2_C0: i32 = -6801 + (1 << 3);

    let frac = LOG2_C0
        + mult16_16_q15(
            n,
            15746
                + mult16_16_q15(
                    n,
                    -5217 + mult16_16_q15(n, 2545 + mult16_16_q15(n, -1401)),
                ),
        );

    shl(i - 13, 10) + shr(frac, 4)
}

/// Computes the fractional part of 2^x for x in Q7 format.
///
/// Returns the result in Q15 format.
fn celt_exp2_frac(x: i32) -> i32 {
    let frac = shl(x, 4);
    16383
        + mult16_16_q15(
            frac,
            22804 + mult16_16_q15(frac, 14819 + mult16_16_q15(10204, frac)),
        )
}

/// Computes 2^x where x is in Q7 format (7 fractional bits).
///
/// Returns the result as a 32-bit integer.
pub fn celt_exp2(x: i32) -> i32 {
    let integer = shr(x, 7);
    if integer < 0 {
        return 0;
    }
    else if integer >= 15 {
        return 0x7fffffff;
    }
    let frac = celt_exp2_frac(x & 0x007f);
    shl(extend32(frac), integer)
}

/// Variable-direction shift right.
/// Shifts right if shift > 0, left if shift < 0.
#[inline]
pub fn vshr32(a: i32, shift: i32) -> i32 {
    if shift > 0 {
        a >> shift
    }
    else {
        a << (-shift)
    }
}

/// Arithmetic shift right.
#[inline]
pub fn shr(a: i32, shift: i32) -> i32 {
    a >> shift
}

/// Logical shift left.
#[inline]
pub fn shl(a: i32, shift: i32) -> i32 {
    a << shift
}

/// Alias for shl (32-bit shift left).
#[inline]
pub fn shl32(a: i32, shift: i32) -> i32 {
    a << shift
}

/// Shift right with rounding.
#[inline]
pub fn pshr(a: i32, shift: i32) -> i32 {
    shr(a + (1 << shift >> 1), shift)
}

/// Sign-extend a 16-bit value to 32 bits.
#[inline]
pub fn extend32(x: i32) -> i32 {
    x
}

/// 16x16 multiplication with Q15 result (divide by 2^15).
#[inline]
pub fn mult16_16_q15(a: i32, b: i32) -> i32 {
    shr((a as i64 * b as i64) as i32, 15)
}

/// Saturate a value to 16-bit range.
#[inline]
pub fn sat16(x: i32) -> i16 {
    x.clamp(i16::MIN as i32, i16::MAX as i32) as i16
}

/// Minimum of two values.
#[inline]
pub fn min32(a: i32, b: i32) -> i32 {
    a.min(b)
}

/// Maximum of two values.
#[inline]
pub fn max32(a: i32, b: i32) -> i32 {
    a.max(b)
}

/// Integer log base 2 for CELT operations.
///
/// Returns the number of bits required to represent x.
/// For x == 0, returns 0.
#[inline]
pub fn celt_ilog2(x: u32) -> u32 {
    if x == 0 {
        0
    } else {
        31 - x.leading_zeros()
    }
}

/// Integer division with rounding towards zero.
#[inline]
pub fn celt_udiv(n: i32, d: i32) -> i32 {
    n / d
}

/// Compute approximate reciprocal square root using fixed-point.
///
/// Input is in Q16 format, output is in Q15 format.
/// Uses Newton-Raphson iteration.
pub fn celt_rsqrt_norm(x: i32) -> i32 {
    // rsqrt lookup table for initial approximation
    static RSQRT_TABLE: [i16; 32] = [
        32767, 31790, 30894, 30070, 29309, 28602, 27945, 27330,
        26755, 26214, 25705, 25225, 24770, 24339, 23930, 23541,
        23170, 22817, 22479, 22155, 21845, 21548, 21263, 20988,
        20724, 20470, 20225, 19988, 19760, 19539, 19326, 19119,
    ];

    debug_assert!(x > 0, "celt_rsqrt_norm() requires positive input");

    // Get initial estimate from table
    let k = (celt_ilog2(x as u32 - 1) >> 1) as usize;
    let k = k.min(31);

    let t = vshr32(x, 2 * k as i32 - 28);
    let t_idx = ((t - 128) >> 3).clamp(0, 31) as usize;
    let mut r = RSQRT_TABLE[t_idx] as i32;

    // Newton-Raphson iteration: r = r * (3 - x * r^2) / 2
    r = mult16_16_q15(r, (49152 - mult16_16_q15(mult16_16_q15(x >> (k as i32), r), r)));
    r = mult16_16_q15(r, (49152 - mult16_16_q15(mult16_16_q15(x >> (k as i32), r), r)));

    // Adjust for the scaling
    if k != 0 {
        r >>= k;
    }

    r.max(1)
}

/// Integer square root approximation.
#[inline]
pub fn celt_sqrt(x: i32) -> i32 {
    if x <= 0 {
        return 0;
    }

    // Newton-Raphson iteration
    let mut guess = 1i32 << (((32 - (x as u32).leading_zeros()) + 1) / 2);
    for _ in 0..4 {
        guess = (guess + x / guess) >> 1;
    }

    guess
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ilog2() {
        assert_eq!(ilog2(1), 0);
        assert_eq!(ilog2(2), 1);
        assert_eq!(ilog2(3), 1);
        assert_eq!(ilog2(4), 2);
        assert_eq!(ilog2(7), 2);
        assert_eq!(ilog2(8), 3);
        assert_eq!(ilog2(255), 7);
        assert_eq!(ilog2(256), 8);
        assert_eq!(ilog2(1024), 10);
    }

    #[test]
    fn test_zlog2() {
        assert_eq!(zlog2(0), 0);
        assert_eq!(zlog2(-5), 0);
        assert_eq!(zlog2(1), 0);
        assert_eq!(zlog2(2), 1);
        assert_eq!(zlog2(8), 3);
    }

    #[test]
    fn test_ec_ilog() {
        assert_eq!(ec_ilog(0), 1);
        assert_eq!(ec_ilog(1), 1);
        assert_eq!(ec_ilog(2), 2);
        assert_eq!(ec_ilog(3), 2);
        assert_eq!(ec_ilog(4), 3);
    }

    #[test]
    fn test_vshr32() {
        assert_eq!(vshr32(16, 2), 4);
        assert_eq!(vshr32(16, -2), 64);
        assert_eq!(vshr32(100, 0), 100);
    }

    #[test]
    fn test_mult16_16_q15() {
        // Q15: 1.0 = 32768
        // 0.5 * 0.5 = 0.25
        let half = 16384; // 0.5 in Q15
        let result = mult16_16_q15(half, half);
        assert!((result - 8192).abs() < 2); // Should be close to 0.25 in Q15
    }

    #[test]
    fn test_sat16() {
        assert_eq!(sat16(100), 100);
        assert_eq!(sat16(40000), 32767);
        assert_eq!(sat16(-40000), -32768);
        assert_eq!(sat16(32767), 32767);
        assert_eq!(sat16(-32768), -32768);
    }
}
