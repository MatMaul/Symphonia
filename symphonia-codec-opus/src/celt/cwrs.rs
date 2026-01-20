// Symphonia
// Copyright (c) 2024 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! Combinatorial Weighted Reed-Solomon (CWRS) codec for PVQ.
//!
//! This module implements the CWRS codec used for encoding and decoding
//! pyramid vector quantization (PVQ) codewords in CELT.

use crate::celt::tables::{CELT_PVQ_U_DATA, CELT_PVQ_U_ROW};
use crate::entropy::RangeDecoder;

/// Look up U(n, k) in the PVQ table.
///
/// U(n, k) represents the number of (n,k) combinations in the PVQ lattice.
#[inline]
pub fn celt_pvq_u(n: usize, k: usize) -> u64 {
    let row = n.min(k);
    let col = n.max(k);
    CELT_PVQ_U_DATA[CELT_PVQ_U_ROW[row] + col] as u64
}

/// Compute V(n, k) = U(n, k) + U(n, k+1).
///
/// V(n, k) is the total number of PVQ codewords for dimension n with k pulses.
#[inline]
pub fn celt_pvq_v(n: usize, k: usize) -> u64 {
    celt_pvq_u(n, k) + celt_pvq_u(n, k + 1)
}

/// Helper to cap a value to u32, treating negative as 0.
#[inline]
fn cap_to_u32(val: i64) -> u64 {
    if val < 0 { 0 } else { val as u64 }
}

/// Decode a PVQ codeword using the CWRS algorithm.
///
/// Given the encoded index `i`, dimension `n`, and number of pulses `k`,
/// this function decodes the pulse vector into `y`.
///
/// # Arguments
/// * `n` - Number of dimensions
/// * `k` - Number of pulses (sum of absolute values)
/// * `i` - Encoded index in range [0, V(n,k))
/// * `y` - Output vector of pulse values (length n)
///
/// # Returns
/// The sum of squared values (Ryy = sum of y[j]^2)
pub fn cwrsi(n: usize, mut k: usize, mut i: u64, y: &mut [i32]) -> i32 {
    debug_assert!(k > 0, "cwrsi() needs at least one pulse");
    debug_assert!(n > 1, "cwrsi() needs at least two dimensions");

    let mut yy = 0i32;
    let mut y_ptr = 0usize;
    let mut n_remaining = n;

    while n_remaining > 2 {
        let p: u64;
        let s: i32;
        let k0: usize;
        let val: i32;

        if k >= n_remaining {
            // k >= n case: use row n
            let row = CELT_PVQ_U_ROW[n_remaining];
            p = CELT_PVQ_U_DATA[row + k + 1] as u64;

            // Determine sign
            s = if i >= p { -1 } else { 0 };
            i -= cap_to_u32((p as i64) & (s as i64));

            k0 = k;
            let q = CELT_PVQ_U_DATA[row + n_remaining] as u64;

            if q > i {
                debug_assert!(p > q);
                k = n_remaining;
                loop {
                    k -= 1;
                    let p_new = CELT_PVQ_U_DATA[CELT_PVQ_U_ROW[k] + n_remaining] as u64;
                    if p_new <= i {
                        i -= p_new;
                        break;
                    }
                }
            }
            else {
                let mut p_search = CELT_PVQ_U_DATA[row + k] as u64;
                while p_search > i {
                    k -= 1;
                    p_search = CELT_PVQ_U_DATA[row + k] as u64;
                }
                i -= p_search;
            }

            val = ((k0 as i32 - k as i32 + s) ^ s);
            y[y_ptr] = val;
            y_ptr += 1;
            yy += val * val;
        }
        else {
            // k < n case
            let p0 = CELT_PVQ_U_DATA[CELT_PVQ_U_ROW[k] + n_remaining] as u64;
            let q = CELT_PVQ_U_DATA[CELT_PVQ_U_ROW[k + 1] + n_remaining] as u64;

            if p0 <= i && i < q {
                // Zero coefficient
                i -= p0;
                y[y_ptr] = 0;
                y_ptr += 1;
            }
            else {
                // Non-zero coefficient
                s = if i >= q { -1 } else { 0 };
                i -= cap_to_u32((q as i64) & (s as i64));

                k0 = k;
                loop {
                    k -= 1;
                    let p_new = CELT_PVQ_U_DATA[CELT_PVQ_U_ROW[k] + n_remaining] as u64;
                    if p_new <= i {
                        i -= p_new;
                        break;
                    }
                }

                val = ((k0 as i32 - k as i32 + s) ^ s);
                y[y_ptr] = val;
                y_ptr += 1;
                yy += val * val;
            }
        }
        n_remaining -= 1;
    }

    // Handle the last two coefficients (n_remaining == 2)
    let p = (2 * k + 1) as u64;
    let s = if i >= p { -1i32 } else { 0 };
    i -= cap_to_u32((p as i64) & (s as i64));

    let k0 = k;
    k = ((i + 1) >> 1) as usize;
    if k != 0 {
        i -= (2 * k - 1) as u64;
    }

    let val = ((k0 as i32 - k as i32 + s) ^ s);
    y[y_ptr] = val;
    y_ptr += 1;
    yy += val * val;

    // Last coefficient
    let s_last = -(i as i32);
    let val_last = ((k as i32 + s_last) ^ s_last);
    y[y_ptr] = val_last;
    yy += val_last * val_last;

    yy
}

/// Decode pulses from the range decoder.
///
/// # Arguments
/// * `y` - Output pulse vector
/// * `n` - Number of dimensions
/// * `k` - Number of pulses
/// * `dec` - Range decoder
///
/// # Returns
/// Sum of squared pulse values (Ryy)
pub fn decode_pulses(y: &mut [i32], n: usize, k: usize, dec: &mut RangeDecoder<'_>) -> i32 {
    let ft = celt_pvq_v(n, k);
    let i = dec.dec_uint(ft as u32) as u64;
    cwrsi(n, k, i, y)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_celt_pvq_u() {
        // Row 0 (k=0): only index 0 is 1, rest are 0
        assert_eq!(celt_pvq_u(0, 0), 1);
        assert_eq!(celt_pvq_u(0, 5), 0); // min=0, max=5 -> row 0, col 5 = 0

        // Row 1 (k=1): U(n,1) = 1 for all n >= 1
        assert_eq!(celt_pvq_u(1, 1), 1);
        assert_eq!(celt_pvq_u(1, 10), 1);

        // Row 2 (k=2): U(n,2) = 2n-1 for n >= 2
        // U(2,2) = 3, U(3,2) = 5, U(4,2) = 7, ...
        assert_eq!(celt_pvq_u(2, 2), 3);
        assert_eq!(celt_pvq_u(3, 2), 5);
        assert_eq!(celt_pvq_u(4, 2), 7);
        assert_eq!(celt_pvq_u(5, 2), 9);
    }

    #[test]
    fn test_celt_pvq_v() {
        // V(n, k) = U(n, k) + U(n, k+1)
        let v = celt_pvq_v(2, 2);
        assert_eq!(v, celt_pvq_u(2, 2) + celt_pvq_u(2, 3));
    }

    #[test]
    fn test_cwrsi_simple() {
        // Simple test with n=2, k=1: there are V(2,1) = U(2,1) + U(2,2) codewords
        // For n=2, k=1, the valid codewords are: [1,0], [-1,0], [0,1], [0,-1]
        // i=0 should give the first one
        let mut y = [0i32; 2];
        let ryy = cwrsi(2, 1, 0, &mut y);

        // Sum of absolute values should be k=1
        let sum: i32 = y.iter().map(|&x| x.abs()).sum();
        assert_eq!(sum, 1);
        assert_eq!(ryy, 1); // 1^2 = 1
    }
}
