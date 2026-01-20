// Symphonia
// Copyright (c) 2024 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! Bit allocation for CELT bands.
//!
//! This module implements the rate control and bit allocation algorithm
//! used to distribute bits among frequency bands.

use crate::celt::constants::BITRES;
use crate::celt::mode::{CeltMode, MAX_FINE_BITS};
use crate::celt::tables::LOG2_FRAC_TABLE;
use crate::entropy::RangeDecoder;

/// Number of interpolation steps for bit allocation.
const ALLOC_STEPS: i32 = 6;

/// Convert pulse count to number of pulses.
///
/// For small values, the mapping is linear. For larger values,
/// it's exponential to allow for large pulse counts.
#[inline]
pub fn get_pulses(i: i32) -> i32 {
    if i < 8 {
        i
    } else {
        (8 + (i & 7)) << ((i >> 3) - 1)
    }
}

/// Convert bits to optimal number of pulses for a band.
///
/// Uses the pulse cache to find the number of pulses that best
/// matches the available bits.
pub fn bits2pulses(m: &CeltMode, band: usize, lm: i32, bits: i32) -> i32 {
    let lm = (lm + 1) as usize;
    let cache = m.cache.bits;
    let cache_ptr = m.cache.index[lm * m.nb_ebands + band] as usize;

    let mut lo = 0i32;
    let mut hi = cache[cache_ptr] as i32;
    let bits = bits - 1;

    // Binary search for the best pulse count
    for _ in 0..crate::celt::constants::LOG_MAX_PSEUDO {
        let mid = (lo + hi + 1) >> 1;
        if cache[cache_ptr + mid as usize] as i32 >= bits {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    let low_val = if lo == 0 { -1 } else { cache[cache_ptr + lo as usize] as i32 };

    if bits - low_val <= cache[cache_ptr + hi as usize] as i32 - bits {
        lo
    } else {
        hi
    }
}

/// Convert pulse count to bits required.
pub fn pulses2bits(m: &CeltMode, band: usize, lm: i32, pulses: i32) -> i32 {
    let lm = (lm + 1) as usize;
    if pulses == 0 {
        return 0;
    }
    m.cache.bits[m.cache.index[lm * m.nb_ebands + band] as usize + pulses as usize] as i32 + 1
}

/// Interpolate bit allocation and convert to pulses.
///
/// This function performs the core bit allocation algorithm, distributing
/// bits among bands based on the static allocation tables and dynamic
/// parameters.
#[allow(clippy::too_many_arguments)]
pub fn interp_bits2pulses(
    m: &CeltMode,
    start: usize,
    end: usize,
    skip_start: usize,
    bits1: &[i32],
    bits2: &[i32],
    thresh: &[i32],
    cap: &[i32],
    total: i32,
    balance: &mut i32,
    skip_rsv: i32,
    intensity: &mut i32,
    intensity_rsv: i32,
    dual_stereo: &mut i32,
    dual_stereo_rsv: i32,
    bits: &mut [i32],
    ebits: &mut [i32],
    fine_priority: &mut [i32],
    channels: usize,
    lm: i32,
    dec: &mut RangeDecoder<'_>,
    prev: i32,
    signal_bandwidth: i32,
) -> usize {
    let c = channels as i32;
    let alloc_floor = c << BITRES;
    let stereo = if c > 1 { 1 } else { 0 };
    let log_m = lm << BITRES;
    let mut total = total;
    let mut skip_rsv = skip_rsv;
    let mut intensity_rsv = intensity_rsv;
    let mut dual_stereo_rsv = dual_stereo_rsv;
    let mut skip_start = skip_start;

    // Binary search for the right allocation
    let mut lo = 0i32;
    let mut hi = 1 << ALLOC_STEPS;

    for _ in 0..ALLOC_STEPS {
        let mid = (lo + hi) >> 1;
        let mut psum = 0;
        let mut done = false;

        for j in (start..end).rev() {
            let tmp = bits1[j] + ((mid * bits2[j]) >> ALLOC_STEPS);
            if tmp >= thresh[j] || done {
                done = true;
                psum += tmp.min(cap[j]);
            } else if tmp >= alloc_floor {
                psum += alloc_floor;
            }
        }

        if psum > total {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    // Compute final allocation with the found interpolation factor
    let mut psum = 0;
    let mut done = false;

    for j in (start..end).rev() {
        let mut tmp = bits1[j] + ((lo * bits2[j]) >> ALLOC_STEPS);
        if tmp < thresh[j] && !done {
            if tmp >= alloc_floor {
                tmp = alloc_floor;
            } else {
                tmp = 0;
            }
        } else {
            done = true;
        }
        tmp = tmp.min(cap[j]);
        bits[j] = tmp;
        psum += tmp;
    }

    // Determine which bands to skip (band folding)
    let mut coded_bands = end;

    loop {
        let j = coded_bands - 1;
        if j <= skip_start {
            total += skip_rsv;
            break;
        }

        let left = total - psum;
        let percoeff = celt_udiv(left, (m.ebands[coded_bands] - m.ebands[start]) as i32);
        let left_after = left - (m.ebands[coded_bands] - m.ebands[start]) as i32 * percoeff;
        let rem = (left_after - (m.ebands[j] - m.ebands[start]) as i32).max(0);
        let band_width = (m.ebands[coded_bands] - m.ebands[j]) as i32;
        let band_bits = bits[j] + percoeff * band_width + rem;

        if band_bits >= thresh[j].max(alloc_floor + (1 << BITRES)) {
            // Decode skip bit
            let threshold = if j < prev as usize { 7 } else { 9 };
            if coded_bands <= start + 2
                || (band_bits > (threshold * band_width << lm << BITRES) >> 4
                    && j as i32 <= signal_bandwidth)
            {
                if dec.dec_bit_logp(1) != 0 {
                    break;
                }
            } else if dec.dec_bit_logp(1) != 0 {
                break;
            }
            psum += 1 << BITRES;
        }

        psum -= bits[j] + intensity_rsv;
        if intensity_rsv > 0 {
            intensity_rsv = LOG2_FRAC_TABLE[j - start] as i32;
        }
        psum += intensity_rsv;

        if band_bits >= alloc_floor {
            psum += alloc_floor;
            bits[j] = alloc_floor;
        } else {
            bits[j] = 0;
        }

        coded_bands -= 1;
    }

    debug_assert!(coded_bands > start);

    // Decode intensity stereo position
    if intensity_rsv > 0 {
        *intensity = start as i32 + dec.dec_uint((coded_bands + 1 - start) as u32) as i32;
    } else {
        *intensity = 0;
    }

    // Decode dual stereo flag
    if *intensity <= start as i32 {
        total += dual_stereo_rsv;
        dual_stereo_rsv = 0;
    }
    if dual_stereo_rsv > 0 {
        *dual_stereo = dec.dec_bit_logp(1);
    } else {
        *dual_stereo = 0;
    }

    // Distribute remaining bits
    let left = total - psum;
    let percoeff = celt_udiv(left, (m.ebands[coded_bands] - m.ebands[start]) as i32);
    let mut left = left - (m.ebands[coded_bands] - m.ebands[start]) as i32 * percoeff;

    for j in start..coded_bands {
        bits[j] += percoeff * (m.ebands[j + 1] - m.ebands[j]) as i32;
    }

    for j in start..coded_bands {
        let tmp = left.min((m.ebands[j + 1] - m.ebands[j]) as i32);
        bits[j] += tmp;
        left -= tmp;
    }

    // Compute fine energy bits
    *balance = 0;

    for j in start..coded_bands {
        debug_assert!(bits[j] >= 0);
        let n0 = (m.ebands[j + 1] - m.ebands[j]) as i32;
        let n = n0 << lm;
        let bit = bits[j] + *balance;

        if n > 1 {
            let excess = (bit - cap[j]).max(0);
            bits[j] = bit - excess;

            let mut den = c * n;
            if c == 2 && n > 2 && *dual_stereo == 0 && j < *intensity as usize {
                den += 1;
            }

            let nc_log_n = den * (m.log_n[j] as i32 + log_m);
            let mut offset = (nc_log_n >> 1) - den * crate::celt::constants::FINE_OFFSET;

            if n == 2 {
                offset += (den << BITRES) >> 2;
            }

            if bits[j] + offset < den * 2 << BITRES {
                offset += nc_log_n >> 2;
            } else if bits[j] + offset < den * 3 << BITRES {
                offset += nc_log_n >> 3;
            }

            ebits[j] = ((bits[j] + offset + (den << (BITRES - 1))).max(0) / den) >> BITRES;

            if c * ebits[j] > bits[j] >> BITRES {
                ebits[j] = bits[j] >> stereo >> BITRES;
            }

            ebits[j] = ebits[j].min(MAX_FINE_BITS);

            if ebits[j] * (den << BITRES) >= bits[j] + offset {
                fine_priority[j] = 1;
            } else {
                fine_priority[j] = 0;
            }

            bits[j] -= c * ebits[j] << BITRES;
        } else {
            let excess = (bit - (c << BITRES)).max(0);
            bits[j] = bit - excess;
            ebits[j] = 0;
            fine_priority[j] = 1;
        }

        // Handle excess bits
        let mut excess = if n > 1 {
            (bit - cap[j]).max(0)
        } else {
            (bit - (c << BITRES)).max(0)
        };

        if excess > 0 {
            let extra_fine = (excess >> (stereo + BITRES)).min(MAX_FINE_BITS - ebits[j]);
            ebits[j] += extra_fine;
            let extra_bits = extra_fine * c << BITRES;
            if extra_bits >= excess - *balance {
                fine_priority[j] = 1;
            } else {
                fine_priority[j] = 0;
            }
            excess -= extra_bits;
        }

        *balance = excess;

        debug_assert!(bits[j] >= 0);
        debug_assert!(ebits[j] >= 0);
    }

    // Handle uncoded bands
    for j in coded_bands..end {
        ebits[j] = bits[j] >> stereo >> BITRES;
        debug_assert!(c * ebits[j] << BITRES == bits[j]);
        bits[j] = 0;
        fine_priority[j] = if ebits[j] < 1 { 1 } else { 0 };
    }

    coded_bands
}

/// Initialize capacity for each band.
pub fn init_caps(m: &CeltMode, cap: &mut [i32], lm: usize, channels: usize) {
    let c = channels as i32;
    for i in 0..m.nb_ebands {
        let n = ((m.ebands[i + 1] - m.ebands[i]) as i32) << lm;
        let cache_idx = m.nb_ebands * (2 * lm + channels - 1) + i;
        cap[i] = ((m.cache.caps[cache_idx] as i32) + 64) * c * n >> 2;
    }
}

/// Compute bit allocation for all bands.
///
/// Returns (coded_bands, intensity, dual_stereo, balance).
#[allow(clippy::too_many_arguments)]
pub fn compute_allocation(
    m: &CeltMode,
    start: usize,
    end: usize,
    offsets: &[i32],
    cap: &[i32],
    alloc_trim: i32,
    total: i32,
    pulses: &mut [i32],
    ebits: &mut [i32],
    fine_priority: &mut [i32],
    channels: usize,
    lm: usize,
    dec: &mut RangeDecoder<'_>,
    prev: bool,
    prev_frames: i32,
    signal_bandwidth: i32,
) -> (usize, usize, bool, i32) {
    let total = total.max(0);
    let len = m.nb_ebands;
    let c = channels as i32;
    let stereo = if c > 1 { 1 } else { 0 };
    let log_m = lm << BITRES;

    // Reserve bits for skip and stereo flags
    let mut skip_start = start;
    let mut skip_rsv = if total >= 1 << BITRES { 1 << BITRES } else { 0 };
    let mut total = total - skip_rsv;

    let mut intensity_rsv = 0;
    let mut dual_stereo_rsv = 0;

    if c == 2 {
        intensity_rsv = LOG2_FRAC_TABLE[end - start] as i32;
        if intensity_rsv > total {
            intensity_rsv = 0;
        } else {
            total -= intensity_rsv;
            dual_stereo_rsv = if total >= 1 << BITRES { 1 << BITRES } else { 0 };
            total -= dual_stereo_rsv;
        }
    }

    // Compute allocation vectors
    let mut bits1 = [0i32; 22];
    let mut bits2 = [0i32; 22];
    let mut thresh = [0i32; 22];
    let mut trim_offset = [0i32; 22];

    for j in start..end {
        thresh[j] = (c << BITRES).max(
            (3 * (m.ebands[j + 1] - m.ebands[j]) as i32) << lm << BITRES >> 4,
        );
        trim_offset[j] = c
            * (m.ebands[j + 1] - m.ebands[j]) as i32
            * (alloc_trim - 5 - lm as i32)
            * (end - j - 1) as i32
            * (1 << (lm + BITRES as usize))
            >> 6;
        if (m.ebands[j + 1] - m.ebands[j]) << lm == 1 {
            trim_offset[j] -= c << BITRES;
        }
    }

    // Binary search for allocation vector
    let mut lo = 1i32;
    let mut hi = m.nb_alloc_vectors as i32 - 1;

    while lo <= hi {
        let mid = (lo + hi) >> 1;
        let mut done = false;
        let mut psum = 0;

        for j in (start..end).rev() {
            let n = (m.ebands[j + 1] - m.ebands[j]) as i32;
            let mut bitsj = ((c * n) * m.alloc_vectors[(mid as usize) * len + j] as i32) << lm >> 2;
            if bitsj > 0 {
                bitsj = (bitsj + trim_offset[j]).max(0);
            }
            bitsj += offsets[j];

            if bitsj >= thresh[j] || done {
                done = true;
                psum += bitsj.min(cap[j]);
            } else if bitsj >= c << BITRES {
                psum += c << BITRES;
            }
        }

        if psum > total {
            hi = mid - 1;
        } else {
            lo = mid + 1;
        }
    }

    let hi = lo;
    let lo = hi - 1;

    // Compute final bits1 and bits2
    for j in start..end {
        let n = (m.ebands[j + 1] - m.ebands[j]) as i32;
        let mut bits1j = (c * n * m.alloc_vectors[lo as usize * len + j] as i32) << lm >> 2;
        let mut bits2j = if hi < m.nb_alloc_vectors as i32 {
            (c * n * m.alloc_vectors[hi as usize * len + j] as i32) << lm >> 2
        } else {
            cap[j]
        };

        if bits1j > 0 {
            bits1j = (bits1j + trim_offset[j]).max(0);
        }
        if bits2j > 0 {
            bits2j = (bits2j + trim_offset[j]).max(0);
        }
        if lo > 0 {
            bits1j += offsets[j];
        }
        bits2j += offsets[j];

        if offsets[j] > 0 {
            skip_start = j;
        }

        bits2j = (bits2j - bits1j).max(0);
        bits1[j] = bits1j;
        bits2[j] = bits2j;
    }

    // Final interpolation
    let mut intensity = 0i32;
    let mut dual_stereo = 0i32;
    let mut balance = 0i32;
    let lm = lm as i32;
    let prev_i32 = if prev { prev_frames } else { 0 };

    let coded_bands = interp_bits2pulses(
        m,
        start,
        end,
        skip_start,
        &bits1,
        &bits2,
        &thresh,
        cap,
        total,
        &mut balance,
        skip_rsv,
        &mut intensity,
        intensity_rsv,
        &mut dual_stereo,
        dual_stereo_rsv,
        pulses,
        ebits,
        fine_priority,
        channels,
        lm,
        dec,
        prev_i32,
        signal_bandwidth,
    );

    (coded_bands, intensity as usize, dual_stereo != 0, balance)
}

/// Integer division that rounds towards zero.
#[inline]
fn celt_udiv(n: i32, d: i32) -> i32 {
    debug_assert!(d > 0);
    n / d
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_pulses() {
        assert_eq!(get_pulses(0), 0);
        assert_eq!(get_pulses(1), 1);
        assert_eq!(get_pulses(7), 7);
        assert_eq!(get_pulses(8), 8); // (8 + 0) << 0 = 8
        assert_eq!(get_pulses(9), 9); // (8 + 1) << 0 = 9
        assert_eq!(get_pulses(16), 16); // (8 + 0) << 1 = 16
    }
}
