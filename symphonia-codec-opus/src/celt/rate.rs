use super::entcode::BITRES;
use super::entdec::EcDec;
use super::intrin::celt_udiv;
use super::modes::CeltMode;
use alloc::vec::Vec;

pub const MAX_PSEUDO: i32 = 40;
pub const LOG_MAX_PSEUDO: i32 = 6;
pub const CELT_MAX_PULSES: i32 = 128;
pub const MAX_FINE_BITS: i32 = 8;
pub const FINE_OFFSET: i32 = 21;
pub const QTHETA_OFFSET: i32 = 4;
pub const QTHETA_OFFSET_TWOPHASE: i32 = 16;

pub const LOG2_FRAC_TABLE: [u8; 24] = [
    0, 8, 13, 16, 19, 21, 23, 24, 26, 27, 28, 29, 30, 31, 32, 32, 33, 34, 34, 35, 36,
    36, 37, 37,
];

const ALLOC_STEPS: i32 = 6;

#[inline]
pub fn get_pulses(i: i32) -> i32 {
    if i < 8 {
        i
    } else {
        (8 + (i & 7)) << ((i >> 3) - 1)
    }
}

pub fn bits2pulses(m: &CeltMode, band: i32, lm: i32, bits: i32) -> i32 {
    let lm = lm + 1;
    let cache_index = m.cache.index[(lm * m.nb_ebands + band) as usize] as usize;
    let cache = &m.cache.bits[cache_index..];
    let mut lo = 0;
    let mut hi = cache[0] as i32;
    let target = bits - 1;
    for _ in 0..LOG_MAX_PSEUDO {
        let mid = (lo + hi + 1) >> 1;
        if cache[mid as usize] as i32 >= target {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let lo_bits = if lo == 0 { -1 } else { cache[lo as usize] as i32 };
    if target - lo_bits <= cache[hi as usize] as i32 - target {
        lo
    } else {
        hi
    }
}

pub fn pulses2bits(m: &CeltMode, band: i32, lm: i32, pulses: i32) -> i32 {
    let lm = lm + 1;
    let cache_index = m.cache.index[(lm * m.nb_ebands + band) as usize] as usize;
    let cache = &m.cache.bits[cache_index..];
    if pulses == 0 {
        0
    } else {
        cache[pulses as usize] as i32 + 1
    }
}

pub fn init_caps(m: &CeltMode, cap: &mut [i32], lm: i32, channels: i32) {
    let nb_ebands = m.nb_ebands as usize;
    let lm = lm as i32;
    for i in 0..nb_ebands {
        let n = (m.ebands[i + 1] as i32 - m.ebands[i] as i32) << lm;
        let idx = (m.nb_ebands * (2 * lm + channels - 1)) as usize + i;
        let base = m.cache.caps[idx] as i32 + 64;
        cap[i] = (base * channels * n) >> 2;
    }
}

fn interp_bits2pulses_decode(
    m: &CeltMode,
    start: i32,
    end: i32,
    skip_start: i32,
    bits1: &[i32],
    bits2: &[i32],
    thresh: &[i32],
    cap: &[i32],
    mut total: i32,
    balance: &mut i32,
    skip_rsv: i32,
    intensity: &mut i32,
    mut intensity_rsv: i32,
    dual_stereo: &mut i32,
    mut dual_stereo_rsv: i32,
    pulses: &mut [i32],
    ebits: &mut [i32],
    fine_priority: &mut [i32],
    channels: i32,
    lm: i32,
    dec: &mut EcDec<'_>,
) -> i32 {
    let alloc_floor = channels << (BITRES as i32);
    let stereo = if channels > 1 { 1 } else { 0 };
    let log_m = lm << (BITRES as i32);

    let mut lo = 0i32;
    let mut hi = 1i32 << ALLOC_STEPS;
    for _ in 0..ALLOC_STEPS {
        let mid = (lo + hi) >> 1;
        let mut psum = 0i32;
        let mut done = false;
        for j in (start..end).rev() {
            let idx = j as usize;
            let tmp = bits1[idx] + ((mid * bits2[idx]) >> ALLOC_STEPS);
            if tmp >= thresh[idx] || done {
                done = true;
                psum += tmp.min(cap[idx]);
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

    let mut psum = 0i32;
    let mut done = false;
    for j in (start..end).rev() {
        let idx = j as usize;
        let mut tmp = bits1[idx] + ((lo * bits2[idx]) >> ALLOC_STEPS);
        if tmp < thresh[idx] && !done {
            tmp = if tmp >= alloc_floor { alloc_floor } else { 0 };
        } else {
            done = true;
        }
        tmp = tmp.min(cap[idx]);
        pulses[idx] = tmp;
        psum += tmp;
    }

    let mut coded_bands = end;
    loop {
        let j = coded_bands - 1;
        if j <= skip_start {
            total += skip_rsv;
            break;
        }
        let left = total - psum;
        let denom = (m.ebands[coded_bands as usize] - m.ebands[start as usize]) as i32;
        let percoeff = celt_udiv(left as u32, denom as u32) as i32;
        let mut left = left - denom * percoeff;
        let rem = (left - (m.ebands[j as usize] - m.ebands[start as usize]) as i32).max(0);
        let band_width = (m.ebands[coded_bands as usize] - m.ebands[j as usize]) as i32;
        let mut band_bits = pulses[j as usize] + percoeff * band_width + rem;
        if band_bits >= thresh[j as usize].max(alloc_floor + (1 << (BITRES as i32))) {
            if dec.dec_bit_logp(1) != 0 {
                break;
            }
            psum += 1 << (BITRES as i32);
            band_bits -= 1 << (BITRES as i32);
        }
        psum -= pulses[j as usize] + intensity_rsv;
        if intensity_rsv > 0 {
            intensity_rsv = LOG2_FRAC_TABLE[(j - start) as usize] as i32;
        }
        psum += intensity_rsv;
        if band_bits >= alloc_floor {
            psum += alloc_floor;
            pulses[j as usize] = alloc_floor;
        } else {
            pulses[j as usize] = 0;
        }
        coded_bands -= 1;
    }

    if intensity_rsv > 0 {
        *intensity = start + dec.dec_uint((coded_bands + 1 - start) as u32) as i32;
    } else {
        *intensity = 0;
    }
    if *intensity <= start {
        total += dual_stereo_rsv;
        dual_stereo_rsv = 0;
    }
    if dual_stereo_rsv > 0 {
        *dual_stereo = dec.dec_bit_logp(1);
    } else {
        *dual_stereo = 0;
    }

    let left = total - psum;
    let denom = (m.ebands[coded_bands as usize] - m.ebands[start as usize]) as i32;
    let percoeff = celt_udiv(left as u32, denom as u32) as i32;
    let mut left = left - denom * percoeff;
    for j in start..coded_bands {
        let idx = j as usize;
        let width = (m.ebands[idx + 1] - m.ebands[idx]) as i32;
        pulses[idx] += percoeff * width;
    }
    for j in start..coded_bands {
        let idx = j as usize;
        let width = (m.ebands[idx + 1] - m.ebands[idx]) as i32;
        let tmp = left.min(width);
        pulses[idx] += tmp;
        left -= tmp;
    }

    let mut running_balance = 0i32;
    for j in start..coded_bands {
        let idx = j as usize;
        let n0 = (m.ebands[idx + 1] - m.ebands[idx]) as i32;
        let n = n0 << lm;
        let bit = pulses[idx] + running_balance;
        let mut excess;
        if n > 1 {
            excess = (bit - cap[idx]).max(0);
            pulses[idx] = bit - excess;

            let den = channels * n
                + if channels == 2 && n > 2 && *dual_stereo == 0 && j < *intensity {
                    1
                } else {
                    0
                };
            let nclogn = den * (m.log_n[idx] as i32 + log_m);
            let mut offset = (nclogn >> 1) - den * FINE_OFFSET;
            if n == 2 {
                offset += (den << (BITRES as i32)) >> 2;
            }
            if pulses[idx] + offset < den * 2 << (BITRES as i32) {
                offset += nclogn >> 2;
            } else if pulses[idx] + offset < den * 3 << (BITRES as i32) {
                offset += nclogn >> 3;
            }

            let mut eb = pulses[idx] + offset + (den << ((BITRES as i32) - 1));
            if eb < 0 {
                eb = 0;
            }
            let mut ebits_j = (celt_udiv(eb as u32, den as u32) as i32) >> (BITRES as i32);
            if channels * ebits_j > (pulses[idx] >> stereo >> (BITRES as i32)) {
                ebits_j = pulses[idx] >> stereo >> (BITRES as i32);
            }
            ebits_j = ebits_j.min(MAX_FINE_BITS);
            fine_priority[idx] = if ebits_j * (den << (BITRES as i32)) >= pulses[idx] + offset {
                1
            } else {
                0
            };
            pulses[idx] -= channels * ebits_j << (BITRES as i32);
            ebits[idx] = ebits_j;
        } else {
            excess = (bit - (channels << (BITRES as i32))).max(0);
            pulses[idx] = bit - excess;
            ebits[idx] = 0;
            fine_priority[idx] = 1;
        }

        if excess > 0 {
            let mut extra_fine = excess >> (stereo + (BITRES as i32));
            extra_fine = extra_fine.min(MAX_FINE_BITS - ebits[idx]);
            ebits[idx] += extra_fine;
            let extra_bits = extra_fine * channels << (BITRES as i32);
            fine_priority[idx] = if extra_bits >= excess - running_balance { 1 } else { 0 };
            excess -= extra_bits;
        }
        running_balance = excess;
    }
    *balance = running_balance;

    for j in coded_bands..end {
        let idx = j as usize;
        ebits[idx] = pulses[idx] >> stereo >> (BITRES as i32);
        pulses[idx] = 0;
        fine_priority[idx] = if ebits[idx] < 1 { 1 } else { 0 };
    }

    coded_bands
}

pub fn clt_compute_allocation(
    m: &CeltMode,
    start: i32,
    end: i32,
    offsets: &[i32],
    cap: &[i32],
    alloc_trim: i32,
    intensity: &mut i32,
    dual_stereo: &mut i32,
    mut total: i32,
    balance: &mut i32,
    pulses: &mut [i32],
    ebits: &mut [i32],
    fine_priority: &mut [i32],
    channels: i32,
    lm: i32,
    dec: &mut EcDec<'_>,
) -> i32 {
    total = total.max(0);
    let len = m.nb_ebands as usize;
    let mut skip_start = start;
    let skip_rsv = if total >= 1 << (BITRES as i32) {
        1 << (BITRES as i32)
    } else {
        0
    };
    total -= skip_rsv;

    let mut intensity_rsv = 0;
    let mut dual_stereo_rsv = 0;
    if channels == 2 {
        intensity_rsv = LOG2_FRAC_TABLE[(end - start) as usize] as i32;
        if intensity_rsv > total {
            intensity_rsv = 0;
        } else {
            total -= intensity_rsv;
            dual_stereo_rsv = if total >= 1 << (BITRES as i32) {
                1 << (BITRES as i32)
            } else {
                0
            };
            total -= dual_stereo_rsv;
        }
    }

    let mut bits1 = vec![0i32; len];
    let mut bits2 = vec![0i32; len];
    let mut thresh = vec![0i32; len];
    let mut trim_offset = vec![0i32; len];

    for j in start..end {
        let idx = j as usize;
        let n = (m.ebands[idx + 1] - m.ebands[idx]) as i32;
        thresh[idx] = (channels << (BITRES as i32))
            .max((3 * n << lm << (BITRES as i32)) >> 4);
        trim_offset[idx] = channels
            * n
            * (alloc_trim - 5 - lm)
            * (end - j - 1)
            * (1 << (lm + (BITRES as i32)))
            >> 6;
        if (n << lm) == 1 {
            trim_offset[idx] -= channels << (BITRES as i32);
        }
    }

    let mut lo = 1i32;
    let mut hi = m.nb_alloc_vectors - 1;
    while lo <= hi {
        let mid = (lo + hi) >> 1;
        let mut done = false;
        let mut psum = 0i32;
        for j in (start..end).rev() {
            let idx = j as usize;
            let n = (m.ebands[idx + 1] - m.ebands[idx]) as i32;
            let mut bitsj = channels
                * n
                * m.alloc_vectors[(mid as usize) * len + idx] as i32
                << lm
                >> 2;
            if bitsj > 0 {
                bitsj = (bitsj + trim_offset[idx]).max(0);
            }
            bitsj += offsets[idx];
            if bitsj >= thresh[idx] || done {
                done = true;
                psum += bitsj.min(cap[idx]);
            } else if bitsj >= channels << (BITRES as i32) {
                psum += channels << (BITRES as i32);
            }
        }
        if psum > total {
            hi = mid - 1;
        } else {
            lo = mid + 1;
        }
    }
    hi = lo;
    lo -= 1;

    for j in start..end {
        let idx = j as usize;
        let n = (m.ebands[idx + 1] - m.ebands[idx]) as i32;
        let mut bits1j = channels
            * n
            * m.alloc_vectors[(lo as usize) * len + idx] as i32
            << lm
            >> 2;
        let mut bits2j = if hi >= m.nb_alloc_vectors {
            cap[idx]
        } else {
            channels
                * n
                * m.alloc_vectors[(hi as usize) * len + idx] as i32
                << lm
                >> 2
        };
        if bits1j > 0 {
            bits1j = (bits1j + trim_offset[idx]).max(0);
        }
        if bits2j > 0 {
            bits2j = (bits2j + trim_offset[idx]).max(0);
        }
        if lo > 0 {
            bits1j += offsets[idx];
        }
        bits2j += offsets[idx];
        if offsets[idx] > 0 {
            skip_start = j;
        }
        bits2j = (bits2j - bits1j).max(0);
        bits1[idx] = bits1j;
        bits2[idx] = bits2j;
    }

    let coded_bands = interp_bits2pulses_decode(
        m,
        start,
        end,
        skip_start,
        &bits1,
        &bits2,
        &thresh,
        cap,
        total,
        balance,
        skip_rsv,
        intensity,
        intensity_rsv,
        dual_stereo,
        dual_stereo_rsv,
        pulses,
        ebits,
        fine_priority,
        channels,
        lm,
        dec,
    );

    coded_bands
}
