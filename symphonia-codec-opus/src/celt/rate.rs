use super::modes::CeltMode;

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
