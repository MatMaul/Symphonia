use super::intrin::ec_ilog;

pub const EC_UINT_BITS: u32 = 8;
pub const BITRES: u32 = 3;

#[inline]
pub fn tell(nbits_total: i32, rng: u32) -> i32 {
    nbits_total - ec_ilog(rng)
}

pub fn tell_frac(nbits_total: i32, rng: u32) -> u32 {
    const CORRECTION: [u32; 8] = [
        35733, 38967, 42495, 46340, 50535, 55109, 60097, 65535,
    ];
    let nbits = (nbits_total as u32) << BITRES;
    let mut l = ec_ilog(rng) as u32;
    debug_assert!(l >= 16);
    let r = rng >> (l - 16);
    let mut b = (r >> 12) - 8;
    if r > CORRECTION[b as usize] {
        b += 1;
    }
    l = (l << 3) + b;
    nbits - l
}
