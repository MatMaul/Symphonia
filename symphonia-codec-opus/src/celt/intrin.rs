#[inline]
pub fn ec_ilog(v: u32) -> i32 {
    if v == 0 {
        0
    } else {
        32 - v.leading_zeros() as i32
    }
}

#[inline]
pub fn ec_mini(a: u32, b: u32) -> u32 {
    if a < b {
        a
    } else {
        b
    }
}

#[inline]
pub fn imul32(a: u32, b: u32) -> u32 {
    a.wrapping_mul(b)
}

#[inline]
pub fn celt_udiv(n: u32, d: u32) -> u32 {
    debug_assert!(d > 0);
    n / d
}

#[inline]
pub fn celt_sudiv(n: i32, d: i32) -> i32 {
    debug_assert!(d > 0);
    if n < 0 {
        -(celt_udiv(n.wrapping_neg() as u32, d as u32) as i32)
    } else {
        celt_udiv(n as u32, d as u32) as i32
    }
}
