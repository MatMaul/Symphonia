use super::intrin::ec_ilog;

#[inline]
fn frac_mul16(a: i32, b: i32) -> i32 {
    (16384 + a * b) >> 15
}

pub fn hysteresis_decision(val: f32, thresholds: &[f32], hysteresis: &[f32], prev: i32) -> i32 {
    let mut i = 0i32;
    while (i as usize) < thresholds.len() {
        if val < thresholds[i as usize] {
            break;
        }
        i += 1;
    }
    if i > prev && val < thresholds[prev as usize] + hysteresis[prev as usize] {
        i = prev;
    }
    if i < prev && prev > 0 && val > thresholds[(prev - 1) as usize] - hysteresis[(prev - 1) as usize] {
        i = prev;
    }
    i
}

pub fn celt_lcg_rand(seed: u32) -> u32 {
    seed.wrapping_mul(1664525).wrapping_add(1013904223)
}

pub fn bitexact_cos(x: i16) -> i16 {
    let tmp = (4096 + (x as i32 * x as i32)) >> 13;
    debug_assert!(tmp <= 32767);
    let mut x2 = tmp;
    x2 = (32767 - x2)
        + frac_mul16(x2, (-7651 + frac_mul16(x2, (8277 + frac_mul16(-626, x2)))));
    debug_assert!(x2 <= 32766);
    (1 + x2) as i16
}

pub fn bitexact_log2tan(isin: i32, icos: i32) -> i32 {
    let lc = ec_ilog(icos as u32);
    let ls = ec_ilog(isin as u32);
    let icos = icos << (15 - lc);
    let isin = isin << (15 - ls);
    (ls - lc) * (1 << 11)
        + frac_mul16(isin, frac_mul16(isin, -2597) + 7932)
        - frac_mul16(icos, frac_mul16(icos, -2597) + 7932)
}
