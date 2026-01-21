use super::intrin::ec_ilog;

#[cfg(feature = "std")]
fn cosf(x: f32) -> f32 {
    x.cos()
}

#[cfg(not(feature = "std"))]
fn cosf(x: f32) -> f32 {
    libm::cosf(x)
}

#[cfg(feature = "std")]
fn sqrtf(x: f32) -> f32 {
    x.sqrt()
}

#[cfg(not(feature = "std"))]
fn sqrtf(x: f32) -> f32 {
    libm::sqrtf(x)
}

pub const PI: f32 = 3.1415926535897931;

pub fn isqrt32(mut val: u32) -> u32 {
    let mut g = 0u32;
    let mut bshift = ((ec_ilog(val) - 1) >> 1) as i32;
    let mut b = 1u32 << bshift;
    while bshift >= 0 {
        let t = (((g << 1) + b) as u64) << bshift;
        if t <= val as u64 {
            g += b;
            val -= t as u32;
        }
        b >>= 1;
        bshift -= 1;
    }
    g
}

pub fn celt_div(a: f32, b: f32) -> f32 {
    a / b
}

pub fn celt_rcp(x: f32) -> f32 {
    1.0 / x
}

pub fn celt_sqrt(x: f32) -> f32 {
    sqrtf(x)
}

pub fn celt_sqrt32(x: f32) -> f32 {
    sqrtf(x)
}

pub fn celt_rsqrt(x: f32) -> f32 {
    1.0 / sqrtf(x)
}

pub fn celt_rsqrt_norm(x: f32) -> f32 {
    celt_rsqrt(x)
}

pub fn celt_rsqrt_norm32(x: f32) -> f32 {
    celt_rsqrt(x)
}

pub fn celt_cos_norm(x: f32) -> f32 {
    cosf(0.5 * PI * x)
}

pub fn celt_atan_norm(x: f32) -> f32 {
    const ATAN2_2_OVER_PI: f32 = 0.636619772367581;
    const A03: f32 = -3.3331659436225891113281250000e-01;
    const A05: f32 = 1.99627041816711425781250000000e-01;
    const A07: f32 = -1.3976582884788513183593750000e-01;
    const A09: f32 = 9.79423448443412780761718750000e-02;
    const A11: f32 = -5.7773590087890625000000000000e-02;
    const A13: f32 = 2.30401363223791122436523437500e-02;
    const A15: f32 = -4.3554059229791164398193359375e-03;

    let x_sq = x * x;
    ATAN2_2_OVER_PI
        * (x
            + x * x_sq
                * (A03
                    + x_sq * (A05 + x_sq * (A07 + x_sq * (A09 + x_sq * (A11 + x_sq * (A13 + x_sq * A15)))))))
}

pub fn celt_atan2p_norm(y: f32, x: f32) -> f32 {
    debug_assert!(x >= 0.0 && y >= 0.0);
    if (x * x + y * y) < 1e-18 {
        return 0.0;
    }
    if y < x {
        celt_atan_norm(y / x)
    } else {
        1.0 - celt_atan_norm(x / y)
    }
}
