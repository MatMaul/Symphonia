use super::entdec::EcDec;

const LAPLACE_LOG_MINP: u32 = 0;
const LAPLACE_MINP: u32 = 1 << LAPLACE_LOG_MINP;
const LAPLACE_NMIN: u32 = 16;

#[inline]
fn imin_u32(a: u32, b: u32) -> u32 {
    if a < b { a } else { b }
}

#[inline]
fn imax_u16(a: u16, b: u16) -> u16 {
    if a > b { a } else { b }
}

#[inline]
fn ec_laplace_get_freq1(fs0: u32, decay: i32) -> u32 {
    let ft = 32768u32 - LAPLACE_MINP * (2 * LAPLACE_NMIN) - fs0;
    ((ft as i64 * (16384 - decay) as i64) >> 15) as u32
}

pub fn ec_laplace_decode(dec: &mut EcDec<'_>, mut fs: u32, decay: i32) -> i32 {
    let mut val = 0i32;
    let fm = dec.decode_bin(15);
    let mut fl = 0u32;
    if fm >= fs {
        val += 1;
        fl = fs;
        fs = ec_laplace_get_freq1(fs, decay) + LAPLACE_MINP;
        while fs > LAPLACE_MINP && fm >= fl + 2 * fs {
            fs *= 2;
            fl += fs;
            fs = (((fs - 2 * LAPLACE_MINP) as i64 * decay as i64) >> 15) as u32;
            fs += LAPLACE_MINP;
            val += 1;
        }
        if fs <= LAPLACE_MINP {
            let di = (fm - fl) >> (LAPLACE_LOG_MINP + 1);
            val += di as i32;
            fl += 2 * di * LAPLACE_MINP;
        }
        if fm < fl + fs {
            val = -val;
        } else {
            fl += fs;
        }
    }
    debug_assert!(fl < 32768);
    debug_assert!(fs > 0);
    debug_assert!(fl <= fm);
    let upper = imin_u32(fl + fs, 32768);
    debug_assert!(fm < upper);
    dec.update(fl, upper, 32768);
    val
}

pub fn ec_laplace_decode_p0(dec: &mut EcDec<'_>, p0: u16, decay: u16) -> i32 {
    let mut sign_icdf = [0u16; 3];
    sign_icdf[0] = 32768 - p0;
    sign_icdf[1] = sign_icdf[0] / 2;
    sign_icdf[2] = 0;
    let mut s = dec.dec_icdf16(&sign_icdf, 15);
    if s == 2 {
        s = -1;
    }
    if s != 0 {
        let mut icdf = [0u16; 8];
        icdf[0] = imax_u16(7, decay);
        for i in 1..7 {
            let decayed = ((icdf[i - 1] as u32 * decay as u32) >> 15) as u16;
            icdf[i] = imax_u16(7 - i as u16, decayed);
        }
        icdf[7] = 0;
        let mut value = 1i32;
        loop {
            let v = dec.dec_icdf16(&icdf, 15);
            value += v;
            if v != 7 {
                break;
            }
        }
        s * value
    } else {
        0
    }
}
