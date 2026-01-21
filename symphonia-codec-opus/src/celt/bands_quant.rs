use super::entcode::BITRES;
use super::entdec::EcDec;
use super::types::{CeltEner, CeltNorm};
use super::modes::CeltMode;
use super::rate::{bits2pulses, get_pulses, pulses2bits, QTHETA_OFFSET, QTHETA_OFFSET_TWOPHASE};
use super::{
    alg_unquant, bitexact_cos, bitexact_log2tan, celt_lcg_rand, compute_qn, isqrt32,
    renormalise_vector, deinterleave_hadamard, interleave_hadamard, haar1, stereo_merge,
};
use super::math::celt_sqrt;
use super::intrin::celt_udiv;

#[derive(Clone, Copy, Default)]
pub struct SplitCtx {
    pub inv: i32,
    pub imid: i32,
    pub iside: i32,
    pub delta: i32,
    pub itheta: i32,
    pub qalloc: i32,
}

pub struct BandCtx<'a> {
    pub m: &'a CeltMode,
    pub i: i32,
    pub intensity: i32,
    pub spread: i32,
    pub tf_change: i32,
    pub remaining_bits: i32,
    pub band_e: Option<&'a [CeltEner]>,
    pub seed: u32,
    pub resynth: bool,
    pub disable_inv: bool,
    pub avoid_split_noise: bool,
}

pub fn quant_band_n1(
    ctx: &mut BandCtx<'_>,
    dec: &mut EcDec<'_>,
    x: &mut [CeltNorm],
    y: Option<&mut [CeltNorm]>,
    lowband_out: Option<&mut [CeltNorm]>,
) -> u32 {
    let bitres = BITRES as u32;
    let bit = 1i32 << bitres;
    let mut sign = 0i32;
    if ctx.remaining_bits >= bit {
        sign = dec.dec_bits(1) as i32;
        ctx.remaining_bits -= bit;
    }
    if ctx.resynth {
        x[0] = if sign != 0 { -1.0 } else { 1.0 };
    }
    if let Some(y) = y {
        let mut sign = 0i32;
        if ctx.remaining_bits >= bit {
            sign = dec.dec_bits(1) as i32;
            ctx.remaining_bits -= bit;
        }
        if ctx.resynth {
            y[0] = if sign != 0 { -1.0 } else { 1.0 };
        }
    }
    if let Some(lowband_out) = lowband_out {
        lowband_out[0] = x[0];
    }
    1
}

#[inline]
fn frac_mul16(a: i32, b: i32) -> i32 {
    (16384 + a * b) >> 15
}

pub fn compute_theta_decode(
    ctx: &mut BandCtx<'_>,
    dec: &mut EcDec<'_>,
    n: i32,
    b: &mut i32,
    b0: i32,
    lm: i32,
    stereo: bool,
    fill: &mut i32,
) -> SplitCtx {
    let mut sctx = SplitCtx::default();
    let pulse_cap = ctx.m.log_n[ctx.i as usize] as i32 + lm * (1 << BITRES);
    let offset = (pulse_cap >> 1)
        - if stereo && n == 2 {
            QTHETA_OFFSET_TWOPHASE
        } else {
            QTHETA_OFFSET
        };
    let mut qn = compute_qn(n, *b, offset, pulse_cap, stereo);
    if stereo && ctx.i >= ctx.intensity {
        qn = 1;
    }
    let tell = dec.tell_frac() as i32;
    let mut itheta = 0;
    let mut inv = 0;

    if qn != 1 {
        if stereo && n > 2 {
            let p0 = 3;
            let x0 = qn / 2;
            let ft = p0 * (x0 + 1) + x0;
            let fs = dec.decode(ft as u32) as i32;
            let x = if fs < (x0 + 1) * p0 {
                fs / p0
            } else {
                x0 + 1 + (fs - (x0 + 1) * p0)
            };
            let fl = if x <= x0 {
                p0 * x
            } else {
                (x - 1 - x0) + (x0 + 1) * p0
            };
            let fh = if x <= x0 {
                p0 * (x + 1)
            } else {
                (x - x0) + (x0 + 1) * p0
            };
            dec.update(fl as u32, fh as u32, ft as u32);
            itheta = x;
        } else if b0 > 1 || stereo {
            itheta = dec.dec_uint((qn + 1) as u32) as i32;
        } else {
            let ft = ((qn >> 1) + 1) * ((qn >> 1) + 1);
            let fm = dec.decode(ft as u32) as i32;
            let (fs, fl);
            if fm < ((qn >> 1) * ((qn >> 1) + 1) >> 1) {
                itheta = (isqrt32((8 * fm + 1) as u32) as i32 - 1) >> 1;
                fs = itheta + 1;
                fl = itheta * (itheta + 1) >> 1;
            } else {
                itheta = (2 * (qn + 1) - isqrt32((8 * (ft - fm - 1) + 1) as u32) as i32)
                    >> 1;
                fs = qn + 1 - itheta;
                fl = ft - ((qn + 1 - itheta) * (qn + 2 - itheta) >> 1);
            }
            dec.update(fl as u32, (fl + fs) as u32, ft as u32);
        }
        itheta = celt_udiv((itheta * 16384) as u32, qn as u32) as i32;
    } else if stereo {
        let bitres = BITRES as u32;
        if *b > (2i32 << bitres) && ctx.remaining_bits > (2i32 << bitres) {
            inv = dec.dec_bit_logp(2);
        }
        if ctx.disable_inv {
            inv = 0;
        }
        itheta = 0;
    }

    let qalloc = dec.tell_frac() as i32 - tell;
    *b -= qalloc;

    let (imid, iside, delta);
    if itheta == 0 {
        imid = 32767;
        iside = 0;
        *fill &= (1i32 << (b0 as u32)) - 1;
        delta = -16384;
    } else if itheta == 16384 {
        imid = 0;
        iside = 32767;
        *fill &= ((1i32 << (b0 as u32)) - 1) << (b0 as u32);
        delta = 16384;
    } else {
        imid = bitexact_cos(itheta as i16) as i32;
        iside = bitexact_cos((16384 - itheta) as i16) as i32;
        delta = frac_mul16((n - 1) << 7, bitexact_log2tan(iside, imid));
    }

    sctx.inv = inv;
    sctx.imid = imid;
    sctx.iside = iside;
    sctx.delta = delta;
    sctx.itheta = itheta;
    sctx.qalloc = qalloc;
    sctx
}

pub fn quant_partition_decode(
    ctx: &mut BandCtx<'_>,
    dec: &mut EcDec<'_>,
    x: &mut [CeltNorm],
    n: i32,
    mut b: i32,
    b_count: i32,
    lowband: Option<&[CeltNorm]>,
    lm: i32,
    gain: f32,
    mut fill: i32,
) -> u32 {
    let b0 = b_count;
    let cache_index = ctx.m.cache.index[(lm + 1) as usize * ctx.m.nb_ebands as usize + ctx.i as usize] as usize;
    let cache = &ctx.m.cache.bits[cache_index..];

    if lm != -1 && b > cache[cache[0] as usize] as i32 + 12 && n > 2 {
        let mut n = n >> 1;
        let (x0, y) = x.split_at_mut(n as usize);
        let mut lm = lm - 1;
        if b_count == 1 {
            fill = (fill & 1) | (fill << 1);
        }
        let b_split = (b_count + 1) >> 1;
        let mut sctx = compute_theta_decode(ctx, dec, n, &mut b, b0, lm, false, &mut fill);
        let mid = sctx.imid as f32 * (1.0 / 32768.0);
        let side = sctx.iside as f32 * (1.0 / 32768.0);
        let mut delta = sctx.delta;

        if b0 > 1 && (sctx.itheta & 0x3fff) != 0 {
            if sctx.itheta > 8192 {
                delta -= delta >> ((4 - lm) as u32);
            } else {
                let adjust = (n << BITRES) >> ((5 - lm) as u32);
                delta = (delta + adjust).min(0);
            }
        }
        let mut mbits = (b - delta) / 2;
        if mbits < 0 {
            mbits = 0;
        }
        if mbits > b {
            mbits = b;
        }
        let mut sbits = b - mbits;
        ctx.remaining_bits -= sctx.qalloc;

        let next_lowband = lowband.map(|lb| &lb[n as usize..]);
        let rebalance = ctx.remaining_bits;
        let mut cm;
        if mbits >= sbits {
            cm = quant_partition_decode(
                ctx,
                dec,
                x0,
                n,
                mbits,
                b_split,
                lowband,
                lm,
                gain * mid,
                fill,
            );
            let rebalance = mbits - (rebalance - ctx.remaining_bits);
            if rebalance > 3 << BITRES && sctx.itheta != 0 {
                sbits += rebalance - (3 << BITRES);
            }
            let cm2 = quant_partition_decode(
                ctx,
                dec,
                y,
                n,
                sbits,
                b_split,
                next_lowband,
                lm,
                gain * side,
                fill >> (b_split as u32),
            );
            cm |= cm2 << ((b0 >> 1) as u32);
        } else {
            cm = quant_partition_decode(
                ctx,
                dec,
                y,
                n,
                sbits,
                b_split,
                next_lowband,
                lm,
                gain * side,
                fill >> (b_split as u32),
            ) << ((b0 >> 1) as u32);
            let rebalance = sbits - (rebalance - ctx.remaining_bits);
            if rebalance > 3 << BITRES && sctx.itheta != 16384 {
                mbits += rebalance - (3 << BITRES);
            }
            cm |= quant_partition_decode(
                ctx,
                dec,
                x0,
                n,
                mbits,
                b_split,
                lowband,
                lm,
                gain * mid,
                fill,
            );
        }
        cm
    } else {
        let mut q = bits2pulses(ctx.m, ctx.i, lm, b);
        let mut curr_bits = pulses2bits(ctx.m, ctx.i, lm, q);
        ctx.remaining_bits -= curr_bits;
        while ctx.remaining_bits < 0 && q > 0 {
            ctx.remaining_bits += curr_bits;
            q -= 1;
            curr_bits = pulses2bits(ctx.m, ctx.i, lm, q);
            ctx.remaining_bits -= curr_bits;
        }
        if q != 0 {
            let k = get_pulses(q);
            alg_unquant(x, n, k, ctx.spread, b_count, dec, gain)
        } else {
            if ctx.resynth {
                let cm_mask = (1u32 << (b_count as u32)) - 1;
                let fill_mask = (fill as u32) & cm_mask;
                if fill_mask == 0 {
                    for v in x.iter_mut().take(n as usize) {
                        *v = 0.0;
                    }
                    return 0;
                }
                if let Some(lowband) = lowband {
                    for j in 0..(n as usize) {
                        ctx.seed = celt_lcg_rand(ctx.seed);
                        let mut tmp = 1.0 / 256.0;
                        if (ctx.seed & 0x8000) == 0 {
                            tmp = -tmp;
                        }
                        x[j] = lowband[j] + tmp;
                    }
                } else {
                    for j in 0..(n as usize) {
                        ctx.seed = celt_lcg_rand(ctx.seed);
                        x[j] = (ctx.seed >> 20) as i32 as f32;
                    }
                }
                renormalise_vector(x, n, gain);
                return fill_mask;
            }
            0
        }
    }
}

const BIT_INTERLEAVE_TABLE: [u8; 16] = [
    0, 1, 1, 1, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3,
];

const BIT_DEINTERLEAVE_TABLE: [u8; 16] = [
    0x00, 0x03, 0x0C, 0x0F, 0x30, 0x33, 0x3C, 0x3F, 0xC0, 0xC3, 0xCC, 0xCF, 0xF0, 0xF3,
    0xFC, 0xFF,
];

pub fn quant_band_decode(
    ctx: &mut BandCtx<'_>,
    dec: &mut EcDec<'_>,
    x: &mut [CeltNorm],
    n: i32,
    b: i32,
    b_count: i32,
    lowband: Option<&[CeltNorm]>,
    lm: i32,
    lowband_out: Option<&mut [CeltNorm]>,
    gain: f32,
    lowband_scratch: Option<&mut [CeltNorm]>,
    mut fill: i32,
) -> u32 {
    let n0 = n;
    let mut n_b = n;
    let b0 = b_count;
    let mut time_divide = 0;
    let mut recombine = 0;
    let long_blocks = b0 == 1;
    let mut tf_change = ctx.tf_change;

    n_b = celt_udiv(n_b as u32, b_count as u32) as i32;

    if n == 1 {
        return quant_band_n1(ctx, dec, x, None, lowband_out);
    }

    if tf_change > 0 {
        recombine = tf_change;
    }

    let mut lowband_buf = lowband_scratch;
    let mut use_scratch = false;
    if let (Some(lb), Some(scratch)) = (lowband, lowband_buf.as_deref_mut()) {
        if recombine != 0 || ((n_b & 1) == 0 && tf_change < 0) || b0 > 1 {
            scratch[..n as usize].copy_from_slice(&lb[..n as usize]);
            use_scratch = true;
        }
    }

    for k in 0..recombine {
        if use_scratch {
            if let Some(scratch) = lowband_buf.as_deref_mut() {
                haar1(scratch, n >> (k as u32), 1 << (k as u32));
            }
        }
        fill = BIT_INTERLEAVE_TABLE[(fill & 0xF) as usize] as i32
            | ((BIT_INTERLEAVE_TABLE[((fill >> 4) & 0xF) as usize] as i32) << 2);
    }
    let b_count = b_count >> (recombine as u32);
    n_b <<= recombine as u32;

    while (n_b & 1) == 0 && tf_change < 0 {
        if use_scratch {
            if let Some(scratch) = lowband_buf.as_deref_mut() {
                haar1(scratch, n_b, b_count);
            }
        }
        fill |= fill << (b_count as u32);
        n_b >>= 1;
        time_divide += 1;
        tf_change += 1;
    }

    let b0 = b_count;
    let n_b0 = n_b;

    if b0 > 1 {
        if use_scratch {
            if let Some(scratch) = lowband_buf.as_deref_mut() {
                deinterleave_hadamard(
                    scratch,
                    n_b >> (recombine as u32),
                    b0 << (recombine as u32),
                    long_blocks,
                );
            }
        }
    }

    let lowband_read = if use_scratch {
        lowband_buf.as_deref().map(|scratch| &scratch[..n as usize])
    } else {
        lowband
    };
    let mut cm = quant_partition_decode(ctx, dec, x, n, b, b_count, lowband_read, lm, gain, fill);

    if ctx.resynth {
        if b0 > 1 {
            interleave_hadamard(
                x,
                n_b >> (recombine as u32),
                b0 << (recombine as u32),
                long_blocks,
            );
        }

        n_b = n_b0;
        let mut b_count = b0;
        for _ in 0..time_divide {
            b_count >>= 1;
            n_b <<= 1;
            cm |= cm >> (b_count as u32);
            haar1(x, n_b, b_count);
        }

        for k in 0..recombine {
            cm = BIT_DEINTERLEAVE_TABLE[(cm & 0xF) as usize] as u32;
            haar1(x, n0 >> (k as u32), 1 << (k as u32));
        }

        let b_count = b_count << (recombine as u32);
        if let Some(lowband_out) = lowband_out {
            let scale = celt_sqrt(n0 as f32);
            for j in 0..(n0 as usize) {
                lowband_out[j] = scale * x[j];
            }
        }
        cm &= (1u32 << (b_count as u32)) - 1;
    }

    cm
}

pub fn quant_band_stereo_decode(
    ctx: &mut BandCtx<'_>,
    dec: &mut EcDec<'_>,
    x: &mut [CeltNorm],
    y: &mut [CeltNorm],
    n: i32,
    mut b: i32,
    b_count: i32,
    lowband: Option<&[CeltNorm]>,
    lm: i32,
    lowband_out: Option<&mut [CeltNorm]>,
    lowband_scratch: Option<&mut [CeltNorm]>,
    mut fill: i32,
) -> u32 {
    if n == 1 {
        return quant_band_n1(ctx, dec, x, Some(y), lowband_out);
    }

    let orig_fill = fill;
    let sctx = compute_theta_decode(ctx, dec, n, &mut b, b_count, lm, true, &mut fill);
    let inv = sctx.inv;
    let itheta = sctx.itheta;
    let delta = sctx.delta;
    let qalloc = sctx.qalloc;

    let mid = sctx.imid as f32 * (1.0 / 32768.0);
    let side = sctx.iside as f32 * (1.0 / 32768.0);

    let mut cm = 0u32;
    if n == 2 {
        let mut sbits = 0;
        if itheta != 0 && itheta != 16384 {
            sbits = 1 << BITRES;
        }
        let mbits = b - sbits;
        let use_y_as_mid = itheta > 8192;
        ctx.remaining_bits -= qalloc + sbits;
        let mut sign = 0i32;
        if sbits != 0 {
            sign = dec.dec_bits(1) as i32;
        }
        let sign = 1 - 2 * sign;
        if use_y_as_mid {
            cm = quant_band_decode(
                ctx,
                dec,
                y,
                n,
                mbits,
                b_count,
                lowband,
                lm,
                lowband_out,
                1.0,
                lowband_scratch,
                orig_fill,
            );
            x[0] = -(sign as f32) * y[1];
            x[1] = (sign as f32) * y[0];
        } else {
            cm = quant_band_decode(
                ctx,
                dec,
                x,
                n,
                mbits,
                b_count,
                lowband,
                lm,
                lowband_out,
                1.0,
                lowband_scratch,
                orig_fill,
            );
            y[0] = -(sign as f32) * x[1];
            y[1] = (sign as f32) * x[0];
        }
        if ctx.resynth {
            x[0] = mid * x[0];
            x[1] = mid * x[1];
            y[0] = side * y[0];
            y[1] = side * y[1];
            let tmp = x[0];
            x[0] = tmp - y[0];
            y[0] = tmp + y[0];
            let tmp = x[1];
            x[1] = tmp - y[1];
            y[1] = tmp + y[1];
        }
    } else {
        let mut mbits = (b - delta) / 2;
        if mbits < 0 {
            mbits = 0;
        }
        if mbits > b {
            mbits = b;
        }
        let mut sbits = b - mbits;
        ctx.remaining_bits -= qalloc;
        let rebalance = ctx.remaining_bits;
        if mbits >= sbits {
            cm = quant_band_decode(
                ctx,
                dec,
                x,
                n,
                mbits,
                b_count,
                lowband,
                lm,
                lowband_out,
                1.0,
                lowband_scratch,
                fill,
            );
            let rebalance = mbits - (rebalance - ctx.remaining_bits);
            if rebalance > 3 << BITRES && itheta != 0 {
                sbits += rebalance - (3 << BITRES);
            }
            let cm2 = quant_band_decode(
                ctx,
                dec,
                y,
                n,
                sbits,
                b_count,
                None,
                lm,
                None,
                side,
                None,
                fill >> (b_count as u32),
            );
            cm |= cm2;
        } else {
            cm = quant_band_decode(
                ctx,
                dec,
                y,
                n,
                sbits,
                b_count,
                None,
                lm,
                None,
                side,
                None,
                fill >> (b_count as u32),
            );
            let rebalance = sbits - (rebalance - ctx.remaining_bits);
            if rebalance > 3 << BITRES && itheta != 16384 {
                mbits += rebalance - (3 << BITRES);
            }
            let cm2 = quant_band_decode(
                ctx,
                dec,
                x,
                n,
                mbits,
                b_count,
                lowband,
                lm,
                lowband_out,
                1.0,
                lowband_scratch,
                fill,
            );
            cm |= cm2;
        }
    }

    if ctx.resynth {
        if n != 2 {
            stereo_merge(x, y, mid, n);
        }
        if inv != 0 {
            for v in y.iter_mut().take(n as usize) {
                *v = -*v;
            }
        }
    }

    cm
}
