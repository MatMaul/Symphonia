use super::entcode::BITRES;
use super::entdec::EcDec;
use super::modes::CeltMode;
use super::rate::{bits2pulses, get_pulses, pulses2bits, QTHETA_OFFSET, QTHETA_OFFSET_TWOPHASE};
use super::types::{CeltEner, CeltNorm};
use super::{
    alg_unquant, bitexact_cos, bitexact_log2tan, celt_lcg_rand, compute_qn, isqrt32,
    renormalise_vector, deinterleave_hadamard, interleave_hadamard, haar1, stereo_merge,
    SPREAD_AGGRESSIVE,
};
use super::math::celt_sqrt;
use super::intrin::{celt_sudiv, celt_udiv};
use alloc::vec;
use alloc::vec::Vec;

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

fn special_hybrid_folding(
    mode: &CeltMode,
    norm: &mut [CeltNorm],
    mut norm2: Option<&mut [CeltNorm]>,
    start: i32,
    m: i32,
    dual_stereo: bool,
) {
    let start = start as usize;
    let nb_ebands = mode.nb_ebands as usize;
    if start + 2 >= nb_ebands {
        return;
    }
    let ebands = mode.ebands;
    let n1 = (m * (ebands[start + 1] as i32 - ebands[start] as i32)) as usize;
    let n2 = (m * (ebands[start + 2] as i32 - ebands[start + 1] as i32)) as usize;
    if n2 <= n1 || n2 > norm.len() {
        return;
    }
    let src_start = 2 * n1 - n2;
    let src_end = src_start + (n2 - n1);
    if src_end > norm.len() {
        return;
    }
    let (norm_prefix, norm_suffix) = norm.split_at_mut(n1);
    norm_suffix[..(n2 - n1)].copy_from_slice(&norm_prefix[src_start..src_end]);
    if dual_stereo {
        if let Some(norm2) = norm2.as_deref_mut() {
            let (norm2_prefix, norm2_suffix) = norm2.split_at_mut(n1);
            norm2_suffix[..(n2 - n1)].copy_from_slice(&norm2_prefix[src_start..src_end]);
        }
    }
}

pub fn quant_all_bands_decode(
    mode: &CeltMode,
    start: i32,
    end: i32,
    x: &mut [CeltNorm],
    mut y: Option<&mut [CeltNorm]>,
    collapse_masks: &mut [u8],
    pulses: &[i32],
    short_blocks: bool,
    spread: i32,
    mut dual_stereo: bool,
    intensity: i32,
    tf_res: &[i32],
    total_bits: i32,
    balance: i32,
    dec: &mut EcDec<'_>,
    lm: i32,
    coded_bands: i32,
    seed: &mut u32,
    disable_inv: bool,
) {
    let m = 1 << lm;
    let b_blocks = if short_blocks { m } else { 1 };
    let channels = if y.is_some() { 2 } else { 1 };
    if channels == 1 {
        dual_stereo = false;
    }
    let ebands = mode.ebands;
    let nb_ebands = mode.nb_ebands as usize;

    debug_assert!(start >= 0 && end <= mode.nb_ebands);
    debug_assert!(pulses.len() >= end as usize);
    debug_assert!(tf_res.len() >= end as usize);
    debug_assert!(collapse_masks.len() >= channels as usize * nb_ebands);

    let norm_offset = m * ebands[start as usize] as i32;
    let last_band_start = m * ebands[nb_ebands - 1] as i32;
    let norm_len = (last_band_start - norm_offset).max(0) as usize;
    let mut norm = vec![0.0f32; norm_len];
    let mut norm2 = if channels == 2 {
        vec![0.0f32; norm_len]
    } else {
        Vec::new()
    };

    let mut max_n = 0i32;
    for i in start..end {
        let width = ebands[i as usize + 1] as i32 - ebands[i as usize] as i32;
        max_n = max_n.max(m * width);
    }
    let scratch_len = max_n.max(0) as usize;
    let mut band_scratch = vec![0.0f32; scratch_len];
    let mut band_scratch_y = if channels == 2 {
        vec![0.0f32; scratch_len]
    } else {
        Vec::new()
    };
    let mut lowband_scratch = vec![0.0f32; scratch_len];

    let mut ctx = BandCtx {
        m: mode,
        i: start,
        intensity,
        spread,
        tf_change: 0,
        remaining_bits: 0,
        band_e: None,
        seed: *seed,
        resynth: true,
        disable_inv,
        avoid_split_noise: b_blocks > 1,
    };

    let mut lowband_offset = 0;
    let mut update_lowband = true;
    let mut balance = balance;

    for i in start..end {
        ctx.i = i;
        let last = i == end - 1;
        let band_start = (m * ebands[i as usize] as i32) as usize;
        let band_end = (m * ebands[i as usize + 1] as i32) as usize;
        let n = (band_end - band_start) as i32;
        debug_assert!(n > 0);

        let tell = dec.tell_frac() as i32;
        if i != start {
            balance -= tell;
        }
        let remaining_bits = total_bits - tell - 1;
        ctx.remaining_bits = remaining_bits;

        let mut b = 0i32;
        if i <= coded_bands - 1 {
            let denom = (coded_bands - i).min(3);
            let curr_balance = celt_sudiv(balance, denom);
            let max_bits = (remaining_bits + 1).min(16383);
            let mut b_tmp = pulses[i as usize] + curr_balance;
            b_tmp = b_tmp.min(max_bits);
            if b_tmp < 0 {
                b_tmp = 0;
            }
            b = b_tmp;
        }

        if ctx.resynth
            && (m * ebands[i as usize] as i32 - n >= m * ebands[start as usize] as i32
                || i == start + 1)
            && (update_lowband || lowband_offset == 0)
        {
            lowband_offset = i;
        }
        if i == start + 1 {
            special_hybrid_folding(
                mode,
                &mut norm,
                if channels == 2 { Some(&mut norm2) } else { None },
                start,
                m,
                dual_stereo,
            );
        }

        ctx.tf_change = tf_res[i as usize];

        let mut effective_lowband = -1;
        let (mut x_cm, mut y_cm);
        if lowband_offset != 0
            && (spread != SPREAD_AGGRESSIVE || b_blocks > 1 || ctx.tf_change < 0)
        {
            effective_lowband =
                (m * ebands[lowband_offset as usize] as i32 - norm_offset - n).max(0);
            let mut fold_start = lowband_offset;
            while fold_start > 0 {
                let prev = fold_start - 1;
                if m * ebands[prev as usize] as i32 <= effective_lowband + norm_offset {
                    break;
                }
                fold_start = prev;
            }
            let mut fold_end = lowband_offset - 1;
            loop {
                fold_end += 1;
                if fold_end >= i {
                    break;
                }
                if m * ebands[fold_end as usize] as i32 >= effective_lowband + norm_offset + n {
                    break;
                }
            }
            x_cm = 0;
            y_cm = 0;
            for fold_i in fold_start..fold_end {
                let base = fold_i as usize * channels as usize;
                x_cm |= collapse_masks[base] as u32;
                y_cm |= collapse_masks[base + channels as usize - 1] as u32;
            }
        } else {
            x_cm = (1u32 << (b_blocks as u32)) - 1;
            y_cm = x_cm;
        }

        if dual_stereo && i == intensity {
            dual_stereo = false;
            if ctx.resynth && !norm.is_empty() {
                let split =
                    (m * ebands[i as usize] as i32 - norm_offset).max(0) as usize;
                let len = split.min(norm.len());
                for j in 0..len {
                    norm[j] = 0.5 * (norm[j] + norm2[j]);
                }
            }
        }

        let use_output = i < mode.eff_ebands;
        let mut x_band: &mut [CeltNorm];
        let mut y_band: Option<&mut [CeltNorm]> = None;
        if use_output {
            debug_assert!(band_end <= x.len());
            x_band = &mut x[band_start..band_end];
            if let Some(y_buf) = y.as_deref_mut() {
                debug_assert!(band_end <= y_buf.len());
                y_band = Some(&mut y_buf[band_start..band_end]);
            }
        } else {
            debug_assert!(scratch_len >= n as usize);
            x_band = &mut band_scratch[..n as usize];
            if channels == 2 {
                y_band = Some(&mut band_scratch_y[..n as usize]);
            }
        }

        let use_lowband_scratch = use_output && !last && scratch_len >= n as usize;
        let lowband_out_offset = (band_start as i32 - norm_offset).max(0) as usize;

        if dual_stereo {
            let x_cm_out = {
                let (lowband, lowband_out) = if last || lowband_out_offset >= norm.len() {
                    let lowband = if effective_lowband >= 0
                        && (effective_lowband as usize) < norm.len()
                    {
                        Some(&norm[effective_lowband as usize..])
                    } else {
                        None
                    };
                    (lowband, None)
                } else {
                    let (norm_prefix, norm_suffix) = norm.split_at_mut(lowband_out_offset);
                    let lowband = if effective_lowband >= 0
                        && (effective_lowband as usize) < norm_prefix.len()
                    {
                        Some(&norm_prefix[effective_lowband as usize..])
                    } else {
                        None
                    };
                    debug_assert!(norm_suffix.len() >= n as usize);
                    (lowband, Some(&mut norm_suffix[..n as usize]))
                };
                let scratch = if use_lowband_scratch {
                    Some(&mut lowband_scratch[..n as usize])
                } else {
                    None
                };
                quant_band_decode(
                    &mut ctx,
                    dec,
                    x_band,
                    n,
                    b / 2,
                    b_blocks,
                    lowband,
                    lm,
                    lowband_out,
                    1.0,
                    scratch,
                    x_cm as i32,
                )
            };
            x_cm = x_cm_out;

            if let Some(y_band) = y_band.as_deref_mut() {
                let y_cm_out = {
                    let (lowband, lowband_out) = if last || lowband_out_offset >= norm2.len() {
                        let lowband = if effective_lowband >= 0
                            && (effective_lowband as usize) < norm2.len()
                        {
                            Some(&norm2[effective_lowband as usize..])
                        } else {
                            None
                        };
                        (lowband, None)
                    } else {
                        let (norm_prefix, norm_suffix) = norm2.split_at_mut(lowband_out_offset);
                        let lowband = if effective_lowband >= 0
                            && (effective_lowband as usize) < norm_prefix.len()
                        {
                            Some(&norm_prefix[effective_lowband as usize..])
                        } else {
                            None
                        };
                        debug_assert!(norm_suffix.len() >= n as usize);
                        (lowband, Some(&mut norm_suffix[..n as usize]))
                    };
                    let scratch = if use_lowband_scratch {
                        Some(&mut lowband_scratch[..n as usize])
                    } else {
                        None
                    };
                    quant_band_decode(
                        &mut ctx,
                        dec,
                        y_band,
                        n,
                        b / 2,
                        b_blocks,
                        lowband,
                        lm,
                        lowband_out,
                        1.0,
                        scratch,
                        y_cm as i32,
                    )
                };
                y_cm = y_cm_out;
            } else {
                y_cm = x_cm;
            }
        } else if let Some(y_band) = y_band.as_deref_mut() {
            let cm = {
                let (lowband, lowband_out) = if last || lowband_out_offset >= norm.len() {
                    let lowband = if effective_lowband >= 0
                        && (effective_lowband as usize) < norm.len()
                    {
                        Some(&norm[effective_lowband as usize..])
                    } else {
                        None
                    };
                    (lowband, None)
                } else {
                    let (norm_prefix, norm_suffix) = norm.split_at_mut(lowband_out_offset);
                    let lowband = if effective_lowband >= 0
                        && (effective_lowband as usize) < norm_prefix.len()
                    {
                        Some(&norm_prefix[effective_lowband as usize..])
                    } else {
                        None
                    };
                    debug_assert!(norm_suffix.len() >= n as usize);
                    (lowband, Some(&mut norm_suffix[..n as usize]))
                };
                let scratch = if use_lowband_scratch {
                    Some(&mut lowband_scratch[..n as usize])
                } else {
                    None
                };
                quant_band_stereo_decode(
                    &mut ctx,
                    dec,
                    x_band,
                    y_band,
                    n,
                    b,
                    b_blocks,
                    lowband,
                    lm,
                    lowband_out,
                    scratch,
                    (x_cm | y_cm) as i32,
                )
            };
            x_cm = cm;
            y_cm = cm;
        } else {
            let cm = {
                let (lowband, lowband_out) = if last || lowband_out_offset >= norm.len() {
                    let lowband = if effective_lowband >= 0
                        && (effective_lowband as usize) < norm.len()
                    {
                        Some(&norm[effective_lowband as usize..])
                    } else {
                        None
                    };
                    (lowband, None)
                } else {
                    let (norm_prefix, norm_suffix) = norm.split_at_mut(lowband_out_offset);
                    let lowband = if effective_lowband >= 0
                        && (effective_lowband as usize) < norm_prefix.len()
                    {
                        Some(&norm_prefix[effective_lowband as usize..])
                    } else {
                        None
                    };
                    debug_assert!(norm_suffix.len() >= n as usize);
                    (lowband, Some(&mut norm_suffix[..n as usize]))
                };
                let scratch = if use_lowband_scratch {
                    Some(&mut lowband_scratch[..n as usize])
                } else {
                    None
                };
                quant_band_decode(
                    &mut ctx,
                    dec,
                    x_band,
                    n,
                    b,
                    b_blocks,
                    lowband,
                    lm,
                    lowband_out,
                    1.0,
                    scratch,
                    (x_cm | y_cm) as i32,
                )
            };
            x_cm = cm;
            y_cm = cm;
        }

        let mask_base = i as usize * channels as usize;
        collapse_masks[mask_base] = x_cm as u8;
        collapse_masks[mask_base + channels as usize - 1] = y_cm as u8;
        balance += pulses[i as usize] + tell;

        update_lowband = b > (n << BITRES);
        ctx.avoid_split_noise = false;
    }

    *seed = ctx.seed;
}
