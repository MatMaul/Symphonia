use super::modes::CeltCoef;
use super::types::CeltSig;

pub const COMBFILTER_MAXPERIOD: i32 = 1024;
pub const COMBFILTER_MINPERIOD: i32 = 15;

const GAINS: [[CeltSig; 3]; 3] = [
    [0.3066406250, 0.2170410156, 0.1296386719],
    [0.4638671875, 0.2680664062, 0.0],
    [0.7998046875, 0.1000976562, 0.0],
];

fn comb_filter_const(
    buf: &mut [CeltSig],
    start: usize,
    t: usize,
    n: usize,
    g10: CeltSig,
    g11: CeltSig,
    g12: CeltSig,
) {
    if n == 0 {
        return;
    }
    debug_assert!(start >= t + 2);
    debug_assert!(start + n <= buf.len());

    let mut x4 = buf[start - t - 2];
    let mut x3 = buf[start - t - 1];
    let mut x2 = buf[start - t];
    let mut x1 = buf[start - t + 1];
    for i in 0..n {
        let x0 = buf[start + i - t + 2];
        let idx = start + i;
        let y = buf[idx] + g10 * x2 + g11 * (x1 + x3) + g12 * (x0 + x4);
        buf[idx] = y;
        x4 = x3;
        x3 = x2;
        x2 = x1;
        x1 = x0;
    }
}

pub fn comb_filter(
    buf: &mut [CeltSig],
    start: usize,
    t0: i32,
    t1: i32,
    n: usize,
    g0: CeltSig,
    g1: CeltSig,
    tapset0: i32,
    tapset1: i32,
    window: &[CeltCoef],
    overlap: usize,
) {
    if g0 == 0.0 && g1 == 0.0 {
        return;
    }
    let tapset0 = tapset0 as usize;
    let tapset1 = tapset1 as usize;
    debug_assert!(tapset0 < GAINS.len() && tapset1 < GAINS.len());

    let t0 = t0.max(COMBFILTER_MINPERIOD) as usize;
    let t1 = t1.max(COMBFILTER_MINPERIOD) as usize;
    let mut overlap = overlap.min(n);
    if g0 == g1 && t0 == t1 && tapset0 == tapset1 {
        overlap = 0;
    }
    debug_assert!(start >= t0 + 2);
    debug_assert!(start >= t1 + 2);
    debug_assert!(start + n <= buf.len());

    let g00 = g0 * GAINS[tapset0][0];
    let g01 = g0 * GAINS[tapset0][1];
    let g02 = g0 * GAINS[tapset0][2];
    let g10 = g1 * GAINS[tapset1][0];
    let g11 = g1 * GAINS[tapset1][1];
    let g12 = g1 * GAINS[tapset1][2];

    let mut x1 = buf[start + 1 - t1];
    let mut x2 = buf[start - t1];
    let mut x3 = buf[start - t1 - 1];
    let mut x4 = buf[start - t1 - 2];
    for i in 0..overlap {
        let x0 = buf[start + i + 2 - t1];
        let f = window[i] * window[i];
        let one_minus_f = 1.0 - f;
        let idx = start + i;
        let y = buf[idx]
            + (one_minus_f * g00) * buf[start + i - t0]
            + (one_minus_f * g01)
                * (buf[start + i - t0 + 1] + buf[start + i - t0 - 1])
            + (one_minus_f * g02)
                * (buf[start + i - t0 + 2] + buf[start + i - t0 - 2])
            + (f * g10) * x2
            + (f * g11) * (x1 + x3)
            + (f * g12) * (x0 + x4);
        buf[idx] = y;
        x4 = x3;
        x3 = x2;
        x2 = x1;
        x1 = x0;
    }
    if g1 == 0.0 {
        return;
    }
    comb_filter_const(buf, start + overlap, t1, n - overlap, g10, g11, g12);
}
