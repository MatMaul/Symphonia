use super::fft::{KissFftCpx, opus_fft_impl};
use super::modes::{CeltCoef, MdctLookup};
use alloc::vec;

pub fn mdct_forward(
    l: &MdctLookup,
    input: &[CeltCoef],
    output: &mut [CeltCoef],
    window: &[CeltCoef],
    overlap: usize,
    shift: usize,
    stride: usize,
) {
    let mut n = l.n as usize;
    let mut trig_offset = 0usize;
    for _ in 0..shift {
        n >>= 1;
        trig_offset += n;
    }
    let trig = &l.trig[trig_offset..];
    let n2 = n >> 1;
    let n4 = n >> 2;
    let st = l.kfft[shift];
    let scale = st.scale;

    debug_assert!(input.len() >= n2 + overlap);
    debug_assert!(window.len() >= overlap);
    debug_assert!(output.len() >= stride * (n2 - 1) + 1);

    let mut f = vec![0.0; n2];
    let mut f2 = vec![KissFftCpx::default(); n4];

    let overlap2 = overlap >> 1;
    let mut xp1 = overlap2;
    let mut xp2 = n2 - 1 + overlap2;
    let mut yp = 0usize;
    let mut wp1 = overlap2;
    let mut wp2 = overlap2 - 1;

    let overlap_quarter = (overlap + 3) >> 2;
    let mut i = 0usize;
    while i < overlap_quarter {
        f[yp] = input[xp1 + n2] * window[wp2] + input[xp2] * window[wp1];
        f[yp + 1] = input[xp1] * window[wp1] - input[xp2 - n2] * window[wp2];
        xp1 += 2;
        xp2 -= 2;
        wp1 += 2;
        wp2 -= 2;
        yp += 2;
        i += 1;
    }

    wp1 = 0;
    wp2 = overlap - 1;
    while i < n4 - overlap_quarter {
        f[yp] = input[xp2];
        f[yp + 1] = input[xp1];
        xp1 += 2;
        xp2 -= 2;
        yp += 2;
        i += 1;
    }

    while i < n4 {
        f[yp] = -input[xp1 - n2] * window[wp1] + input[xp2] * window[wp2];
        f[yp + 1] = input[xp1] * window[wp2] + input[xp2 + n2] * window[wp1];
        xp1 += 2;
        xp2 -= 2;
        wp1 += 2;
        wp2 -= 2;
        yp += 2;
        i += 1;
    }

    for i in 0..n4 {
        let t0 = trig[i];
        let t1 = trig[n4 + i];
        let re = f[2 * i];
        let im = f[2 * i + 1];
        let yr = re * t0 - im * t1;
        let yi = im * t0 + re * t1;
        let idx = st.bitrev[i] as usize;
        f2[idx] = KissFftCpx { r: yr * scale, i: yi * scale };
    }

    opus_fft_impl(st, &mut f2);

    let mut yp1 = 0usize;
    let mut yp2 = stride * (n2 - 1);
    for i in 0..n4 {
        let t0 = trig[i];
        let t1 = trig[n4 + i];
        let yr = f2[i].i * t1 - f2[i].r * t0;
        let yi = f2[i].r * t1 + f2[i].i * t0;
        output[yp1] = yr;
        output[yp2] = yi;
        yp1 += 2 * stride;
        yp2 -= 2 * stride;
    }
}

pub fn mdct_backward(
    l: &MdctLookup,
    input: &[CeltCoef],
    output: &mut [CeltCoef],
    window: &[CeltCoef],
    overlap: usize,
    shift: usize,
    stride: usize,
) {
    let mut n = l.n as usize;
    let mut trig_offset = 0usize;
    for _ in 0..shift {
        n >>= 1;
        trig_offset += n;
    }
    let trig = &l.trig[trig_offset..];
    let n2 = n >> 1;
    let n4 = n >> 2;
    let st = l.kfft[shift];

    debug_assert!(input.len() >= stride * (n2 - 1) + 1);
    debug_assert!(output.len() >= n2 + overlap);
    debug_assert!(window.len() >= overlap);

    let mut f2 = vec![KissFftCpx::default(); n4];

    let mut xp1 = 0usize;
    let mut xp2 = stride * (n2 - 1);
    for i in 0..n4 {
        let t0 = trig[i];
        let t1 = trig[n4 + i];
        let x1 = input[xp1];
        let x2 = input[xp2];
        let yr = x2 * t0 + x1 * t1;
        let yi = x1 * t0 - x2 * t1;
        let rev = st.bitrev[i] as usize;
        f2[rev] = KissFftCpx { r: yi, i: yr };
        xp1 += 2 * stride;
        if i + 1 < n4 {
            xp2 -= 2 * stride;
        }
    }

    opus_fft_impl(st, &mut f2);

    let mut tmp = vec![0.0; n2];
    for i in 0..n4 {
        tmp[2 * i] = f2[i].r;
        tmp[2 * i + 1] = f2[i].i;
    }

    let mut yp0 = 0usize;
    let mut yp1 = n2 - 2;
    for i in 0..((n4 + 1) >> 1) {
        let re = tmp[yp0 + 1];
        let im = tmp[yp0];
        let t0 = trig[i];
        let t1 = trig[n4 + i];
        let yr = re * t0 + im * t1;
        let yi = re * t1 - im * t0;

        let re2 = tmp[yp1 + 1];
        let im2 = tmp[yp1];
        tmp[yp0] = yr;
        tmp[yp1 + 1] = yi;

        let t0b = trig[n4 - i - 1];
        let t1b = trig[n2 - i - 1];
        let yr2 = re2 * t0b + im2 * t1b;
        let yi2 = re2 * t1b - im2 * t0b;
        tmp[yp1] = yr2;
        tmp[yp0 + 1] = yi2;

        yp0 += 2;
        yp1 -= 2;
    }

    let out_offset = overlap >> 1;
    output[out_offset..out_offset + n2].copy_from_slice(&tmp[..n2]);

    let mut xp = overlap - 1;
    let mut yp = 0usize;
    let mut wp1 = 0usize;
    let mut wp2 = overlap - 1;
    for _ in 0..overlap / 2 {
        let x1 = output[xp];
        let x2 = output[yp];
        output[yp] = x2 * window[wp2] - x1 * window[wp1];
        output[xp] = x2 * window[wp1] + x1 * window[wp2];
        xp -= 1;
        yp += 1;
        wp1 += 1;
        wp2 -= 1;
    }
}
