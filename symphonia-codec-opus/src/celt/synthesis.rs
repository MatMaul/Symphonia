use alloc::vec;

use super::bands::denormalise_bands;
use super::mdct::mdct_backward;
use super::modes::CeltMode;
use super::types::{CeltGlog, CeltNorm, CeltSig};

pub fn celt_synthesis(
    mode: &CeltMode,
    x: &[CeltNorm],
    out_syn: &mut [&mut [CeltSig]],
    old_band_e: &[CeltGlog],
    start: i32,
    eff_end: i32,
    channels: i32,
    cc: i32,
    is_transient: bool,
    lm: i32,
    downsample: i32,
    silence: bool,
) {
    let overlap = mode.overlap as usize;
    let nb_ebands = mode.nb_ebands as usize;
    let n = (mode.short_mdct_size << lm) as usize;
    let m = 1 << lm;

    let (b, nb, shift) = if is_transient {
        (m, mode.short_mdct_size as usize, mode.max_lm as usize)
    } else {
        (1, (mode.short_mdct_size << lm) as usize, (mode.max_lm - lm) as usize)
    };

    let mut freq = vec![0.0f32; n];

    if cc == 2 && channels == 1 {
        let mut freq2 = vec![0.0f32; n];
        denormalise_bands(
            mode,
            x,
            &mut freq,
            old_band_e,
            start,
            eff_end,
            m,
            downsample,
            silence,
        );
        freq2.copy_from_slice(&freq);
        for block in 0..(b as usize) {
            let offset = nb * block;
            mdct_backward(
                &mode.mdct,
                &freq2[block..],
                &mut out_syn[0][offset..],
                mode.window,
                overlap,
                shift,
                b as usize,
            );
        }
        for block in 0..(b as usize) {
            let offset = nb * block;
            mdct_backward(
                &mode.mdct,
                &freq[block..],
                &mut out_syn[1][offset..],
                mode.window,
                overlap,
                shift,
                b as usize,
            );
        }
    } else if cc == 1 && channels == 2 {
        let mut freq2 = vec![0.0f32; n];
        denormalise_bands(
            mode,
            x,
            &mut freq,
            old_band_e,
            start,
            eff_end,
            m,
            downsample,
            silence,
        );
        denormalise_bands(
            mode,
            &x[n..],
            &mut freq2,
            &old_band_e[nb_ebands..],
            start,
            eff_end,
            m,
            downsample,
            silence,
        );
        for i in 0..n {
            freq[i] = 0.5 * (freq[i] + freq2[i]);
        }
        for block in 0..(b as usize) {
            let offset = nb * block;
            mdct_backward(
                &mode.mdct,
                &freq[block..],
                &mut out_syn[0][offset..],
                mode.window,
                overlap,
                shift,
                b as usize,
            );
        }
    } else {
        for ch in 0..(cc as usize) {
            denormalise_bands(
                mode,
                &x[ch * n..],
                &mut freq,
                &old_band_e[ch * nb_ebands..],
                start,
                eff_end,
                m,
                downsample,
                silence,
            );
            for block in 0..(b as usize) {
                let offset = nb * block;
                mdct_backward(
                    &mode.mdct,
                    &freq[block..],
                    &mut out_syn[ch][offset..],
                    mode.window,
                    overlap,
                    shift,
                    b as usize,
                );
            }
        }
    }
}
