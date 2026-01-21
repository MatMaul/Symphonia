use super::bands::anti_collapse;
use super::bands_quant::quant_all_bands_decode;
use super::entcode::BITRES;
use super::entdec::EcDec;
use super::modes::CeltMode;
use super::quant_bands::{unquant_coarse_energy, unquant_energy_finalise, unquant_fine_energy};
use super::rate::{clt_compute_allocation, init_caps};
use super::synthesis::celt_synthesis;
use super::pitch::{comb_filter, COMBFILTER_MINPERIOD};
use super::types::{CeltGlog, CeltSig, VERY_SMALL};
use super::vq::SPREAD_NORMAL;
use alloc::vec::Vec;

pub const TRIM_ICDF: [u8; 11] = [126, 124, 119, 109, 87, 41, 19, 9, 4, 2, 0];
pub const SPREAD_ICDF: [u8; 4] = [25, 23, 2, 0];
const TAPSET_ICDF: [u8; 3] = [2, 1, 0];
const DECODE_BUFFER_SIZE: usize = 2048;
const POSTFILTER_GAIN_STEP: f32 = 0.09375;

const TF_SELECT_TABLE: [[i8; 8]; 4] = [
    [0, -1, 0, -1, 0, -1, 0, -1],
    [0, -1, 0, -2, 1, 0, 1, -1],
    [0, -2, 0, -3, 2, 0, 1, -1],
    [0, -2, 0, -3, 3, 0, 1, -1],
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CeltDecodeError {
    InvalidFrameSize,
    BufferTooSmall,
    InvalidPacket,
}

pub struct CeltDecoder {
    mode: &'static CeltMode,
    channels: i32,
    stream_channels: i32,
    downsample: i32,
    start: i32,
    end: i32,
    disable_inv: bool,
    rng: u32,
    old_band_e: Vec<CeltGlog>,
    old_log_e: Vec<CeltGlog>,
    old_log_e2: Vec<CeltGlog>,
    preemph_mem: [CeltSig; 2],
    decode_mem: Vec<CeltSig>,
    postfilter_period: i32,
    postfilter_period_old: i32,
    postfilter_gain: CeltSig,
    postfilter_gain_old: CeltSig,
    postfilter_tapset: i32,
    postfilter_tapset_old: i32,
}

impl CeltDecoder {
    pub fn new(mode: &'static CeltMode, channels: i32) -> Result<Self, CeltDecodeError> {
        if channels != 1 && channels != 2 {
            return Err(CeltDecodeError::InvalidPacket);
        }
        let nb_ebands = mode.nb_ebands as usize;
        let overlap = mode.overlap as usize;
        let decode_mem_len = channels as usize * (DECODE_BUFFER_SIZE + overlap);
        let mut decoder = Self {
            mode,
            channels,
            stream_channels: channels,
            downsample: 1,
            start: 0,
            end: mode.nb_ebands,
            disable_inv: false,
            rng: 0,
            old_band_e: vec![0.0f32; 2 * nb_ebands],
            old_log_e: vec![-28.0f32; 2 * nb_ebands],
            old_log_e2: vec![-28.0f32; 2 * nb_ebands],
            preemph_mem: [0.0f32; 2],
            decode_mem: vec![0.0f32; decode_mem_len],
            postfilter_period: 0,
            postfilter_period_old: 0,
            postfilter_gain: 0.0,
            postfilter_gain_old: 0.0,
            postfilter_tapset: 0,
            postfilter_tapset_old: 0,
        };
        decoder.reset();
        Ok(decoder)
    }

    pub fn reset(&mut self) {
        let nb_ebands = self.mode.nb_ebands as usize;
        self.rng = 0;
        self.old_band_e.fill(0.0);
        self.old_log_e[..2 * nb_ebands].fill(-28.0);
        self.old_log_e2[..2 * nb_ebands].fill(-28.0);
        self.preemph_mem = [0.0; 2];
        self.decode_mem.fill(0.0);
        self.postfilter_period = 0;
        self.postfilter_period_old = 0;
        self.postfilter_gain = 0.0;
        self.postfilter_gain_old = 0.0;
        self.postfilter_tapset = 0;
        self.postfilter_tapset_old = 0;
    }

    pub fn set_stream_channels(&mut self, channels: i32) {
        let channels = channels.clamp(1, self.channels);
        self.stream_channels = channels;
    }

    pub fn set_end_band(&mut self, end: i32) {
        let max_end = self.mode.nb_ebands;
        let min_end = self.start + 1;
        self.end = end.clamp(min_end, max_end);
    }

    pub fn decode_frame(
        &mut self,
        data: &[u8],
        out: &mut [&mut [CeltSig]],
        frame_size: i32,
    ) -> Result<i32, CeltDecodeError> {
        if data.is_empty() {
            return Err(CeltDecodeError::InvalidPacket);
        }
        if out.len() < self.channels as usize {
            return Err(CeltDecodeError::BufferTooSmall);
        }

        let mut lm = 0i32;
        while lm <= self.mode.max_lm {
            if (self.mode.short_mdct_size << (lm as u32)) == frame_size {
                break;
            }
            lm += 1;
        }
        if lm > self.mode.max_lm {
            return Err(CeltDecodeError::InvalidFrameSize);
        }

        let n = (self.mode.short_mdct_size << (lm as u32)) as usize;
        if n > DECODE_BUFFER_SIZE {
            return Err(CeltDecodeError::InvalidFrameSize);
        }
        let overlap = self.mode.overlap as usize;
        let decode_stride = DECODE_BUFFER_SIZE + overlap;
        let frame_start = DECODE_BUFFER_SIZE - n;
        for ch in 0..(self.channels as usize) {
            if out[ch].len() < n {
                return Err(CeltDecodeError::BufferTooSmall);
            }
        }

        let start = self.start;
        let end = self.end;
        let eff_end = end.min(self.mode.eff_ebands);
        let nb_ebands = self.mode.nb_ebands as usize;
        let mut dec = EcDec::new(data);
        let bitres = BITRES as u32;
        let total_bits = data.len() as i32 * 8;
        let mut tell = dec.tell();

        if self.stream_channels == 1 {
            for i in 0..nb_ebands {
                let idx = i + nb_ebands;
                if self.old_band_e[i] < self.old_band_e[idx] {
                    self.old_band_e[i] = self.old_band_e[idx];
                }
            }
        }

        let mut silence = false;
        if tell >= total_bits {
            silence = true;
        } else if tell == 1 {
            silence = dec.dec_bit_logp(15) != 0;
        }
        if silence {
            let remaining = (total_bits - dec.tell()).max(0) as u32;
            if remaining > 0 {
                dec.dec_bits(remaining);
            }
            tell = total_bits;
        }

        let mut postfilter_gain = 0.0f32;
        let mut postfilter_pitch = 0i32;
        let mut postfilter_tapset = 0i32;
        if start == 0 && tell + 16 <= total_bits {
            if dec.dec_bit_logp(1) != 0 {
                let octave = dec.dec_uint(6) as u32;
                postfilter_pitch =
                    ((16u32 << octave) + dec.dec_bits(4 + octave)).saturating_sub(1) as i32;
                let qg = dec.dec_bits(3) as i32;
                if dec.tell() + 2 <= total_bits {
                    postfilter_tapset = dec.dec_icdf(&TAPSET_ICDF, 2);
                }
                postfilter_gain = POSTFILTER_GAIN_STEP * (qg as f32 + 1.0);
            }
            tell = dec.tell();
        }

        let mut is_transient = 0;
        if lm > 0 && tell + 3 <= total_bits {
            is_transient = dec.dec_bit_logp(3);
            tell = dec.tell();
        }
        let short_blocks = is_transient != 0;

        let intra_ener = if tell + 3 <= total_bits {
            dec.dec_bit_logp(3) != 0
        } else {
            false
        };

        unquant_coarse_energy(
            self.mode,
            start,
            end,
            &mut self.old_band_e,
            intra_ener,
            &mut dec,
            self.stream_channels,
            lm,
        );

        let mut tf_res = vec![0i32; nb_ebands];
        tf_decode(start, end, is_transient != 0, &mut tf_res, lm, &mut dec);

        let mut spread_decision = SPREAD_NORMAL;
        if dec.tell() + 4 <= total_bits {
            spread_decision = dec.dec_icdf(&SPREAD_ICDF, 5);
        }

        let mut cap = vec![0i32; nb_ebands];
        init_caps(self.mode, &mut cap, lm, self.stream_channels);

        let mut offsets = vec![0i32; nb_ebands];
        let mut dynalloc_logp = 6i32;
        let mut total_bits_q3 = total_bits << bitres;
        let mut tell_frac = dec.tell_frac() as i32;
        for i in start..end {
            let idx = i as usize;
            let width = self.stream_channels
                * (self.mode.ebands[idx + 1] as i32 - self.mode.ebands[idx] as i32)
                << (lm as u32);
            let quanta = (width << bitres).min((6 << bitres).max(width));
            let mut dynalloc_loop_logp = dynalloc_logp;
            let mut boost = 0i32;
            while tell_frac + (dynalloc_loop_logp << bitres) < total_bits_q3
                && boost < cap[idx]
            {
                let flag = dec.dec_bit_logp(dynalloc_loop_logp as u32);
                tell_frac = dec.tell_frac() as i32;
                if flag == 0 {
                    break;
                }
                boost += quanta;
                total_bits_q3 -= quanta;
                dynalloc_loop_logp = 1;
            }
            offsets[idx] = boost;
            if boost > 0 {
                dynalloc_logp = (dynalloc_logp - 1).max(2);
            }
        }

        let alloc_trim = if tell_frac + (6 << bitres) <= total_bits_q3 {
            dec.dec_icdf(&TRIM_ICDF, 7)
        } else {
            5
        };

        let mut bits =
            ((data.len() as i32 * 8) << bitres) - dec.tell_frac() as i32 - 1;
        let anti_collapse_rsv =
            if short_blocks && lm >= 2 && bits >= ((lm + 2) << bitres) {
                1 << bitres
            } else {
                0
            };
        bits -= anti_collapse_rsv;

        let mut pulses = vec![0i32; nb_ebands];
        let mut fine_quant = vec![0i32; nb_ebands];
        let mut fine_priority = vec![0i32; nb_ebands];
        let mut intensity = 0;
        let mut dual_stereo = 0;
        let mut balance = 0i32;
        let coded_bands = clt_compute_allocation(
            self.mode,
            start,
            end,
            &offsets,
            &cap,
            alloc_trim,
            &mut intensity,
            &mut dual_stereo,
            bits,
            &mut balance,
            &mut pulses,
            &mut fine_quant,
            &mut fine_priority,
            self.stream_channels,
            lm,
            &mut dec,
        );

        unquant_fine_energy(
            self.mode,
            start,
            end,
            &mut self.old_band_e,
            None,
            &fine_quant,
            &mut dec,
            self.stream_channels,
        );

        self.shift_decode_mem(n);

        let mut x = vec![0.0f32; self.stream_channels as usize * n];
        let (x0, x1) = if self.stream_channels == 2 {
            let (left, right) = x.split_at_mut(n);
            (left, Some(right))
        } else {
            (x.as_mut_slice(), None)
        };
        let mut collapse_masks = vec![0u8; self.stream_channels as usize * nb_ebands];
        quant_all_bands_decode(
            self.mode,
            start,
            end,
            x0,
            x1,
            &mut collapse_masks,
            &pulses,
            short_blocks,
            spread_decision,
            dual_stereo != 0,
            intensity,
            &tf_res,
            (data.len() as i32 * 8 << bitres) - anti_collapse_rsv,
            balance,
            &mut dec,
            lm,
            coded_bands,
            &mut self.rng,
            self.disable_inv,
        );

        let mut anti_collapse_on = false;
        if anti_collapse_rsv > 0 {
            anti_collapse_on = dec.dec_bits(1) != 0;
        }

        unquant_energy_finalise(
            self.mode,
            start,
            end,
            Some(&mut self.old_band_e),
            &fine_quant,
            &fine_priority,
            data.len() as i32 * 8 - dec.tell(),
            &mut dec,
            self.stream_channels,
        );

        if anti_collapse_on {
            anti_collapse(
                self.mode,
                &mut x,
                &collapse_masks,
                lm,
                self.stream_channels,
                n as i32,
                start,
                end,
                &self.old_band_e,
                &self.old_log_e,
                &self.old_log_e2,
                &pulses,
                self.rng,
                false,
            );
        }

        if silence {
            for v in &mut self.old_band_e {
                *v = -28.0;
            }
        }

        {
            let mut out_refs: Vec<&mut [CeltSig]> = self
                .decode_mem
                .chunks_exact_mut(decode_stride)
                .take(self.channels as usize)
                .map(|buf| &mut buf[frame_start..frame_start + n])
                .collect();
            celt_synthesis(
                self.mode,
                &x,
                &mut out_refs,
                &self.old_band_e,
                start,
                eff_end,
                self.stream_channels,
                self.channels,
                short_blocks,
                lm,
                self.downsample,
                silence,
            );
        }

        self.apply_postfilter(
            n,
            frame_start,
            postfilter_pitch,
            postfilter_gain,
            postfilter_tapset,
            lm,
        );

        self.update_energy_state(short_blocks, start, end);
        let mut out_refs: Vec<&mut [CeltSig]> =
            out.iter_mut().take(self.channels as usize).map(|plane| &mut plane[..n]).collect();
        for ch in 0..(self.channels as usize) {
            let base = ch * decode_stride;
            let src = &self.decode_mem[base + frame_start..base + frame_start + n];
            out_refs[ch].copy_from_slice(src);
        }
        self.deemphasis_simple(out_refs.as_mut_slice(), n);

        Ok(n as i32)
    }

    fn shift_decode_mem(&mut self, n: usize) {
        debug_assert!(n <= DECODE_BUFFER_SIZE);
        let overlap = self.mode.overlap as usize;
        let stride = DECODE_BUFFER_SIZE + overlap;
        for ch in 0..(self.channels as usize) {
            let base = ch * stride;
            let buf = &mut self.decode_mem[base..base + stride];
            buf.copy_within(n..stride, 0);
        }
    }

    fn apply_postfilter(
        &mut self,
        n: usize,
        frame_start: usize,
        postfilter_pitch: i32,
        postfilter_gain: CeltSig,
        postfilter_tapset: i32,
        lm: i32,
    ) {
        self.postfilter_period = self.postfilter_period.max(COMBFILTER_MINPERIOD);
        self.postfilter_period_old = self.postfilter_period_old.max(COMBFILTER_MINPERIOD);

        let overlap = self.mode.overlap as usize;
        let stride = DECODE_BUFFER_SIZE + overlap;
        let short_mdct_size = self.mode.short_mdct_size as usize;

        for ch in 0..(self.channels as usize) {
            let base = ch * stride;
            let buf = &mut self.decode_mem[base..base + stride];
            comb_filter(
                buf,
                frame_start,
                self.postfilter_period_old,
                self.postfilter_period,
                short_mdct_size,
                self.postfilter_gain_old,
                self.postfilter_gain,
                self.postfilter_tapset_old,
                self.postfilter_tapset,
                self.mode.window,
                overlap,
            );
            if lm != 0 {
                comb_filter(
                    buf,
                    frame_start + short_mdct_size,
                    self.postfilter_period,
                    postfilter_pitch,
                    n - short_mdct_size,
                    self.postfilter_gain,
                    postfilter_gain,
                    self.postfilter_tapset,
                    postfilter_tapset,
                    self.mode.window,
                    overlap,
                );
            }
        }

        self.postfilter_period_old = self.postfilter_period;
        self.postfilter_gain_old = self.postfilter_gain;
        self.postfilter_tapset_old = self.postfilter_tapset;
        self.postfilter_period = postfilter_pitch;
        self.postfilter_gain = postfilter_gain;
        self.postfilter_tapset = postfilter_tapset;
        if lm != 0 {
            self.postfilter_period_old = self.postfilter_period;
            self.postfilter_gain_old = self.postfilter_gain;
            self.postfilter_tapset_old = self.postfilter_tapset;
        }
    }

    fn update_energy_state(&mut self, is_transient: bool, start: i32, end: i32) {
        let nb_ebands = self.mode.nb_ebands as usize;
        if self.stream_channels == 1 {
            let offset = nb_ebands;
            let (left, right) = self.old_band_e.split_at_mut(offset);
            right[..nb_ebands].copy_from_slice(&left[..nb_ebands]);
        }

        if !is_transient {
            self.old_log_e2.copy_from_slice(&self.old_log_e);
            self.old_log_e.copy_from_slice(&self.old_band_e);
        } else {
            for i in 0..self.old_log_e.len() {
                if self.old_log_e[i] > self.old_band_e[i] {
                    self.old_log_e[i] = self.old_band_e[i];
                }
            }
        }

        for c in 0..2 {
            let base = c * nb_ebands;
            for i in 0..(start as usize) {
                self.old_band_e[base + i] = 0.0;
                self.old_log_e[base + i] = -28.0;
                self.old_log_e2[base + i] = -28.0;
            }
            for i in (end as usize)..nb_ebands {
                self.old_band_e[base + i] = 0.0;
                self.old_log_e[base + i] = -28.0;
                self.old_log_e2[base + i] = -28.0;
            }
        }
    }

    fn deemphasis_simple(&mut self, out: &mut [&mut [CeltSig]], n: usize) {
        let coef0 = self.mode.preemph[0];
        for ch in 0..(self.channels as usize) {
            let mut m = self.preemph_mem[ch];
            for sample in out[ch].iter_mut().take(n) {
                let tmp = *sample + VERY_SMALL + m;
                m = coef0 * tmp;
                *sample = tmp;
            }
            self.preemph_mem[ch] = m;
        }
    }
}

pub fn tf_decode(
    start: i32,
    end: i32,
    is_transient: bool,
    tf_res: &mut [i32],
    lm: i32,
    dec: &mut EcDec<'_>,
) {
    let mut budget = dec.storage() as i32 * 8;
    let mut tell = dec.tell();
    let mut logp = if is_transient { 2 } else { 4 };
    let tf_select_rsv = lm > 0 && tell + logp + 1 <= budget;
    if tf_select_rsv {
        budget -= 1;
    }

    let mut tf_changed = 0;
    let mut curr = 0;
    for i in start..end {
        if tell + logp <= budget {
            curr ^= dec.dec_bit_logp(logp as u32);
            tell = dec.tell();
            tf_changed |= curr;
        }
        tf_res[i as usize] = curr;
        logp = if is_transient { 4 } else { 5 };
    }

    let mut tf_select = 0;
    if tf_select_rsv {
        let base = 4 * (is_transient as i32);
        let idx0 = (base + tf_changed) as usize;
        let idx1 = (base + 2 + tf_changed) as usize;
        if TF_SELECT_TABLE[lm as usize][idx0] != TF_SELECT_TABLE[lm as usize][idx1] {
            tf_select = dec.dec_bit_logp(1);
        }
    }

    let base = 4 * (is_transient as i32);
    for i in start..end {
        let idx = (base + 2 * tf_select + tf_res[i as usize]) as usize;
        tf_res[i as usize] = TF_SELECT_TABLE[lm as usize][idx] as i32;
    }
}
