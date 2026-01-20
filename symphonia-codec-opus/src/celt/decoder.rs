// Symphonia
// Copyright (c) 2024 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! CELT decoder implementation.
//!
//! This module implements the main CELT decoder state machine that orchestrates
//! all the decoding stages: entropy decoding, band energy decoding, vector
//! quantization, IMDCT, and output synthesis.

#![allow(dead_code)]

#[cfg(not(feature = "std"))]
use alloc::{vec, vec::Vec};
#[cfg(feature = "std")]
use std::{vec, vec::Vec};

use crate::celt::bands::{anti_collapse, celt_lcg_rand, quant_all_bands};
use crate::celt::constants::{
    BITRES, COMBFILTER_MINPERIOD, DB_SHIFT, DECODE_BUFFER_SIZE, LPC_ORDER, MAX_FRAME_SIZE, Q15ONE,
    SPREAD_NORMAL,
};
use crate::celt::mode::CeltMode;
use crate::celt::quant_bands::{
    unquant_coarse_energy, unquant_energy_finalise, unquant_fine_energy,
};
use crate::celt::rate::{compute_allocation, init_caps};
use crate::celt::synthesis::{celt_synthesis, comb_filter, deemphasis, tf_decode};
use crate::celt::tables::{SPREAD_ICDF, TAPSET_ICDF, TRIM_ICDF};
use crate::entropy::RangeDecoder;

/// CELT decoder state.
pub struct CeltDecoder {
    /// Mode configuration (48kHz parameters).
    mode: CeltMode,

    /// Overlap size.
    overlap: usize,

    /// Number of channels in output.
    channels: usize,

    /// Number of channels in stream (may differ during stereo intensity coding).
    stream_channels: usize,

    /// Downsample factor (1, 2, 3, 4, or 6 for 48/24/16/12/8 kHz output).
    downsample: usize,

    /// Start band (usually 0).
    start: usize,

    /// End band (usually 21 for fullband).
    end: usize,

    /// Whether signalling bits are present.
    signalling: bool,

    /// Random seed for noise generation.
    rng: u32,

    /// Error flag.
    error: bool,

    /// Last pitch index (for PLC).
    last_pitch_index: usize,

    /// Consecutive lost frame count.
    loss_count: usize,

    /// Post-filter period (current frame).
    postfilter_period: usize,

    /// Post-filter period (previous frame).
    postfilter_period_old: usize,

    /// Post-filter gain (current frame, Q15).
    postfilter_gain: i32,

    /// Post-filter gain (previous frame, Q15).
    postfilter_gain_old: i32,

    /// Post-filter tapset (current frame).
    postfilter_tapset: usize,

    /// Post-filter tapset (previous frame).
    postfilter_tapset_old: usize,

    /// De-emphasis memory (one per channel).
    preemph_mem: Vec<i32>,

    /// Decode buffer (stores previous frame's tail for overlap-add).
    decode_mem: Vec<Vec<i32>>,

    /// LPC coefficients (for PLC).
    lpc: Vec<Vec<i32>>,

    /// Old band energies (for prediction).
    old_ebands: Vec<i32>,

    /// Old log energies (for anti-collapse).
    old_log_e: Vec<i32>,

    /// Old log energies (second generation, for anti-collapse).
    old_log_e2: Vec<i32>,

    /// Background log energies (for noise floor).
    background_log_e: Vec<i32>,
}

impl CeltDecoder {
    /// Create a new CELT decoder.
    ///
    /// # Arguments
    /// * `sampling_rate` - Target output sampling rate (8000, 12000, 16000, 24000, or 48000)
    /// * `channels` - Number of output channels (1 or 2)
    pub fn new(sampling_rate: u32, channels: usize) -> Result<Self, &'static str> {
        if channels < 1 || channels > 2 {
            return Err("channels must be 1 or 2");
        }

        let downsample = resampling_factor(sampling_rate)?;
        let mode = CeltMode::new();
        let nb_ebands = mode.nb_ebands;

        // Initialize state
        let mut decoder = Self {
            mode,
            overlap: 120,
            channels,
            stream_channels: channels,
            downsample,
            start: 0,
            end: nb_ebands,
            signalling: true,
            rng: 0,
            error: false,
            last_pitch_index: 0,
            loss_count: 0,
            postfilter_period: 0,
            postfilter_period_old: 0,
            postfilter_gain: 0,
            postfilter_gain_old: 0,
            postfilter_tapset: 0,
            postfilter_tapset_old: 0,
            preemph_mem: vec![0; channels],
            decode_mem: vec![vec![0; DECODE_BUFFER_SIZE + 120]; channels],
            lpc: vec![vec![0; LPC_ORDER]; channels],
            old_ebands: vec![0; 2 * nb_ebands],
            old_log_e: vec![minus_28db(); 2 * nb_ebands],
            old_log_e2: vec![minus_28db(); 2 * nb_ebands],
            background_log_e: vec![minus_28db(); 2 * nb_ebands],
        };

        Ok(decoder)
    }

    /// Reset the decoder state (e.g., after a seek).
    pub fn reset(&mut self) {
        let nb_ebands = self.mode.nb_ebands;

        self.rng = 0;
        self.error = false;
        self.last_pitch_index = 0;
        self.loss_count = 0;
        self.postfilter_period = 0;
        self.postfilter_period_old = 0;
        self.postfilter_gain = 0;
        self.postfilter_gain_old = 0;
        self.postfilter_tapset = 0;
        self.postfilter_tapset_old = 0;

        for mem in &mut self.preemph_mem {
            *mem = 0;
        }

        for buf in &mut self.decode_mem {
            buf.fill(0);
        }

        for lpc in &mut self.lpc {
            lpc.fill(0);
        }

        self.old_ebands.fill(0);
        self.old_log_e.fill(minus_28db());
        self.old_log_e2.fill(minus_28db());
        self.background_log_e.fill(minus_28db());
    }

    /// Decode a CELT frame.
    ///
    /// # Arguments
    /// * `data` - Encoded CELT packet data (or None for PLC)
    /// * `pcm` - Output buffer for decoded samples (interleaved f32)
    /// * `frame_size` - Number of samples per channel to decode
    ///
    /// # Returns
    /// Number of samples per channel actually decoded
    pub fn decode(
        &mut self,
        data: Option<&[u8]>,
        pcm: &mut [f32],
        frame_size: usize,
    ) -> Result<usize, &'static str> {
        let frame_size = frame_size * self.downsample;

        // Validate frame size
        let lm = self.find_lm(frame_size)?;
        let m = 1usize << lm;
        let n = m * self.mode.short_mdct_size;

        if n > MAX_FRAME_SIZE {
            return Err("frame size too large");
        }

        // Set up output pointers
        let mut out_syn = vec![vec![0i32; DECODE_BUFFER_SIZE]; self.channels];
        let out_syn_offsets: Vec<usize> = vec![DECODE_BUFFER_SIZE - n; self.channels];

        let eff_end = self.end.min(self.mode.eff_ebands);

        // Handle packet loss or missing data
        if data.is_none() || data.map(|d| d.len()).unwrap_or(0) <= 1 {
            self.decode_lost(n, lm);
            self.apply_postfilter_and_output(&mut out_syn, &out_syn_offsets, pcm, n)?;
            return Ok(frame_size / self.downsample);
        }

        let data = data.unwrap();
        let length = data.len();

        if length > 1275 {
            return Err("packet too large");
        }

        // Initialize entropy decoder
        let mut dec = RangeDecoder::init(data);

        // If mono stream, copy energy to both channels
        let c = self.stream_channels;
        let nb_ebands = self.mode.nb_ebands;

        if c == 1 {
            for i in 0..nb_ebands {
                self.old_ebands[i] = self.old_ebands[i].max(self.old_ebands[nb_ebands + i]);
            }
        }

        let total_bits = length as i32 * 8;
        let mut tell = dec.tell();

        // Decode silence flag
        let silence = if tell >= total_bits {
            true
        }
        else if tell == 1 {
            dec.dec_bit_logp(15) != 0
        }
        else {
            false
        };

        if silence {
            // Pretend we used all the bits
            tell = total_bits;
        }

        // Decode post-filter parameters
        let (postfilter_pitch, postfilter_gain, postfilter_tapset) =
            self.decode_postfilter(&mut dec, total_bits, tell)?;
        tell = dec.tell();

        // Decode transient flag
        let is_transient =
            if lm > 0 && tell + 3 <= total_bits { dec.dec_bit_logp(3) != 0 } else { false };
        tell = dec.tell();

        let short_blocks = if is_transient { m } else { 0 };

        // Decode intra energy flag
        let intra_ener = if tell + 3 <= total_bits { dec.dec_bit_logp(3) != 0 } else { false };

        // Decode coarse energy
        unquant_coarse_energy(
            &self.mode,
            self.start,
            self.end,
            &mut self.old_ebands,
            intra_ener,
            &mut dec,
            c,
            lm,
        );

        // Decode TF resolution
        let mut tf_res = vec![0i32; nb_ebands];
        tf_decode(self.start, self.end, is_transient, &mut tf_res, lm, &mut dec);

        tell = dec.tell();

        // Decode spread decision
        let spread_decision =
            if tell + 4 <= total_bits { dec.dec_icdf(SPREAD_ICDF, 5) } else { SPREAD_NORMAL };

        // Initialize caps
        let mut cap = vec![0i32; nb_ebands];
        init_caps(&self.mode, &mut cap, lm, c);

        // Decode dynamic allocation
        let mut offsets = vec![0i32; nb_ebands];
        let mut dynalloc_logp = 6;
        let total_bits_scaled = total_bits << BITRES;
        tell = dec.tell_frac();

        for i in self.start..self.end {
            let width = c as i32 * ((self.mode.ebands[i + 1] - self.mode.ebands[i]) as i32) << lm;
            let quanta = (width << BITRES).min((6 << BITRES).max(width));

            let mut boost = 0;
            while tell + (dynalloc_logp << BITRES) < total_bits_scaled && boost < cap[i] {
                if dec.dec_bit_logp(dynalloc_logp as u32) == 0 {
                    break;
                }
                boost += quanta;
                dynalloc_logp = 1;
                tell = dec.tell_frac();
            }

            offsets[i] = boost;
            if boost > 0 {
                dynalloc_logp = dynalloc_logp.max(2) - 1;
            }
        }

        // Decode trim
        let alloc_trim =
            if tell + (6 << BITRES) <= total_bits_scaled { dec.dec_icdf(TRIM_ICDF, 7) } else { 5 };

        // Calculate remaining bits and anti-collapse reserve
        let mut bits = (length as i32 * 8 << BITRES) - dec.tell_frac() - 1;
        let anti_collapse_rsv = if is_transient && lm >= 2 && bits >= ((lm as i32 + 2) << BITRES) {
            1 << BITRES
        }
        else {
            0
        };
        bits -= anti_collapse_rsv;

        // Decode fine energy and pulses allocation
        let mut fine_quant = vec![0i32; nb_ebands];
        let mut pulses = vec![0i32; nb_ebands];
        let mut fine_priority = vec![0i32; nb_ebands];

        let (coded_bands, mut intensity, mut dual_stereo, balance) = compute_allocation(
            &self.mode,
            self.start,
            self.end,
            &offsets,
            &cap,
            alloc_trim,
            bits,
            &mut pulses,
            &mut fine_quant,
            &mut fine_priority,
            c,
            lm,
            &mut dec,
            false,
            0,
            0,
        );

        // Decode fine energy
        unquant_fine_energy(
            &self.mode,
            self.start,
            self.end,
            &mut self.old_ebands,
            &fine_quant,
            &mut dec,
            c,
        );

        // Shift decode memory
        for ch in 0..self.channels {
            self.decode_mem[ch].copy_within(n..DECODE_BUFFER_SIZE - n + self.overlap / 2, 0);
        }

        // Allocate coefficients
        let mut x = vec![vec![0i32; n]; c];
        let mut collapse_masks = vec![0i16; c * nb_ebands];

        // Decode bands - need to split x to satisfy borrow checker
        {
            let (x0, rest) = x.split_at_mut(1);
            let y = if c == 2 { Some(rest[0].as_mut_slice()) } else { None };
            quant_all_bands(
                &self.mode,
                self.start,
                self.end,
                &mut x0[0],
                y,
                &mut collapse_masks,
                &pulses,
                short_blocks != 0,
                spread_decision,
                dual_stereo,
                intensity,
                &tf_res,
                total_bits_scaled - anti_collapse_rsv,
                balance,
                &mut dec,
                lm,
                coded_bands,
                &mut self.rng,
            );
        }

        // Decode anti-collapse flag
        let anti_collapse_on = if anti_collapse_rsv > 0 { dec.dec_bits(1) != 0 } else { false };

        // Decode final fine energy bits
        unquant_energy_finalise(
            &self.mode,
            self.start,
            self.end,
            &mut self.old_ebands,
            &fine_quant,
            &fine_priority,
            (length as i32 * 8 - dec.tell()) as i32,
            &mut dec,
            c,
        );

        // Apply anti-collapse
        if anti_collapse_on {
            anti_collapse(
                &self.mode,
                &mut x,
                &collapse_masks,
                lm,
                c,
                n,
                self.start,
                self.end,
                &self.old_ebands,
                &self.old_log_e,
                &self.old_log_e2,
                &pulses,
                &mut self.rng,
            );
        }

        // Set silence energy if needed
        if silence {
            for e in self.old_ebands.iter_mut().take(c * nb_ebands) {
                *e = minus_28db();
            }
        }

        // Synthesis: IMDCT + overlap-add
        celt_synthesis(
            &self.mode,
            &x,
            &mut out_syn,
            &self.old_ebands,
            self.start,
            eff_end,
            c,
            self.channels,
            is_transient,
            lm,
            self.downsample,
            silence,
        );

        // Capture dec state before methods that borrow self mutably
        let dec_tell = dec.tell();
        let dec_has_error = dec.has_error();

        // Apply post-filter
        self.apply_postfilter(
            &mut out_syn,
            &out_syn_offsets,
            n,
            lm,
            postfilter_pitch,
            postfilter_gain,
            postfilter_tapset,
        );

        // Update state
        self.update_state(
            c,
            is_transient,
            postfilter_pitch,
            postfilter_gain,
            postfilter_tapset,
            lm,
        );

        // De-emphasis and output
        self.apply_deemphasis(&out_syn, &out_syn_offsets, pcm, n)?;

        self.loss_count = 0;

        if dec_tell > length as i32 * 8 {
            return Err("decoder overread");
        }

        if dec_has_error {
            self.error = true;
        }

        Ok(frame_size / self.downsample)
    }

    /// Find LM (log2 of frame size multiplier) from frame size.
    fn find_lm(&self, frame_size: usize) -> Result<usize, &'static str> {
        for lm in 0..=self.mode.max_lm {
            if self.mode.short_mdct_size << lm == frame_size {
                return Ok(lm);
            }
        }
        Err("invalid frame size")
    }

    /// Decode post-filter parameters.
    fn decode_postfilter(
        &self,
        dec: &mut RangeDecoder<'_>,
        total_bits: i32,
        mut tell: i32,
    ) -> Result<(usize, i32, usize), &'static str> {
        let mut postfilter_pitch = 0usize;
        let mut postfilter_gain = 0i32;
        let mut postfilter_tapset = 0usize;

        if self.start == 0 && tell + 16 <= total_bits {
            if dec.dec_bit_logp(1) != 0 {
                let octave = dec.dec_uint(6) as usize;
                postfilter_pitch = (16 << octave) + dec.dec_bits(4 + octave as i32) as usize - 1;

                let qg = dec.dec_bits(3);
                tell = dec.tell();

                if tell + 2 <= total_bits {
                    postfilter_tapset = dec.dec_icdf(TAPSET_ICDF, 2) as usize;
                }

                // 0.09375 in Q15 ≈ 3072
                postfilter_gain = 3072 * (qg + 1);
            }
        }

        Ok((postfilter_pitch, postfilter_gain, postfilter_tapset))
    }

    /// Decode lost frame (packet loss concealment).
    fn decode_lost(&mut self, _n: usize, _lm: usize) {
        // MVP: Simple noise-based PLC
        // For production, implement proper PLC with pitch prediction

        for c in 0..self.channels {
            // Decay energy
            for e in self.old_ebands.iter_mut() {
                *e = (*e - 512).max(minus_28db()); // ~0.5 dB decay
            }

            // Generate noise
            for sample in self.decode_mem[c].iter_mut() {
                self.rng = celt_lcg_rand(self.rng);
                *sample = (self.rng as i32) >> 20;
            }
        }

        self.loss_count += 1;
    }

    /// Apply post-filter to synthesis output.
    fn apply_postfilter(
        &mut self,
        out_syn: &mut [Vec<i32>],
        out_syn_offsets: &[usize],
        n: usize,
        lm: usize,
        postfilter_pitch: usize,
        postfilter_gain: i32,
        postfilter_tapset: usize,
    ) {
        for c in 0..self.channels {
            self.postfilter_period = self.postfilter_period.max(COMBFILTER_MINPERIOD);
            self.postfilter_period_old = self.postfilter_period_old.max(COMBFILTER_MINPERIOD);

            // Clone before mutable borrow to satisfy borrow checker
            // TODO check if we can avoid clone
            let src = out_syn[c].clone();
            comb_filter(
                &mut out_syn[c],
                out_syn_offsets[c],
                &src,
                out_syn_offsets[c],
                self.postfilter_period_old,
                self.postfilter_period,
                self.mode.short_mdct_size,
                self.postfilter_gain_old,
                self.postfilter_gain,
                self.postfilter_tapset_old,
                self.postfilter_tapset,
                self.mode.window,
                self.mode.overlap,
            );

            if lm != 0 {
                let src = out_syn[c].clone();
                comb_filter(
                    &mut out_syn[c],
                    out_syn_offsets[c] + self.mode.short_mdct_size,
                    &src,
                    out_syn_offsets[c] + self.mode.short_mdct_size,
                    self.postfilter_period,
                    postfilter_pitch,
                    n - self.mode.short_mdct_size,
                    self.postfilter_gain,
                    postfilter_gain,
                    self.postfilter_tapset,
                    postfilter_tapset,
                    self.mode.window,
                    self.mode.overlap,
                );
            }
        }
    }

    /// Update decoder state for next frame.
    fn update_state(
        &mut self,
        c: usize,
        is_transient: bool,
        postfilter_pitch: usize,
        postfilter_gain: i32,
        postfilter_tapset: usize,
        lm: usize,
    ) {
        let nb_ebands = self.mode.nb_ebands;

        // Save post-filter state
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

        // Copy mono to stereo if needed
        if c == 1 {
            self.old_ebands.copy_within(0..nb_ebands, nb_ebands);
        }

        // Update log energy history
        if !is_transient {
            self.old_log_e2.copy_from_slice(&self.old_log_e);
            self.old_log_e.copy_from_slice(&self.old_ebands);

            // Update background energy
            let max_bg_increase = if self.loss_count < 10 { 32 } else { 1024 };
            for i in 0..2 * nb_ebands {
                self.background_log_e[i] =
                    (self.background_log_e[i] + max_bg_increase).min(self.old_ebands[i]);
            }
        }
        else {
            for i in 0..2 * nb_ebands {
                self.old_log_e[i] = self.old_log_e[i].min(self.old_ebands[i]);
            }
        }

        // Zero out unused bands
        for ch in 0..2 {
            for i in 0..self.start {
                self.old_ebands[ch * nb_ebands + i] = 0;
                self.old_log_e[ch * nb_ebands + i] = minus_28db();
                self.old_log_e2[ch * nb_ebands + i] = minus_28db();
            }
            for i in self.end..nb_ebands {
                self.old_ebands[ch * nb_ebands + i] = 0;
                self.old_log_e[ch * nb_ebands + i] = minus_28db();
                self.old_log_e2[ch * nb_ebands + i] = minus_28db();
            }
        }

        self.rng = self.rng; // Keep RNG state
    }

    /// Apply de-emphasis and write output.
    fn apply_deemphasis(
        &mut self,
        out_syn: &[Vec<i32>],
        out_syn_offsets: &[usize],
        pcm: &mut [f32],
        n: usize,
    ) -> Result<(), &'static str> {
        deemphasis(
            out_syn,
            out_syn_offsets,
            pcm,
            n,
            self.channels,
            self.downsample,
            &self.mode.preemph,
            &mut self.preemph_mem,
            false,
        );
        Ok(())
    }

    /// Apply post-filter and write output (used in PLC).
    fn apply_postfilter_and_output(
        &mut self,
        out_syn: &mut [Vec<i32>],
        out_syn_offsets: &[usize],
        pcm: &mut [f32],
        n: usize,
    ) -> Result<(), &'static str> {
        // Apply post-filter (with old parameters)
        for c in 0..self.channels {
            // Clone before mutable borrow to satisfy borrow checker
            // TODO check if we can avoir clone
            let src = out_syn[c].clone();
            comb_filter(
                &mut out_syn[c],
                out_syn_offsets[c],
                &src,
                out_syn_offsets[c],
                self.postfilter_period_old,
                self.postfilter_period,
                self.mode.short_mdct_size,
                self.postfilter_gain_old,
                self.postfilter_gain,
                self.postfilter_tapset_old,
                self.postfilter_tapset,
                self.mode.window,
                self.mode.overlap,
            );
        }

        self.apply_deemphasis(out_syn, out_syn_offsets, pcm, n)
    }

    /// Set the number of stream channels.
    pub fn set_channels(&mut self, channels: usize) -> Result<(), &'static str> {
        if channels < 1 || channels > 2 {
            return Err("channels must be 1 or 2");
        }
        self.stream_channels = channels;
        Ok(())
    }

    /// Set the start band.
    pub fn set_start_band(&mut self, band: usize) -> Result<(), &'static str> {
        if band >= self.mode.nb_ebands {
            return Err("start band out of range");
        }
        self.start = band;
        Ok(())
    }

    /// Set the end band.
    pub fn set_end_band(&mut self, band: usize) -> Result<(), &'static str> {
        if band < 1 || band > self.mode.nb_ebands {
            return Err("end band out of range");
        }
        self.end = band;
        Ok(())
    }

    /// Get the decoder's lookahead (overlap samples).
    pub fn get_lookahead(&self) -> usize {
        self.overlap / self.downsample
    }

    /// Get the final range value (for testing).
    pub fn get_final_range(&self) -> u32 {
        self.rng
    }
}

/// Get the resampling factor for a given sampling rate.
fn resampling_factor(rate: u32) -> Result<usize, &'static str> {
    match rate {
        48000 => Ok(1),
        24000 => Ok(2),
        16000 => Ok(3),
        12000 => Ok(4),
        8000 => Ok(6),
        _ => Err("unsupported sampling rate"),
    }
}

/// Return -28 dB in Q10 fixed-point.
#[inline]
fn minus_28db() -> i32 {
    -(((28.0 * (1 << DB_SHIFT) as f64) + 0.5) as i32)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decoder_creation() {
        let dec = CeltDecoder::new(48000, 2);
        assert!(dec.is_ok());

        let dec = CeltDecoder::new(48000, 3);
        assert!(dec.is_err());

        let dec = CeltDecoder::new(44100, 2);
        assert!(dec.is_err());
    }

    #[test]
    fn test_resampling_factor() {
        assert_eq!(resampling_factor(48000), Ok(1));
        assert_eq!(resampling_factor(24000), Ok(2));
        assert_eq!(resampling_factor(16000), Ok(3));
        assert_eq!(resampling_factor(12000), Ok(4));
        assert_eq!(resampling_factor(8000), Ok(6));
        assert!(resampling_factor(44100).is_err());
    }

    #[test]
    fn test_decoder_reset() {
        let mut dec = CeltDecoder::new(48000, 2).unwrap();
        dec.rng = 12345;
        dec.loss_count = 5;

        dec.reset();

        assert_eq!(dec.rng, 0);
        assert_eq!(dec.loss_count, 0);
    }
}
