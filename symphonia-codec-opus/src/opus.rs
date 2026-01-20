// Symphonia
// Copyright (c) 2024 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! Opus decoder implementation.
//!
//! This module implements the main OpusDecoder that wraps the CELT decoder
//! and implements the symphonia AudioDecoder trait.

#[cfg(not(feature = "std"))]
use alloc::{boxed::Box, vec, vec::Vec};
#[cfg(feature = "std")]
use std::{vec, vec::Vec};

use symphonia_core::audio::{
    AsGenericAudioBufferRef, AudioBuffer, AudioMut, AudioSpec, GenericAudioBufferRef,
};
use symphonia_core::codecs::CodecInfo;
use symphonia_core::codecs::audio::well_known::CODEC_ID_OPUS;
use symphonia_core::codecs::audio::{AudioCodecParameters, AudioDecoderOptions};
use symphonia_core::codecs::audio::{AudioDecoder, FinalizeResult};
use symphonia_core::codecs::registry::{RegisterableAudioDecoder, SupportedAudioCodec};
use symphonia_core::errors::{decode_error, unsupported_error, Result};
use symphonia_core::formats::Packet;
use symphonia_core::support_audio_codec;

use crate::celt::decoder::CeltDecoder;

/// Maximum frame duration in ms.
const MAX_FRAME_DURATION_MS: usize = 120;

/// Maximum sample rate.
const MAX_SAMPLE_RATE: u32 = 48000;

/// Maximum samples per frame (120ms @ 48kHz).
const MAX_FRAME_SIZE: usize = MAX_FRAME_DURATION_MS * MAX_SAMPLE_RATE as usize / 1000;

/// Opus decoder.
pub struct OpusDecoder {
    /// Decoder parameters.
    params: AudioCodecParameters,

    /// CELT decoder (for CELT-only and hybrid modes).
    celt_decoder: CeltDecoder,

    /// Output sample rate.
    sample_rate: u32,

    /// Number of channels.
    channels: usize,

    /// Output audio buffer.
    buf: AudioBuffer<f32>,

    /// Temporary decode buffer.
    decode_buf: Vec<f32>,
}

impl OpusDecoder {
    /// Create a new Opus decoder.
    pub fn try_new(params: &AudioCodecParameters, _opts: &AudioDecoderOptions) -> Result<Self> {
        // This decoder only supports Opus.
        if params.codec != CODEC_ID_OPUS {
            return unsupported_error("opus: invalid codec");
        }

        // Get channel count from parameters.
        let channels = params.channels.as_ref().map(|c| c.count()).unwrap_or(2);

        if channels < 1 || channels > 2 {
            return unsupported_error("opus: only mono and stereo are supported");
        }

        // Get sample rate (default to 48kHz).
        let sample_rate = params.sample_rate.unwrap_or(48000);

        // Validate sample rate.
        if !matches!(sample_rate, 8000 | 12000 | 16000 | 24000 | 48000) {
            return unsupported_error("opus: unsupported sample rate");
        }

        // Create CELT decoder.
        let celt_decoder = CeltDecoder::new(sample_rate, channels)
            .map_err(|e| symphonia_core::errors::Error::DecodeError(e))?;

        // Create output buffer.
        let spec = AudioSpec::new(sample_rate, params.channels.clone().unwrap_or_default());
        let buf = AudioBuffer::new(spec, MAX_FRAME_SIZE);

        Ok(Self {
            params: params.clone(),
            celt_decoder,
            sample_rate,
            channels,
            buf,
            decode_buf: vec![0.0; MAX_FRAME_SIZE * channels],
        })
    }

    /// Parse the TOC byte to extract mode, bandwidth, and frame count info.
    fn parse_toc(toc: u8) -> (OpusMode, Bandwidth, FrameConfig) {
        let config = (toc >> 3) & 0x1f;
        let stereo = (toc >> 2) & 0x01 != 0;
        let frame_code = toc & 0x03;

        let (mode, bandwidth) = match config {
            0..=3 => (OpusMode::Silk, Bandwidth::Narrowband),
            4..=7 => (OpusMode::Silk, Bandwidth::Mediumband),
            8..=11 => (OpusMode::Silk, Bandwidth::Wideband),
            12..=13 => (OpusMode::Hybrid, Bandwidth::Superwideband),
            14..=15 => (OpusMode::Hybrid, Bandwidth::Fullband),
            16..=19 => (OpusMode::Celt, Bandwidth::Narrowband),
            20..=23 => (OpusMode::Celt, Bandwidth::Wideband),
            24..=27 => (OpusMode::Celt, Bandwidth::Superwideband),
            28..=31 => (OpusMode::Celt, Bandwidth::Fullband),
            _ => unreachable!(),
        };

        // Frame duration depends on mode (RFC 6716 Section 3.1):
        // - SILK/Hybrid modes (configs 0-15): 10ms, 20ms, 40ms, 60ms
        // - CELT-only modes (configs 16-31): 2.5ms, 5ms, 10ms, 20ms
        let frame_duration = match config {
            // SILK and Hybrid modes
            0 | 4 | 8 | 12 | 14 => FrameDuration::Ms10,
            1 | 5 | 9 | 13 | 15 => FrameDuration::Ms20,
            2 | 6 | 10 => FrameDuration::Ms40,
            3 | 7 | 11 => FrameDuration::Ms60,
            // CELT-only modes
            16 | 20 | 24 | 28 => FrameDuration::Ms2_5,
            17 | 21 | 25 | 29 => FrameDuration::Ms5,
            18 | 22 | 26 | 30 => FrameDuration::Ms10,
            19 | 23 | 27 | 31 => FrameDuration::Ms20,
            _ => unreachable!(),
        };

        let frame_config = FrameConfig {
            duration: frame_duration,
            stereo,
            code: frame_code,
        };

        (mode, bandwidth, frame_config)
    }

    /// Get frame size in samples for a given duration at 48kHz.
    fn frame_size_48k(duration: FrameDuration) -> usize {
        match duration {
            FrameDuration::Ms2_5 => 120,
            FrameDuration::Ms5 => 240,
            FrameDuration::Ms10 => 480,
            FrameDuration::Ms20 => 960,
            FrameDuration::Ms40 => 1920,
            FrameDuration::Ms60 => 2880,
        }
    }

    /// Get frame size in samples for the configured sample rate.
    fn frame_size(&self, duration: FrameDuration) -> usize {
        let size_48k = Self::frame_size_48k(duration);
        size_48k * self.sample_rate as usize / 48000
    }

    fn decode_inner(&mut self, packet: &Packet) -> Result<()> {
        let data = packet.buf();

        if data.is_empty() {
            return decode_error("opus: empty packet");
        }

        // Parse TOC byte.
        let toc = data[0];
        let (mode, _bandwidth, frame_config) = Self::parse_toc(toc);

        // For MVP, only support CELT mode.
        if mode != OpusMode::Celt {
            return unsupported_error("opus: only CELT mode is currently supported");
        }

        // Determine number of frames and their sizes.
        let (frame_count, frame_sizes) = self.parse_frame_lengths(&data[1..], frame_config.code)?;

        // Calculate total frame size in samples.
        let frame_samples = self.frame_size(frame_config.duration);
        let total_samples = frame_samples * frame_count;

        if total_samples > MAX_FRAME_SIZE {
            return decode_error("opus: frame too large");
        }

        // Clear output buffer and render space.
        self.buf.clear();
        self.buf.render_uninit(Some(total_samples));

        // Decode each frame.
        let mut data_offset = 1 + self.toc_extra_bytes(frame_config.code, frame_count);
        let mut sample_offset = 0;

        for i in 0..frame_count {
            let frame_len = frame_sizes[i];
            let frame_data = if frame_len > 0 {
                Some(&data[data_offset..data_offset + frame_len])
            } else {
                None
            };

            // Decode CELT frame.
            let decode_buf = &mut self.decode_buf[..frame_samples * self.channels];
            let decoded = self.celt_decoder
                .decode(frame_data, decode_buf, frame_samples)
                .map_err(|e| symphonia_core::errors::Error::DecodeError(e))?;

            // Copy to output buffer (de-interleave).
            for ch in 0..self.channels {
                let plane = self.buf.plane_mut(ch).unwrap();
                for s in 0..decoded {
                    plane[sample_offset + s] = decode_buf[s * self.channels + ch];
                }
            }

            sample_offset += decoded;
            data_offset += frame_len;
        }

        // Trim buffer to actual samples decoded.
        self.buf.trim(0, total_samples.saturating_sub(sample_offset));

        Ok(())
    }

    /// Parse frame lengths from packet data.
    fn parse_frame_lengths(
        &self,
        data: &[u8],
        code: u8,
    ) -> Result<(usize, Vec<usize>)> {
        match code {
            0 => {
                // One frame.
                Ok((1, vec![data.len()]))
            }
            1 => {
                // Two equal-sized frames.
                if data.len() % 2 != 0 {
                    return decode_error("opus: invalid packet for code 1");
                }
                let frame_len = data.len() / 2;
                Ok((2, vec![frame_len, frame_len]))
            }
            2 => {
                // Two frames with different sizes.
                if data.is_empty() {
                    return decode_error("opus: invalid packet for code 2");
                }

                let (len1, bytes_read) = self.parse_frame_length(data)?;
                let len2 = data.len() - bytes_read - len1;

                Ok((2, vec![len1, len2]))
            }
            3 => {
                // Multiple frames (VBR or CBR).
                if data.is_empty() {
                    return decode_error("opus: invalid packet for code 3");
                }

                let count_byte = data[0];
                let vbr = (count_byte & 0x80) != 0;
                let padding = (count_byte & 0x40) != 0;
                let frame_count = (count_byte & 0x3f) as usize;

                if frame_count == 0 {
                    return decode_error("opus: zero frames in code 3 packet");
                }

                let mut offset = 1;

                // Skip padding bytes.
                if padding {
                    let mut pad_len = 0usize;
                    loop {
                        if offset >= data.len() {
                            return decode_error("opus: padding overflow");
                        }
                        let p = data[offset] as usize;
                        offset += 1;
                        pad_len += p;
                        if p < 255 {
                            break;
                        }
                    }
                    // Padding is at the end, so we just note its length.
                    let remaining = data.len() - offset;
                    if pad_len > remaining {
                        return decode_error("opus: padding exceeds packet");
                    }
                }

                let data_remaining = &data[offset..];

                if vbr {
                    // Variable bitrate: each frame has its own length.
                    let mut sizes = Vec::with_capacity(frame_count);
                    let mut len_offset = 0;
                    let mut total_len = 0;

                    for _ in 0..frame_count - 1 {
                        let (len, bytes) = self.parse_frame_length(&data_remaining[len_offset..])?;
                        sizes.push(len);
                        total_len += len;
                        len_offset += bytes;
                    }

                    // Last frame gets remaining bytes.
                    let last_len = data_remaining.len() - len_offset - total_len;
                    sizes.push(last_len);

                    Ok((frame_count, sizes))
                } else {
                    // Constant bitrate.
                    if data_remaining.len() % frame_count != 0 {
                        return decode_error("opus: CBR packet size not divisible");
                    }
                    let frame_len = data_remaining.len() / frame_count;
                    Ok((frame_count, vec![frame_len; frame_count]))
                }
            }
            _ => decode_error("opus: invalid frame code"),
        }
    }

    /// Parse a variable-length frame size.
    fn parse_frame_length(&self, data: &[u8]) -> Result<(usize, usize)> {
        if data.is_empty() {
            return decode_error("opus: unexpected end of packet");
        }

        let first = data[0] as usize;
        if first < 252 {
            Ok((first, 1))
        } else {
            if data.len() < 2 {
                return decode_error("opus: unexpected end of packet");
            }
            let second = data[1] as usize;
            Ok((second * 4 + first, 2))
        }
    }

    /// Calculate extra bytes after TOC for frame length encoding.
    fn toc_extra_bytes(&self, code: u8, frame_count: usize) -> usize {
        match code {
            0 | 1 => 0,
            2 => 1, // At least 1 byte for first frame length.
            3 => 1 + frame_count - 1, // Count byte + (n-1) length bytes (minimum).
            _ => 0,
        }
    }
}

impl AudioDecoder for OpusDecoder {
    fn reset(&mut self) {
        self.celt_decoder.reset();
        self.buf.clear();
    }

    fn codec_info(&self) -> &CodecInfo {
        // Only one codec is supported.
        &Self::supported_codecs().first().unwrap().info
    }

    fn codec_params(&self) -> &AudioCodecParameters {
        &self.params
    }

    fn decode(&mut self, packet: &Packet) -> Result<GenericAudioBufferRef<'_>> {
        if let Err(e) = self.decode_inner(packet) {
            self.buf.clear();
            Err(e)
        } else {
            Ok(self.buf.as_generic_audio_buffer_ref())
        }
    }

    fn finalize(&mut self) -> FinalizeResult {
        Default::default()
    }

    fn last_decoded(&self) -> GenericAudioBufferRef<'_> {
        self.buf.as_generic_audio_buffer_ref()
    }
}

impl RegisterableAudioDecoder for OpusDecoder {
    fn try_registry_new(
        params: &AudioCodecParameters,
        opts: &AudioDecoderOptions,
    ) -> Result<Box<dyn AudioDecoder>>
    where
        Self: Sized,
    {
        Ok(Box::new(OpusDecoder::try_new(params, opts)?))
    }

    fn supported_codecs() -> &'static [SupportedAudioCodec] {
        &[support_audio_codec!(CODEC_ID_OPUS, "opus", "Opus")]
    }
}

/// Opus coding mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OpusMode {
    /// SILK mode (speech).
    Silk,
    /// Hybrid mode (SILK + CELT).
    Hybrid,
    /// CELT mode (audio).
    Celt,
}

/// Opus bandwidth.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
enum Bandwidth {
    Narrowband,
    Mediumband,
    Wideband,
    Superwideband,
    Fullband,
}

/// Frame duration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FrameDuration {
    Ms2_5,
    Ms5,
    Ms10,
    Ms20,
    Ms40,
    Ms60,
}

/// Frame configuration from TOC byte.
#[derive(Debug, Clone, Copy)]
struct FrameConfig {
    duration: FrameDuration,
    stereo: bool,
    code: u8,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_toc() {
        // CELT fullband, mono, 20ms, code 0
        let (mode, bw, config) = OpusDecoder::parse_toc(0b11101_0_00);
        assert_eq!(mode, OpusMode::Celt);
        assert_eq!(bw, Bandwidth::Fullband);
        assert!(!config.stereo);
        assert_eq!(config.code, 0);
    }

    #[test]
    fn test_frame_size() {
        assert_eq!(OpusDecoder::frame_size_48k(FrameDuration::Ms10), 480);
        assert_eq!(OpusDecoder::frame_size_48k(FrameDuration::Ms20), 960);
    }
}
