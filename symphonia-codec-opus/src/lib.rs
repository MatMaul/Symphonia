// Symphonia
// Copyright (c) 2026 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#![cfg_attr(not(feature = "std"), no_std)]
#![warn(rust_2018_idioms)]
#![forbid(unsafe_code)]
// The following lints are allowed in all Symphonia crates. Please see clippy.toml for their
// justification.
#![allow(clippy::comparison_chain)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::identity_op)]
#![allow(clippy::manual_range_contains)]

extern crate alloc;

mod celt;

use symphonia_core::audio::{AsGenericAudioBufferRef, AudioBuffer, AudioMut, AudioSpec, Channels};
use symphonia_core::codecs::CodecInfo;
use symphonia_core::codecs::audio::{AudioCodecParameters, AudioDecoder, AudioDecoderOptions, FinalizeResult};
use symphonia_core::codecs::audio::well_known::CODEC_ID_OPUS;
use symphonia_core::codecs::registry::{RegisterableAudioDecoder, SupportedAudioCodec};
use symphonia_core::errors::{Error, Result, decode_error, unsupported_error};
use symphonia_core::formats::Packet;
use symphonia_core::io::{BufReader, ReadBytes};
use symphonia_core::support_audio_codec;

use celt::{CeltDecodeError, CeltDecoder};

/// Opus decoder (CELT-only for now).
pub struct OpusDecoder {
    params: AudioCodecParameters,
    celt: CeltDecoder,
    buf: AudioBuffer<f32>,
    channels: usize,
}

impl OpusDecoder {
    pub fn try_new(params: &AudioCodecParameters, _opts: &AudioDecoderOptions) -> Result<Self> {
        if params.codec != CODEC_ID_OPUS {
            return unsupported_error("opus: invalid codec");
        }
        let extra_data = match params.extra_data.as_ref() {
            Some(buf) => buf,
            None => return unsupported_error("opus: missing extra data"),
        };

        let (channels, _pre_skip) = parse_opus_head(extra_data)?;
        if channels > 2 {
            return unsupported_error("opus: multichannel not supported");
        }
        let mode = match celt::mode_from_static(48_000, 960) {
            Some(mode) => mode,
            None => return unsupported_error("opus: unsupported static mode"),
        };
        let celt = CeltDecoder::new(mode, channels as i32)
            .map_err(|_| Error::DecodeError("opus: decoder init failed"))?;

        let spec = match params.channels.as_ref() {
            Some(ch) => AudioSpec::new(48_000, ch.clone()),
            None => AudioSpec::new(48_000, Channels::Discrete(channels as u16)),
        };

        Ok(Self {
            params: params.clone(),
            celt,
            buf: AudioBuffer::new(spec, 0),
            channels: channels as usize,
        })
    }

    fn decode_inner(&mut self, packet: &Packet) -> Result<()> {
        let parsed = parse_opus_packet(packet.buf())?;
        if parsed.config < 16 {
            return unsupported_error("opus: silk/hybrid not supported");
        }

        let total_frames = parsed.frames.len();
        if total_frames == 0 {
            return decode_error("opus: empty packet");
        }

        let frame_size = parsed.frame_size as usize;
        let total_samples = frame_size * total_frames;

        self.buf.clear();
        self.buf.grow_capacity(total_samples);
        self.buf.render_silence(Some(total_samples));

        for (idx, frame) in parsed.frames.iter().enumerate() {
            let offset = idx * frame_size;
            let mut out_refs: Vec<&mut [f32]> = self
                .buf
                .iter_planes_mut()
                .take(self.channels)
                .map(|plane| &mut plane[offset..offset + frame_size])
                .collect();
            match self.celt.decode_frame(frame, &mut out_refs, parsed.frame_size as i32) {
                Ok(_) => {}
                Err(CeltDecodeError::BufferTooSmall) => {
                    return decode_error("opus: output buffer too small");
                }
                Err(_) => return decode_error("opus: CELT decode failed"),
            }
        }

        Ok(())
    }
}

impl AudioDecoder for OpusDecoder {
    fn reset(&mut self) {
        self.celt.reset();
    }

    fn codec_info(&self) -> &CodecInfo {
        &Self::supported_codecs().first().unwrap().info
    }

    fn codec_params(&self) -> &AudioCodecParameters {
        &self.params
    }

    fn decode(&mut self, packet: &Packet) -> Result<symphonia_core::audio::GenericAudioBufferRef<'_>> {
        if let Err(err) = self.decode_inner(packet) {
            self.buf.clear();
            Err(err)
        } else {
            Ok(self.buf.as_generic_audio_buffer_ref())
        }
    }

    fn finalize(&mut self) -> FinalizeResult {
        Default::default()
    }

    fn last_decoded(&self) -> symphonia_core::audio::GenericAudioBufferRef<'_> {
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
        &[support_audio_codec!(CODEC_ID_OPUS, "opus", "Opus (CELT-only)")]
    }
}

struct ParsedOpusPacket<'a> {
    toc: u8,
    config: u8,
    frame_size: usize,
    frames: Vec<&'a [u8]>,
}

fn parse_opus_head(buf: &[u8]) -> Result<(u8, u16)> {
    const OPUS_HEAD_MAGIC: &[u8; 8] = b"OpusHead";
    if buf.len() < 19 {
        return decode_error("opus: identification header too short");
    }
    let mut reader = BufReader::new(buf);
    let mut magic = [0u8; 8];
    reader.read_buf_exact(&mut magic)?;
    if magic != *OPUS_HEAD_MAGIC {
        return decode_error("opus: invalid identification header");
    }
    let _version = reader.read_u8()?;
    let channels = reader.read_u8()?;
    if channels == 0 {
        return decode_error("opus: invalid channel count");
    }
    let pre_skip = reader.read_u16()?;
    Ok((channels, pre_skip))
}

fn opus_packet_get_samples_per_frame(toc: u8, sample_rate: u32) -> usize {
    if (toc & 0x80) != 0 {
        let audiosize = ((toc >> 3) & 0x3) as u32;
        (sample_rate << audiosize) as usize / 400
    } else if (toc & 0x60) == 0x60 {
        if (toc & 0x08) != 0 {
            (sample_rate / 50) as usize
        } else {
            (sample_rate / 100) as usize
        }
    } else {
        let audiosize = ((toc >> 3) & 0x3) as u32;
        if audiosize == 3 {
            (sample_rate * 60 / 1000) as usize
        } else {
            (sample_rate << audiosize) as usize / 100
        }
    }
}

fn parse_size(data: &[u8]) -> Result<(usize, usize)> {
    if data.is_empty() {
        return decode_error("opus: truncated size field");
    }
    if data[0] < 252 {
        Ok((data[0] as usize, 1))
    } else if data.len() < 2 {
        decode_error("opus: truncated size field")
    } else {
        Ok(((4 * data[1] as usize) + data[0] as usize, 2))
    }
}

fn parse_opus_packet(buf: &[u8]) -> Result<ParsedOpusPacket<'_>> {
    if buf.is_empty() {
        return decode_error("opus: empty packet");
    }
    let toc = buf[0];
    let config = toc >> 3;
    let frame_size = opus_packet_get_samples_per_frame(toc, 48_000);

    let mut idx = 1usize;
    let mut remaining = buf.len() - 1;
    let mut sizes: Vec<usize> = Vec::new();

    match toc & 0x3 {
        0 => {
            sizes.push(remaining);
        }
        1 => {
            if remaining & 0x1 != 0 {
                return decode_error("opus: invalid CBR frame sizes");
            }
            let size = remaining / 2;
            sizes.push(size);
            sizes.push(size);
        }
        2 => {
            let (size0, bytes) = parse_size(&buf[idx..])?;
            idx += bytes;
            if remaining < bytes {
                return decode_error("opus: invalid size field");
            }
            remaining -= bytes;
            if size0 > remaining {
                return decode_error("opus: invalid VBR frame size");
            }
            sizes.push(size0);
            sizes.push(remaining - size0);
        }
        _ => {
            if remaining < 1 {
                return decode_error("opus: invalid packet header");
            }
            let ch = buf[idx];
            idx += 1;
            remaining -= 1;
            let count = (ch & 0x3f) as usize;
            if count == 0 || frame_size * count > 5760 {
                return decode_error("opus: invalid frame count");
            }
            if (ch & 0x40) != 0 {
                loop {
                    if remaining == 0 {
                        return decode_error("opus: truncated padding");
                    }
                    let p = buf[idx];
                    idx += 1;
                    remaining -= 1;
                    let tmp = if p == 255 { 254 } else { p } as usize;
                    if remaining < tmp {
                        return decode_error("opus: invalid padding length");
                    }
                    remaining -= tmp;
                    if p != 255 {
                        break;
                    }
                }
            }

            let vbr = (ch & 0x80) != 0;
            if vbr {
                let mut remaining_payload = remaining;
                for _ in 0..(count - 1) {
                    let (size, bytes) = parse_size(&buf[idx..])?;
                    idx += bytes;
                    if remaining_payload < bytes {
                        return decode_error("opus: invalid size field");
                    }
                    remaining_payload -= bytes;
                    if size > remaining_payload {
                        return decode_error("opus: invalid VBR frame size");
                    }
                    sizes.push(size);
                    remaining_payload -= size;
                }
                if remaining_payload > 1275 {
                    return decode_error("opus: invalid packet size");
                }
                sizes.push(remaining_payload);
            } else {
                if remaining % count != 0 {
                    return decode_error("opus: invalid CBR packet size");
                }
                let size = remaining / count;
                if size > 1275 {
                    return decode_error("opus: invalid packet size");
                }
                sizes.resize(count, size);
            }
        }
    }

    let mut frames = Vec::with_capacity(sizes.len());
    if sizes.iter().any(|&size| size > 1275) {
        return decode_error("opus: invalid frame size");
    }
    for size in sizes {
        let end = idx + size;
        if end > buf.len() {
            return decode_error("opus: truncated packet");
        }
        frames.push(&buf[idx..end]);
        idx = end;
    }

    Ok(ParsedOpusPacket { toc, config, frame_size, frames })
}
