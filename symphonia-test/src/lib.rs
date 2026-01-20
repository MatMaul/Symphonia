// Symphonia Test Framework
// Copyright (c) 2019-2022 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! A generic test framework for comparing Symphonia's decoding output against FFmpeg.
//!
//! This crate provides utilities for verifying that Symphonia produces correct decoded audio
//! by comparing its output sample-by-sample against a reference decoder (FFmpeg by default).
//!
//! # Example
//!
//! ```no_run
//! use symphonia_test::{TestConfig, test_decode_files};
//!
//! let config = TestConfig::default();
//! let files = vec!["test1.flac", "test2.flac"];
//! let results = test_decode_files(&files, &config);
//!
//! for result in results {
//!     assert!(result.is_ok(), "Failed: {:?}", result);
//! }
//! ```

use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::process::{Child, Command, Stdio};

use symphonia::core::audio::GenericAudioBufferRef;
use symphonia::core::codecs::CodecParameters;
use symphonia::core::codecs::audio::{AudioDecoder, AudioDecoderOptions};
use symphonia::core::errors::Error as SymphoniaError;
use symphonia::core::formats::probe::Hint;
use symphonia::core::formats::{FormatOptions, FormatReader, TrackType};
use symphonia::core::io::{MediaSourceStream, ReadOnlySource};
use symphonia::core::meta::MetadataOptions;

/// Error type for the test framework.
#[derive(Debug, thiserror::Error)]
pub enum TestError {
    #[error("Symphonia error: {0}")]
    Symphonia(#[from] SymphoniaError),

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error(
        "Sample mismatch at packet {packet}, sample {sample}: symphonia={symphonia}, ffmpeg={ffmpeg}, delta={delta}"
    )]
    SampleMismatch { packet: u64, sample: u64, symphonia: f32, ffmpeg: f32, delta: f32 },

    #[error("Decoder configuration mismatch: {0}")]
    ConfigMismatch(String),

    #[error("Reference decoder (FFmpeg) failed to start: {0}")]
    FfmpegSpawnFailed(String),

    #[error("No audio track found in file")]
    NoAudioTrack,

    #[error(
        "Remaining samples mismatch: symphonia has {symphonia} extra, ffmpeg has {ffmpeg} extra"
    )]
    RemainingSamplesMismatch { symphonia: u64, ffmpeg: u64 },
}

pub type TestResult<T> = Result<T, TestError>;

/// Reference decoder to use for comparison.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum RefDecoder {
    #[default]
    Ffmpeg,
}

/// Configuration for the test framework.
#[derive(Debug, Clone)]
pub struct TestConfig {
    /// Reference decoder to use.
    pub ref_decoder: RefDecoder,

    /// Maximum allowable sample delta (absolute value).
    /// Default is ~2^-17 (-102.4dB).
    pub max_sample_delta: f32,

    /// Whether to enable gapless decoding.
    pub gapless: bool,

    /// Whether to continue after decode errors.
    pub keep_going: bool,

    /// Stop after the first failed sample.
    pub stop_on_first_failure: bool,

    /// Maximum number of sample mismatches to report before stopping.
    pub max_failures_to_report: usize,
}

impl Default for TestConfig {
    fn default() -> Self {
        Self {
            ref_decoder: RefDecoder::default(),
            max_sample_delta: 0.00001, // ~2^-17 (-102.4dB)
            gapless: true,
            keep_going: false,
            stop_on_first_failure: false,
            max_failures_to_report: 10,
        }
    }
}

/// Statistics from a comparison test.
#[derive(Debug, Default, Clone)]
pub struct TestStats {
    /// Total number of packets decoded.
    pub packets: u64,
    /// Total number of samples compared.
    pub samples: u64,
    /// Number of failed samples.
    pub failed_samples: u64,
    /// Number of packets with at least one failed sample.
    pub failed_packets: u64,
    /// Maximum absolute sample delta observed.
    pub max_delta: f32,
    /// Number of unchecked samples remaining in Symphonia decoder.
    pub symphonia_remaining: u64,
    /// Number of unchecked samples remaining in FFmpeg decoder.
    pub ffmpeg_remaining: u64,
}

impl TestStats {
    /// Returns true if the test passed (no failed samples and no remaining samples mismatch).
    pub fn passed(&self) -> bool {
        self.failed_samples == 0 && self.symphonia_remaining == 0 && self.ffmpeg_remaining == 0
    }
}

/// Wrapper around a reference decoder process.
struct RefProcess {
    child: Child,
}

impl RefProcess {
    fn spawn_ffmpeg(path: &str, gapless: bool) -> TestResult<Self> {
        let mut cmd = Command::new("ffmpeg");

        // Gapless argument must come before everything else.
        if !gapless {
            cmd.arg("-flags2").arg("skip_manual");
        }

        cmd.arg("-nostats")
            .arg("-hide_banner")
            .arg("-i")
            .arg(path)
            .arg("-map")
            .arg("0:a:0")
            .arg("-c:a")
            .arg("pcm_f32le")
            .arg("-f")
            .arg("wav")
            .arg("-")
            .stdout(Stdio::piped())
            .stderr(Stdio::null());

        let child = cmd.spawn().map_err(|e| TestError::FfmpegSpawnFailed(e.to_string()))?;

        Ok(Self { child })
    }

    fn spawn(decoder: RefDecoder, path: &str, gapless: bool) -> TestResult<Self> {
        match decoder {
            RefDecoder::Ffmpeg => Self::spawn_ffmpeg(path, gapless),
        }
    }
}

/// Internal decoder instance wrapper.
struct DecoderInstance {
    format: Box<dyn FormatReader>,
    decoder: Box<dyn AudioDecoder>,
    track_id: u32,
}

impl DecoderInstance {
    fn open(mss: MediaSourceStream<'static>, fmt_opts: FormatOptions) -> TestResult<Self> {
        let meta_opts = MetadataOptions::default();
        let dec_opts = AudioDecoderOptions::default();
        let hint = Hint::new();

        let format = symphonia::default::get_probe().probe(&hint, mss, fmt_opts, meta_opts)?;

        let track = format.default_track(TrackType::Audio).ok_or(TestError::NoAudioTrack)?;

        let codec_params = match &track.codec_params {
            Some(CodecParameters::Audio(params)) => params,
            _ => return Err(TestError::NoAudioTrack),
        };

        let decoder =
            symphonia::default::get_codecs().make_audio_decoder(codec_params, &dec_opts)?;

        let track_id = track.id;

        Ok(Self { format, decoder, track_id })
    }

    fn samples_per_frame(&self) -> Option<u64> {
        self.decoder.codec_params().channels.as_ref().map(|ch| ch.count() as u64)
    }

    fn next_audio_buf(
        &mut self,
        keep_going: bool,
    ) -> TestResult<Option<GenericAudioBufferRef<'_>>> {
        loop {
            let packet = match self.format.next_packet() {
                Ok(Some(packet)) => packet,
                Ok(None) => return Ok(None),
                Err(SymphoniaError::IoError(err))
                    if err.kind() == std::io::ErrorKind::UnexpectedEof =>
                {
                    // WavReader will always return an UnexpectedEof when it ends because the
                    // reference decoder is piping the decoded audio and cannot write out the
                    // actual length of the media.
                    return Ok(None);
                }
                Err(err) => return Err(err.into()),
            };

            if packet.track_id() != self.track_id {
                continue;
            }

            match self.decoder.decode(&packet) {
                Ok(_) => break,
                Err(SymphoniaError::DecodeError(_)) if keep_going => continue,
                Err(err) => return Err(err.into()),
            }
        }

        Ok(Some(self.decoder.last_decoded()))
    }

    fn flush(&mut self, keep_going: bool) -> TestResult<u64> {
        let mut samples = 0u64;

        while let Some(buf) = self.next_audio_buf(keep_going)? {
            samples += buf.samples_interleaved() as u64;
        }

        Ok(samples)
    }
}

/// Compare Symphonia's decoding output against FFmpeg for a single file.
///
/// Returns `TestStats` containing detailed comparison statistics.
pub fn compare_decode<P: AsRef<Path>>(path: P, config: &TestConfig) -> TestResult<TestStats> {
    let path = path.as_ref();
    let path_str = path.to_string_lossy();

    // Spawn reference decoder
    let mut ref_process = RefProcess::spawn(config.ref_decoder, &path_str, config.gapless)?;

    // Open reference decoder output as a Symphonia stream (WAV format)
    let ref_ms =
        Box::new(ReadOnlySource::new(BufReader::new(ref_process.child.stdout.take().unwrap())));
    let ref_mss = MediaSourceStream::new(ref_ms, Default::default());
    let mut ref_inst = DecoderInstance::open(ref_mss, FormatOptions::default())?;

    // Open target file with Symphonia
    let tgt_ms = Box::new(File::open(path)?);
    let tgt_mss = MediaSourceStream::new(tgt_ms, Default::default());
    let tgt_fmt_opts = FormatOptions { enable_gapless: config.gapless, ..Default::default() };
    let mut tgt_inst = DecoderInstance::open(tgt_mss, tgt_fmt_opts)?;

    // Verify samples per frame match
    let samples_per_frame = tgt_inst.samples_per_frame().unwrap_or(1);
    if samples_per_frame != ref_inst.samples_per_frame().unwrap_or(1) {
        return Err(TestError::ConfigMismatch(
            "samples per frame mismatch between target and reference".to_string(),
        ));
    }

    let mut stats = TestStats::default();

    // Buffers for samples
    let mut ref_sample_buf: Vec<f32> = Vec::new();
    let mut ref_sample_pos = 0;
    let mut ref_sample_cnt = 0;

    let mut tgt_sample_buf: Vec<f32> = Vec::new();
    let mut tgt_sample_pos = 0;
    let mut tgt_sample_cnt = 0;

    let mut early_exit = false;

    'outer: loop {
        // Decode next target buffer
        match tgt_inst.next_audio_buf(config.keep_going)? {
            Some(buf) => buf.copy_to_vec_interleaved(&mut tgt_sample_buf),
            None => break,
        }

        tgt_sample_cnt = tgt_sample_buf.len();
        tgt_sample_pos = 0;

        let mut pkt_failed_samples = 0u64;

        while tgt_sample_pos < tgt_sample_cnt {
            // Need more reference samples
            if ref_sample_pos == ref_sample_cnt {
                match ref_inst.next_audio_buf(true)? {
                    Some(buf) => buf.copy_to_vec_interleaved(&mut ref_sample_buf),
                    None => break 'outer,
                }
                ref_sample_cnt = ref_sample_buf.len();
                ref_sample_pos = 0;
            }

            let ref_samples = &ref_sample_buf[ref_sample_pos..];
            let tgt_samples = &tgt_sample_buf[tgt_sample_pos..];

            let n_test = std::cmp::min(ref_samples.len(), tgt_samples.len());

            for (&t, &r) in tgt_samples[..n_test].iter().zip(&ref_samples[..n_test]) {
                let delta = t.clamp(-1.0, 1.0) - r.clamp(-1.0, 1.0);

                if delta.abs() > config.max_sample_delta {
                    pkt_failed_samples += 1;

                    if config.stop_on_first_failure {
                        return Err(TestError::SampleMismatch {
                            packet: stats.packets,
                            sample: stats.samples,
                            symphonia: t,
                            ffmpeg: r,
                            delta,
                        });
                    }
                }

                stats.max_delta = stats.max_delta.max(delta.abs());
                stats.samples += 1;
            }

            ref_sample_pos += n_test;
            tgt_sample_pos += n_test;
        }

        stats.failed_samples += pkt_failed_samples;
        stats.failed_packets += u64::from(pkt_failed_samples > 0);
        stats.packets += 1;

        if stats.failed_samples > 0 && config.stop_on_first_failure {
            early_exit = true;
            break;
        }
    }

    // Count remaining samples if we didn't exit early
    if !early_exit {
        let tgt_remaining = (tgt_sample_cnt - tgt_sample_pos) as u64 + tgt_inst.flush(true)?;
        let ref_remaining = (ref_sample_cnt - ref_sample_pos) as u64 + ref_inst.flush(true)?;

        stats.symphonia_remaining = tgt_remaining;
        stats.ffmpeg_remaining = ref_remaining;
    }

    Ok(stats)
}

/// Test a single file and return Ok if it passes, Err otherwise.
pub fn test_decode<P: AsRef<Path>>(path: P, config: &TestConfig) -> TestResult<TestStats> {
    let stats = compare_decode(&path, config)?;

    if stats.failed_samples > 0 {
        return Err(TestError::SampleMismatch {
            packet: 0,
            sample: 0,
            symphonia: 0.0,
            ffmpeg: 0.0,
            delta: stats.max_delta,
        });
    }

    if stats.symphonia_remaining != 0 || stats.ffmpeg_remaining != 0 {
        return Err(TestError::RemainingSamplesMismatch {
            symphonia: stats.symphonia_remaining,
            ffmpeg: stats.ffmpeg_remaining,
        });
    }

    Ok(stats)
}

/// Test multiple files and return results for each.
pub fn test_decode_files<P: AsRef<Path>>(
    paths: &[P],
    config: &TestConfig,
) -> Vec<(String, TestResult<TestStats>)> {
    paths
        .iter()
        .map(|p| {
            let path_str = p.as_ref().to_string_lossy().to_string();
            let result = test_decode(p, config);
            (path_str, result)
        })
        .collect()
}

/// Collect all files matching a glob pattern in a directory.
pub fn collect_test_files<P: AsRef<Path>>(
    dir: P,
    extension: &str,
) -> std::io::Result<Vec<std::path::PathBuf>> {
    let dir = dir.as_ref();
    let mut files = Vec::new();

    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();

        if path.is_file() {
            if let Some(ext) = path.extension() {
                if ext.eq_ignore_ascii_case(extension) {
                    files.push(path);
                }
            }
        }
    }

    files.sort();
    Ok(files)
}

/// Encode a WAV file to a target format using the provided encoder function,
/// then compare Symphonia's decoding output against FFmpeg's.
///
/// The encoder function receives the WAV file path and a destination path for the
/// encoded file. It should return `Ok(())` on success.
///
/// # Arguments
///
/// * `wav_path` - Path to the input WAV file
/// * `encoded_extension` - File extension for the encoded file (e.g., "opus", "mp3")
/// * `encoder` - Function that encodes the WAV file to the target format
/// * `config` - Test configuration
///
/// Returns `TestStats` containing detailed comparison statistics.
pub fn compare_decode_from_wav<P, F>(
    wav_path: P,
    encoded_extension: &str,
    encoder: F,
    config: &TestConfig,
) -> TestResult<TestStats>
where
    P: AsRef<Path>,
    F: FnOnce(&Path, &Path) -> std::io::Result<()>,
{
    let wav_path = wav_path.as_ref();

    // Create a temporary file for the encoded output
    let encoded_path = std::env::temp_dir().join(format!(
        "symphonia_test_{}.{}",
        std::process::id(),
        encoded_extension
    ));

    // Encode the WAV file
    encoder(wav_path, &encoded_path)
        .map_err(|e| TestError::FfmpegSpawnFailed(format!("Encoder failed: {}", e)))?;

    // Run the comparison test on the encoded file
    let result = compare_decode(&encoded_path, config);

    // Clean up the temporary file
    let _ = std::fs::remove_file(&encoded_path);

    result
}

/// Test a WAV file by encoding it and comparing decode results.
///
/// This is a convenience wrapper around `compare_decode_from_wav` that validates
/// the results and returns an error if the test fails.
pub fn test_decode_from_wav<P, F>(
    wav_path: P,
    encoded_extension: &str,
    encoder: F,
    config: &TestConfig,
) -> TestResult<TestStats>
where
    P: AsRef<Path>,
    F: FnOnce(&Path, &Path) -> std::io::Result<()>,
{
    let stats = compare_decode_from_wav(wav_path, encoded_extension, encoder, config)?;

    if stats.failed_samples > 0 {
        return Err(TestError::SampleMismatch {
            packet: 0,
            sample: 0,
            symphonia: 0.0,
            ffmpeg: 0.0,
            delta: stats.max_delta,
        });
    }

    if stats.symphonia_remaining != 0 || stats.ffmpeg_remaining != 0 {
        return Err(TestError::RemainingSamplesMismatch {
            symphonia: stats.symphonia_remaining,
            ffmpeg: stats.ffmpeg_remaining,
        });
    }

    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_config_default() {
        let config = TestConfig::default();
        assert_eq!(config.ref_decoder, RefDecoder::Ffmpeg);
        assert!(config.gapless);
        assert!(!config.keep_going);
    }

    #[test]
    fn test_stats_passed() {
        let mut stats = TestStats::default();
        assert!(stats.passed());

        stats.failed_samples = 1;
        assert!(!stats.passed());

        stats.failed_samples = 0;
        stats.symphonia_remaining = 1;
        assert!(!stats.passed());
    }
}
