use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use symphonia_test::{compare_decode_from_wav, TestConfig};

/// Get the path to the Opus test files directory.
fn test_dir() -> PathBuf {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    PathBuf::from(manifest_dir).join("tests/files")
}

/// Check if FFmpeg is available on the system.
fn ffmpeg_available() -> bool {
    Command::new("ffmpeg")
        .arg("-version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok()
}

/// Check if opusenc is available on the system.
fn opusenc_available() -> bool {
    Command::new("opusenc")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok()
}

/// Encode a WAV file to Opus using opusenc with the specified arguments.
fn encode_opus(wav_path: &Path, opus_path: &Path, args: &[&str]) -> std::io::Result<()> {
    let status = Command::new("opusenc")
        .args(args)
        .arg(wav_path)
        .arg(opus_path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()?;

    if status.success() {
        Ok(())
    } else {
        Err(std::io::Error::new(
            std::io::ErrorKind::Other,
            format!("opusenc failed with exit code: {:?}", status.code()),
        ))
    }
}

/// Run a comparison test for a WAV file encoded to CELT-only Opus.
///
/// Uses opusenc with `--music --vbr --bitrate 128` to produce a CELT-only stream.
pub fn test_opus_celt_file(filename: &str) -> Result<(), String> {
    if !ffmpeg_available() {
        return Err("FFmpeg not available".to_string());
    }

    if !opusenc_available() {
        return Err("opusenc not available".to_string());
    }

    let wav_path = test_dir().join(filename);
    if !wav_path.exists() {
        return Err(format!("file not found: {}", wav_path.display()));
    }

    let config = TestConfig::default();

    let encoder = |wav: &Path, opus: &Path| {
        encode_opus(wav, opus, &["--music", "--vbr", "--bitrate", "128"])
    };

    match compare_decode_from_wav(&wav_path, "opus", encoder, &config) {
        Ok(stats) => {
            if stats.passed() {
                Ok(())
            } else {
                Err(format!(
                    "Test failed for {}: {} failed samples out of {}, max delta: {}, \
                     symphonia remaining: {}, ffmpeg remaining: {}",
                    filename,
                    stats.failed_samples,
                    stats.samples,
                    stats.max_delta,
                    stats.symphonia_remaining,
                    stats.ffmpeg_remaining
                ))
            }
        }
        Err(e) => Err(format!("Test error for {}: {}", filename, e)),
    }
}
