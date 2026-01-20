use std::path::PathBuf;

use symphonia_test::{TestConfig, compare_decode};

/// Get the path to the FLAC test files directory.
fn test_dir() -> PathBuf {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    PathBuf::from(manifest_dir).join("tests/files")
}

/// Check if FFmpeg is available on the system.
fn ffmpeg_available() -> bool {
    std::process::Command::new("ffmpeg")
        .arg("-version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .is_ok()
}

/// Run a comparison test for a single file, returning detailed error info on failure.
pub fn test_flac_file(filename: &str) -> Result<(), String> {
    if !ffmpeg_available() {
        return Err("FFmpeg not available".to_string());
    }

    let path = test_dir().join(filename);
    if !path.exists() {
        return Err(format!("file not found: {}, please initialize submodule", path.display()));
    }

    let config = TestConfig::default();

    match compare_decode(&path, &config) {
        Ok(stats) => {
            if stats.passed() {
                Ok(())
            }
            else {
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
