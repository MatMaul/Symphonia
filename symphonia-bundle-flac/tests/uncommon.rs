// =============================================================================
// Uncommon FLAC Files - Edge Cases and Unusual Configurations
// =============================================================================

use crate::util::test_flac_file;

mod util;

#[test]
#[ignore]
// TODO fix
fn flac_01_changing_samplerate() {
    test_flac_file("uncommon/01 - changing samplerate.flac").unwrap();
}

#[test]
#[ignore]
// TODO fix
fn flac_02_increasing_number_of_channels() {
    test_flac_file("uncommon/02 - increasing number of channels.flac").unwrap();
}

#[test]
#[ignore]
// TODO fix
fn flac_03_decreasing_number_of_channels() {
    test_flac_file("uncommon/03 - decreasing number of channels.flac").unwrap();
}

#[test]
fn flac_04_changing_bitdepth() {
    test_flac_file("uncommon/04 - changing bitdepth.flac").unwrap();
}

#[test]
#[ignore]
// TODO fix
fn flac_05_32bps_audio() {
    test_flac_file("uncommon/05 - 32bps audio.flac").unwrap();
}

#[test]
#[ignore]
// TODO fix
fn flac_06_samplerate_768khz() {
    test_flac_file("uncommon/06 - samplerate 768kHz.flac").unwrap();
}

#[test]
fn flac_07_15_bit_per_sample() {
    test_flac_file("uncommon/07 - 15 bit per sample.flac").unwrap();
}

#[test]
fn flac_08_blocksize_65535() {
    test_flac_file("uncommon/08 - blocksize 65535.flac").unwrap();
}

#[test]
fn flac_09_rice_partition_order_15() {
    test_flac_file("uncommon/09 - Rice partition order 15.flac").unwrap();
}

#[test]
#[ignore]
// TODO fix
fn flac_10_file_starting_at_frame_header() {
    test_flac_file("uncommon/10 - file starting at frame header.flac").unwrap();
}

#[test]
#[ignore]
// TODO fix
fn flac_11_file_starting_with_unparsable_data() {
    test_flac_file("uncommon/11 - file starting with unparsable data.flac").unwrap();
}
