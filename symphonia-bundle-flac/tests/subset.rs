// =============================================================================
// Block Size Tests (Files 01-10)
// =============================================================================

use crate::util::test_flac_file;

mod util;

#[test]
fn flac_01_blocksize_4096() {
    test_flac_file("subset/01 - blocksize 4096.flac").unwrap();
}

#[test]
fn flac_02_blocksize_4608() {
    test_flac_file("subset/02 - blocksize 4608.flac").unwrap();
}

#[test]
fn flac_03_blocksize_16() {
    test_flac_file("subset/03 - blocksize 16.flac").unwrap();
}

#[test]
fn flac_04_blocksize_192() {
    test_flac_file("subset/04 - blocksize 192.flac").unwrap();
}

#[test]
fn flac_05_blocksize_254() {
    test_flac_file("subset/05 - blocksize 254.flac").unwrap();
}

#[test]
fn flac_06_blocksize_512() {
    test_flac_file("subset/06 - blocksize 512.flac").unwrap();
}

#[test]
fn flac_07_blocksize_725() {
    test_flac_file("subset/07 - blocksize 725.flac").unwrap();
}

#[test]
fn flac_08_blocksize_1000() {
    test_flac_file("subset/08 - blocksize 1000.flac").unwrap();
}

#[test]
fn flac_09_blocksize_1937() {
    test_flac_file("subset/09 - blocksize 1937.flac").unwrap();
}

#[test]
fn flac_10_blocksize_2304() {
    test_flac_file("subset/10 - blocksize 2304.flac").unwrap();
}

// =============================================================================
// Advanced Encoding Features (Files 11-23)
// =============================================================================

#[test]
fn flac_11_partition_order_8() {
    test_flac_file("subset/11 - partition order 8.flac").unwrap();
}

#[test]
fn flac_12_qlp_precision_15_bit() {
    test_flac_file("subset/12 - qlp precision 15 bit.flac").unwrap();
}

#[test]
fn flac_13_qlp_precision_2_bit() {
    test_flac_file("subset/13 - qlp precision 2 bit.flac").unwrap();
}

#[test]
fn flac_14_wasted_bits() {
    test_flac_file("subset/14 - wasted bits.flac").unwrap();
}

#[test]
fn flac_15_only_verbatim_subframes() {
    test_flac_file("subset/15 - only verbatim subframes.flac").unwrap();
}

#[test]
fn flac_16_partition_order_8_escaped() {
    test_flac_file("subset/16 - partition order 8 containing escaped partitions.flac").unwrap();
}

#[test]
fn flac_17_all_fixed_orders() {
    test_flac_file("subset/17 - all fixed orders.flac").unwrap();
}

#[test]
fn flac_18_precision_search() {
    test_flac_file("subset/18 - precision search.flac").unwrap();
}

#[test]
fn flac_19_samplerate_35467hz() {
    test_flac_file("subset/19 - samplerate 35467Hz.flac").unwrap();
}

#[test]
fn flac_20_samplerate_39khz() {
    test_flac_file("subset/20 - samplerate 39kHz.flac").unwrap();
}

#[test]
fn flac_21_samplerate_22050hz() {
    test_flac_file("subset/21 - samplerate 22050Hz.flac").unwrap();
}

#[test]
fn flac_22_12_bit_per_sample() {
    test_flac_file("subset/22 - 12 bit per sample.flac").unwrap();
}

#[test]
fn flac_23_8_bit_per_sample() {
    test_flac_file("subset/23 - 8 bit per sample.flac").unwrap();
}

// =============================================================================
// Variable Block Size (Files 24-27)
// =============================================================================

#[test]
fn flac_24_variable_blocksize_flake_264() {
    test_flac_file("subset/24 - variable blocksize file created with flake revision 264.flac").unwrap();
}

#[test]
fn flac_25_variable_blocksize_flake_264_smaller() {
    test_flac_file("subset/25 - variable blocksize file created with flake revision 264, modified to create smaller blocks.flac").unwrap();
}

#[test]
fn flac_26_variable_blocksize_cuetools() {
    test_flac_file("subset/26 - variable blocksize file created with CUETools.Flake 2.1.6.flac").unwrap();
}

#[test]
fn flac_27_variable_blocksize_flake_011() {
    test_flac_file("subset/27 - old format variable blocksize file created with Flake 0.11.flac").unwrap();
}

// =============================================================================
// High Resolution Audio (Files 28-37)
// =============================================================================

#[test]
fn flac_28_high_resolution_default() {
    test_flac_file("subset/28 - high resolution audio, default settings.flac").unwrap();
}

#[test]
fn flac_29_high_resolution_blocksize_16384() {
    test_flac_file("subset/29 - high resolution audio, blocksize 16384.flac").unwrap();
}

#[test]
fn flac_30_high_resolution_blocksize_13456() {
    test_flac_file("subset/30 - high resolution audio, blocksize 13456.flac").unwrap();
}

#[test]
fn flac_31_high_resolution_32nd_order() {
    test_flac_file("subset/31 - high resolution audio, using only 32nd order predictors.flac").unwrap();
}

#[test]
fn flac_32_high_resolution_escaped_partitions() {
    test_flac_file("subset/32 - high resolution audio, partition order 8 containing escaped partitions.flac").unwrap();
}

#[test]
fn flac_33_samplerate_192khz() {
    test_flac_file("subset/33 - samplerate 192kHz.flac").unwrap();
}

#[test]
#[ignore]
// Quite slow
fn flac_34_samplerate_192khz_32nd_order() {
    test_flac_file("subset/34 - samplerate 192kHz, using only 32nd order predictors.flac").unwrap();
}

#[test]
fn flac_35_samplerate_134560hz() {
    test_flac_file("subset/35 - samplerate 134560Hz.flac").unwrap();
}

#[test]
fn flac_36_samplerate_384khz() {
    test_flac_file("subset/36 - samplerate 384kHz.flac").unwrap();
}

#[test]
fn flac_37_20_bit_per_sample() {
    test_flac_file("subset/37 - 20 bit per sample.flac").unwrap();
}

// =============================================================================
// Multi-Channel Audio (Files 38-44)
// =============================================================================

#[test]
fn flac_38_3_channels() {
    test_flac_file("subset/38 - 3 channels (3.0).flac").unwrap();
}

#[test]
fn flac_39_4_channels() {
    test_flac_file("subset/39 - 4 channels (4.0).flac").unwrap();
}

#[test]
fn flac_40_5_channels() {
    test_flac_file("subset/40 - 5 channels (5.0).flac").unwrap();
}

#[test]
fn flac_41_6_channels() {
    test_flac_file("subset/41 - 6 channels (5.1).flac").unwrap();
}

#[test]
fn flac_42_7_channels() {
    test_flac_file("subset/42 - 7 channels (6.1).flac").unwrap();
}

#[test]
fn flac_43_8_channels() {
    test_flac_file("subset/43 - 8 channels (7.1).flac").unwrap();
}

#[test]
fn flac_44_8_channel_192khz_24bit_32nd_order() {
    test_flac_file("subset/44 - 8-channel surround, 192kHz, 24 bit, using only 32nd order predictors.flac").unwrap();
}

// =============================================================================
// Metadata Edge Cases (Files 45-59)
// =============================================================================

#[test]
fn flac_45_no_total_samples() {
    test_flac_file("subset/45 - no total number of samples set.flac").unwrap();
}

#[test]
fn flac_46_no_min_max_framesize() {
    test_flac_file("subset/46 - no min-max framesize set.flac").unwrap();
}

#[test]
fn flac_47_only_streaminfo() {
    test_flac_file("subset/47 - only STREAMINFO.flac").unwrap();
}

#[test]
fn flac_48_extremely_large_seektable() {
    test_flac_file("subset/48 - Extremely large SEEKTABLE.flac").unwrap();
}

#[test]
fn flac_49_extremely_large_padding() {
    test_flac_file("subset/49 - Extremely large PADDING.flac").unwrap();
}

#[test]
fn flac_50_extremely_large_picture() {
    test_flac_file("subset/50 - Extremely large PICTURE.flac").unwrap();
}

#[test]
#[ignore]
// TODO fix
fn flac_51_extremely_large_vorbiscomment() {
    test_flac_file("subset/51 - Extremely large VORBISCOMMENT.flac").unwrap();
}

#[test]
fn flac_52_extremely_large_application() {
    test_flac_file("subset/52 - Extremely large APPLICATION.flac").unwrap();
}

#[test]
fn flac_53_cuesheet_many_indexes() {
    test_flac_file("subset/53 - CUESHEET with very many indexes.flac").unwrap();
}

#[test]
#[ignore]
// TODO fix
fn flac_54_1000x_vorbiscomment() {
    test_flac_file("subset/54 - 1000x repeating VORBISCOMMENT.flac").unwrap();
}

#[test]
#[ignore]
// TODO fix
fn flac_55_combined_48_53() {
    test_flac_file("subset/55 - file 48-53 combined.flac").unwrap();
}

#[test]
fn flac_56_jpg_picture() {
    test_flac_file("subset/56 - JPG PICTURE.flac").unwrap();
}

#[test]
fn flac_57_png_picture() {
    test_flac_file("subset/57 - PNG PICTURE.flac").unwrap();
}

#[test]
fn flac_58_gif_picture() {
    test_flac_file("subset/58 - GIF PICTURE.flac").unwrap();
}

#[test]
fn flac_59_avif_picture() {
    test_flac_file("subset/59 - AVIF PICTURE.flac").unwrap();
}

// =============================================================================
// Miscellaneous (Files 60-64)
// =============================================================================

#[test]
fn flac_60_mono_audio() {
    test_flac_file("subset/60 - mono audio.flac").unwrap();
}

#[test]
fn flac_61_predictor_overflow_16bit() {
    test_flac_file("subset/61 - predictor overflow check, 16-bit.flac").unwrap();
}

#[test]
fn flac_62_predictor_overflow_20bit() {
    test_flac_file("subset/62 - predictor overflow check, 20-bit.flac").unwrap();
}

#[test]
fn flac_63_predictor_overflow_24bit() {
    test_flac_file("subset/63 - predictor overflow check, 24-bit.flac").unwrap();
}

#[test]
fn flac_64_rice_escape_code_zero() {
    test_flac_file("subset/64 - rice partitions with escape code zero.flac").unwrap();
}
