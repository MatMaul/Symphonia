use super::modes::{CeltCoef, CeltMode, KissFftState, KissTwiddleCpx, MdctLookup, PulseCache};

include!("generated/static_modes_tables.rs");

pub static FFT_STATE_48000_960_0: KissFftState = KissFftState {
    nfft: 480,
    scale: 0.0020833334,
    shift: -1,
    factors: [5, 96, 3, 32, 4, 8, 2, 4, 4, 1, 0, 0, 0, 0, 0, 0],
    bitrev: &FFT_BITREV_480,
    twiddles: &FFT_TWIDDLES_48000_960,
};

pub static FFT_STATE_48000_960_1: KissFftState = KissFftState {
    nfft: 240,
    scale: 0.0041666669,
    shift: 1,
    factors: [5, 48, 3, 16, 4, 4, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    bitrev: &FFT_BITREV_240,
    twiddles: &FFT_TWIDDLES_48000_960,
};

pub static FFT_STATE_48000_960_2: KissFftState = KissFftState {
    nfft: 120,
    scale: 0.0083333338,
    shift: 2,
    factors: [5, 24, 3, 8, 2, 4, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    bitrev: &FFT_BITREV_120,
    twiddles: &FFT_TWIDDLES_48000_960,
};

pub static FFT_STATE_48000_960_3: KissFftState = KissFftState {
    nfft: 60,
    scale: 0.016666668,
    shift: 3,
    factors: [5, 12, 3, 4, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    bitrev: &FFT_BITREV_60,
    twiddles: &FFT_TWIDDLES_48000_960,
};

pub static MODE_48000_960_120: CeltMode = CeltMode {
    fs: 48000,
    overlap: 120,
    nb_ebands: 21,
    eff_ebands: 21,
    preemph: [0.85000610, 0.0000000, 1.0000000, 1.0000000],
    ebands: &EBAND_5MS,
    max_lm: 3,
    nb_short_mdcts: 8,
    short_mdct_size: 120,
    nb_alloc_vectors: 11,
    alloc_vectors: &BAND_ALLOCATION,
    log_n: &LOG_N_400,
    window: &WINDOW_120,
    mdct: MdctLookup {
        n: 1920,
        maxshift: 3,
        kfft: [&FFT_STATE_48000_960_0, &FFT_STATE_48000_960_1, &FFT_STATE_48000_960_2, &FFT_STATE_48000_960_3],
        trig: &MDCT_TWIDDLES_960,
    },
    cache: PulseCache {
        size: 392,
        index: &CACHE_INDEX_50,
        bits: &CACHE_BITS_50,
        caps: &CACHE_CAPS_50,
    },
};

pub static STATIC_MODES: [&CeltMode; 1] = [&MODE_48000_960_120];

pub fn mode_from_static(fs: i32, frame_size: i32) -> Option<&'static CeltMode> {
    if fs == 48000 && frame_size == 960 {
        Some(&MODE_48000_960_120)
    } else {
        None
    }
}