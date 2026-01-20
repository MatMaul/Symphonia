// Symphonia
// Copyright (c) 2024 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! CELT decoder constants.

/// Q15 representation of 1.0
pub const Q15ONE: i32 = 32767;

/// Signal scaling factor
pub const CELT_SIG_SCALE: f32 = 32768.0;

/// Signal shift for integer representation
pub const SIG_SHIFT: i32 = 12;

/// Normalization scaling
pub const NORM_SCALING: i32 = 16384;

/// Decibel shift for energy quantization (Q10 format)
pub const DB_SHIFT: i32 = 10;

/// Bit resolution for rate calculations
pub const BITRES: i32 = 3;

/// Small epsilon value
pub const EPSILON: i32 = 1;

/// Very small value (for comparisons)
pub const VERY_SMALL: i32 = 0;

/// Very large 16-bit value
pub const VERY_LARGE16: i16 = 32767;

/// Q15 one (16-bit)
pub const Q15_ONE: i16 = 32767;

/// Maximum comb filter period
pub const COMBFILTER_MAXPERIOD: usize = 1024;

/// Minimum comb filter period
pub const COMBFILTER_MINPERIOD: usize = 15;

/// Size of decode buffer (must hold maximum frame + overlap)
pub const DECODE_BUFFER_SIZE: usize = 2048;

/// Bit allocation table size
pub const BITALLOC_SIZE: usize = 11;

/// Maximum pitch period
pub const MAX_PERIOD: usize = 1024;

/// Total number of modes (for CELT-only, this is 1)
pub const TOTAL_MODES: usize = 1;

/// Maximum pseudo-random value
pub const MAX_PSEUDO: i32 = 40;

/// Log of maximum pseudo-random value
pub const LOG_MAX_PSEUDO: i32 = 6;

/// Maximum pulses in PVQ
pub const CELT_MAX_PULSES: i32 = 128;

/// Maximum fine energy bits
pub const MAX_FINE_BITS: i32 = 8;

/// Offset for fine energy bits
pub const FINE_OFFSET: i32 = 21;

/// Offset for theta quantization
pub const QTHETA_OFFSET: i32 = 4;

/// Offset for theta quantization (two-phase)
pub const QTHETA_OFFSET_TWOPHASE: i32 = 16;

/// Maximum pitch lag for PLC
pub const PLC_PITCH_LAG_MAX: i32 = 720;

/// Minimum pitch lag for PLC
pub const PLC_PITCH_LAG_MIN: i32 = 100;

/// LPC filter order
pub const LPC_ORDER: usize = 24;

/// Maximum number of bands
pub const CELT_MAX_BANDS: usize = 21;

/// Spread values for PVQ
pub const SPREAD_NONE: i32 = 0;
pub const SPREAD_LIGHT: i32 = 1;
pub const SPREAD_NORMAL: i32 = 2;
pub const SPREAD_AGGRESSIVE: i32 = 3;

/// Overlap for 48kHz (120 samples = 2.5ms)
pub const OVERLAP: usize = 120;

/// Short MDCT size for 48kHz (120 samples = 2.5ms)
pub const SHORT_MDCT_SIZE: usize = 120;

/// Maximum LM (log2 of multiplier: 2.5ms * 8 = 20ms)
pub const MAX_LM: usize = 3;

/// Maximum frame size in samples (20ms @ 48kHz = 960)
pub const MAX_FRAME_SIZE: usize = 960;

/// De-emphasis filter coefficients (alpha ~= 0.85)
pub const PREEMPH_COEF: i32 = 27853;

/// Sampling rate
pub const CELT_SAMPLING_RATE: i32 = 48000;
