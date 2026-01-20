// Symphonia
// Copyright (c) 2024 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! CELT mode configuration.
//!
//! Defines the mode-dependent parameters for CELT decoding.

use crate::celt::tables;

/// Number of bands in CELT.
pub const NB_EBANDS: usize = 21;

/// Number of allocation vectors.
pub const NB_ALLOC_VECTORS: usize = 11;

/// Maximum fine quantization bits.
pub const MAX_FINE_BITS: i32 = 8;

/// Cache for pulse-to-bits mapping.
#[derive(Clone)]
pub struct PulseCache {
    /// Cache bits lookup
    pub bits: &'static [i16],
    /// Cache index
    pub index: &'static [i16],
    /// Maximum bits capacity
    pub caps: &'static [i16],
}

impl Default for PulseCache {
    fn default() -> Self {
        Self::new()
    }
}

impl PulseCache {
    /// Create a new pulse cache for the standard mode.
    pub const fn new() -> Self {
        Self {
            bits: tables::CACHE_BITS50,
            index: tables::CACHE_INDEX50,
            caps: tables::CACHE_CAPS50,
        }
    }
}

/// CELT mode configuration.
///
/// This structure contains all the mode-dependent parameters needed
/// for encoding and decoding CELT frames.
#[derive(Clone)]
pub struct CeltMode {
    /// Number of frequency bands.
    pub nb_ebands: usize,

    /// Number of effective bands (depends on bandwidth).
    pub eff_ebands: usize,

    /// Band boundaries (in MDCT bins for short frames).
    pub ebands: &'static [i16],

    /// Log of band sizes.
    pub log_n: &'static [i16],

    /// Number of allocation vectors.
    pub nb_alloc_vectors: usize,

    /// Static allocation table.
    pub alloc_vectors: &'static [i16],

    /// Pulse cache for bits-to-pulses mapping.
    pub cache: PulseCache,

    /// Short MDCT size (2.5ms @ 48kHz = 120).
    pub short_mdct_size: usize,
}

impl Default for CeltMode {
    fn default() -> Self {
        Self::new()
    }
}

impl CeltMode {
    /// Create the standard CELT mode for 48kHz.
    pub const fn new() -> Self {
        Self {
            nb_ebands: NB_EBANDS,
            eff_ebands: NB_EBANDS,
            ebands: tables::EBAND5MS,
            log_n: tables::LOG_N400,
            nb_alloc_vectors: NB_ALLOC_VECTORS,
            alloc_vectors: tables::BAND_ALLOCATION,
            cache: PulseCache::new(),
            short_mdct_size: 120,
        }
    }

    /// Get the number of MDCT bins for a given band at a given LM.
    #[inline]
    pub fn band_size(&self, band: usize, lm: usize) -> usize {
        let m = 1 << lm;
        m * (self.ebands[band + 1] - self.ebands[band]) as usize
    }

    /// Get the starting bin for a given band at a given LM.
    #[inline]
    pub fn band_start(&self, band: usize, lm: usize) -> usize {
        let m = 1 << lm;
        m * self.ebands[band] as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_celt_mode() {
        let mode = CeltMode::new();
        assert_eq!(mode.nb_ebands, 21);
        assert_eq!(mode.ebands.len(), 22); // 21 bands + 1 endpoint

        // Check band sizes for LM=0 (2.5ms)
        assert_eq!(mode.band_size(0, 0), 1); // Band 0: 1 bin
        assert_eq!(mode.band_size(8, 0), 2); // Band 8: 2 bins (10-8)

        // Check band sizes for LM=3 (20ms)
        assert_eq!(mode.band_size(0, 3), 8); // Band 0: 8 bins
    }
}
