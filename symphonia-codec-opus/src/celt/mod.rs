// Symphonia
// Copyright (c) 2024 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! CELT (Constrained Energy Lapped Transform) decoder.
//!
//! This module implements the CELT codec, which is used for general audio
//! and high bitrates in Opus. CELT uses MDCT for frequency domain coding
//! combined with vector quantization.

pub mod bands;
pub mod constants;
pub mod cwrs;
pub mod decoder;
pub mod laplace;
pub mod mode;
pub mod quant_bands;
pub mod rate;
pub mod synthesis;
pub mod tables;
pub mod vq;
