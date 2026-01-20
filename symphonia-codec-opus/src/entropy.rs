// Symphonia
// Copyright (c) 2024 The Project Symphonia Developers.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! Range decoder for entropy coding.
//!
//! This module implements a range coder for decoding bitstreams in Opus/CELT packets.
//! The range coder is used for all variable-length coding operations.

// Constants for the range coder
const EC_WINDOW_SIZE: usize = 32;
const EC_SYM_BITS: i32 = 8;
const EC_CODE_BITS: i32 = 32;
const EC_SYM_MAX: u32 = 0x0000_00FF;
const EC_CODE_SHIFT: i32 = 23;
const EC_CODE_TOP: u32 = 0x8000_0000;
const EC_CODE_BOT: u32 = 0x0080_0000;
const EC_CODE_EXTRA: i32 = 7;

/// Range decoder for entropy-coded bitstreams.
///
/// This implements the range coding algorithm used in Opus/CELT for
/// efficient variable-length coding.
pub struct RangeDecoder<'a> {
    /// Input buffer
    buf: &'a [u8],
    /// Total buffer size
    storage: usize,
    /// Read offset from start
    offs: usize,
    /// Read offset from end (for raw bits)
    end_offs: usize,
    /// Bits read from end
    end_window: u32,
    /// Number of valid bits in end_window
    nend_bits: i32,
    /// Total bits read
    nbits_total: i32,
    /// Range
    rng: u32,
    /// Value
    val: u32,
    /// Temporary storage for normalization
    rem: u8,
    /// Error flag
    error: bool,
}

impl<'a> RangeDecoder<'a> {
    /// Initialize the range decoder with a data buffer.
    pub fn init(buf: &'a [u8]) -> Self {
        let storage = buf.len();

        let mut dec = RangeDecoder {
            buf,
            storage,
            offs: 0,
            end_offs: 0,
            end_window: 0,
            nend_bits: 0,
            nbits_total: EC_CODE_BITS + 1
                - ((EC_CODE_BITS - EC_CODE_EXTRA) / EC_SYM_BITS) * EC_SYM_BITS,
            rng: 1 << EC_CODE_EXTRA,
            val: 0,
            rem: 0,
            error: false,
        };

        dec.rem = dec.read_byte();
        dec.val = (dec.rng - 1 - ((dec.rem as u32) >> (EC_SYM_BITS - EC_CODE_EXTRA))) as u32;
        dec.dec_normalize();
        dec
    }

    /// Read a byte from the start of the buffer.
    #[inline]
    fn read_byte(&mut self) -> u8 {
        if self.offs < self.storage {
            let val = self.buf[self.offs];
            self.offs += 1;
            val
        }
        else {
            0
        }
    }

    /// Read a byte from the end of the buffer (for raw bits).
    #[inline]
    fn read_byte_from_end(&mut self) -> u8 {
        if self.end_offs < self.storage {
            self.end_offs += 1;
            self.buf[self.storage - self.end_offs]
        }
        else {
            0
        }
    }

    /// Normalize the range to keep it in valid bounds.
    #[inline]
    fn dec_normalize(&mut self) {
        while self.rng <= EC_CODE_BOT {
            self.nbits_total += EC_SYM_BITS;
            self.rng <<= EC_SYM_BITS;

            // Use up the remaining bits from our last symbol
            let mut sym = self.rem as u32;

            // Read the next value from the input
            self.rem = self.read_byte();

            // Take the rest of the bits we need from this new symbol
            sym = ((sym << EC_SYM_BITS) | (self.rem as u32)) >> (EC_SYM_BITS - EC_CODE_EXTRA);

            // And subtract them from val, capped to be less than EC_CODE_TOP
            self.val = (((self.val as u64) << EC_SYM_BITS)
                + ((EC_SYM_MAX & !sym) as u64)) as u32
                & (EC_CODE_TOP - 1);
        }
    }

    /// Decode a bit with a given log probability.
    ///
    /// Returns 1 with probability 1/(2^logp), 0 otherwise.
    #[inline]
    pub fn decode_bit_logp(&mut self, logp: u32) -> bool {
        let r = self.rng;
        let d = self.val;
        let s = r >> logp;
        let ret = d < s;

        if !ret {
            self.val = d - s;
            self.rng = r - s;
        }
        else {
            self.rng = s;
        }

        self.dec_normalize();
        ret
    }

    /// Decode using an inverse cumulative distribution function.
    ///
    /// Returns the symbol index decoded from the ICDF table.
    pub fn decode_icdf(&mut self, icdf: &[u16], ftb: u32) -> i32 {
        let mut s = self.rng;
        let d = self.val;
        let r = s >> ftb;
        let mut ret = -1i32;
        let mut t;

        loop {
            ret += 1;
            t = s;
            s = r * (icdf[ret as usize] as u32);
            if d >= s {
                break;
            }
        }

        self.val = d - s;
        self.rng = t - s;
        self.dec_normalize();
        ret
    }

    /// Decode using an ICDF table with an offset.
    pub fn decode_icdf_offset(&mut self, icdf: &[u16], offset: usize, ftb: u32) -> i32 {
        let mut s = self.rng;
        let d = self.val;
        let r = s >> ftb;
        let mut ret = (offset as i32) - 1;
        let mut t;

        loop {
            ret += 1;
            t = s;
            s = r * (icdf[ret as usize] as u32);
            if d >= s {
                break;
            }
        }

        self.val = d - s;
        self.rng = t - s;
        self.dec_normalize();
        ret - (offset as i32)
    }

    /// Decode an unsigned integer uniformly distributed in [0, ft).
    pub fn decode_uint(&mut self, mut ft: u32) -> u32 {
        debug_assert!(ft > 1, "decode_uint: ft must be > 1");

        ft -= 1;
        let ftb = ec_ilog(ft);

        if ftb > EC_SYM_BITS as u32 {
            let ftb_red = ftb - (EC_SYM_BITS as u32);
            let ft_red = (ft >> ftb_red) + 1;

            let s = self.decode_value(ft_red);
            self.dec_update(s, s + 1, ft_red);

            let t = (s << ftb_red) | self.decode_bits(ftb_red as i32);
            if t <= ft {
                return t;
            }
            self.error = true;
            ft
        }
        else {
            ft += 1;
            let s = self.decode_value(ft);
            self.dec_update(s, s + 1, ft);
            s
        }
    }

    /// Decode raw bits from the end of the stream.
    pub fn decode_bits(&mut self, bits: i32) -> u32 {
        let mut window = self.end_window;
        let mut available = self.nend_bits;

        while available < bits {
            window |= (self.read_byte_from_end() as u32) << available;
            available += EC_SYM_BITS;
        }

        let ret = window & ((1 << bits) - 1);
        window >>= bits;
        available -= bits;

        self.end_window = window;
        self.nend_bits = available;
        ret
    }

    /// Helper: Decode a value from the range.
    #[inline]
    fn decode_value(&mut self, ft: u32) -> u32 {
        let ext = self.rng / ft;
        let s = self.val / ext;
        ft - s.min(ft) - 1
    }


    /// Get the number of bits read so far.
    pub fn tell(&self) -> i32 {
        self.nbits_total - ec_ilog(self.rng) as i32
    }

    /// Get the number of bits read with fractional precision.
    pub fn tell_frac(&self) -> i32 {
        let correction = [35733, 38967, 42495, 46340, 50535, 55109, 60097, 65535];
        let mut nbits = self.nbits_total << 3;
        let l = ec_ilog(self.rng);
        let r = (self.rng >> (l - 16)) as usize;
        nbits -= (l << 3) as i32;
        nbits += 3;
        nbits += (correction[(r & 7) as usize] as i32) >> (l - 9);
        nbits
    }

    /// Check if an error occurred during decoding.
    pub fn has_error(&self) -> bool {
        self.error
    }

    /// Get the storage size in bytes.
    pub fn storage(&self) -> usize {
        self.storage
    }

    /// Decode a value from [0, ft) using ftb bits of precision.
    ///
    /// Used internally for Laplace decoding and other operations.
    pub fn decode_bin(&mut self, ftb: u32) -> u32 {
        self.rng >>= ftb;
        let val = self.val / self.rng;
        let ft = 1u32 << ftb;
        ft.saturating_sub(val.min(ft) + 1)
    }

    /// Update the decoder state after decoding a symbol.
    ///
    /// # Arguments
    /// * `fl` - Lower bound of the symbol interval
    /// * `fh` - Upper bound of the symbol interval
    /// * `ft` - Total probability (must be a power of 2 for dec_update_bin)
    pub fn dec_update(&mut self, fl: u32, fh: u32, ft: u32) {
        let s = self.rng / ft * (ft - fh);
        self.val -= s;
        if fl > 0 {
            self.rng = self.rng / ft * (fh - fl);
        } else {
            self.rng -= s;
        }
        self.dec_normalize();
    }

    /// Decode a value from [0, ft) and return it.
    ///
    /// The caller must call dec_update with the symbol bounds.
    pub fn decode(&mut self, ft: u32) -> u32 {
        let ext = self.rng / ft;
        let s = self.val / ext;
        ft - s.min(ft) - 1
    }

    /// Decode a symbol with probability 1/2^logp.
    ///
    /// Returns 1 with probability 1/(2^logp), 0 otherwise.
    pub fn dec_bit_logp(&mut self, logp: u32) -> i32 {
        let r = self.rng;
        let d = self.val;
        let s = r >> logp;

        if d < s {
            self.rng = s;
            self.dec_normalize();
            1
        } else {
            self.val = d - s;
            self.rng = r - s;
            self.dec_normalize();
            0
        }
    }

    /// Decode a symbol using an ICDF table, returning the raw index.
    ///
    /// This variant uses i16 ICDF tables which is the common format.
    pub fn dec_icdf(&mut self, icdf: &[i16], ftb: u32) -> i32 {
        let mut s = self.rng;
        let d = self.val;
        let r = s >> ftb;
        let mut ret = -1i32;
        let mut t;

        loop {
            ret += 1;
            t = s;
            s = r * (icdf[ret as usize] as u32);
            if d >= s {
                break;
            }
        }

        self.val = d - s;
        self.rng = t - s;
        self.dec_normalize();
        ret
    }

    /// Decode raw bits (like dec_bits but returns i32).
    pub fn dec_bits(&mut self, bits: i32) -> i32 {
        self.decode_bits(bits) as i32
    }

    /// Decode an unsigned integer uniformly from [0, ft).
    pub fn dec_uint(&mut self, ft: u32) -> u32 {
        self.decode_uint(ft)
    }
}

/// Compute the integer log base 2 (number of bits required).
#[inline]
fn ec_ilog(x: u32) -> u32 {
    if x == 0 {
        1
    }
    else {
        32 - x.leading_zeros()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ec_ilog() {
        assert_eq!(ec_ilog(0), 1);
        assert_eq!(ec_ilog(1), 1);
        assert_eq!(ec_ilog(2), 2);
        assert_eq!(ec_ilog(3), 2);
        assert_eq!(ec_ilog(4), 3);
        assert_eq!(ec_ilog(255), 8);
        assert_eq!(ec_ilog(256), 9);
    }

    #[test]
    fn test_range_decoder_init() {
        let data = vec![0x12, 0x34, 0x56, 0x78];
        let _dec = RangeDecoder::init(&data);
        // Just check it initializes without panicking
    }

    #[test]
    fn test_decode_bits() {
        // Simple test: decode some raw bits
        let data = vec![0xFF, 0xFF, 0xFF, 0xFF];
        let mut dec = RangeDecoder::init(&data);

        let val = dec.decode_bits(4);
        assert_eq!(val, 0xF); // Should get 4 bits of 1s
    }
}
