#![allow(dead_code)]

mod bands;
mod cwrs;
mod entcode;
mod entdec;
mod fft;
mod intrin;
mod laplace;
mod math;
mod mfrngcod;
mod mdct;
mod modes;
mod rate;
mod static_modes;

pub(crate) use entdec::EcDec;
pub(crate) use fft::{KissFftCpx, opus_fft, opus_fft_impl, opus_ifft};
pub(crate) use cwrs::{decode_pulses, get_required_bits};
pub(crate) use laplace::{ec_laplace_decode, ec_laplace_decode_p0};
pub(crate) use math::{
    celt_atan2p_norm, celt_atan_norm, celt_cos_norm, celt_div, celt_rcp, celt_rsqrt,
    celt_rsqrt_norm, celt_rsqrt_norm32, celt_sqrt, celt_sqrt32, isqrt32, PI,
};
pub(crate) use rate::{
    bits2pulses, get_pulses, pulses2bits, CELT_MAX_PULSES, FINE_OFFSET, LOG2_FRAC_TABLE,
    LOG_MAX_PSEUDO, MAX_FINE_BITS, MAX_PSEUDO, QTHETA_OFFSET, QTHETA_OFFSET_TWOPHASE,
};
pub(crate) use mdct::{mdct_backward, mdct_forward};
pub(crate) use modes::CeltMode;
pub(crate) use static_modes::mode_from_static;
pub(crate) use bands::{bitexact_cos, bitexact_log2tan, celt_lcg_rand, hysteresis_decision};
