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
mod quant_bands;
mod types;
mod vq;
mod synthesis;
mod bands_quant;

pub(crate) use entdec::EcDec;
pub(crate) use fft::{KissFftCpx, opus_fft, opus_fft_impl, opus_ifft};
pub(crate) use cwrs::{decode_pulses, get_required_bits};
pub(crate) use laplace::{ec_laplace_decode, ec_laplace_decode_p0};
pub(crate) use math::{
    celt_atan2p_norm, celt_atan_norm, celt_cos_norm, celt_div, celt_exp2, celt_exp2_db,
    celt_log2, celt_log2_db, celt_rcp, celt_rsqrt, celt_rsqrt_norm, celt_rsqrt_norm32,
    celt_sqrt, celt_sqrt32, isqrt32, PI,
};
pub(crate) use quant_bands::{unquant_coarse_energy, unquant_fine_energy, unquant_energy_finalise, E_MEANS};
pub(crate) use rate::{
    bits2pulses, get_pulses, pulses2bits, CELT_MAX_PULSES, FINE_OFFSET, LOG2_FRAC_TABLE,
    LOG_MAX_PSEUDO, MAX_FINE_BITS, MAX_PSEUDO, QTHETA_OFFSET, QTHETA_OFFSET_TWOPHASE,
};
pub(crate) use mdct::{mdct_backward, mdct_forward};
pub(crate) use modes::CeltMode;
pub(crate) use static_modes::mode_from_static;
pub(crate) use bands::{
    anti_collapse, bitexact_cos, bitexact_log2tan, celt_lcg_rand, denormalise_bands,
    compute_qn, deinterleave_hadamard, haar1, hysteresis_decision, intensity_stereo,
    interleave_hadamard, stereo_merge, stereo_split,
};
pub(crate) use types::{CeltEner, CeltGlog, CeltNorm, CeltRes, CeltSig, EPSILON, Q15_ONE, Q31_ONE, VERY_SMALL};
pub(crate) use vq::{alg_unquant, exp_rotation, renormalise_vector, SPREAD_AGGRESSIVE, SPREAD_LIGHT, SPREAD_NONE, SPREAD_NORMAL};
pub(crate) use synthesis::celt_synthesis;
pub(crate) use bands_quant::{
    BandCtx, SplitCtx, compute_theta_decode, quant_band_decode, quant_band_n1,
    quant_band_stereo_decode, quant_partition_decode,
};
