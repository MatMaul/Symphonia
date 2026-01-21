## Symphonia Opus (CELT) port notes

Purpose: quick orientation for future sessions, plus a running status of what
is already ported. Keep this file updated as work progresses.

### Repo structure (Opus-specific)
- `symphonia-codec-opus/src/lib.rs`: crate root; currently only exposes the `celt`
  module (no full decoder wired yet).
- `symphonia-codec-opus/src/celt/`: CELT port work-in-progress.
- `symphonia-codec-opus/opus-1.6.1/`: reference C implementation (Opus 1.6.1).

### Implemented CELT pieces (Rust)
- `symphonia-codec-opus/src/celt/entdec.rs`: range decoder (ec_dec).
- `symphonia-codec-opus/src/celt/entcode.rs`: range encoder helpers (partial).
- `symphonia-codec-opus/src/celt/mfrngcod.rs`: range coder core constants/logic.
- `symphonia-codec-opus/src/celt/static_modes.rs`: static mode tables + windows.
- `symphonia-codec-opus/src/celt/modes.rs`: `CeltMode` structs and lookup data.
- `symphonia-codec-opus/src/celt/fft.rs`: FFT core (Kiss FFT port).
- `symphonia-codec-opus/src/celt/mdct.rs`: MDCT forward/backward (core).
- `symphonia-codec-opus/src/celt/laplace.rs`: Laplace decoding.
- `symphonia-codec-opus/src/celt/cwrs.rs`: PVQ pulse tables + pulse decode.
- `symphonia-codec-opus/src/celt/math.rs`: math helpers (float path).
- `symphonia-codec-opus/src/celt/rate.rs`: rate utilities (pulses/bits).
- `symphonia-codec-opus/src/celt/bands.rs`: band helper utilities (hysteresis,
  bitexact cos/log2tan, LCG, denormalise, anti_collapse).
- `symphonia-codec-opus/src/celt/quant_bands.rs`: energy unquantization tables
  + decode helpers (coarse/fine/finalise).
- `symphonia-codec-opus/src/celt/vq.rs`: PVQ decode helpers (exp_rotation,
  alg_unquant, renormalise).
- `symphonia-codec-opus/src/celt/synthesis.rs`: CELT synthesis wrapper
  (denormalise + MDCT).
- `symphonia-codec-opus/src/celt/types.rs`: shared CELT float type aliases and
  constants.

### Not yet implemented (next targets)
- Decoder-side `quant_all_bands` flow from `bands.c`.
- Finish wiring CELT decoder with `celt_synthesis` and EcDec/static modes.
- Full decoder integration with Symphonia (`AudioDecoder` + packet decode).

### Architectural notes
- Optimized/arch-specific code is not ported yet. Keep opt functions isolated
  so future SIMD/asm ports can swap in cleanly.
- Both float and fixed-point exist in C; only float is ported so far.
  A cargo feature for fixed-point is still pending.

### Reference entry points (C)
- `opus-1.6.1/celt/celt_decoder.c`: `celt_decode_with_ec` decode flow.
- `opus-1.6.1/celt/bands.c`: `quant_all_bands`, helpers, and stereo decisions.
- `opus-1.6.1/celt/quant_bands.c`: energy quant/unquant.
- `opus-1.6.1/celt/vq.c`: PVQ decoding + spread/rotation helpers.
