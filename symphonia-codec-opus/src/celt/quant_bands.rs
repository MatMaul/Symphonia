use super::entdec::EcDec;
use super::laplace::ec_laplace_decode;
use super::modes::CeltMode;
use super::rate::MAX_FINE_BITS;
use super::types::CeltGlog;

include!("generated/quant_bands_tables.rs");

pub fn unquant_coarse_energy(
    m: &CeltMode,
    start: i32,
    end: i32,
    old_ebands: &mut [CeltGlog],
    intra: bool,
    dec: &mut EcDec<'_>,
    channels: i32,
    lm: i32,
) {
    let prob_model = &E_PROB_MODEL[lm as usize][if intra { 1 } else { 0 }];
    let mut prev = [0.0f32; 2];
    let (coef, beta) = if intra {
        (0.0, BETA_INTRA)
    } else {
        (PRED_COEF[lm as usize], BETA_COEF[lm as usize])
    };
    let budget = (dec.storage() * 8) as i32;
    let nb_ebands = m.nb_ebands as usize;

    for i in start..end {
        for c in 0..channels {
            let tell = dec.tell();
            let mut qi = if budget - tell >= 15 {
                let pi = 2 * i.min(20);
                ec_laplace_decode(
                    dec,
                    (prob_model[pi as usize] as u32) << 7,
                    (prob_model[pi as usize + 1] as i32) << 6,
                )
            } else if budget - tell >= 2 {
                let q = dec.dec_icdf(&SMALL_ENERGY_ICDF, 2);
                (q >> 1) ^ -((q & 1) as i32)
            } else if budget - tell >= 1 {
                -dec.dec_bit_logp(1)
            } else {
                -1
            };
            let q = qi as f32;
            let idx = (i as usize) + (c as usize) * nb_ebands;
            let old_e = old_ebands[idx].max(-9.0);
            let tmp = coef * old_e + prev[c as usize] + q;
            old_ebands[idx] = tmp;
            prev[c as usize] = prev[c as usize] + q - beta * q;
        }
    }
}

pub fn unquant_fine_energy(
    m: &CeltMode,
    start: i32,
    end: i32,
    old_ebands: &mut [CeltGlog],
    prev_quant: Option<&[i32]>,
    extra_quant: &[i32],
    dec: &mut EcDec<'_>,
    channels: i32,
) {
    let nb_ebands = m.nb_ebands as usize;
    let budget = (dec.storage() * 8) as i32;

    for i in start..end {
        let extra = extra_quant[i as usize];
        if extra <= 0 {
            continue;
        }
        if dec.tell() + channels * extra > budget {
            continue;
        }
        let prev = prev_quant.map(|q| q[i as usize]).unwrap_or(0);
        let extra_shift = 1.0 / (1u32 << extra) as f32;
        let prev_shift = 1.0 / (1u32 << prev) as f32;
        for c in 0..channels {
            let q2 = dec.dec_bits(extra as u32) as f32;
            let offset = ((q2 + 0.5) * extra_shift - 0.5) * prev_shift;
            let idx = (i as usize) + (c as usize) * nb_ebands;
            old_ebands[idx] += offset;
        }
    }
}

pub fn unquant_energy_finalise(
    m: &CeltMode,
    start: i32,
    end: i32,
    old_ebands: Option<&mut [CeltGlog]>,
    fine_quant: &[i32],
    fine_priority: &[i32],
    mut bits_left: i32,
    dec: &mut EcDec<'_>,
    channels: i32,
) {
    let nb_ebands = m.nb_ebands as usize;
    if let Some(old_ebands) = old_ebands {
        for prio in 0..2 {
            for i in start..end {
                if bits_left < channels {
                    return;
                }
                let fine = fine_quant[i as usize];
                if fine >= MAX_FINE_BITS || fine_priority[i as usize] != prio {
                    continue;
                }
                let fine_shift = 1.0 / (1u32 << (fine + 1)) as f32;
                for c in 0..channels {
                    let q2 = dec.dec_bits(1) as f32;
                    let offset = (q2 - 0.5) * fine_shift;
                    let idx = (i as usize) + (c as usize) * nb_ebands;
                    old_ebands[idx] += offset;
                    bits_left -= 1;
                }
            }
        }
    } else {
        for prio in 0..2 {
            for i in start..end {
                if bits_left < channels {
                    return;
                }
                let fine = fine_quant[i as usize];
                if fine >= MAX_FINE_BITS || fine_priority[i as usize] != prio {
                    continue;
                }
                for _ in 0..channels {
                    dec.dec_bits(1);
                    bits_left -= 1;
                }
            }
        }
    }
}
