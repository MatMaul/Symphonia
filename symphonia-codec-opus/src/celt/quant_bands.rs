use super::entdec::EcDec;
use super::laplace::ec_laplace_decode;
use super::modes::CeltMode;
use super::rate::MAX_FINE_BITS;
use super::types::CeltGlog;

pub const E_MEANS: [CeltGlog; 25] = [
    6.437500, 6.250000, 5.750000, 5.312500, 5.062500, 4.812500, 4.500000, 4.375000, 4.875000,
    4.687500, 4.562500, 4.437500, 4.875000, 4.625000, 4.312500, 4.500000, 4.375000, 4.625000,
    4.750000, 4.437500, 3.750000, 3.750000, 3.750000, 3.750000, 3.750000,
];

const PRED_COEF: [f32; 4] = [
    29440.0 / 32768.0,
    26112.0 / 32768.0,
    21248.0 / 32768.0,
    16384.0 / 32768.0,
];

const BETA_COEF: [f32; 4] = [
    30147.0 / 32768.0,
    22282.0 / 32768.0,
    12124.0 / 32768.0,
    6554.0 / 32768.0,
];

const BETA_INTRA: f32 = 4915.0 / 32768.0;

const E_PROB_MODEL: [[[u8; 42]; 2]; 4] = [
    [
        [
            72, 127, 65, 129, 66, 128, 65, 128, 64, 128, 62, 128, 64, 128, 64, 128, 92, 78,
            92, 79, 92, 78, 90, 79, 116, 41, 115, 40, 114, 40, 132, 26, 132, 26, 145, 17,
            161, 12, 176, 10, 177, 11,
        ],
        [
            24, 179, 48, 138, 54, 135, 54, 132, 53, 134, 56, 133, 55, 132, 55, 132, 61, 114,
            70, 96, 74, 88, 75, 88, 87, 74, 89, 66, 91, 67, 100, 59, 108, 50, 120, 40,
            122, 37, 97, 43, 78, 50,
        ],
    ],
    [
        [
            83, 78, 84, 81, 88, 75, 86, 74, 87, 71, 90, 73, 93, 74, 93, 74, 109, 40, 114, 36,
            117, 34, 117, 34, 143, 17, 145, 18, 146, 19, 162, 12, 165, 10, 178, 7, 189, 6,
            190, 8, 177, 9,
        ],
        [
            23, 178, 54, 115, 63, 102, 66, 98, 69, 99, 74, 89, 71, 91, 73, 91, 78, 89, 86, 80,
            92, 66, 93, 64, 102, 59, 103, 60, 104, 60, 117, 52, 123, 44, 138, 35, 133, 31,
            97, 38, 77, 45,
        ],
    ],
    [
        [
            61, 90, 93, 60, 105, 42, 107, 41, 110, 45, 116, 38, 113, 38, 112, 38, 124, 26,
            132, 27, 136, 19, 140, 20, 155, 14, 159, 16, 158, 18, 170, 13, 177, 10, 187, 8,
            192, 6, 175, 9, 159, 10,
        ],
        [
            21, 178, 59, 110, 71, 86, 75, 85, 84, 83, 91, 66, 88, 73, 87, 72, 92, 75, 98, 72,
            105, 58, 107, 54, 115, 52, 114, 55, 112, 56, 129, 51, 132, 40, 150, 33, 140, 29,
            98, 35, 77, 42,
        ],
    ],
    [
        [
            42, 121, 96, 66, 108, 43, 111, 40, 117, 44, 123, 32, 120, 36, 119, 33, 127, 33,
            134, 34, 139, 21, 147, 23, 152, 20, 158, 25, 154, 26, 166, 21, 173, 16, 184, 13,
            184, 10, 150, 13, 139, 15,
        ],
        [
            22, 178, 63, 114, 74, 82, 84, 83, 92, 82, 103, 62, 96, 72, 96, 67, 101, 73, 107, 72,
            113, 55, 118, 52, 125, 52, 118, 52, 117, 55, 135, 49, 137, 39, 157, 32, 145, 29,
            97, 33, 77, 40,
        ],
    ],
];

const SMALL_ENERGY_ICDF: [u8; 3] = [2, 1, 0];

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
                    (prob_model[pi as usize + 1] as u32) << 6,
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
