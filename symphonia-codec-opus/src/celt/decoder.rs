use super::entdec::EcDec;

pub const TRIM_ICDF: [u8; 11] = [126, 124, 119, 109, 87, 41, 19, 9, 4, 2, 0];
pub const SPREAD_ICDF: [u8; 4] = [25, 23, 2, 0];

const TF_SELECT_TABLE: [[i8; 8]; 4] = [
    [0, -1, 0, -1, 0, -1, 0, -1],
    [0, -1, 0, -2, 1, 0, 1, -1],
    [0, -2, 0, -3, 2, 0, 1, -1],
    [0, -2, 0, -3, 3, 0, 1, -1],
];

pub fn tf_decode(
    start: i32,
    end: i32,
    is_transient: bool,
    tf_res: &mut [i32],
    lm: i32,
    dec: &mut EcDec<'_>,
) {
    let mut budget = dec.storage() as i32 * 8;
    let mut tell = dec.tell();
    let mut logp = if is_transient { 2 } else { 4 };
    let tf_select_rsv = lm > 0 && tell + logp + 1 <= budget;
    if tf_select_rsv {
        budget -= 1;
    }

    let mut tf_changed = 0;
    let mut curr = 0;
    for i in start..end {
        if tell + logp <= budget {
            curr ^= dec.dec_bit_logp(logp as u32);
            tell = dec.tell();
            tf_changed |= curr;
        }
        tf_res[i as usize] = curr;
        logp = if is_transient { 4 } else { 5 };
    }

    let mut tf_select = 0;
    if tf_select_rsv {
        let base = 4 * (is_transient as i32);
        let idx0 = (base + tf_changed) as usize;
        let idx1 = (base + 2 + tf_changed) as usize;
        if TF_SELECT_TABLE[lm as usize][idx0] != TF_SELECT_TABLE[lm as usize][idx1] {
            tf_select = dec.dec_bit_logp(1);
        }
    }

    let base = 4 * (is_transient as i32);
    for i in start..end {
        let idx = (base + 2 * tf_select + tf_res[i as usize]) as usize;
        tf_res[i as usize] = TF_SELECT_TABLE[lm as usize][idx] as i32;
    }
}
