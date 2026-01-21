use super::entcode::{tell, tell_frac, EC_UINT_BITS};
use super::intrin::{celt_udiv, ec_mini, ec_ilog, imul32};
use super::mfrngcod::{
    EC_CODE_BOT, EC_CODE_EXTRA, EC_CODE_TOP, EC_SYM_BITS, EC_SYM_MAX, EC_WINDOW_SIZE,
};

pub struct EcDec<'a> {
    buf: &'a [u8],
    storage: u32,
    end_offs: u32,
    end_window: u32,
    nend_bits: i32,
    nbits_total: i32,
    offs: u32,
    rng: u32,
    val: u32,
    ext: u32,
    rem: u32,
    error: i32,
}

impl<'a> EcDec<'a> {
    pub fn new(buf: &'a [u8]) -> Self {
        let mut dec = Self {
            buf,
            storage: buf.len().min(u32::MAX as usize) as u32,
            end_offs: 0,
            end_window: 0,
            nend_bits: 0,
            nbits_total: 0,
            offs: 0,
            rng: 0,
            val: 0,
            ext: 0,
            rem: 0,
            error: 0,
        };
        dec.init();
        dec
    }

    pub fn range_bytes(&self) -> u32 {
        self.offs
    }

    pub fn storage(&self) -> u32 {
        self.storage
    }

    pub fn error(&self) -> i32 {
        self.error
    }

    pub fn tell(&self) -> i32 {
        tell(self.nbits_total, self.rng)
    }

    pub fn tell_frac(&self) -> u32 {
        tell_frac(self.nbits_total, self.rng)
    }

    pub fn decode(&mut self, ft: u32) -> u32 {
        let s;
        self.ext = celt_udiv(self.rng, ft);
        s = self.val / self.ext;
        ft - ec_mini(s + 1, ft)
    }

    pub fn decode_bin(&mut self, bits: u32) -> u32 {
        let s;
        self.ext = self.rng >> bits;
        s = self.val / self.ext;
        (1u32 << bits) - ec_mini(s + 1, 1u32 << bits)
    }

    pub fn update(&mut self, fl: u32, fh: u32, ft: u32) {
        let s = imul32(self.ext, ft - fh);
        self.val = self.val.wrapping_sub(s);
        self.rng = if fl > 0 {
            imul32(self.ext, fh - fl)
        } else {
            self.rng.wrapping_sub(s)
        };
        self.normalize();
    }

    pub fn dec_bit_logp(&mut self, logp: u32) -> i32 {
        let r = self.rng;
        let d = self.val;
        let s = r >> logp;
        let ret = d < s;
        if !ret {
            self.val = d - s;
        }
        self.rng = if ret { s } else { r - s };
        self.normalize();
        ret as i32
    }

    pub fn dec_icdf(&mut self, icdf: &[u8], ftb: u32) -> i32 {
        let mut s = self.rng;
        let d = self.val;
        let r = s >> ftb;
        let mut ret: i32 = -1;
        let mut t;
        loop {
            t = s;
            ret += 1;
            s = imul32(r, icdf[ret as usize] as u32);
            if d >= s {
                break;
            }
        }
        self.val = d - s;
        self.rng = t - s;
        self.normalize();
        ret
    }

    pub fn dec_icdf16(&mut self, icdf: &[u16], ftb: u32) -> i32 {
        let mut s = self.rng;
        let d = self.val;
        let r = s >> ftb;
        let mut ret: i32 = -1;
        let mut t;
        loop {
            t = s;
            ret += 1;
            s = imul32(r, icdf[ret as usize] as u32);
            if d >= s {
                break;
            }
        }
        self.val = d - s;
        self.rng = t - s;
        self.normalize();
        ret
    }

    pub fn dec_uint(&mut self, mut ft: u32) -> u32 {
        debug_assert!(ft > 1);
        ft -= 1;
        let mut ftb = ec_ilog(ft) as i32;
        if ftb > EC_UINT_BITS as i32 {
            ftb -= EC_UINT_BITS as i32;
            let ft_top = (ft >> (ftb as u32)) + 1;
            let s = self.decode(ft_top);
            self.update(s, s + 1, ft_top);
            let t = (s << (ftb as u32)) | self.dec_bits(ftb as u32);
            if t <= ft {
                return t;
            }
            self.error = 1;
            ft
        } else {
            ft += 1;
            let s = self.decode(ft);
            self.update(s, s + 1, ft);
            s
        }
    }

    pub fn dec_bits(&mut self, bits: u32) -> u32 {
        let mut window = self.end_window;
        let mut available = self.nend_bits;
        if (available as u32) < bits {
            loop {
                window |= self.read_byte_from_end() << (available as u32);
                available += EC_SYM_BITS as i32;
                if available > (EC_WINDOW_SIZE - EC_SYM_BITS) as i32 {
                    break;
                }
            }
        }
        let ret = window & ((1u32 << bits) - 1);
        window >>= bits;
        available -= bits as i32;
        self.end_window = window;
        self.nend_bits = available;
        self.nbits_total += bits as i32;
        ret
    }

    fn init(&mut self) {
        use super::mfrngcod::{EC_CODE_BITS, EC_CODE_EXTRA, EC_SYM_BITS};

        self.end_offs = 0;
        self.end_window = 0;
        self.nend_bits = 0;
        self.nbits_total = EC_CODE_BITS as i32 + 1
            - (((EC_CODE_BITS - EC_CODE_EXTRA) / EC_SYM_BITS) * EC_SYM_BITS) as i32;
        self.offs = 0;
        self.rng = 1u32 << EC_CODE_EXTRA;
        self.rem = self.read_byte();
        self.val = self
            .rng
            .wrapping_sub(1)
            .wrapping_sub(self.rem >> (EC_SYM_BITS - EC_CODE_EXTRA));
        self.error = 0;
        self.normalize();
    }

    fn read_byte(&mut self) -> u32 {
        if self.offs < self.storage {
            let b = self.buf[self.offs as usize] as u32;
            self.offs += 1;
            b
        } else {
            0
        }
    }

    fn read_byte_from_end(&mut self) -> u32 {
        if self.end_offs < self.storage {
            self.end_offs += 1;
            let index = self.storage - self.end_offs;
            self.buf[index as usize] as u32
        } else {
            0
        }
    }

    fn normalize(&mut self) {
        while self.rng <= EC_CODE_BOT {
            let mut sym = self.rem;
            self.nbits_total += EC_SYM_BITS as i32;
            self.rng <<= EC_SYM_BITS;
            self.rem = self.read_byte();
            sym = ((sym << EC_SYM_BITS) | self.rem) >> (EC_SYM_BITS - EC_CODE_EXTRA);
            let sym_mask = EC_SYM_MAX & !sym;
            self.val = (self.val << EC_SYM_BITS).wrapping_add(sym_mask) & (EC_CODE_TOP - 1);
        }
    }
}
