use super::modes::{CeltCoef, KissFftState, KissTwiddleCpx, MAX_FACTORS};

#[derive(Clone, Copy, Default)]
pub struct KissFftCpx {
    pub r: CeltCoef,
    pub i: CeltCoef,
}

impl KissFftCpx {
    #[inline]
    fn add(self, other: Self) -> Self {
        Self { r: self.r + other.r, i: self.i + other.i }
    }

    #[inline]
    fn sub(self, other: Self) -> Self {
        Self { r: self.r - other.r, i: self.i - other.i }
    }

    #[inline]
    fn add_assign(&mut self, other: Self) {
        self.r += other.r;
        self.i += other.i;
    }

    #[inline]
    fn sub_assign(&mut self, other: Self) {
        self.r -= other.r;
        self.i -= other.i;
    }

    #[inline]
    fn mul_twiddle(self, tw: KissTwiddleCpx) -> Self {
        Self { r: self.r * tw.r - self.i * tw.i, i: self.r * tw.i + self.i * tw.r }
    }

    #[inline]
    fn mul_scalar(self, s: CeltCoef) -> Self {
        Self { r: self.r * s, i: self.i * s }
    }
}

fn kf_bfly2(fout: &mut [KissFftCpx], m: usize, n: usize) {
    if m == 1 {
        let mut idx = 0;
        for _ in 0..n {
            let t = fout[idx + 1];
            let f0 = fout[idx];
            fout[idx + 1] = f0.sub(t);
            fout[idx] = f0.add(t);
            idx += 2;
        }
    } else {
        const TW: CeltCoef = 0.7071067812;
        let mut idx = 0;
        for _ in 0..n {
            let base = idx;
            let base2 = idx + 4;

            let t0 = fout[base2];
            let f0 = fout[base];
            fout[base2] = f0.sub(t0);
            fout[base] = f0.add(t0);

            let t1 = KissFftCpx {
                r: (fout[base2 + 1].r + fout[base2 + 1].i) * TW,
                i: (fout[base2 + 1].i - fout[base2 + 1].r) * TW,
            };
            let f1 = fout[base + 1];
            fout[base2 + 1] = f1.sub(t1);
            fout[base + 1] = f1.add(t1);

            let t2 = KissFftCpx { r: fout[base2 + 2].i, i: -fout[base2 + 2].r };
            let f2 = fout[base + 2];
            fout[base2 + 2] = f2.sub(t2);
            fout[base + 2] = f2.add(t2);

            let t3 = KissFftCpx {
                r: (fout[base2 + 3].i - fout[base2 + 3].r) * TW,
                i: (-(fout[base2 + 3].i + fout[base2 + 3].r)) * TW,
            };
            let f3 = fout[base + 3];
            fout[base2 + 3] = f3.sub(t3);
            fout[base + 3] = f3.add(t3);

            idx += 8;
        }
    }
}

fn kf_bfly4(
    fout: &mut [KissFftCpx],
    fstride: usize,
    st: &KissFftState,
    m: usize,
    n: usize,
    mm: usize,
) {
    if m == 1 {
        for i in 0..n {
            let base = i * mm;
            let f0 = fout[base];
            let f1 = fout[base + 1];
            let f2 = fout[base + 2];
            let f3 = fout[base + 3];

            let scratch0 = f0.sub(f2);
            let mut f0_new = f0.add(f2);
            let scratch1 = f1.add(f3);
            let scratch1b = f1.sub(f3);

            fout[base + 2] = f0_new.sub(scratch1);
            f0_new.add_assign(scratch1);
            fout[base] = f0_new;

            fout[base + 1] = KissFftCpx { r: scratch0.r + scratch1b.i, i: scratch0.i - scratch1b.r };
            fout[base + 3] = KissFftCpx { r: scratch0.r - scratch1b.i, i: scratch0.i + scratch1b.r };
        }
    } else {
        let m2 = 2 * m;
        let m3 = 3 * m;
        for i in 0..n {
            let base = i * mm;
            let mut tw1 = 0usize;
            let mut tw2 = 0usize;
            let mut tw3 = 0usize;
            for j in 0..m {
                let idx = base + j;
                let s0 = fout[idx + m].mul_twiddle(st.twiddles[tw1]);
                let s1 = fout[idx + m2].mul_twiddle(st.twiddles[tw2]);
                let s2 = fout[idx + m3].mul_twiddle(st.twiddles[tw3]);

                let scratch5 = fout[idx].sub(s1);
                let mut f0_new = fout[idx].add(s1);
                let scratch3 = s0.add(s2);
                let scratch4 = s0.sub(s2);

                fout[idx + m2] = f0_new.sub(scratch3);
                f0_new.add_assign(scratch3);
                fout[idx] = f0_new;

                fout[idx + m] = KissFftCpx { r: scratch5.r + scratch4.i, i: scratch5.i - scratch4.r };
                fout[idx + m3] = KissFftCpx { r: scratch5.r - scratch4.i, i: scratch5.i + scratch4.r };

                tw1 += fstride;
                tw2 += fstride * 2;
                tw3 += fstride * 3;
            }
        }
    }
}

fn kf_bfly3(
    fout: &mut [KissFftCpx],
    fstride: usize,
    st: &KissFftState,
    m: usize,
    n: usize,
    mm: usize,
) {
    let epi3 = st.twiddles[fstride * m];
    let m2 = 2 * m;
    for i in 0..n {
        let base = i * mm;
        let mut tw1 = 0usize;
        let mut tw2 = 0usize;
        for j in 0..m {
            let idx = base + j;
            let s1 = fout[idx + m].mul_twiddle(st.twiddles[tw1]);
            let s2 = fout[idx + m2].mul_twiddle(st.twiddles[tw2]);

            let s3 = s1.add(s2);
            let mut s0 = s1.sub(s2);

            fout[idx + m] = KissFftCpx {
                r: fout[idx].r - 0.5 * s3.r,
                i: fout[idx].i - 0.5 * s3.i,
            };

            s0 = s0.mul_scalar(epi3.i);
            fout[idx].add_assign(s3);

            fout[idx + m2] = KissFftCpx {
                r: fout[idx + m].r + s0.i,
                i: fout[idx + m].i - s0.r,
            };
            fout[idx + m] = KissFftCpx {
                r: fout[idx + m].r - s0.i,
                i: fout[idx + m].i + s0.r,
            };

            tw1 += fstride;
            tw2 += fstride * 2;
        }
    }
}

fn kf_bfly5(
    fout: &mut [KissFftCpx],
    fstride: usize,
    st: &KissFftState,
    m: usize,
    n: usize,
    mm: usize,
) {
    let ya = st.twiddles[fstride * m];
    let yb = st.twiddles[fstride * 2 * m];
    for i in 0..n {
        let base = i * mm;
        for u in 0..m {
            let f0 = base + u;
            let f1 = f0 + m;
            let f2 = f1 + m;
            let f3 = f2 + m;
            let f4 = f3 + m;

            let scratch0 = fout[f0];
            let scratch1 = fout[f1].mul_twiddle(st.twiddles[u * fstride]);
            let scratch2 = fout[f2].mul_twiddle(st.twiddles[2 * u * fstride]);
            let scratch3 = fout[f3].mul_twiddle(st.twiddles[3 * u * fstride]);
            let scratch4 = fout[f4].mul_twiddle(st.twiddles[4 * u * fstride]);

            let scratch7 = scratch1.add(scratch4);
            let scratch10 = scratch1.sub(scratch4);
            let scratch8 = scratch2.add(scratch3);
            let scratch9 = scratch2.sub(scratch3);

            fout[f0].r += scratch7.r + scratch8.r;
            fout[f0].i += scratch7.i + scratch8.i;

            let scratch5 = KissFftCpx {
                r: scratch0.r + scratch7.r * ya.r + scratch8.r * yb.r,
                i: scratch0.i + scratch7.i * ya.r + scratch8.i * yb.r,
            };

            let scratch6 = KissFftCpx {
                r: scratch10.i * ya.i + scratch9.i * yb.i,
                i: -(scratch10.r * ya.i + scratch9.r * yb.i),
            };

            fout[f1] = scratch5.sub(scratch6);
            fout[f4] = scratch5.add(scratch6);

            let scratch11 = KissFftCpx {
                r: scratch0.r + scratch7.r * yb.r + scratch8.r * ya.r,
                i: scratch0.i + scratch7.i * yb.r + scratch8.i * ya.r,
            };
            let scratch12 = KissFftCpx {
                r: scratch9.i * ya.i - scratch10.i * yb.i,
                i: scratch10.r * yb.i - scratch9.r * ya.i,
            };

            fout[f2] = scratch11.add(scratch12);
            fout[f3] = scratch11.sub(scratch12);
        }
    }
}

pub fn opus_fft_impl(st: &KissFftState, fout: &mut [KissFftCpx]) {
    let mut fstride = [0usize; MAX_FACTORS + 1];
    let mut l = 0usize;
    let mut p = st.factors[0] as usize;
    let mut m = st.factors[1] as usize;
    let shift = if st.shift > 0 { st.shift as usize } else { 0 };

    fstride[0] = 1;
    loop {
        fstride[l + 1] = fstride[l] * p;
        l += 1;
        if m == 1 {
            break;
        }
        p = st.factors[2 * l] as usize;
        m = st.factors[2 * l + 1] as usize;
    }

    m = st.factors[2 * l - 1] as usize;
    for i in (0..l).rev() {
        let m2 = if i != 0 { st.factors[2 * i - 1] as usize } else { 1 };
        let radix = st.factors[2 * i] as usize;
        match radix {
            2 => {
                kf_bfly2(fout, m, fstride[i]);
            }
            4 => {
                kf_bfly4(fout, fstride[i] << shift, st, m, fstride[i], m2);
            }
            3 => {
                kf_bfly3(fout, fstride[i] << shift, st, m, fstride[i], m2);
            }
            5 => {
                kf_bfly5(fout, fstride[i] << shift, st, m, fstride[i], m2);
            }
            _ => {}
        }
        m = m2;
    }
}

pub fn opus_fft(st: &KissFftState, fin: &[KissFftCpx], fout: &mut [KissFftCpx]) {
    debug_assert_ne!(fin.as_ptr(), fout.as_ptr());
    debug_assert_eq!(fin.len(), st.nfft as usize);
    debug_assert_eq!(fout.len(), st.nfft as usize);

    for i in 0..st.nfft as usize {
        let x = fin[i];
        let idx = st.bitrev[i] as usize;
        fout[idx] = KissFftCpx { r: x.r * st.scale, i: x.i * st.scale };
    }

    opus_fft_impl(st, fout);
}

pub fn opus_ifft(st: &KissFftState, fin: &[KissFftCpx], fout: &mut [KissFftCpx]) {
    debug_assert_ne!(fin.as_ptr(), fout.as_ptr());
    debug_assert_eq!(fin.len(), st.nfft as usize);
    debug_assert_eq!(fout.len(), st.nfft as usize);

    for i in 0..st.nfft as usize {
        let idx = st.bitrev[i] as usize;
        fout[idx] = fin[i];
    }
    for v in fout.iter_mut() {
        v.i = -v.i;
    }
    opus_fft_impl(st, fout);
    for v in fout.iter_mut() {
        v.i = -v.i;
    }
}
