pub type CeltCoef = f32;

pub const MAX_FACTORS: usize = 8;

#[derive(Clone, Copy)]
pub struct PulseCache {
    pub size: i32,
    pub index: &'static [i16],
    pub bits: &'static [u8],
    pub caps: &'static [u8],
}

#[derive(Clone, Copy)]
pub struct KissTwiddleCpx {
    pub r: CeltCoef,
    pub i: CeltCoef,
}

#[derive(Clone, Copy)]
pub struct KissFftState {
    pub nfft: i32,
    pub scale: CeltCoef,
    pub shift: i32,
    pub factors: [i16; 2 * MAX_FACTORS],
    pub bitrev: &'static [i16],
    pub twiddles: &'static [KissTwiddleCpx],
}

#[derive(Clone, Copy)]
pub struct MdctLookup {
    pub n: i32,
    pub maxshift: i32,
    pub kfft: [&'static KissFftState; 4],
    pub trig: &'static [CeltCoef],
}

#[derive(Clone, Copy)]
pub struct CeltMode {
    pub fs: i32,
    pub overlap: i32,
    pub nb_ebands: i32,
    pub eff_ebands: i32,
    pub preemph: [CeltCoef; 4],
    pub ebands: &'static [i16],
    pub max_lm: i32,
    pub nb_short_mdcts: i32,
    pub short_mdct_size: i32,
    pub nb_alloc_vectors: i32,
    pub alloc_vectors: &'static [u8],
    pub log_n: &'static [i16],
    pub window: &'static [CeltCoef],
    pub mdct: MdctLookup,
    pub cache: PulseCache,
}
