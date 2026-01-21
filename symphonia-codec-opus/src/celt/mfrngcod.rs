pub const EC_SYM_BITS: u32 = 8;
pub const EC_CODE_BITS: u32 = 32;
pub const EC_SYM_MAX: u32 = (1u32 << EC_SYM_BITS) - 1;
pub const EC_CODE_SHIFT: u32 = EC_CODE_BITS - EC_SYM_BITS - 1;
pub const EC_CODE_TOP: u32 = 1u32 << (EC_CODE_BITS - 1);
pub const EC_CODE_BOT: u32 = EC_CODE_TOP >> EC_SYM_BITS;
pub const EC_CODE_EXTRA: u32 = ((EC_CODE_BITS - 2) % EC_SYM_BITS) + 1;
pub const EC_WINDOW_SIZE: u32 = 32;
