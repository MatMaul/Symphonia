use crate::util::test_ogg_vorbis_file;

mod util;

#[test]
fn vorbis_1_0_test() {
    test_ogg_vorbis_file("1.0-test.ogg").unwrap();
}
