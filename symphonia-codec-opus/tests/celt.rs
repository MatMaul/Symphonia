mod util;

#[test]
fn test_celt_128kbps() {
    if let Err(e) = util::test_opus_celt_file("test01.wav") {
        panic!("{}", e);
    }
}
