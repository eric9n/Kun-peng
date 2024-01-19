use kr2r::mmscanner::MinimizerScanner;

fn main() {
    let seq: Vec<u8> = b"ACGATCGACGACG".to_vec();
    let mut scanner = MinimizerScanner::default(seq, 10, 5);
    let m1 = scanner.next_minimizer();
    println!("m1 {:?}", m1);
    assert_eq!(m1, Some(728));
    let m2 = scanner.next_minimizer();
    println!("m2 {:?}", m2);
    assert_eq!(m2, Some(536));
}
