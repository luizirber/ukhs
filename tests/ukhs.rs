use proptest::prelude::*;

use ukhs::UKHS;

use std::collections::HashSet;
use std::iter::FromIterator;

proptest! {
    #[test]
    fn oracle_check(seq in "[ACGT]{20,}") {
        // TODO: implement for other k,w combinations
        let k = 7;
        let w = 20;
        let ukhs = UKHS::new(k, w).unwrap();

        let it = ukhs.iter_sequence(seq.as_bytes());
        let mut unikmers: Vec<String> = it.map(|(_, x)| x).collect();

        let it = ukhs.hash_iter_sequence(seq.as_bytes()).unwrap();
        let ukhs_hash: Vec<(u64, u64)> = it.collect();

        assert!(
            ukhs_hash.len() >= seq.len() - w + 1,
            "iter len {}, should be at least {}",
            ukhs_hash.len(),
            seq.len() - w + 1
        );

        let mut ukhs_unhash: Vec<String> = ukhs_hash
            .iter()
            .map(|(_, hash)| ukhs.kmer_for_ukhs_hash(*hash).unwrap())
            .collect();

        assert_eq!(unikmers, ukhs_unhash);
    }
}
