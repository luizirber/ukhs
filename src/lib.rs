#![allow(clippy::unreadable_literal)]
#![allow(clippy::len_without_is_empty)]

pub mod errors;

use std::collections::{HashMap, VecDeque};
use std::iter::Peekable;
use std::str;

use bbhash::MPHF;
use nthash::{ntf64, NtHashForwardIterator};
use once_cell::sync::Lazy;

pub use crate::errors::UKHSError as Error;

static UKHS_HASHES: Lazy<HashMap<(usize, usize), &'static str>> = Lazy::new(|| {
    let mut map = HashMap::new();
    map.insert((7, 20), include_str!("../data/res_7_20_4_0.txt"));
    map.insert((9, 20), include_str!("../data/res_9_20_4_0.txt"));
    map.insert((9, 30), include_str!("../data/res_9_30_4_0.txt"));
    // TODO: add others
    map
});

#[derive(Clone)]
pub struct UKHS {
    k: usize,
    w: usize,
    mphf: MPHF,
    revmap: Vec<u64>,
    kmers: Vec<String>,
    kmers_hashes: Vec<u64>,
}

impl<'a> UKHS {
    pub fn new(k: usize, w: usize) -> Result<UKHS, Error> {
        if k > w {
            return Err(Error::KSizeOutOfWRange { ksize: k, wsize: w });
        }

        let w_round = (w / 10) * 10;
        let entries = UKHS_HASHES[&(k, w_round)];
        let mut kmers: Vec<String> = entries
            .split('\n')
            .filter_map(|s| {
                if s.len() == k {
                    Some(String::from(s))
                } else {
                    None
                }
            })
            .collect();

        // TODO: is the order relevant for interoperability?
        // for now this is necessary to make binary_search work
        kmers.sort_unstable();

        let kmers_hashes: Vec<u64> = kmers.iter().map(|h| ntf64(h.as_bytes(), 0, k)).collect();

        let mphf = MPHF::new(kmers_hashes.clone(), 1, 1.0); // TODO: any way to avoid this clone?
        let mut revmap = vec![0; kmers_hashes.len()];
        for hash in &kmers_hashes {
            revmap[mphf.lookup(*hash).unwrap() as usize] = *hash;
        }

        Ok(UKHS {
            k,
            w,
            mphf,
            revmap,
            kmers,
            kmers_hashes,
        })
    }

    pub fn len(&self) -> usize {
        self.kmers_hashes.len()
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn w(&self) -> usize {
        self.w
    }

    pub fn query_bucket(&self, hash: u64) -> Option<usize> {
        if let Some(pos) = self.mphf.lookup(hash) {
            if self.revmap[pos as usize] == hash {
                return Some(pos as usize);
            }
        }
        None
    }

    /// Creates a new UKHSIterator with internal state properly initialized.
    pub fn iter_sequence(&'a self, seq: &'a [u8]) -> UKHSIterator<'a> {
        let mut max_idx = seq.len() - self.k + 1;

        if self.k > seq.len() || self.w > seq.len() {
            // In these cases the sequence is too small for having any k-mers or
            // w-mers, so the iterator will return None right away.
            max_idx = 0;
        }

        UKHSIterator {
            seq,
            ukhs: self,
            current_k_idx: 0,
            current_w_idx: 0,
            max_idx,
        }
    }

    /// Creates a new UKHSHashIterator with internal state properly initialized.
    pub fn hash_iter_sequence(&'a self, seq: &'a [u8]) -> Result<UKHSHashIterator<'a>, Error> {
        if self.k > seq.len() {
            return Err(Error::KSizeOutOfRange {
                ksize: self.k,
                sequence: String::from_utf8(seq.to_vec()).unwrap(),
            });
        }

        if self.w > seq.len() {
            return Err(Error::WSizeOutOfRange {
                wsize: self.w,
                sequence: String::from_utf8(seq.to_vec()).unwrap(),
            });
        }

        let mut nthash_k_iter = NtHashForwardIterator::new(seq, self.k)?.peekable();
        let mut nthash_w_iter = NtHashForwardIterator::new(seq, self.w)?.peekable();

        let current_w_hash = nthash_w_iter.next().unwrap();
        let mut current_unikmers = VecDeque::with_capacity(self.k);

        for i in 0..=self.w - self.k {
            let k_hash = nthash_k_iter.next().unwrap();

            if self.contains(k_hash) {
                current_unikmers.push_back((i, k_hash));
            }
        }

        let current_w_idx = 0;
        let current_k_idx = 0;
        let max_k_pos = seq.len() - self.k + 1;

        Ok(UKHSHashIterator {
            nthash_k_iter,
            nthash_w_iter,
            ukhs: self,
            current_w_hash,
            current_w_idx,
            current_k_idx,
            max_k_pos,
            current_unikmers,
        })
    }

    pub fn contains(&self, hash: u64) -> bool {
        if let Some(pos) = self.mphf.lookup(hash) {
            if self.revmap[pos as usize] == hash {
                return true;
            }
        }
        false
    }

    pub fn contains_kmer(&self, kmer: &str) -> bool {
        self.kmers.binary_search(&kmer.into()).is_ok()
    }

    pub fn kmer_for_ukhs_hash(&self, hash: u64) -> Option<String> {
        if let Some(pos) = self.kmers_hashes.iter().position(|&x| x == hash) {
            Some(self.kmers[pos].clone())
        } else {
            None
        }
    }
}

/// An iterator for finding universal hitting k-mers in a sequence.
///
/// ```
///     use ukhs::UKHS;
///
///     # fn main() -> Result<(), ukhs::Error> {
///     let seq = b"ACACCGTAGCCTCCAGATGC";
///     let ukhs = UKHS::new(7, 20)?;
///
///     let it = ukhs.iter_sequence(seq);
///     let ukhs: Vec<(String, String)> = it.collect();
///     assert_eq!(ukhs,
///                [
///                    ("ACACCGTAGCCTCCAGATGC".into(), "ACACCGT".into()),
///                    ("ACACCGTAGCCTCCAGATGC".into(), "CCGTAGC".into()),
///                    ("ACACCGTAGCCTCCAGATGC".into(), "AGCCTCC".into()),
///                    ("ACACCGTAGCCTCCAGATGC".into(), "GCCTCCA".into())
///                ]);
///     # Ok(())
///     # }
/// ```
pub struct UKHSIterator<'a> {
    seq: &'a [u8],
    ukhs: &'a UKHS,
    current_k_idx: usize,
    current_w_idx: usize,
    max_idx: usize,
}

impl<'a> Iterator for UKHSIterator<'a> {
    type Item = (String, String);

    fn next(&mut self) -> Option<Self::Item> {
        while self.current_k_idx != self.max_idx {
            let last_k_pos = self.current_w_idx + self.ukhs.w() - self.ukhs.k() + 1;

            if self.current_k_idx == last_k_pos {
                self.current_w_idx += 1;
                self.current_k_idx = self.current_w_idx;
            }

            let current_kmer =
                str::from_utf8(&self.seq[self.current_k_idx..self.current_k_idx + self.ukhs.k])
                    .unwrap();

            self.current_k_idx += 1;

            if self.ukhs.contains_kmer(current_kmer) {
                let wmer_start = self.current_w_idx;
                let current_wmer =
                    str::from_utf8(&self.seq[wmer_start..wmer_start + self.ukhs.w]).unwrap();

                return Some((current_wmer.into(), current_kmer.into()));
            };
        }
        None
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.max_idx, Some(self.max_idx))
    }
}

impl<'a> ExactSizeIterator for UKHSIterator<'a> {}

/// An iterator for finding universal hitting k-mers in a sequence.
/// It uses ntHash for efficient k-mer hashing.
///
/// ```
///     use ukhs::UKHS;
///
///     # fn main() -> Result<(), ukhs::Error> {
///     let seq = b"ACACCGTAGCCTCCAGATGC";
///     let ukhs = UKHS::new(7, 20)?;
///
///     let it = ukhs.hash_iter_sequence(seq)?;
///     let matches: Vec<(u64, u64)> = it.collect();
///     assert_eq!(matches,
///                [
///                    (0x37137c91412bb512, 0xfbd9591aa929c685),
///                    (0x37137c91412bb512, 0x9cd9a1bcb779d6ad),
///                    (0x37137c91412bb512, 0x46fa47d28c0ffba5),
///                    (0x37137c91412bb512, 0xf482addc6edbc920)
///                ]);
///     # Ok(())
///     # }
/// ```
///
/// If you're only interested in the UKHS hashes (and not which w-mer it comes
/// from) you can also collect only the UKHS hashes:
///
/// ```
///     use ukhs::UKHS;
///
///     # fn main() -> Result<(), ukhs::Error> {
///     let ukhs = UKHS::new(7, 20)?;
///     let seq = b"ACACCGTAGCCTCCAGATGC";
///     let it = ukhs.hash_iter_sequence(seq)?;
///     let hashes: Vec<u64> = it.map(|(_, hash)| hash).collect();
///     assert_eq!(hashes, [0xfbd9591aa929c685, 0x9cd9a1bcb779d6ad, 0x46fa47d28c0ffba5, 0xf482addc6edbc920]);
///     # Ok(())
///     # }
/// ```
pub struct UKHSHashIterator<'a> {
    ukhs: &'a UKHS,
    nthash_k_iter: Peekable<NtHashForwardIterator<'a>>,
    nthash_w_iter: Peekable<NtHashForwardIterator<'a>>,
    current_w_hash: u64,
    current_w_idx: usize,
    current_k_idx: usize,
    max_k_pos: usize,
    current_unikmers: VecDeque<(usize, u64)>,
}

impl<'a> Iterator for UKHSHashIterator<'a> {
    type Item = (u64, u64);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // We're past the last possible k-mer; stop.
            if self.current_k_idx >= self.max_k_pos {
                return None;
            }

            let last_k_pos = self.current_w_idx + self.ukhs.w() - self.ukhs.k() + 1;

            // First state: draining queue
            if self.current_k_idx < last_k_pos {
                self.current_k_idx += 1;

                if let Some((pos, k_hash)) = self
                    .current_unikmers
                    .iter()
                    .find(|(p, _)| *p >= self.current_k_idx - 1)
                {
                    self.current_k_idx = *pos + 1;
                    return Some((self.current_w_hash, *k_hash));
                } else {
                    // In this case we went through all the current_unikmers;
                    // set current_k_idx to last_k_pos
                    self.current_k_idx = last_k_pos;
                }
            } else {
                // Second state:
                // - if front of current_unikmers < current_idx, pop it

                self.current_w_idx += 1;

                if let Some((pos, _)) = self.current_unikmers.front() {
                    if *pos < self.current_w_idx {
                        self.current_unikmers.pop_front();
                    }
                }

                // - push_back next k-mer (if it is inside next w-mer)
                let new_k_hash = self.nthash_k_iter.next().unwrap();

                if self.ukhs.contains(new_k_hash) {
                    self.current_unikmers.push_back((last_k_pos, new_k_hash));
                }

                // - advance to next w-mer
                self.current_w_hash = self.nthash_w_iter.next().unwrap();

                self.current_k_idx = self.current_w_idx;
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        // TODO: max_idx is the is the maximum,
        // but very unlikely (every k-mer would have to be a UKHS k-mer.
        self.nthash_k_iter.size_hint()
    }
}

impl<'a> ExactSizeIterator for UKHSHashIterator<'a> {}

#[cfg(test)]
mod test {
    use super::*;

    use std::collections::HashSet;
    use std::iter::FromIterator;

    #[test]
    fn basic_check() {
        let seq = b"ACACCGTAGCCTCCAGATGC";
        let w = 20;
        let k = 7;
        let ukhs = UKHS::new(k, w).unwrap();

        let it = ukhs.iter_sequence(seq);
        let mut unikmers: Vec<String> = it.map(|(_, x)| x).collect();
        unikmers.sort_unstable();

        assert_eq!(unikmers, ["ACACCGT", "AGCCTCC", "CCGTAGC", "GCCTCCA"]);

        let it = ukhs.hash_iter_sequence(seq).unwrap();
        let ukhs_hash: Vec<(u64, u64)> = it.collect();
        assert!(ukhs_hash.len() > seq.len() - w);

        let ukhs_unhash_set: HashSet<String> = ukhs_hash
            .iter()
            .map(|(_, hash)| ukhs.kmer_for_ukhs_hash(*hash).unwrap())
            .collect();
        let mut ukhs_unhash = Vec::<String>::from_iter(ukhs_unhash_set.into_iter());
        ukhs_unhash.sort_unstable();

        assert_eq!(unikmers, ukhs_unhash);
        assert_eq!(
            ukhs_hash,
            [
                (0x37137c91412bb512, 0xfbd9591aa929c685),
                (0x37137c91412bb512, 0x9cd9a1bcb779d6ad),
                (0x37137c91412bb512, 0x46fa47d28c0ffba5),
                (0x37137c91412bb512, 0xf482addc6edbc920)
            ]
        );
    }

    #[test]
    fn longer_check() {
        let seq = b"ACACCGTAGCCTCCAGATGCGTAG";
        /*
                 ACACCGTAGCCTCCAGATGCGTAG
                 *  *   **       *

        w-mer 0: ACACCGTAGCCTCCAGATGC
                 *  *   **
        w-mer 1: CACCGTAGCCTCCAGATGCG
                   *   **
        w-mer 2: ACCGTAGCCTCCAGATGCGT
                  *   **
        w-mer 3: CCGTAGCCTCCAGATGCGTA
                 *   **       *
        w-mer 4: CGTAGCCTCCAGATGCGTAG
                    **       *
        */
        let k = 7;
        let w = 20;
        let ukhs = UKHS::new(k, w).unwrap();

        let it = ukhs.iter_sequence(seq);
        let unikmers: Vec<String> = it.map(|(_, x)| x).collect();

        assert_eq!(
            unikmers,
            [
                "ACACCGT", "CCGTAGC", "AGCCTCC", "GCCTCCA", // w-mer 0
                "CCGTAGC", "AGCCTCC", "GCCTCCA", // w-mer 1
                "CCGTAGC", "AGCCTCC", "GCCTCCA", // w-mer 2
                "CCGTAGC", "AGCCTCC", "GCCTCCA", "ATGCGTA", // w-mer 3
                "AGCCTCC", "GCCTCCA", "ATGCGTA" // w-mer 4
            ]
        );

        let it = ukhs.hash_iter_sequence(seq).unwrap();
        let ukhs_hash: Vec<(u64, u64)> = it.collect();

        assert!(
            ukhs_hash.len() > seq.len() - w,
            "iter len {}, should be at least {}",
            ukhs_hash.len(),
            seq.len() - w + 1
        );

        let ukhs_unhash: Vec<String> = ukhs_hash
            .iter()
            .map(|(_, hash)| ukhs.kmer_for_ukhs_hash(*hash).unwrap())
            .collect();

        assert_eq!(unikmers, ukhs_unhash);

        assert_eq!(
            ukhs_hash,
            [
                // w-mer 0
                // (ACACCGTAGCCTCCAGATGC, ACACCGT)
                (0x37137c91412bb512, 0xfbd9591aa929c685),
                // (ACACCGTAGCCTCCAGATGC, CCGTAGC)
                (0x37137c91412bb512, 0x9cd9a1bcb779d6ad),
                // (ACACCGTAGCCTCCAGATGC, AGCCTCC)
                (0x37137c91412bb512, 0x46fa47d28c0ffba5),
                // (ACACCGTAGCCTCCAGATGC, GCCTCCA)
                (0x37137c91412bb512, 0xf482addc6edbc920),
                // w-mer 1
                // (CACCGTAGCCTCCAGATGCG, CCGTAGC)
                (0xf52d9b92474381bf, 0x9cd9a1bcb779d6ad),
                // (CACCGTAGCCTCCAGATGCG, AGCCTCC)
                (0xf52d9b92474381bf, 0x46fa47d28c0ffba5),
                // (CACCGTAGCCTCCAGATGCG, GCCTCCA)
                (0xf52d9b92474381bf, 0xf482addc6edbc920),
                // w-mer 2
                // (ACCGTAGCCTCCAGATGCGT, CCGTAGC)
                (0xdb5854d371a65e15, 0x9cd9a1bcb779d6ad),
                // (ACCGTAGCCTCCAGATGCGT, AGCCTCC)
                (0xdb5854d371a65e15, 0x46fa47d28c0ffba5),
                // (ACCGTAGCCTCCAGATGCGT, GCCTCCA)
                (0xdb5854d371a65e15, 0xf482addc6edbc920),
                // w-mer 3
                // (CCGTAGCCTCCAGATGCGTA, CCGTAGC)
                (0x31020e7531c970e0, 0x9cd9a1bcb779d6ad),
                // (CCGTAGCCTCCAGATGCGTA, AGCCTCC)
                (0x31020e7531c970e0, 0x46fa47d28c0ffba5),
                // (CCGTAGCCTCCAGATGCGTA, GCCTCCA)
                (0x31020e7531c970e0, 0xf482addc6edbc920),
                // (CCGTAGCCTCCAGATGCGTA, ATGCGTA)
                (0x31020e7531c970e0, 0x6903a07436e4ffa1),
                // w-mer 4
                // (CCGTAGCCTCCAGATGCGTA, CCGTAGC)
                // (CGTAGCCTCCAGATGCGTAG, AGCCTCC)
                (0x5a6008385506dbd8, 0x46fa47d28c0ffba5),
                // (CGTAGCCTCCAGATGCGTAG, GCCTCCA)
                (0x5a6008385506dbd8, 0xf482addc6edbc920),
                // (CGTAGCCTCCAGATGCGTAG, TGCGTAG)
                (0x5a6008385506dbd8, 0x6903a07436e4ffa1)
            ]
        );
    }
}
