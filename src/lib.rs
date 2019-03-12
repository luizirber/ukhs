#![allow(clippy::unreadable_literal)]

pub mod errors;

use std::collections::HashMap;
use std::str;

use boomphf::MPHF;
use failure::{Error, SyncFailure};
use lazy_static::lazy_static;
use nthash::{ntf64, NtHashForwardIterator};

use crate::errors::UKHSError;

lazy_static! {
    static ref UKHS_HASHES: HashMap<(usize, usize), &'static str> = {
        let mut map = HashMap::new();
        map.insert((7, 20), include_str!("../data/res_7_20_4_0.txt"));
        // TODO: add others
        map
    };
}

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
            return Err(UKHSError::KSizeOutOfWRange { ksize: k, wsize: w }.into());
        }

        let entries = UKHS_HASHES[&(k, w)];
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

    /// Creates a new UKHSIterator with internal state properly initialized.
    pub fn iter_sequence(&'a self, seq: &'a [u8]) -> Result<UKHSIterator<'a>, Error> {
        if self.k > seq.len() {
            return Err(UKHSError::KSizeOutOfRange {
                ksize: self.k,
                sequence: String::from_utf8(seq.to_vec()).unwrap(),
            }
            .into());
        }

        if self.w > seq.len() {
            return Err(UKHSError::WSizeOutOfRange {
                wsize: self.w,
                sequence: String::from_utf8(seq.to_vec()).unwrap(),
            }
            .into());
        }

        let current_idx = 0;
        let max_idx = seq.len() - self.k + 1;

        Ok(UKHSIterator {
            seq,
            ukhs: self,
            current_idx,
            max_idx,
        })
    }

    /// Creates a new UKHSHashIterator with internal state properly initialized.
    pub fn hash_iter_sequence(&'a self, seq: &'a [u8]) -> Result<UKHSHashIterator<'a>, Error> {
        if self.k > seq.len() {
            return Err(UKHSError::KSizeOutOfRange {
                ksize: self.k,
                sequence: String::from_utf8(seq.to_vec()).unwrap(),
            }
            .into());
        }

        if self.w > seq.len() {
            return Err(UKHSError::WSizeOutOfRange {
                wsize: self.w,
                sequence: String::from_utf8(seq.to_vec()).unwrap(),
            }
            .into());
        }

        let max_idx = seq.len() - self.k + 1;
        let current_idx = 0;

        let nthash_k_iter = NtHashForwardIterator::new(seq, self.k).map_err(SyncFailure::new)?;
        let mut nthash_w_iter =
            NtHashForwardIterator::new(seq, self.w).map_err(SyncFailure::new)?;

        let current_w_hash = nthash_w_iter.next().unwrap();

        Ok(UKHSHashIterator {
            nthash_k_iter,
            nthash_w_iter,
            ukhs: self,
            current_w_hash,
            current_idx,
            max_idx,
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
///     # use failure::Error;
///     use ukhs::UKHS;
///
///     # fn main() -> Result<(), Error> {
///     let seq = b"ACACCGTAGCCTCCAGATGC";
///     let ukhs = UKHS::new(7, 20)?;
///
///     let it = ukhs.iter_sequence(seq)?;
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
    current_idx: usize,
    max_idx: usize,
}

impl<'a> Iterator for UKHSIterator<'a> {
    type Item = (String, String);

    fn next(&mut self) -> Option<Self::Item> {
        while self.current_idx != self.max_idx {
            let current_kmer =
                str::from_utf8(&self.seq[self.current_idx..self.current_idx + self.ukhs.k])
                    .unwrap();

            let wmer_start = self.current_idx / self.ukhs.w;
            let current_wmer =
                str::from_utf8(&self.seq[wmer_start..wmer_start + self.ukhs.w]).unwrap();

            self.current_idx += 1;

            if self.ukhs.contains_kmer(current_kmer) {
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
///     # use failure::Error;
///     use ukhs::UKHS;
///
///     # fn main() -> Result<(), Error> {
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
///     # use failure::Error;
///     use ukhs::UKHS;
///
///     # fn main() -> Result<(), Error> {
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
    nthash_k_iter: NtHashForwardIterator<'a>,
    nthash_w_iter: NtHashForwardIterator<'a>,
    current_w_hash: u64,
    current_idx: usize,
    max_idx: usize,
}

impl<'a> UKHSHashIterator<'a> {}

impl<'a> Iterator for UKHSHashIterator<'a> {
    type Item = (u64, u64);

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(k_hash) = self.nthash_k_iter.next() {
            self.current_idx += 1;
            if self.current_idx % self.ukhs.w == 0 {
                self.current_w_hash = self.nthash_w_iter.next().unwrap();
            }

            if self.ukhs.contains(k_hash) {
                return Some((self.current_w_hash, k_hash));
            };
        }
        None
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.max_idx, Some(self.max_idx))
    }
}

impl<'a> ExactSizeIterator for UKHSHashIterator<'a> {}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn basic_check() {
        let seq = b"ACACCGTAGCCTCCAGATGC";
        let ukhs = UKHS::new(7, 20).unwrap();

        let it = ukhs.iter_sequence(seq).unwrap();
        let mut unikmers: Vec<String> = it.map(|(_, x)| x).collect();
        unikmers.sort_unstable();

        assert_eq!(unikmers, ["ACACCGT", "AGCCTCC", "CCGTAGC", "GCCTCCA"]);

        let it = ukhs.hash_iter_sequence(seq).unwrap();
        let ukhs_hash: Vec<(u64, u64)> = it.collect();
        let mut ukhs_unhash: Vec<String> = ukhs_hash
            .iter()
            .map(|(_, hash)| ukhs.kmer_for_ukhs_hash(*hash).unwrap())
            .collect();
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
}
