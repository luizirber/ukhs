#![allow(clippy::unreadable_literal)]

pub mod errors;

use std::str;

use failure::{Error, SyncFailure};
use nthash::NtHashForwardIterator;
use phf;

use crate::errors::UKHSError;

include!(concat!(env!("OUT_DIR"), "/codegen.rs"));

/// An iterator for finding universal hitting k-mers in a sequence.
///
/// ```
///     # use failure::Error;
///     use ukhs::UKHSIterator;
///
///     # fn main() -> Result<(), Error> {
///     let seq = b"ACACCGTAGCCTCCAGATGC";
///     let it = UKHSIterator::new(seq, 7, 20)?;
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
#[derive(Debug)]
pub struct UKHSIterator<'a> {
    seq: &'a [u8],
    k: usize,
    w: usize,
    current_idx: usize,
    max_idx: usize,
}

impl<'a> UKHSIterator<'a> {
    /// Creates a new UKHSIterator with internal state properly initialized.
    pub fn new(seq: &'a [u8], k: usize, w: usize) -> Result<UKHSIterator<'a>, Error> {
        if k > seq.len() {
            return Err(UKHSError::KSizeOutOfRange {
                ksize: k,
                sequence: String::from_utf8(seq.to_vec()).unwrap(),
            }
            .into());
        }

        if k > w {
            return Err(UKHSError::KSizeOutOfWRange { ksize: k, wsize: w }.into());
        }

        if w > seq.len() {
            return Err(UKHSError::WSizeOutOfRange {
                wsize: w,
                sequence: String::from_utf8(seq.to_vec()).unwrap(),
            }
            .into());
        }

        let current_idx = 0;
        let max_idx = seq.len() - k + 1;

        Ok(UKHSIterator {
            seq,
            k,
            w,
            current_idx,
            max_idx,
        })
    }
}

impl<'a> Iterator for UKHSIterator<'a> {
    type Item = (String, String);

    fn next(&mut self) -> Option<Self::Item> {
        while self.current_idx != self.max_idx {
            let current_kmer =
                str::from_utf8(&self.seq[self.current_idx..self.current_idx + self.k]).unwrap();

            let wmer_start = self.current_idx / self.w;
            let current_wmer = str::from_utf8(&self.seq[wmer_start..wmer_start + self.w]).unwrap();

            self.current_idx += 1;

            if UKHS_HASHES.contains(current_kmer) {
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
///     use ukhs::UKHSHashIterator;
///
///     # fn main() -> Result<(), Error> {
///     let seq = b"ACACCGTAGCCTCCAGATGC";
///     let it = UKHSHashIterator::new(seq, 7, 20)?;
///     let ukhs: Vec<(u64, u64)> = it.collect();
///     assert_eq!(ukhs,
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
///     use ukhs::UKHSHashIterator;
///
///     # fn main() -> Result<(), Error> {
///     let seq = b"ACACCGTAGCCTCCAGATGC";
///     let it = UKHSHashIterator::new(seq, 7, 20)?;
///     let ukhs: Vec<u64> = it.map(|(_, hash)| hash).collect();
///     assert_eq!(ukhs, [0xfbd9591aa929c685, 0x9cd9a1bcb779d6ad, 0x46fa47d28c0ffba5, 0xf482addc6edbc920]);
///     # Ok(())
///     # }
/// ```
#[derive(Debug)]
pub struct UKHSHashIterator<'a> {
    nthash_k_iter: NtHashForwardIterator<'a>,
    nthash_w_iter: NtHashForwardIterator<'a>,
    k: usize,
    w: usize,
    current_w_hash: u64,
    current_idx: usize,
    max_idx: usize,
}

impl<'a> UKHSHashIterator<'a> {
    /// Creates a new UKHSHashIterator with internal state properly initialized.
    pub fn new(seq: &'a [u8], k: usize, w: usize) -> Result<UKHSHashIterator<'a>, Error> {
        if k > seq.len() {
            return Err(UKHSError::KSizeOutOfRange {
                ksize: k,
                sequence: String::from_utf8(seq.to_vec()).unwrap(),
            }
            .into());
        }

        if k > w {
            return Err(UKHSError::KSizeOutOfWRange { ksize: k, wsize: w }.into());
        }

        if w > seq.len() {
            return Err(UKHSError::WSizeOutOfRange {
                wsize: w,
                sequence: String::from_utf8(seq.to_vec()).unwrap(),
            }
            .into());
        }

        let max_idx = seq.len() - k + 1;
        let current_idx = 0;

        let nthash_k_iter = NtHashForwardIterator::new(seq, k).map_err(SyncFailure::new)?;
        let mut nthash_w_iter = NtHashForwardIterator::new(seq, w).map_err(SyncFailure::new)?;

        let current_w_hash = nthash_w_iter.next().unwrap();

        Ok(UKHSHashIterator {
            nthash_k_iter,
            nthash_w_iter,
            k,
            w,
            current_w_hash,
            current_idx,
            max_idx,
        })
    }

    pub fn kmer_for_ukhs_hash(hash: u64) -> Option<String> {
        if let Some(kmer) = UKHS_NTHASHES.get(&hash) {
            Some((*kmer).to_string())
        } else {
            None
        }
    }
}

impl<'a> Iterator for UKHSHashIterator<'a> {
    type Item = (u64, u64);

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(k_hash) = self.nthash_k_iter.next() {
            self.current_idx += 1;
            if self.current_idx % self.w == 0 {
                self.current_w_hash = self.nthash_w_iter.next().unwrap();
            }

            if UKHS_NTHASHES.get(&k_hash).is_some() {
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
        let it = UKHSIterator::new(seq, 7, 20).unwrap();
        let mut ukhs: Vec<String> = it.map(|(_, x)| x).collect();
        ukhs.sort_unstable();

        assert_eq!(ukhs, ["ACACCGT", "AGCCTCC", "CCGTAGC", "GCCTCCA"]);

        let it = UKHSHashIterator::new(seq, 7, 20).unwrap();
        let ukhs_hash: Vec<(u64, u64)> = it.collect();
        let mut ukhs_unhash: Vec<String> = ukhs_hash
            .iter()
            .map(|(_, hash)| UKHSHashIterator::kmer_for_ukhs_hash(*hash).unwrap())
            .collect();
        ukhs_unhash.sort_unstable();

        assert_eq!(ukhs, ukhs_unhash);
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
