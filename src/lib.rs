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
///     let ukhs: Vec<String> = it.collect();
///     assert_eq!(ukhs, ["ACACCGT", "CCGTAGC", "AGCCTCC", "GCCTCCA"]);
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
    type Item = String;

    fn next(&mut self) -> Option<Self::Item> {
        while self.current_idx != self.max_idx {
            let current_kmer =
                str::from_utf8(&self.seq[self.current_idx..self.current_idx + self.k]).unwrap();

            self.current_idx += 1;

            if UKHS_HASHES.contains(current_kmer) {
                return Some(current_kmer.into());
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
///     let ukhs: Vec<u64> = it.collect();
///     assert_eq!(ukhs, [0xfbd9591aa929c685, 0x9cd9a1bcb779d6ad, 0x46fa47d28c0ffba5, 0xf482addc6edbc920]);
///     # Ok(())
///     # }
/// ```
#[derive(Debug)]
pub struct UKHSHashIterator<'a> {
    nthash_iter: NtHashForwardIterator<'a>,
    k: usize,
    w: usize,
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

        let nthash_iter = NtHashForwardIterator::new(seq, k).map_err(SyncFailure::new)?;

        Ok(UKHSHashIterator {
            nthash_iter,
            k,
            w,
            max_idx,
        })
    }
}

impl<'a> Iterator for UKHSHashIterator<'a> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(current_hashed) = self.nthash_iter.next() {
            if let Some(_) = UKHS_NTHASHES.get(&current_hashed) {
                return Some(current_hashed);
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
        let ukhs: Vec<String> = it.collect();

        assert_eq!(ukhs, ["ACACCGT", "CCGTAGC", "AGCCTCC", "GCCTCCA"]);

        let it = UKHSHashIterator::new(seq, 7, 20).unwrap();
        let ukhs_hash: Vec<u64> = it.collect();
        let ukhs_unhash: Vec<String> = ukhs_hash
            .iter()
            .map(|x| UKHS_NTHASHES.get(x).unwrap().to_string())
            .collect();

        assert_eq!(ukhs, ukhs_unhash);
    }
}
