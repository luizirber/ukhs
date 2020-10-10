use thiserror::Error;

#[derive(Debug, Error)]
pub enum UKHSError {
    #[error("K size {ksize} is out of range for sequence {sequence:?}")]
    KSizeOutOfRange { ksize: usize, sequence: String },

    #[error("K size {ksize} is out of range for window range {wsize}")]
    KSizeOutOfWRange { ksize: usize, wsize: usize },

    #[error("Window size {wsize} is out of range for sequence {sequence:?}")]
    WSizeOutOfRange { wsize: usize, sequence: String },

    #[error(transparent)]
    NtHashError(#[from] nthash::result::Error),
}
