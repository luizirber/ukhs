use failure::Fail;

#[derive(Debug, Fail)]
pub enum UKHSError {
    #[fail(display = "K size {} is out of range for sequence {}", ksize, sequence)]
    KSizeOutOfRange { ksize: usize, sequence: String },

    #[fail(
        display = "K size {} is out of range for window range {}",
        ksize, wsize
    )]
    KSizeOutOfWRange { ksize: usize, wsize: usize },

    #[fail(
        display = "Window size {} is out of range for sequence {}",
        wsize, sequence
    )]
    WSizeOutOfRange { wsize: usize, sequence: String },
}
