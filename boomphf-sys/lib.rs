pub mod boomphf {
    use std::os::raw::{c_float, c_int, c_ulonglong, c_void};

    pub type Mphf = c_void;
    pub type MphfMutPtr = *mut Mphf;

    #[link(name = "bbhashadapter", kind = "static")]
    extern "C" {
        pub fn new_mphf(
            nelem: c_ulonglong,
            kmers: *const c_ulonglong,
            num_thread: c_int,
            gamma: c_float,
        ) -> MphfMutPtr;
    }
}

#[cfg(test)]
mod tests {
    use crate::boomphf::new_mphf;

    #[test]
    fn it_works() {
        unsafe {
            let mphf = new_mphf(3, vec![1, 2, 3].as_ptr(), 1, 0.0);
        }
        assert_eq!(2 + 2, 4);
    }
}
