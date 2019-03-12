pub mod boomphf {
    use std::os::raw::{c_float, c_int, c_uchar, c_ulonglong, c_void};

    pub type Mphf = c_void;
    pub type MphfMutPtr = *mut Mphf;

    #[link(name = "bbhashadapter", kind = "static")]
    extern "C" {
        #[link_name = "mphf_new"]
        pub fn new_mphf(
            nelem: c_ulonglong,
            kmers: *const c_ulonglong,
            num_thread: c_int,
            gamma: c_float,
        ) -> MphfMutPtr;

        #[link_name = "mphf_lookup"]
        pub fn lookup(this: MphfMutPtr, element: c_ulonglong) -> c_ulonglong;

        #[link_name = "mphf_save"]
        pub fn save(this: MphfMutPtr, path: *const c_uchar, path_size: c_ulonglong);

        #[link_name = "mphf_load"]
        pub fn load(path: *const c_uchar, path_size: c_ulonglong) -> MphfMutPtr;
    }

}

#[cfg(test)]
mod tests {
    use crate::boomphf::{load, lookup, new_mphf, save};
    use tempfile::NamedTempFile;

    #[test]
    fn it_works() {
        unsafe {
            let mphf = new_mphf(10, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9].as_ptr(), 1, 1.0);

            assert_eq!(lookup(mphf, 9), 8);

            let file = NamedTempFile::new().unwrap();
            let path = file.path().to_str().unwrap();

            save(mphf, path.as_ptr(), path.len() as u64);

            let loaded_mphf = load(path.as_ptr(), path.len() as u64);
            assert_eq!(lookup(loaded_mphf, 9), 8);
        }
    }
}
