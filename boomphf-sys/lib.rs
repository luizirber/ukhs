pub mod boomphf {
    use std::os::raw::c_void;

    pub type Mphf = c_void;
    pub type MphfMutPtr = *mut Mphf;

    #[link(name = "bbhashadapter", kind = "static")]
    extern "C" {
        #[link_name = "mphf_new"]
        pub fn new_mphf(nelem: u64, kmers: *const u64, num_thread: i32, gamma: f64) -> MphfMutPtr;

        #[link_name = "mphf_lookup"]
        pub fn lookup(this: MphfMutPtr, element: u64) -> u64;

        #[link_name = "mphf_save"]
        pub fn save(this: MphfMutPtr, path: *const u8, path_size: u64);

        #[link_name = "mphf_load"]
        pub fn load(path: *const u8, path_size: u64) -> MphfMutPtr;
    }

}

#[cfg(test)]
mod tests {
    use tempfile::NamedTempFile;

    use crate::boomphf::{load, lookup, new_mphf, save};

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
