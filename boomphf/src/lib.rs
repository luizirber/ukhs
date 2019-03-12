use std::path::Path;

use boomphf_sys::boomphf;

pub struct MPHF {
    inner: boomphf::MphfMutPtr,
}

impl MPHF {
    pub fn new(elements: Vec<u64>, num_threads: i32, gamma: f64) -> MPHF {
        let inner = unsafe {
            boomphf::new_mphf(elements.len() as u64, elements.as_ptr(), num_threads, gamma)
        };

        MPHF { inner }
    }

    pub fn lookup(&self, elem: u64) -> Option<u64> {
        let pos = unsafe { boomphf::lookup(self.inner, elem) };
        if pos == u64::max_value() {
            None
        } else {
            Some(pos)
        }
    }

    pub fn save<P: AsRef<Path>>(&self, path: P) {
        let filename = path.as_ref().to_str().unwrap();

        unsafe {
            boomphf::save(self.inner, filename.as_ptr(), filename.len() as u64);
        };
    }

    pub fn load<P: AsRef<Path>>(path: P) -> MPHF {
        let filename = path.as_ref().to_str().unwrap();

        let inner = unsafe { boomphf::load(filename.as_ptr(), filename.len() as u64) };
        MPHF { inner }
    }
}

#[cfg(test)]
mod tests {
    use tempfile::NamedTempFile;

    use crate::MPHF;

    #[test]
    fn it_works() {
        let mphf = MPHF::new(vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 1, 1.0);

        assert_eq!(mphf.lookup(9).unwrap(), 8);
        assert_eq!(mphf.lookup(25), None);

        let file = NamedTempFile::new().unwrap();
        let path = file.path();

        mphf.save(path);

        let loaded_mphf = MPHF::load(path);
        assert_eq!(loaded_mphf.lookup(9).unwrap(), 8);
        assert_eq!(loaded_mphf.lookup(25), None);
    }
}
