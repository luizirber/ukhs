use std::path::Path;

use bbhash_sys::boomphf;

pub struct MPHF {
    inner: boomphf::MphfMutPtr,
}

// TODO: this should be OK, because once it's build you can't modified it
// BUT IS IT?!?!
impl Clone for MPHF {
    fn clone(&self) -> Self {
        MPHF { inner: self.inner }
    }
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

// TODO: need to make sure this is true!
unsafe impl Send for MPHF {}
unsafe impl Sync for MPHF {}

#[cfg(test)]
mod tests {
    use boomphf;
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

    #[test]
    fn boomphf_oracle_smoketest() {
        // sample set of objects
        let possible_objects = vec![1, 10, 1000, 23, 457, 856, 845, 124, 912];
        let n = possible_objects.len();

        // generate a minimal perfect hash function of these items
        let boo = boomphf::Mphf::new(1.7, &possible_objects.clone());
        let bb = MPHF::new(possible_objects.clone(), 1, 1.7);

        // Get hash value of all objects
        let mut boo_hashes = Vec::new();
        let mut bb_hashes = Vec::new();
        for v in possible_objects {
            let boo_hash = boo.hash(&v);
            boo_hashes.push(boo_hash);
            let bb_hash = bb.lookup(v).unwrap();
            bb_hashes.push(bb_hash);

            dbg!((boo_hash, bb_hash));
            //assert_eq!(boo_hash, bb_hash);
        }
        boo_hashes.sort();
        bb_hashes.sort();

        // Expected hash output is set of all integers from 0..n
        let expected_hashes: Vec<u64> = (0..n as u64).collect();
        assert!(boo_hashes == expected_hashes);
        assert!(bb_hashes == expected_hashes)
    }

}
