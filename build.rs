use phf_codegen;

use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use nthash::ntf64;

fn main() {
    let path = Path::new(&env::var("OUT_DIR").unwrap()).join("codegen.rs");
    let mut file = BufWriter::new(File::create(&path).unwrap());

    write!(
        &mut file,
        "static UKHS_HASHES: phf::Set<&'static str> =
          "
    )
    .unwrap();

    let entries = include_str!("data/res_7_20_4_0.txt");

    let mut builder = phf_codegen::Set::new();
    for key in entries.split('\n') {
        builder.entry(key);
    }
    builder.build(&mut file).unwrap();

    write!(
        &mut file,
        ";
        
        static UKHS_NTHASHES: phf::Map<u64, &'static str> =
          "
    )
    .unwrap();

    let mut builder = phf_codegen::Map::new();
    for value in entries.split('\n') {
        if value.len() == 7 {
            let key = ntf64(value.as_bytes(), 0, 7);
            builder.entry(key, format!("\"{}\"", value).as_str());
        }
    }
    builder.build(&mut file).unwrap();

    write!(&mut file, ";\n").unwrap();
}