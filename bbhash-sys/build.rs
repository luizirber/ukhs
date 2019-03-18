fn main() {
    let mut c = cc::Build::new();
    c.cpp(true);
    c.warnings(false);
    c.flag("-fno-rtti");
    c.flag("-std=c++11");
    //c.flag("-fno-exceptions");
    c.include("BBHash");
    c.file("ffi.cpp");
    c.compile("libbbhashadapter.a");
}
