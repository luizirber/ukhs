#include <BooPHF.h>
#include <fstream>

typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
typedef boomphf::mphf<uint64_t, hasher_t> mphf_t;

extern "C" {
  mphf_t* mphf_new(unsigned long long int n, uint64_t *input_range, int num_thread, double gamma) {
    std::vector<uint64_t> vec;
    vec.assign(input_range, input_range + n);
    return new mphf_t(n, vec, num_thread, gamma, false, false);
  }

  uint64_t mphf_lookup(mphf_t* self, uint64_t elem) {
    return self->lookup(elem);
  }

  void mphf_save(mphf_t* self, const char *path, uint64_t size) {
    std::string cpppath(path, size);
    std::ofstream output(cpppath, std::ios::binary);
    self->save(output);
  }

  mphf_t* mphf_load(const char *path, uint64_t size) {
    std::string cpppath(path, size);
    std::ifstream input(cpppath, std::ios::binary);
    std::vector<uint64_t> vec;

    mphf_t* obj = new mphf_t(0, vec, 1, 1, false, false);
    obj->load(input);
    return obj;
  }
}
