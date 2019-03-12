#include <BooPHF.h>

typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
typedef boomphf::mphf<uint64_t, hasher_t> mphf_t;

extern "C" {
  mphf_t* new_mphf(unsigned long long int n, uint64_t *input_range, int num_thread, double gamma) {
    std::vector<uint64_t> vec;
    vec.assign(input_range, input_range + n);
    return new mphf_t(n, vec, num_thread, gamma);
  }
}
