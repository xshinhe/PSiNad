#include "Kernel_Random.h"


namespace kids {

Kernel_Random::rng_t Kernel_Random::rand_rng;
std::uniform_int_distribution<int> Kernel_Random::rand_uid{0, 1};         ///< catalog distribution
std::uniform_real_distribution<kids_real> Kernel_Random::rand_udd{0, 1};  ///< uniform distribution
std::normal_distribution<kids_real> Kernel_Random::rand_nd{0.0, 1.0};     ///< normal distribution
std::poisson_distribution<int> Kernel_Random::rand_pd{1.0};               ///< possion distribution

};  // namespace kids
