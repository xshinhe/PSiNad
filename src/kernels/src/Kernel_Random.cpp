#include "kids/Kernel_Random.h"

#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"

namespace PROJECT_NS {

Kernel_Random::rng_t                      Kernel_Random::rand_rng;
std::uniform_int_distribution<int>        Kernel_Random::rand_uid{0, 1};     ///< catalog distribution
std::uniform_real_distribution<kids_real> Kernel_Random::rand_udd{0, 1};     ///< uniform distribution
std::normal_distribution<kids_real>       Kernel_Random::rand_nd{0.0, 1.0};  ///< normal distribution
std::poisson_distribution<int>            Kernel_Random::rand_pd{1.0};       ///< possion distribution

const std::string Kernel_Random::getName() { return "Kernel_Random"; }

int Kernel_Random::getType() const { return utils::hash(FUNCTION_NAME); }

int Kernel_Random::rand_catalog(int* res_arr, int N, bool reset, int begin, int end) {
    if (reset) rand_uid.param(uid_range{begin, end});
    for (int i = 0; i < N; ++i) { res_arr[i] = rand_uid(rand_rng); }
    return 0;
}

int Kernel_Random::rand_uniform(kids_real* res_arr, int N, kids_real sigma) {
    for (int i = 0; i < N; ++i) { res_arr[i] = sigma * rand_udd(rand_rng); }
    return 0;
}

int Kernel_Random::rand_gaussian(kids_real* res_arr, int N, kids_real sigma,
                                 kids_real mu) {  // parse mu after sigma !!!
    for (int i = 0; i < N; ++i) { res_arr[i] = mu + sigma * rand_nd(rand_rng); }
    return 0;
}

int Kernel_Random::rand_exponent(kids_real* res_arr, int N) {
    kids_real randu;
    for (int i = 0; i < N; ++i) {
        rand_uniform(&randu);
        res_arr[i] = -std::log(1.0 - randu);
    }
    return 0;
}

int Kernel_Random::rand_poisson(int* res_arr, int N, kids_real lambda) {
    // when lambda > 10, it can be approximated by rand_gaussian(res_arr, N, sqrt(lambda), lambda);
    rand_pd.param(pd_range{lambda});
    for (int i = 0; i < N; ++i) res_arr[i] = rand_pd(rand_rng);
    return 0;
}

int Kernel_Random::rand_simplex(kids_real* res_arr, int N, kids_real constr) {
    for (int i = 1; i < N; ++i) res_arr[i] = rand_udd(rand_rng);
    res_arr[0] = 0;
    std::qsort(res_arr, N, sizeof(*res_arr), [](const void* a, const void* b) {
        kids_real arg1 = *static_cast<const kids_real*>(a);
        kids_real arg2 = *static_cast<const kids_real*>(b);

        if (arg1 < arg2) return -1;
        if (arg1 > arg2) return 1;
        return 0;
    });
    for (int i = 0; i < N - 1; ++i) { res_arr[i] = res_arr[i + 1] - res_arr[i]; }
    res_arr[N - 1] = 1.0f - res_arr[N - 1];

    if (constr != 1.0f)
        for (int i = 1; i < N; ++i) res_arr[i] *= constr;
    return 0;
}

int Kernel_Random::rand_sphere(kids_real* res_arr, int N, kids_real constr) {
    kids_real norm = 0.0f;
    for (int i = 0; i < N; ++i) {
        res_arr[i] = rand_nd(rand_rng);
        norm += res_arr[i] * res_arr[i];
    }
    norm = std::sqrt(norm);
    for (int i = 0; i < N; ++i) { res_arr[i] /= norm; }

    if (constr != 1.0f)
        for (int i = 0; i < N; ++i) res_arr[i] *= constr;
    return 0;
}

void Kernel_Random::setInputDataSet_impl(std::shared_ptr<DataSet>& DS) {
    seed = DS->def_int("random.seed", rng_t::state_size);  //
}

Status& Kernel_Random::initializeKernel_impl(Status& stat) {
    if (count_calc == 0) {
        std::random_device source;

        if (!restart) {
            for (int i = 0; i < rng_t::state_size; ++i) seed[i] = source();
        }
        std::seed_seq sseq(seed, seed + rng_t::state_size);
        rand_rng = rng_t(sseq);
    }
    return stat;
}

};  // namespace PROJECT_NS
