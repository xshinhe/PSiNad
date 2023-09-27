#include "random_utils.h"

#include <algorithm>   //=> std::generate, as syntactic sugar
#include <functional>  //=> std::ref, std::bind (functional coding)
#include <iterator>    //=> std::begin, std::end
#include <random>

/*=================================================================
=            specify random generators & distributions            =
=================================================================*/

static bool if_has_initialized = false;

rng_t rand_rng;

rng_t::result_type rng_seeds[rng_t::state_size];

rng_t SeededEngine(rng_t::result_type* init_seeds) {
    if (init_seeds == nullptr) {
        // default: generated seeds by std::random_device
        std::random_device source;
        std::generate(std::begin(rng_seeds), std::end(rng_seeds), std::ref(source));
    } else {
        for (int i = 0; i < rng_t::state_size; ++i) rng_seeds[i] = init_seeds[i];
    }
    if_has_initialized = true;
    std::seed_seq sseq(std::begin(rng_seeds), std::end(rng_seeds));
    return rng_t(sseq);
}

/*----------  random distributions (localy defined) -------------*/

using uid_range = std::uniform_int_distribution<int>::param_type;
using pd_range  = std::poisson_distribution<int>::param_type;

static std::uniform_int_distribution<int> rand_uid{0, 1};      ///< catalog distribution
static std::uniform_real_distribution<double> rand_udd(0, 1);  ///< uniform distribution
static std::normal_distribution<double> rand_nd(0.0, 1.0);     ///< normal distribution
static std::poisson_distribution<int> rand_pd(1.0);


/*=====  End of specify random generators & distributions  ======*/


/*===============================================================
=            functions call random number for arrays            =
===============================================================*/


int rand_catalog(int* res_arr, const int& N, const bool& reset, const int& begin, const int& end) {
    if (reset) rand_uid.param(uid_range{begin, end});
    for (int i = 0; i < N; ++i) { res_arr[i] = rand_uid(rand_rng); }
    return 0;
}

int rand_uniform(double* res_arr, const int& N, const double& sigma) {
    for (int i = 0; i < N; ++i) { res_arr[i] = sigma * rand_udd(rand_rng); }
    return 0;
}

int rand_gaussian(double* res_arr, const int& N, const double& sigma,
                  const double& mu) {  // parse mu after sigma !!!
    for (int i = 0; i < N; ++i) { res_arr[i] = mu + sigma * rand_nd(rand_rng); }
    return 0;
}

int rand_exponent(double* res_arr, const int& N) {
    double randu;
    for (int i = 0; i < N; ++i) {
        rand_uniform(&randu);
        res_arr[i] = -std::log(1.0 - randu);
    }
    return 0;
}

int rand_poisson(int* res_arr, const int& N, const double& lambda) {
    // when lambda > 10, it can be approximated by rand_gaussian(res_arr, N, sqrt(lambda), lambda);
    rand_pd.param(pd_range{lambda});
    for (int i = 0; i < N; ++i) res_arr[i] = rand_pd(rand_rng);

    /** @deprecated following shows "knuth algorithm", but it's slow
     *  double p, L, randu;
     *  p = 1.0;
     *  L = std::exp(-lambda);
     *  for (int i = 0; i < N; ++i) {
     *      res_arr[i] = 0;
     *      while (p > L) {
     *          rand_uniform(&randu);
     *          p *= randu;
     *          res_arr[i] += 1;
     *      }
     *      res_arr[i] -= 1;
     *  }
     */
    return 0;
}

int rand_simplex(double* res_arr, const int& N, const double& constr) {
    for (int i = 1; i < N; ++i) res_arr[i] = rand_udd(rand_rng);
    res_arr[0] = 0;
    std::qsort(res_arr, N, sizeof(*res_arr), [](const void* a, const void* b) {
        double arg1 = *static_cast<const double*>(a);
        double arg2 = *static_cast<const double*>(b);

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

int rand_sphere(double* res_arr, const int& N, const double& constr) {
    double norm = 0.0f;
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
