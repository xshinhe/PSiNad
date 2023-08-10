#ifndef Kernel_Random_H
#define Kernel_Random_H

#include <random>

#include "../core/Kernel.h"


namespace PROJECT_NS {
/**
 * @brief Kernel_Random manipulation of random engine and numbers
 */
class Kernel_Random : public Kernel {
   public:
    using rng_t     = std::mt19937;  ///< random number generator type
    using uid_range = std::uniform_int_distribution<int>::param_type;
    using pd_range  = std::poisson_distribution<int>::param_type;

    static rng_t rand_rng;

    static std::uniform_int_distribution<int> rand_uid;        ///< catalog distribution
    static std::uniform_real_distribution<num_real> rand_udd;  ///< uniform distribution
    static std::normal_distribution<num_real> rand_nd;         ///< normal distribution
    static std::poisson_distribution<int> rand_pd;             ///< possion distribution

    inline virtual const std::string name() { return "Kernel_Random"; }

    static inline int rand_catalog(int* res_arr, int N = 1, bool reset = false, int begin = 0, int end = 1) {
        if (reset) rand_uid.param(uid_range{begin, end});
        for (int i = 0; i < N; ++i) { res_arr[i] = rand_uid(rand_rng); }
        return 0;
    }

    static inline int rand_uniform(num_real* res_arr, int N = 1, num_real sigma = 1.0) {
        for (int i = 0; i < N; ++i) { res_arr[i] = sigma * rand_udd(rand_rng); }
        return 0;
    }

    static inline int rand_gaussian(num_real* res_arr, int N = 1, num_real sigma = 1.0,
                                    num_real mu = 0.0) {  // parse mu after sigma !!!
        for (int i = 0; i < N; ++i) { res_arr[i] = mu + sigma * rand_nd(rand_rng); }
        return 0;
    }

    static inline int rand_exponent(num_real* res_arr, int N = 1) {
        num_real randu;
        for (int i = 0; i < N; ++i) {
            rand_uniform(&randu);
            res_arr[i] = -std::log(1.0 - randu);
        }
        return 0;
    }

    static inline int rand_poisson(int* res_arr, int N = 1, num_real lambda = 1.0f) {
        // when lambda > 10, it can be approximated by rand_gaussian(res_arr, N, sqrt(lambda), lambda);
        rand_pd.param(pd_range{lambda});
        for (int i = 0; i < N; ++i) res_arr[i] = rand_pd(rand_rng);
        return 0;
    }

    static inline int rand_simplex(num_real* res_arr, int N = 1, num_real constr = 1.0f) {
        for (int i = 1; i < N; ++i) res_arr[i] = rand_udd(rand_rng);
        res_arr[0] = 0;
        std::qsort(res_arr, N, sizeof(*res_arr), [](const void* a, const void* b) {
            num_real arg1 = *static_cast<const num_real*>(a);
            num_real arg2 = *static_cast<const num_real*>(b);

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

    static inline int rand_sphere(num_real* res_arr, int N = 1, num_real constr = 1.0f) {
        num_real norm = 0.0f;
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

   private:
    int* seed;
    bool restart;

    virtual void init_data_impl(DataSet* S) {
        seed = S->reg<int>("random.seed", rng_t::state_size);  //
    }

    virtual void init_calc_impl(int stat = -1) {
        if (count_calc == 0) {
            std::random_device source;

            if (!restart) {
                for (int i = 0; i < rng_t::state_size; ++i) seed[i] = source();
            }
            std::seed_seq sseq(seed, seed + rng_t::state_size);
            rand_rng = rng_t(sseq);
        }
    }
};

};  // namespace PROJECT_NS


#endif  // Kernel_Random_H
