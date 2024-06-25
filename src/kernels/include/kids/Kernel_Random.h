#ifndef Kernel_Random_H
#define Kernel_Random_H

#include <random>

#include "kids/Kernel.h"

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

    static std::uniform_int_distribution<int>        rand_uid;  ///< catalog distribution
    static std::uniform_real_distribution<kids_real> rand_udd;  ///< uniform distribution
    static std::normal_distribution<kids_real>       rand_nd;   ///< normal distribution
    static std::poisson_distribution<int>            rand_pd;   ///< possion distribution

    virtual const std::string getName();

    virtual int getType() const;

    static int rand_catalog(int* res_arr, int N = 1, bool reset = false, int begin = 0, int end = 1);

    static int rand_uniform(kids_real* res_arr, int N = 1, kids_real sigma = 1.0);

    static int rand_gaussian(kids_real* res_arr, int N = 1, kids_real sigma = 1.0, kids_real mu = 0.0);

    static int rand_exponent(kids_real* res_arr, int N = 1);

    static int rand_poisson(int* res_arr, int N = 1, kids_real lambda = 1.0f);

    static int rand_simplex(kids_real* res_arr, int N = 1, kids_real constr = 1.0f);

    static int rand_sphere(kids_real* res_arr, int N = 1, kids_real constr = 1.0f);

   private:
    int* seed;
    bool restart;

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& initializeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Random_H
