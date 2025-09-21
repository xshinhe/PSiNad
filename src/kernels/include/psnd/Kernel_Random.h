/**@file        Kernel_Conserve.h
 * @brief       this file provides Kernel_Conserve class enabling energy tracing
 *              and conservation.
 *
 * @author      Xin He
 * @date        2024-03
 * @version     1.0
 * @copyright   GNU Lesser General Public License (LGPL)
 *
 *              Copyright (c) 2024 Xin He, Liu-Group
 *
 *  This software is a product of academic research conducted by Professor Liu's
 *  Group at the College of Chemistry and Molecular Engineering, Peking University.
 *  All rights are reserved by Peking University.
 *  You should have received a copy of the GNU Lesser General Public License along
 *  with this software. If not, see <https://www.gnu.org/licenses/lgpl-3.0.en.html>
 **********************************************************************************
 * @par revision:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-06     <td> Updated.
 * </table>
 **********************************************************************************
 */

#ifndef Kernel_Random_H
#define Kernel_Random_H

#include <random>

#include "psnd/Kernel.h"

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
    static std::uniform_real_distribution<psnd_real> rand_udd;  ///< uniform distribution
    static std::normal_distribution<psnd_real>       rand_nd;   ///< normal distribution
    static std::poisson_distribution<int>            rand_pd;   ///< possion distribution

    virtual const std::string getName();

    virtual int getType() const;

    static void setSeed(int* seed);

    static int rand_catalog(int* res_arr, int N = 1, bool reset = false, int begin = 0, int end = 1);

    static int rand_uniform(psnd_real* res_arr, int N = 1, psnd_real sigma = 1.0);

    static int rand_gaussian(psnd_real* res_arr, int N = 1, psnd_real sigma = 1.0, psnd_real mu = 0.0);

    static int rand_exponent(psnd_real* res_arr, int N = 1);

    static int rand_poisson(int* res_arr, int N = 1, psnd_real lambda = 1.0f);

    static int rand_simplex(psnd_real* res_arr, int N = 1, psnd_real constr = 1.0f);

    static int rand_sphere(psnd_real* res_arr, int N = 1, psnd_real constr = 1.0f);

   private:
    span<psnd_int> seed;

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& initializeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Random_H
