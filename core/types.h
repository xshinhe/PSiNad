/**
 * @file types.h
 * @author xshinhe
 * @version 1.0
 * @date 2019-01-01
 * @brief specific precision of numerical types
 *
 */

// #include "generate/version.h"

#ifndef NUM_TYPES_H
#define NUM_TYPES_H

#include <complex> /* using C++ std::complex */

namespace PROJECT_NS {

using num_int     = int;                   ///< precisions for int type
using num_real    = double;                ///< precisions for real type
using num_complex = std::complex<double>;  ///< precisions for complex type

};  // namespace PROJECT_NS

#endif  // NUM_TYPES_H
