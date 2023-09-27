/**
 * @file array_macro.h
 * @author xshinhe
 * @version 1.0
 * @date 2019-01-01
 * @brief utils for manipulation of array
 * @details
 *
 * @see https://en.cppreference.com/w/cpp/algorithm/generate
 *      https://en.cppreference.com/w/cpp/header/functional
 *      https://en.cppreference.com/w/cpp/iterator
 *      https://en.cppreference.com/w/cpp/header/random
 */

#ifndef Array_Utils_H
#define Array_Utils_H

/*======================================================
=            Declare the implement of array            =
======================================================*/

/**
 *
 * Declare implement of array. Options: Eigen, xtensor and macro. See more:
 *
 *      @test: performance tests are on $DIR/benckmark/array_benckmark.md
 *
 * @note: all arrays are actually saved in 1-d C-like array style, for better balance
 *      with other language like Fortran.
 *
 * @brief: realize basic interface for array manipulations:
 *
 *  interface to Eigen/xtensor: base on inline functions and templated functions. Here
 *      xtensor(or xtensor-blas) further interface to mkl functions by `flens` project.
 *
 *  interface to macro: base on macro definitions. Easy to read and use, but little low
 *      efficiency.
 *
 *  we recomment interface to Eigen by default.
 */

#define ARRAY_IMPLEMENT_MACRO
#define ARRAY_IMPLEMENT_EIEGN
#define ARRAY_IMPLEMENT_XTENSOR

/*=====  End of Declare the implement of array  ======*/



/*----------  for debugging array information ----------*/

#include "io_utils.h"

#define ARRAY_SHOW(_A, _n1, _n2)                                                     \
    ({                                                                               \
        std::cout << "Show Array <" << #_A << ">\n";                                 \
        int _idxA = 0;                                                               \
        for (int _i = 0; _i < (_n1); ++_i) {                                         \
            for (int _j = 0; _j < (_n2); ++_j) std::cout << FMT(4) << (_A)[_idxA++]; \
            std::cout << std::endl;                                                  \
        }                                                                            \
    })

/*----------  common array information ----------*/
#include <cmath>
template <class T>
bool ARRAY_ISFINITE(T* A, size_t n) {
    bool is_finite = true;
    for (int i = 0; is_finite && i < n; ++i) is_finite = std::isfinite(std::abs(A[i]));
    return is_finite;
}

// #include "array_macro.h"
#include "array_eigen.h"
#include "array_mkl.h"
// #include "array_xtensor.h"

// #ifdef ARRAY_IMPLEMENT_EIEGN
// #elif defined(ARRAY_IMPLEMENT_XTENSOR)
// #elif defined(ARRAY_IMPLEMENT_MACRO)
// #else
// #error "Unknown Array Implement"
// #endif // ARRAY_IMPLEMENT_EIEGN

#endif  // Array_Utils_H
