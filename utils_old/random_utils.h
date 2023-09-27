/**
 * @file random_utils.h
 * @author xshinhe
 * @version 1.0
 * @date 2019-01-01
 * @brief generate random numbers
 * @details
 *
 *      https://en.cppreference.com/w/cpp/header/random
 */

#ifndef RANDOM_UTILS_H
#define RANDOM_UTILS_H

#include <random>

/*=================================================================
=            specify random generators & distributions            =
=================================================================*/


typedef std::mt19937 rng_t;
// typedef std::mt19937_64 rng_t;

/**
 * @brief global random engine (for mpi process, each sub-process has one)
 */
extern rng_t rand_rng;  ///< @note: defined in random_utils.cpp

/**
 * @brief rng_seeds can be accessible
 */
extern rng_t::result_type rng_seeds[rng_t::state_size];  ///< @note: defined in random_utils.cpp

/**
 * @brief rng_seeds can be resettable
 * @details
 *      1) if you first run a new simulation, you should call with:
 *      ```
 *          rand_rng = SeededEngine(nullptr);
 *      ```
 *      2) if you want to recover a result with known seeds, you should first initialize rng_seeds,
 *      and pass it into the functon,
 *      ```
 *          rand_rng = SeededEngine(rng_seeds);
 *      ```
 *
 * @param init_seeds [description]
 * @return [description]
 */
rng_t SeededEngine(rng_t::result_type* init_seeds);  //

/* (local) mersenne_twister_engine */


/*=====  End of specify random generators & distributions  ======*/


/*===============================================================
=            functions call random number for arrays            =
===============================================================*/


int rand_catalog(int* res_arr, const int& N = 1, const bool& reset = false, const int& begin = 0, const int& end = 1);

int rand_uniform(double* res_arr, const int& N = 1, const double& sigma = 1.0);

int rand_gaussian(double* res_arr, const int& N = 1, const double& sigma = 1.0, const double& mu = 0.0);

int rand_exponent(double* res_arr, const int& N = 1);

int rand_poisson(int* res_arr, const int& N = 1, const double& lambda = 1.0f);

int rand_poisson2(int* res_arr, const int& N = 1, const double& lambda = 1.0f);

int rand_simplex(double* res_arr, const int& N = 1, const double& constr = 1.0f);

int rand_sphere(double* res_arr, const int& N = 1, const double& constr = 1.0f);

#endif  // RANDOM_UTILS_H
