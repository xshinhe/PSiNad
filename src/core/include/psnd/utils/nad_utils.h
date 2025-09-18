/**
 * @file nad_utils.h
 * @author xshinhe
 * @version 1.0
 * @date 2019-01-01
 * @brief utils for non-adiabatic dynamics
 */

#ifndef NAD_UTILS_H
#define NAD_UTILS_H

#include <cmath>
#include <map>

#include "types.h"

#ifndef GAMMA_WIGNER
#define GAMMA_WIGNER(_F) ((sqrt((double) (_F) + 1) - 1) / (double) (_F))
#endif  //  GAMMA_WIGNER


namespace representation {
enum _enum { diabatic, adiabatic, onthefly, force, density };
const std::map<std::string, _enum> _dict = {
    {"#dia", diabatic},
    {"#adia", adiabatic},
    {"#otf", onthefly},
    {"#den", density},
};
};  // namespace representation


int eac_mvc(num_complex* eac, num_real* mvc, int fdim);

int mvc_eac(num_real* mvc, num_complex* eac, int fdim);

int rho_eac(num_complex* rho, num_complex* eac, int fdim);

int samp_mvc_focus(num_real* mvc, int fdim);

int samp_mvc_gauss(num_real* mvc, num_real variance, int fdim);

int samp_mvc_sphere(num_real* mvc, num_real Rc2, int fdim);

int solve_transform(                                                 ///< solve transform problem
    num_complex* H, num_complex* dH, num_complex* ddH,               // backups
    num_complex* S, num_real* L, num_complex* dL, num_complex* ddL,  // eigenproblem in representation::adiabatic
    num_real* T, num_real* E, num_real* dE, num_real* ddE,           // eigenproblem in representation::diabatic
    num_real* V, num_real* dV, num_real* ddV,                        // original representation::diabatic
    num_real* nr, num_real* np, num_real* nm,                        // phase space
    int rep_type, int level,                                         // flags
    int rdim, int fdim,                                              // shapes
    num_real* workr, num_complex* workc,                             // temperoray array
    bool refered);

/*
    [ref]: J Phys Chem Lett. 2020 Oct 1;11(19):8283-8291. doi: 10.1021/acs.jpclett.0c02533.
    this transform solves eigen solution for (phase corrected) effective Hamiltonian.
*/
int solve_transform_correctphase(                                    ///< solve transform problem with pc
    num_complex* H, num_complex* rho,                                // additional variables
    num_complex* S, num_real* L, num_complex* dL, num_complex* ddL,  // eigenproblem in representation::adiabatic
    num_real* T, num_real* E, num_real* dE, num_real* ddE,           // eigenproblem in representation::diabatic
    num_real* V, num_real* dV, num_real* ddV,                        // original representation::diabatic
    num_real* nr, num_real* np, num_real* nm,                        // phase space
    int rep_type, int level,                                         // flags
    int rdim, int fdim,                                              // shapes
    num_real* workr, num_complex* workc                              // temperoray array
);

int solve_Ut(num_complex* U,                      ///< time-evolve propagator
             num_complex* S, num_real* L,         // eigen solution in representation::adiabatic
             num_real* T, num_real* E,            // eigen solution in representation::diabatic
             num_real dtime,                      // time step
             int rep_type,                        // represation
             int rdim, int fdim,                  // shapes
             num_real* workr, num_complex* workc  // temperoray array
);

int update_eac(num_complex* eac, num_complex* U, int rdim, int fdim, num_real* workr, num_complex* workc);


// !> evolution of gamma matrix during dt
int update_rho(num_complex* rho, num_complex* U, int rdim, int fdim, num_real* workr, num_complex* workc);

int update_drho(num_complex* drho, num_complex* rho, num_complex* U, int rdim, int fdim, num_real* workr,
                num_complex* workc);

#endif  // NAD_UTILS_H
