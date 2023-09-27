#include "thermostat.h"

#include "../model.h"

Thermostat::Thermostat(const Param& iparm) : Model(iparm) {
    std::string thermo_flag = Param_GetT(std::string, parm, "thermo_flag", "#lang");
    thermo_type             = ThermostatPolicy::_dict.at(thermo_flag);

    Param_GetV(gammal, iparm, 0.1);  // friction coefficeint

    // parse temperature
    num_real temp = Param_GetQ(phys::temperature_d, parm, "temp", 300.0f);
    beta          = Param_GetT(num_real, parm, "beta", -1.0f);
    if (beta < 0) beta = 1 / (phys::au::k * temp);

    Param_GetV(nthstp, iparm, 0);   // each thermo step (set it zero to be NVE)
    Param_GetV(nchain, iparm, 1);   // chain length of nhc
    Param_GetV(nrespa, iparm, 10);  // for mult-timescale
}


int Thermostat::init_alloc(const int& N) {
    if (thermo_type == ThermostatPolicy::NHC) {  // allow you to re-allocate it
        ndofs = N;
        ALLOCATE_PTR_TO_VECTOR(nhc_r, nchain * ndofs);
        ALLOCATE_PTR_TO_VECTOR(nhc_p, nchain * ndofs);
        ALLOCATE_PTR_TO_VECTOR(nhc_Q, nchain);
        // optimal gammal = 20*dt
        for (int j = 0; j < nchain * ndofs; ++j) nhc_r[j] = 0.0f, nhc_p[j] = 0.0f;
        for (int j = 0; j < nchain; ++j) nhc_Q[j] = gammal * gammal / beta;
        nhc_Q[0] *= ndofs;
        has_allocated = true;
    }
    return 0;
}

Thermostat::~Thermostat() {
    if (thermo_type == ThermostatPolicy::NHC && has_allocated) { delete[] nhc_r, delete[] nhc_p, delete[] nhc_Q; }
}


bool Thermostat::dothermo(const int& istep) { return (nthstp > 0) && (istep % nthstp == 0); }

int Thermostat::evolve(num_real* nr, num_real* np, num_real* nm, const num_real& dt, const int& N, const int& start,
                       num_real gammal_in) {
    if (gammal_in < 0) gammal_in = gammal;
    num_real c1  = exp(-gammal_in * dt);
    num_real c2p = sqrt(1.0 - c1 * c1);
    num_real randu;
    switch (thermo_type) {
        case ThermostatPolicy::None:
            break;
        case ThermostatPolicy::Langevin:
            for (int j = 0; j < N; ++j) {
                rand_gaussian(&randu);
                np[j] = c1 * np[j] + c2p * sqrt(nm[j] / beta) * randu;
            }
            break;
        case ThermostatPolicy::Andersen:
            for (int j = 0; j < N; ++j) {
                rand_uniform(&randu);
                if (randu < 1 - c1) {
                    rand_gaussian(&randu);
                    np[j] = sqrt(nm[j] / beta) * randu;
                }
            }
            break;
        case ThermostatPolicy::NHC:
            run_NHC(nr, np, nm, dt, N, start, gammal_in);
            break;
        default:
            LOG(FATAL);
    }
    return 0;
}

int Thermostat::run_NHC(num_real* nr, num_real* np, num_real* nm, const num_real& dt, const int& N, const int& start,
                        num_real gammal_in) {
    LOG_IF(FATAL, !has_allocated) << "need auxilirary space";
    LOG_IF(FATAL, start + N > ndofs) << "exceed auxilirary space";

    num_real smalldt  = dt / nrespa;
    num_real hsmalldt = 0.5f * smalldt;
    num_real qsmalldt = 0.25f * smalldt;
    num_real Te       = 1.0f / (phys::au::k * beta);
    num_real* rx      = nhc_r + start;  // thermostat with (start, start + N) DOFs
    num_real* px      = nhc_p + start;  // thermostat with (start, start + N) DOFs
    for (int i = 0, ix = 0; i < N; ++i, rx += nchain, px += nchain) {
        // sweep nchain
        for (int k = 0; k < nrespa; ++k) {
            for (int s = 0; s < size_sy; ++s) {
                // update nhc_p from nchain tail to head
                if (nchain > 1)
                    px[nchain - 1] += (px[nchain - 2] * px[nchain - 2] / nhc_Q[nchain - 2] - Te) * wgt_sy[s] * hsmalldt;
                for (int h = nchain - 2; h > 0; --h) {
                    px[h] *= exp(-px[h + 1] / nhc_Q[h + 1] * wgt_sy[s] * qsmalldt);
                    px[h] += (px[h - 1] * px[h - 1] / nhc_Q[h - 1] - Te) * wgt_sy[s] * hsmalldt;
                    px[h] *= exp(-px[h + 1] / nhc_Q[h + 1] * wgt_sy[s] * qsmalldt);
                }
                if (nchain > 1) px[0] -= px[1] / nhc_Q[1] * wgt_sy[s] * qsmalldt;
                px[0] += (np[i] * np[i] / nm[i] - Te) * wgt_sy[s] * hsmalldt;
                if (nchain > 1) px[0] -= px[1] / nhc_Q[1] * wgt_sy[s] * qsmalldt;

                // update nhc_x and np
                for (int h = 0; h < nchain; ++h) rx[h] += px[h] / nhc_Q[h] * wgt_sy[s] * smalldt;
                np[i] *= exp(-px[0] / nhc_Q[0] * wgt_sy[s] * smalldt);

                // update nhc_p from head to tail
                if (nchain > 1) px[0] -= px[1] / nhc_Q[1] * wgt_sy[s] * qsmalldt;
                px[0] += (np[i] * np[i] / nm[i] - Te) * wgt_sy[s] * hsmalldt;
                if (nchain > 1) px[0] -= px[1] / nhc_Q[1] * wgt_sy[s] * qsmalldt;
                for (int h = 1; h < nchain - 1; ++h) {
                    px[h] *= exp(-px[h + 1] / nhc_Q[h + 1] * wgt_sy[s] * qsmalldt);
                    px[h] += (px[h - 1] * px[h - 1] / nhc_Q[h - 1] - Te) * wgt_sy[s] * hsmalldt;
                    px[h] *= exp(-px[h + 1] / nhc_Q[h + 1] * wgt_sy[s] * qsmalldt);
                }
                if (nchain > 1)
                    px[nchain - 1] += (px[nchain - 2] * px[nchain - 2] / nhc_Q[nchain - 2] - Te) * wgt_sy[s] * hsmalldt;
            }
        }
    }
    return 0;
}
