#ifndef Thermostat_H
#define Thermostat_H
#include <map>

#include "../core/Kernel.h"

#endif  // Thermostat_H

/// Suzuki-Yoshida decomposition framework
const int ndofs_sy               = 7;
const kids_real wgt_sy[ndofs_sy] = {0.784513610477560f, 0.235573213359357f, -1.17767998417887f, 1.3151863206839063f,
                                    -1.17767998417887f, 0.235573213359357f, 0.784513610477560f};
namespace thermo_policy {
enum _enum { None, Langevin, Andersen, NHC };
const std::map<std::string, _enum> _dict = {{"#none", None}, {"#lang", Langevin}, {"#and", Andersen}, {"#nhc", NHC}};
};  // namespace thermo_policy

class Kernel_Update_T final : public Kernel {
   public:
    Kernel_Update_T(double scale) : Kernel(), scale{scale} {};

    inline virtual const std::string name() { return "Kernel_Update_T"; }

   private:
    Option thermo_option;

    int ndofs;

    double scale;
    double dt;
    double beta;
    double gammal;
    double randu;

    double *p, *m;                          ///< saved in "integrator"
    double *c1, *c2p;                       ///< auxiliary
    double *nhc_x, *nhc_p, *nhc_G, *nhc_Q;  ///< auxiliary

    virtual void read_param_impl(Param* PM) {
        ndofs  = PM->get<int>("N", LOC());
        dt     = PM->get<double>("dt", LOC());
        gammal = PM->get<double>("gammal", LOC(), 0.1);
        dt *= scale;

        thermo_option.flag = PM->get<std::string>("thermo_flag", LOC(), "#lang");
        thermo_option.type = thermo_policy::_dict.at(thermo_option.flag);

        // read temperature
        double temp = PM->get<double>("temp", LOC(), phys::temperature_d, 1.0f);
        beta        = 1.0f / (phys::au::k * temp);  // never ignore k_Boltzman

        nthstp = PM->get<int>("nthstp", LOC(), 0);
        nchain = PM->get<int>("nchain", LOC(), 1);
        nrespa = PM->get<int>("nrespa", LOC(), 10);
    }

    virtual void init_data_impl(DataSet* DS) {
        m = DS->def<double>("integrator.m", ndofs);
        p = DS->def<double>("integrator.p", ndofs);


        // if Langevin dynamics, set optimal c1 & c2p
        c1  = DS->def<double>("integrator.param.c1", ndofs);
        c2p = DS->def<double>("integrator.param.c2p", ndofs);
        for (int i = 0; i < ndofs; ++i) {
            c1[i]  = exp(-gammal * dt);
            c2p[i] = sqrt(1.0 - c1[i] * c1[i]);
        }

        // for NHC
        if (thermo_option.type == thermo_policy::NHC) {  // allow you to re-allocate it
            nhc_x = DS->def<double>("integrator.nhc.x", nchain * ndofs);
            nhc_p = DS->def<double>("integrator.nhc.p", nchain * ndofs);
            nhc_G = DS->def<double>("integrator.nhc.G", nchain * ndofs);  ///< p^2/m - Te
            nhc_Q = DS->def<double>("integrator.nhc.Q", nchain);

            // optimal gammal = 20*dt
            for (int i = 0; i < nchain; ++i) nhc_Q[i] = gammal * gammal / beta;
            nhc_Q[0] *= ndofs;
        }
    }

    virtual int exec_kernel_impl(int stat = -1) {
        switch (thermo_option.type) {
            case thermo_policy::Langevin:
                for (int i = 0; i < ndofs; ++i) {
                    Kernel_Random::rand_gaussian(&randu);
                    p[i] = c1[i] * p[i] + c2p[i] * sqrt(m[i] / beta) * randu;
                }
                break;
            case thermo_policy::Andersen:
                for (int i = 0; i < ndofs; ++i) {
                    Kernel_Random::rand_uniform(&randu);
                    if (randu < 1 - c1[i]) {
                        Kernel_Random::rand_gaussian(&randu);
                        p[i] = sqrt(m[i] / beta) * randu;
                    }
                }
                break;
            case thermo_policy::NHC: {
                kids_real smalldt  = dt / nrespa;
                kids_real hsmalldt = 0.5f * smalldt;
                kids_real qsmalldt = 0.25f * smalldt;
                kids_real Te       = 1.0f / (phys::au::k * beta);
                kids_real* xi      = nhc_x;
                kids_real* pi      = nhc_p;
                for (int i = 0; i < ndofs; ++i, xi += nchain, pi += nchain) {
                    // sweep nchain
                    for (int k = 0; k < nrespa; ++k) {
                        for (int s = 0; s < size_sy; ++s) {
                            // update nhc_p from nchain tail to head
                            if (nchain > 1)
                                pi[nchain - 1] +=
                                    (pi[nchain - 2] * pi[nchain - 2] / nhc_Q[nchain - 2] - Te) * wgt_sy[s] * hsmalldt;
                            for (int h = nchain - 2; h > 0; --h) {
                                pi[h] *= exp(-pi[h + 1] / nhc_Q[h + 1] * wgt_sy[s] * qsmalldt);
                                pi[h] += (pi[h - 1] * pi[h - 1] / nhc_Q[h - 1] - Te) * wgt_sy[s] * hsmalldt;
                                pi[h] *= exp(-pi[h + 1] / nhc_Q[h + 1] * wgt_sy[s] * qsmalldt);
                            }
                            if (nchain > 1) pi[0] -= pi[1] / nhc_Q[1] * wgt_sy[s] * qsmalldt;
                            pi[0] += (p[i] * p[i] / m[i] - Te) * wgt_sy[s] * hsmalldt;
                            if (nchain > 1) pi[0] -= pi[1] / nhc_Q[1] * wgt_sy[s] * qsmalldt;

                            // update nhc_x and p
                            for (int h = 0; h < nchain; ++h) xi[h] += pi[h] / nhc_Q[h] * wgt_sy[s] * smalldt;
                            p[i] *= exp(-pi[0] / nhc_Q[0] * wgt_sy[s] * smalldt);

                            // update nhc_p from head to tail
                            if (nchain > 1) pi[0] -= pi[1] / nhc_Q[1] * wgt_sy[s] * qsmalldt;
                            pi[0] += (p[i] * p[i] / m[i] - Te) * wgt_sy[s] * hsmalldt;
                            if (nchain > 1) pi[0] -= pi[1] / nhc_Q[1] * wgt_sy[s] * qsmalldt;
                            for (int h = 1; h < nchain - 1; ++h) {
                                pi[h] *= exp(-pi[h + 1] / nhc_Q[h + 1] * wgt_sy[s] * qsmalldt);
                                pi[h] += (pi[h - 1] * pi[h - 1] / nhc_Q[h - 1] - Te) * wgt_sy[s] * hsmalldt;
                                pi[h] *= exp(-pi[h + 1] / nhc_Q[h + 1] * wgt_sy[s] * qsmalldt);
                            }
                            if (nchain > 1)
                                pi[nchain - 1] +=
                                    (pi[nchain - 2] * pi[nchain - 2] / nhc_Q[nchain - 2] - Te) * wgt_sy[s] * hsmalldt;
                        }
                    }
                }
                break;
            }
            case thermo_policy::None:
            default:
                break;
        }
        return 0;
    }
};
