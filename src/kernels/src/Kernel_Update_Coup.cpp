#include "kids/Kernel_Update_Coup.h"

#include <algorithm>

#include "kids/Kernel_Elec_Utils.h"
#include "kids/Kernel_Random.h"
#include "kids/Kernel_Representation.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_Update_Coup::getName() { return "Kernel_Update_Coup"; }

int Kernel_Update_Coup::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Update_Coup::setInputParam_impl(std::shared_ptr<Param> PM) {
    swarm_type = SwarmCoupPolicy::_from(_param->get_string({"solver.swarm_flag"}, LOC(), "ccps"));
    dt         = _param->get_real({"model.dt", "solver.dt"}, LOC(), phys::time_d);
    sigma_nuc  = _param->get_real({"solver.sigma_nuc"}, LOC(), -1.0e0);  // negative means use variance
    sigma_ele  = _param->get_real({"solver.sigma_ele"}, LOC(), 0.2e0);

    gamma1 = _param->get_real({"solver.gamma"}, LOC(), elec_utils::gamma_wigner(Dimension::F));
    if (gamma1 < -1.5) gamma1 = elec_utils::gamma_opt(Dimension::F);
    if (gamma1 < -0.5) gamma1 = elec_utils::gamma_wigner(Dimension::F);
    gamma2 = (1 - gamma1) / (1.0f + Dimension::F * gamma1);
    xi1    = (1 + Dimension::F * gamma1);
    xi2    = (1 + Dimension::F * gamma2);
}

void Kernel_Update_Coup::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    relwgt      = DS->def(DATA::integrator::COUP::relwgt);
    gf_x        = DS->def(DATA::integrator::COUP::gf_x);
    gf_p        = DS->def(DATA::integrator::COUP::gf_p);
    gf_c        = DS->def(DATA::integrator::COUP::gf_c);
    gf_all      = DS->def(DATA::integrator::COUP::gf_all);
    avgx        = DS->def(DATA::integrator::COUP::avgx);
    varx        = DS->def(DATA::integrator::COUP::varx);
    avgp        = DS->def(DATA::integrator::COUP::avgp);
    varp        = DS->def(DATA::integrator::COUP::varp);
    avgxf       = DS->def(DATA::integrator::COUP::avgxf);
    varxf       = DS->def(DATA::integrator::COUP::varxf);
    xintercept  = DS->def(DATA::integrator::COUP::xintercept);
    xinterceptf = DS->def(DATA::integrator::COUP::xinterceptf);
    xslope      = DS->def(DATA::integrator::COUP::xslope);
    term_1      = DS->def(DATA::integrator::COUP::term_1);
    term_2      = DS->def(DATA::integrator::COUP::term_2);
    fadiat      = DS->def(DATA::integrator::COUP::fadiat);
    pb          = DS->def(DATA::integrator::COUP::pb);  // quantum momentum

    // pf_cross      = DS->def<kids_bool>(DATA::integrator::pf_cross);
    U = DS->def(DATA::integrator::U);
    // Udt     = DS->def(DATA::integrator::Udt);
    Ucdt    = DS->def(DATA::integrator::Ucdt);
    dV      = DS->def(DATA::model::dV);
    dE      = DS->def(DATA::model::rep::dE);
    x       = DS->def(DATA::integrator::x);
    p       = DS->def(DATA::integrator::p);
    m       = DS->def(DATA::integrator::m);
    f       = DS->def(DATA::integrator::f);
    c       = DS->def(DATA::integrator::c);
    rho_ele = DS->def(DATA::integrator::rho_ele);
    K1      = DS->def(DATA::integrator::K1);
    K2      = DS->def(DATA::integrator::K2);
    rhored  = DS->def(DATA::integrator::rhored);  // for saving K1K2
    T       = DS->def(DATA::model::rep::T);
    T_init  = DS->def(DATA::init::T);
}

Status& Kernel_Update_Coup::initializeKernel_impl(Status& stat) {
    // set zero for fadiat
    for (int i = 0; i < Dimension::P * Dimension::NF; ++i) fadiat[i] = 0.0e0;
    return stat;
}

Status& Kernel_Update_Coup::executeKernel_impl(Status& stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {  // transform ele_repr_type into
        auto c       = this->c.subspan(iP * Dimension::F, Dimension::F);
        auto rho_ele = this->rho_ele.subspan(iP * Dimension::FF, Dimension::FF);
        auto T       = this->T.subspan(iP * Dimension::FF, Dimension::FF);
        Kernel_Representation::transform(c.data(), T.data(), Dimension::F,      //
                                         Kernel_Representation::inp_repr_type,  //
                                         Kernel_Representation::nuc_repr_type,  //
                                         SpacePolicy::H);
        Kernel_Representation::transform(rho_ele.data(), T.data(), Dimension::F,  //
                                         Kernel_Representation::inp_repr_type,    //
                                         Kernel_Representation::nuc_repr_type,    //
                                         SpacePolicy::L);
    }

    // common block (overlap property)
    // 2) calculate some position average & variance
    for (int j = 0, jk = 0; j < Dimension::N; ++j) {
        // for total average & variance
        avgx[j] = 0.0e0;
        avgp[j] = 0.0e0;
        for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {  //
            avgx[j] += this->x[iP1 * Dimension::N + j];
            avgp[j] += this->p[iP1 * Dimension::N + j];
        }
        avgx[j] /= Dimension::P;
        avgp[j] /= Dimension::P;
        varx[j] = 0.0e0;
        varp[j] = 0.0e0;
        for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {
            varx[j] += (this->x[iP1 * Dimension::N + j] - avgx[j]) * (this->x[iP1 * Dimension::N + j] - avgx[j]);
            varp[j] += (this->p[iP1 * Dimension::N + j] - avgp[j]) * (this->p[iP1 * Dimension::N + j] - avgp[j]);
        }
        varx[j] /= Dimension::P;
        varp[j] /= Dimension::P;
        continue;
        // state-specific position average & variance (not used now, but in Gross's JCTC)
        // remember that here always tr[rho_ele]=1 so it is just physical density, rho_ele should be in adiabatic rep.
        for (int k = 0; k < Dimension::F; ++k, ++jk) {
            avgxf[jk]   = 0.0e0;
            double norm = 0.0e0;
            for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {
                double wgthere = std::real(rho_ele[iP1 * Dimension::FF + k * Dimension::Fadd1]);
                avgxf[jk] += this->x[iP1 * Dimension::N + j];
                norm += wgthere;
            }
            avgxf[jk] /= norm;
            varxf[jk] = 0.0e0;
            norm      = 0.0e0;
            for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {
                double wgthere = std::real(rho_ele[iP1 * Dimension::FF + k * Dimension::Fadd1]);
                varxf[jk] += (x[iP1 * Dimension::N + j] - avgxf[jk]) *  //
                             (x[iP1 * Dimension::N + j] - avgxf[jk]) * wgthere;
                norm += wgthere;
            }
            varxf[jk] /= norm;
        }
    }
    // PRINT_ARRAY(avgx, 1, Dimension::N);
    // PRINT_ARRAY(avgp, 1, Dimension::N);
    // PRINT_ARRAY(varx, 1, Dimension::N);
    // PRINT_ARRAY(varp, 1, Dimension::N);
    // PRINT_ARRAY(varx, Dimension::N, 1);
    // PRINT_ARRAY(varxf, Dimension::N, Dimension::F);

    // // overlap with width
    // 1)) constant scheme read from Param
    // 2)) proportional to variance
    // double norm_nuc = 1.0e0 / std::sqrt(2.0e0 * phys::math::pi) / sigma_nuc;
    // double norm_ele = 1.0e0 / std::sqrt(2.0e0 * phys::math::pi) / sigma_ele;
    for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {
        span<kids_real>    x1 = this->x.subspan(iP1 * Dimension::N, Dimension::N);
        span<kids_real>    p1 = this->p.subspan(iP1 * Dimension::N, Dimension::N);
        span<kids_complex> c1 = this->c.subspan(iP1 * Dimension::F, Dimension::F);

        for (int iP2 = 0; iP2 < Dimension::P; ++iP2) {
            span<kids_real>    x2 = this->x.subspan(iP2 * Dimension::N, Dimension::N);
            span<kids_real>    p2 = this->p.subspan(iP2 * Dimension::N, Dimension::N);
            span<kids_complex> c2 = this->c.subspan(iP2 * Dimension::F, Dimension::F);

            std::size_t P1P2 = iP1 * Dimension::P + iP2;
            gf_x[P1P2]       = 1.0e0;
            gf_p[P1P2]       = 1.0e0;
            gf_c[P1P2]       = 1.0e0;
            for (int j = 0; j < Dimension::N; ++j) {
                double vareff;
                vareff = (sigma_nuc < 0.0) ? varx[j] * Dimension::N : sigma_nuc * sigma_nuc;
                gf_x[P1P2] *= std::exp(-0.5e0 * (x1[j] - x2[j]) * (x1[j] - x2[j]) / vareff);
                vareff = (sigma_nuc < 0.0) ? varp[j] * Dimension::N : 0.25e0 / (sigma_nuc * sigma_nuc);
                gf_p[P1P2] *= std::exp(-0.5e0 * (p1[j] - p2[j]) * (p1[j] - p2[j]) / vareff);
            }
            for (int k = 0; k < Dimension::F; ++k) {
                gf_c[P1P2] *=                                                              //
                    std::exp(-0.5e0 * std::abs(c1[k] - c2[k]) * std::abs(c1[k] - c2[k]) /  //
                             (sigma_ele * sigma_ele) / (double) Dimension::F);
            }
        }
    }
    for (int i = 0; i < Dimension::P * Dimension::P; ++i) gf_all[i] = gf_x[i] * gf_p[i] * gf_c[i];

    // force matrix in nuc_repr_type
    auto& Force = (Kernel_Representation::nuc_repr_type == RepresentationPolicy::Adiabatic) ? this->dE : this->dV;

    // different procedure for hcps & ccps
    if (swarm_type == SwarmCoupPolicy::hcps) {
        // 1) the gradient integrated with time (the adiabatic representation must be aviable!)
        for (int iP1 = 0, P1jk = 0; iP1 < Dimension::P; ++iP1) {
            for (int j = 0; j < Dimension::N; ++j) {
                for (int k = 0; k < Dimension::F; ++k, ++P1jk) {
                    fadiat[P1jk] += Force[iP1 * Dimension::NFF + j * Dimension::FF + k * Dimension::Fadd1] * dt;
                }
            }
        }
        // PRINT_ARRAY(fadiat, Dimension::P, Dimension::N);

        // the quantum momentum: 1. state-specific intercept (not the swarm specific)
        for (int j = 0, jik = 0; j < Dimension::N; ++j) {
            for (int i = 0; i < Dimension::F; ++i) {
                for (int k = 0; k < Dimension::F; ++k, ++jik) {
                    xinterceptf[jik] = 0.0e0;
                    if (i == k) continue;
                    double norm = 0.0e0;
                    for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {
                        double w1 = std::abs(c[iP1 * Dimension::F + i] * c[iP1 * Dimension::F + i]);
                        double w2 = std::abs(c[iP1 * Dimension::F + k] * c[iP1 * Dimension::F + k]);
                        double w3 = fadiat[iP1 * Dimension::NF + j * Dimension::F + k] -
                                    fadiat[iP1 * Dimension::NF + j * Dimension::F + i];  // (i,k) > k - i
                        norm += w1 * w2 * w3;
                        xinterceptf[jik] += w1 * w2 * w3 * this->x[iP1 * Dimension::N + j];
                    }
                    xinterceptf[jik] /= norm;
                }
            }
        }
        // PRINT_ARRAY(xinterceptf, Dimension::N, Dimension::FF);  // should be small
        // the quantum momentum: 2. slope type is state dependent in JCTC
        // the quantum momentum: 2. slope type is traj dependent in JPCL
        for (int iP1 = 0, P1j = 0; iP1 < Dimension::P; ++iP1) {
            for (int j = 0; j < Dimension::N; ++j, ++P1j) {
                xslope[P1j] = 0.0;
                double norm = 0.0e0;
                for (int iP2 = 0; iP2 < Dimension::P; ++iP2) {
                    norm += gf_x[iP1 * Dimension::P + iP2];
                    xslope[P1j] += 0.5e0 / varx[j] * gf_x[iP1 * Dimension::P + iP2];
                }
                xslope[P1j] /= norm;
            }
        }
        // PRINT_ARRAY(xslope, Dimension::P, Dimension::N);  // should be small
        // continue;
        // return stat;

        // HCPS EOM
        for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {
            for (int i = 0; i < Dimension::F; ++i) {
                for (int k = 0; k < Dimension::F; ++k) {
                    double val = 0.0e0;
                    for (int j = 0; j < Dimension::N; ++j) {
                        val += 2.0e0 / m[j] * xslope[iP1 * Dimension::N + j] *
                               (x[iP1 * Dimension::N + j] - xinterceptf[j * Dimension::FF + i * Dimension::F + k]) *
                               fadiat[iP1 * Dimension::NF + j * Dimension::F + i];
                    }
                    val *= std::abs(c[iP1 * Dimension::F + i] * c[iP1 * Dimension::F + i]);
                    val *= std::abs(c[iP1 * Dimension::F + k] * c[iP1 * Dimension::F + k]);
                    term_1[iP1 * Dimension::FF + i * Dimension::F + k] = val;
                }
            }
            // HCPS (step 1: update momentum)
            for (int j = 0; j < Dimension::N; ++j) {
                double fcoup = 0.0e0;
                for (int i = 0; i < Dimension::F; ++i) {
                    for (int k = 0; k < Dimension::F; ++k) {
                        if (i == k) continue;
                        fcoup += term_1[iP1 * Dimension::FF + i * Dimension::F + k] *
                                 (fadiat[iP1 * Dimension::NF + j * Dimension::F + k] -
                                  fadiat[iP1 * Dimension::NF + j * Dimension::F + i]);
                    }
                }
                p[iP1 * Dimension::N + j] += fcoup * dt;  // note fadiat is gradient than force
            }

            // HCPS (step 1: decoherence in adiabatic)
            auto Ucdt1 = Ucdt.subspan(iP1 * Dimension::FF, Dimension::FF);  // [adia]
            auto c1    = c.subspan(iP1 * Dimension::F, Dimension::F);       // [adia]
            for (int ik = 0; ik < Dimension::FF; ++ik) Ucdt1[ik] = 0.0e0;
            double norm = 0.0e0;
            for (int i = 0; i < Dimension::F; ++i) {
                double ucoup = 0.0;
                for (int k = 0; k < Dimension::F; ++k) {
                    if (i == k) continue;
                    for (int j = 0; j < Dimension::N; ++j) {
                        ucoup += 1.0e0 / m[j] * xslope[iP1 * Dimension::N + j] *
                                 (x[iP1 * Dimension::N + j] - xinterceptf[j * Dimension::FF + i * Dimension::F + k]) *
                                 std::abs(c[iP1 * Dimension::F + k] * c[iP1 * Dimension::F + k]) *
                                 (fadiat[iP1 * Dimension::NF + j * Dimension::F + k] -
                                  fadiat[iP1 * Dimension::NF + j * Dimension::F + i]);
                    }
                }
                // c = c + ucoup * c * dt = exp(x*dt) * c ==> x = ln(1+ucoup*dt) / dt
                Ucdt1[i * Dimension::Fadd1] = 1.0e0 + ucoup * dt;
                norm += (1.0e0 + ucoup * dt) * (1.0e0 + ucoup * dt) * std::abs(c1[i] * c1[i]);
            }
            norm = std::sqrt(norm);
            for (int ik = 0; ik < Dimension::FF; ++ik) Ucdt1[ik] /= norm;

            // convert Ucdt from nuc_repr_type to ele_repr_type & merge in U in ele_repr_type
            auto T1 = T.subspan(iP1 * Dimension::FF, Dimension::FF);
            Kernel_Representation::transform(Ucdt1.data(), T1.data(), Dimension::F,  //
                                             Kernel_Representation::nuc_repr_type,   //
                                             Kernel_Representation::ele_repr_type,   //
                                             SpacePolicy::L);
            auto U1 = U.subspan(iP1 * Dimension::FF, Dimension::FF);
            ARRAY_MATMUL(U1.data(), Ucdt1.data(), U1.data(), Dimension::F, Dimension::F, Dimension::F);
        }
        // PRINT_ARRAY(U, Dimension::P, Dimension::FF);
        // exit(0);
    } else if (swarm_type == SwarmCoupPolicy::ccps) {
        // calc relative weight
        for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {
            for (int j = 0; j < Dimension::N; ++j) {
                relwgt[iP1 * Dimension::N + j] = 0.0e0;
                for (int iP2 = 0; iP2 < Dimension::P; ++iP2) {
                    relwgt[iP1 * Dimension::N + j] +=
                        (p[iP1 * Dimension::N + j] - p[iP2 * Dimension::N + j]) * gf_all[iP1 * Dimension::P + iP2];
                }
            }
        }
        // CCPS EOM (better evaluated in adiabatic!)
        for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {
            span<kids_real>    Force1 = Force.subspan(iP1 * Dimension::NFF, Dimension::NFF);
            span<kids_complex> K11    = this->K1.subspan(iP1 * Dimension::FF, Dimension::FF);
            span<kids_complex> rhoe1  = this->rho_ele.subspan(iP1 * Dimension::FF, Dimension::FF);
            elec_utils::ker_from_rho(K1.data(), rhoe1.data(), xi1, gamma1, Dimension::F);
            for (int j1 = 0; j1 < Dimension::N; ++j1) {
                double fcoup = 0.0e0;
                double norm  = 0.0e0;
                for (int iP2 = 0; iP2 < Dimension::P; ++iP2) {
                    double wgthere = gf_x[iP1 * Dimension::P + iP2] * gf_p[iP1 * Dimension::P + iP2];
                    norm += wgthere;
                    span<kids_real>    Force2 = Force.subspan(iP2 * Dimension::NFF, Dimension::NFF);
                    span<kids_complex> K22    = this->K2.subspan(iP2 * Dimension::FF, Dimension::FF);
                    span<kids_complex> rhoe2  = this->rho_ele.subspan(iP2 * Dimension::FF, Dimension::FF);
                    elec_utils::ker_from_rho(K2.data(), rhoe2.data(), xi2, gamma2, Dimension::F);

                    ARRAY_MATMUL(rhored.data(), K11.data(), K22.data(),  //
                                 Dimension::F, Dimension::F, Dimension::F);

                    double* ftmp1 = Force1.data() + j1 * Dimension::FF;
                    double* ftmp2 = Force2.data() + j1 * Dimension::FF;
                    double  f1j1  = std::real(ARRAY_TRACE2(ftmp1, rhored.data(), Dimension::F, Dimension::F));
                    // don't use averaged force!! bad result!!
                    // f1j1 = 0.5e0 * f1j1 + std::real(ARRAY_TRACE2(ftmp2, rhored.data(), Dimension::F, Dimension::F));
                    double sign21 =
                        std::copysign(1.0, relwgt[iP2 * Dimension::N + j1] * relwgt[iP1 * Dimension::N + j1]);
                    fcoup += wgthere * sign21 *  //
                             std::log(std::abs(relwgt[iP2 * Dimension::N + j1] / relwgt[iP1 * Dimension::N + j1])) *
                             f1j1;
                }
                fcoup *= Dimension::F / norm;
                // if (std::abs(f[iP1 * Dimension::N + j1]) < std::abs(fcoup)) {
                //     std::cout << f[iP1 * Dimension::N + j1] << " ? " << fcoup << std::endl;
                // }
                this->p[iP1 * Dimension::N + j1] -= fcoup * dt;  // because gradient
            }
        }
        // exit(0);
    }

    for (int iP = 0; iP < Dimension::P; ++iP) {  // transform ele_repr_type into
        auto c       = this->c.subspan(iP * Dimension::F, Dimension::F);
        auto rho_ele = this->rho_ele.subspan(iP * Dimension::FF, Dimension::FF);
        auto T       = this->T.subspan(iP * Dimension::FF, Dimension::FF);
        Kernel_Representation::transform(c.data(), T.data(), Dimension::F,      //
                                         Kernel_Representation::nuc_repr_type,  //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::H);
        Kernel_Representation::transform(rho_ele.data(), T.data(), Dimension::F,  //
                                         Kernel_Representation::nuc_repr_type,    //
                                         Kernel_Representation::inp_repr_type,    //
                                         SpacePolicy::L);
    }
    return stat;
}

};  // namespace PROJECT_NS
