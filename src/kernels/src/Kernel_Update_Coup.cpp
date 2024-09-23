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
    dt         = _param->get_real({"model.dt", "solver.dt"}, LOC(), phys::time_d);  //
    sigma_nuc  = _param->get_real({"solver.sigma_nuc"}, LOC(), 0.5e0);
    sigma_ele  = _param->get_real({"solver.sigma_ele"}, LOC(), 0.2e0);
}

void Kernel_Update_Coup::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    relwgt      = DS->def(DATA::integrator::COUP::relwgt);
    gf_x        = DS->def(DATA::integrator::COUP::gf_x);
    gf_p        = DS->def(DATA::integrator::COUP::gf_p);
    gf_c        = DS->def(DATA::integrator::COUP::gf_c);
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
    U       = DS->def(DATA::integrator::U);
    Udt     = DS->def(DATA::integrator::Udt);
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
    if (swarm_type == SwarmCoupPolicy::hcps) {
        // for hcps approach, the ele_repr_type must be Adiabatic! U & Udt are also adiabatic!

        // 1) integrate fadia with time
        for (int iP1 = 0, P1jk = 0; iP1 < Dimension::P; ++iP1) {
            for (int j = 0; j < Dimension::N; ++j) {
                for (int k = 0; k < Dimension::F; ++k, ++P1jk) {
                    fadiat[P1jk] += dE[iP1 * Dimension::NFF + j * Dimension::FF + k * Dimension::Fadd1] * dt;
                }
            }
        }
        // 2) calculate some position average & variance
        for (int j = 0, jk = 0; j < Dimension::N; ++j) {
            // for total average & variance
            avgx[j] = 0.0e0;
            for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {  //
                avgx[j] += this->x[iP1 * Dimension::N + j];
            }
            avgx[j] /= Dimension::P;
            varx[j] = 0.0e0;
            for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {
                varx[j] += (x[iP1 * Dimension::N + j] - avgx[j]) * (x[iP1 * Dimension::N + j] - avgx[j]);
            }
            varx[j] /= Dimension::P;
            // state-specific position average & variance
            // remember that here always tr[rho_ele]=1 and it is just physical density
            for (int k = 0; k < Dimension::F; ++k, ++jk) {
                avgxf[jk]   = 0.0e0;
                double norm = 0.0e0;
                for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {
                    double wgthere = std::real(rho_ele[iP1 * Dimension::FF + k * Dimension::Fadd1]);
                    avgxf[jk] += x[iP1 * Dimension::N + j];
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
        // state-specific intercept (not the swarm specific)
        for (int j = 0, jik = 0; j < Dimension::N; ++j) {
            for (int i = 0; i < Dimension::F; ++i) {
                for (int k = 0; k < Dimension::F; ++k, ++jik) {
                    xinterceptf[jik] = 0.0e0;
                    double norm;
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
        // slope type 1: sum over all state (cannot understand!): JCTC
        // slope type 2: sum over all traj: JPCL
        for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {
            span<kids_real> x1 = this->x.subspan(iP1 * Dimension::N, Dimension::N);
            for (int iP2 = 0; iP2 < Dimension::P; ++iP2) {
                span<kids_real> x2   = this->x.subspan(iP2 * Dimension::N, Dimension::N);
                std::size_t     P1P2 = iP1 * Dimension::P + iP2;
                gf_x[P1P2]           = 1.0e0;
                for (int j = 0; j < Dimension::N; ++j) {
                    double norm_nuc = 1.0e0 / std::sqrt(phys::math::twopi * varx[j]);
                    gf_x[P1P2] *= norm_nuc * std::exp(-0.5e0 * (x1[j] - x2[j]) * (x1[j] - x2[j]) / varx[j]);
                }
            }
        }
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
        // HCPS EOM
        for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {
            for (int i = 0; i < Dimension::F; ++i) {
                for (int k = 0; k < Dimension::F; ++k) {
                    double val = 0.0e0;
                    for (int j = 0; j < Dimension::N; ++j) {
                        val += 2.0e0 / m[j] *
                               (xslope[iP1 * Dimension::N + j] * x[iP1 * Dimension::N + j] -
                                xinterceptf[j * Dimension::FF + i * Dimension::F + k]) *
                               fadiat[iP1 * Dimension::NF + j * Dimension::F + i];
                    }
                    val *= std::abs(c[iP1 * Dimension::F + i] * c[iP1 * Dimension::F + i]);
                    val *= std::abs(c[iP1 * Dimension::F + k] * c[iP1 * Dimension::F + k]);
                    term_1[iP1 * Dimension::FF + i * Dimension::F + k] = val;
                }
            }
            // HCPS (step 1: momentum)
            for (int j = 0; j < Dimension::N; ++j) {
                double fcoup = 0.0e0;
                for (int i = 0; i < Dimension::F; ++i) {
                    for (int k = 0; k < Dimension::F; ++k) {
                        fcoup += term_1[iP1 * Dimension::FF + i * Dimension::F + k] *
                                 (fadiat[iP1 * Dimension::NF + j * Dimension::F + k] -
                                  fadiat[iP1 * Dimension::NF + j * Dimension::F + i]);
                    }
                }
                p[iP1 * Dimension::N + j] += fcoup * dt;  // because fadiat is gradient than force
            }
            // HCPS (step 1: decoherence)
            for (int i = 0; i < Dimension::F; ++i) {
                double ucoup = 0.0;
                for (int k = 0; k < Dimension::F; ++k) {
                    for (int j = 0; j < Dimension::N; ++j) {
                        ucoup += 1.0e0 / m[j] *
                                 (xslope[iP1 * Dimension::N + j] * x[iP1 * Dimension::N + j] -
                                  xinterceptf[j * Dimension::FF + i * Dimension::F + k]) *
                                 std::abs(c[iP1 * Dimension::F + k] * c[iP1 * Dimension::F + k]) *
                                 (fadiat[iP1 * Dimension::NF + j * Dimension::F + k] -
                                  fadiat[iP1 * Dimension::NF + j * Dimension::F + i]);
                    }
                }
                ucoup = std::exp(ucoup * dt);  // cause fadiat is gradient
                U[iP1 * Dimension::FF + i * Dimension::Fadd1] *= ucoup;
            }
        }
    } else if (swarm_type == SwarmCoupPolicy::ccps) {
        // PRINT_ARRAY(x, Dimension::P, Dimension::N);
        for (int j = 0, jk = 0; j < Dimension::N; ++j) {
            // for total average & variance
            avgx[j] = 0.0e0;
            avgp[j] = 0.0e0;
            for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {  //
                avgx[j] += x[iP1 * Dimension::N + j];
                avgp[j] += p[iP1 * Dimension::N + j];
            }
            avgx[j] /= Dimension::P;
            avgp[j] /= Dimension::P;
            varx[j] = 0.0e0;
            varp[j] = 0.0e0;
            for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {
                varx[j] += (x[iP1 * Dimension::N + j] - avgx[j]) * (x[iP1 * Dimension::N + j] - avgx[j]);
                varp[j] += (p[iP1 * Dimension::N + j] - avgp[j]) * (p[iP1 * Dimension::N + j] - avgp[j]);
            }
            varx[j] /= Dimension::P;
            varp[j] /= Dimension::P;
        }
        // PRINT_ARRAY(avgx, 1, Dimension::N);
        // PRINT_ARRAY(varx, 1, Dimension::N);
        // PRINT_ARRAY(varp, 1, Dimension::N);
        // constant width
        double norm_nuc = 1.0e0 / std::sqrt(2.0e0 * phys::math::pi) / sigma_nuc;
        double norm_ele = 1.0e0 / std::sqrt(2.0e0 * phys::math::pi) / sigma_ele;
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
                    gf_x[P1P2] *= std::exp(-0.5e0 * (x1[j] - x2[j]) * (x1[j] - x2[j]) / varx[j]);
                    gf_p[P1P2] *= std::exp(-0.5e0 * (p1[j] - p2[j]) * (p1[j] - p2[j]) / varp[j]);
                }
                for (int k = 0; k < Dimension::F; ++k) {
                    gf_c[P1P2] *=                                                              //
                        std::exp(-0.5e0 * std::abs(c1[k] - c2[k]) * std::abs(c1[k] - c2[k]) /  //
                                 (sigma_ele * sigma_ele));
                }
            }
        }
        // PRINT_ARRAY(gf_x, Dimension::P, Dimension::P);
        // PRINT_ARRAY(gf_p, Dimension::P, Dimension::P);
        // PRINT_ARRAY(gf_c, Dimension::P, Dimension::P);
        // PRINT_ARRAY(p, Dimension::P, Dimension::N);

        // calc relative weight
        for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {
            for (int j = 0; j < Dimension::N; ++j) {
                relwgt[iP1 * Dimension::N + j] = 0.0e0;
                for (int iP2 = 0; iP2 < Dimension::P; ++iP2) {
                    relwgt[iP1 * Dimension::N + j] += (p[iP1 * Dimension::N + j] - p[iP2 * Dimension::N + j]) *
                                                      gf_x[iP1 * Dimension::P + iP2] *
                                                      gf_p[iP1 * Dimension::P + iP2] *  //
                                                      gf_c[iP1 * Dimension::P + iP2];
                }
            }
        }
        // PRINT_ARRAY(relwgt, Dimension::P, Dimension::N);
        // PRINT_ARRAY(K1.data(), Dimension::P, Dimension::FF);
        // PRINT_ARRAY(K2.data(), Dimension::P, Dimension::FF);
        for (int iP = 0; iP < Dimension::P; ++iP) {
            auto K1 = this->K1.subspan(iP * Dimension::FF, Dimension::FF);
            auto K2 = this->K2.subspan(iP * Dimension::FF, Dimension::FF);
            auto T  = this->T.subspan(iP * Dimension::FF, Dimension::FF);
            Kernel_Representation::transform(K1.data(), T.data(), Dimension::F,     //
                                             Kernel_Representation::inp_repr_type,  //
                                             RepresentationPolicy::Adiabatic,       //
                                             SpacePolicy::L);
            Kernel_Representation::transform(K2.data(), T.data(), Dimension::F,     //
                                             Kernel_Representation::inp_repr_type,  //
                                             RepresentationPolicy::Adiabatic,       //
                                             SpacePolicy::L);
        }

        // CCPS EOM (evaluated in adiabatic!)
        for (int iP1 = 0; iP1 < Dimension::P; ++iP1) {
            span<kids_real>    Force1 = this->dE.subspan(iP1 * Dimension::NFF, Dimension::NFF);
            span<kids_complex> K11    = this->K1.subspan(iP1 * Dimension::FF, Dimension::FF);
            for (int j1 = 0; j1 < Dimension::N; ++j1) {
                double fcoup = 0.0e0;
                double norm  = 0.0e0;
                for (int iP2 = 0; iP2 < Dimension::P; ++iP2) {
                    double wgthere = gf_x[iP1 * Dimension::P + iP2] * gf_p[iP1 * Dimension::P + iP2];
                    norm += wgthere;
                    span<kids_complex> K22 = this->K2.subspan(iP2 * Dimension::FF, Dimension::FF);  // inverse kernel
                    ARRAY_MATMUL(rhored.data(), K11.data(), K22.data(), Dimension::F, Dimension::F, Dimension::F);

                    double* ftmp = Force1.data() + j1 * Dimension::FF;
                    double  f1j1 = std::real(ARRAY_TRACE2(ftmp, rhored.data(), Dimension::F, Dimension::F));
                    fcoup += wgthere *  //
                             (relwgt[iP2 * Dimension::N + j1] / relwgt[iP1 * Dimension::N + j1] - 1.0e0) * f1j1;
                }
                fcoup *= Dimension::F / norm;
                if (std::abs(f[iP1 * Dimension::N + j1]) < std::abs(fcoup)) {
                    std::cout << f[iP1 * Dimension::N + j1] << " ? " << fcoup << std::endl;
                }
                this->p[iP1 * Dimension::N + j1] -= fcoup * dt;  // because gradient
            }
        }

        for (int iP = 0; iP < Dimension::P; ++iP) {
            auto K1 = this->K1.subspan(iP * Dimension::FF, Dimension::FF);
            auto K2 = this->K2.subspan(iP * Dimension::FF, Dimension::FF);
            auto T  = this->T.subspan(iP * Dimension::FF, Dimension::FF);
            Kernel_Representation::transform(K1.data(), T.data(), Dimension::F,     //
                                             RepresentationPolicy::Adiabatic,       //
                                             Kernel_Representation::inp_repr_type,  //
                                             SpacePolicy::L);
            Kernel_Representation::transform(K2.data(), T.data(), Dimension::F,     //
                                             RepresentationPolicy::Adiabatic,       //
                                             Kernel_Representation::inp_repr_type,  //
                                             SpacePolicy::L);
        }
        // exit(0);
    }
    return stat;
}

};  // namespace PROJECT_NS
