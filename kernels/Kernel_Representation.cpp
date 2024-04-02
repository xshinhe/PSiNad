#include "Kernel_Representation.h"

#include "../core/linalg.h"
#include "../core/vars_list.h"
#include "Kernel_NADForce.h"

namespace PROJECT_NS {

int Kernel_Representation::transform(kids_complex* A, kids_real* T, int fdim,  //
                                     RepresentationPolicy::_type from, RepresentationPolicy::_type to,
                                     SpacePolicy::_type Stype) {
    if (from == to) return 0;
    int lda = (Stype == SpacePolicy::L) ? fdim : 1;
    if (from == RepresentationPolicy::Diabatic && to == RepresentationPolicy::Adiabatic) {
        ARRAY_MATMUL_TRANS1(A, T, A, fdim, fdim, lda);
        if (Stype == SpacePolicy::L) ARRAY_MATMUL(A, A, T, fdim, fdim, fdim);
    }
    if (from == RepresentationPolicy::Adiabatic && to == RepresentationPolicy::Diabatic) {
        ARRAY_MATMUL(A, T, A, fdim, fdim, lda);
        if (Stype == SpacePolicy::L) ARRAY_MATMUL_TRANS2(A, A, T, fdim, fdim, fdim);
    }
    return 0;
}

void Kernel_Representation::read_param_impl(Param* PM) {
    std::string rep_string = PM->get<std::string>("representation_flag", LOC(), "Diabatic");
    representation_type    = RepresentationPolicy::_from(rep_string);
    inp_repr_type          = RepresentationPolicy::_from(PM->get<std::string>("inp_repr_flag", LOC(), rep_string));
    ele_repr_type          = RepresentationPolicy::_from(PM->get<std::string>("ele_repr_flag", LOC(), rep_string));
    nuc_repr_type          = RepresentationPolicy::_from(PM->get<std::string>("nuc_repr_flag", LOC(), rep_string));
    tcf_repr_type          = RepresentationPolicy::_from(PM->get<std::string>("tcf_repr_flag", LOC(), rep_string));
    phase_correction       = PM->get<bool>("phase_correction", LOC(), false);
    basis_switch           = PM->get<bool>("basis_switch", LOC(), false);
}

void Kernel_Representation::init_data_impl(DataSet* DS) {
    V  = DS->def(DATA::model::V);
    dV = DS->def(DATA::model::dV);
    // ddV = DS->def(DATA::model::ddV);
    E_copy = DS->def(DATA::integrator::E);
    E      = DS->def(DATA::model::rep::E);
    T      = DS->def(DATA::model::rep::T);
    Told   = DS->def(DATA::model::rep::Told);
    dE     = DS->def(DATA::model::rep::dE);
    // ddE = DS->def(DATA::model::rep::ddE);
    L = DS->def(DATA::model::rep::L);
    R = DS->def(DATA::model::rep::R);
    H = DS->def(DATA::model::rep::H);

    m       = DS->def(DATA::integrator::m);
    p       = DS->def(DATA::integrator::p);
    occ_nuc = DS->def(DATA::integrator::occ_nuc);
    rho_ele = DS->def(DATA::integrator::rho_ele);

    TtTold = DS->def(DATA::integrator::tmp::TtTold);
    ve     = DS->def(DATA::integrator::tmp::ve);
    vedE   = DS->def(DATA::integrator::tmp::vedE);
}

void Kernel_Representation::init_calc_impl(int stat) {
    do_refer = false;
    exec_kernel(stat);
    do_refer = true;
    _DataSet->def("init.T", T, Dimension::PFF);
}

int Kernel_Representation::exec_kernel_impl(int stat) {
    if (Dimension::F <= 1) return 0;

    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real*    V    = this->V + iP * Dimension::FF;
        kids_real*    E    = this->E + iP * Dimension::F;
        kids_real*    T    = this->T + iP * Dimension::FF;
        kids_real*    Told = this->Told + iP * Dimension::FF;
        kids_real*    dV   = this->dV + iP * Dimension::NFF;
        kids_real*    dE   = this->dE + iP * Dimension::NFF;
        kids_real*    L    = this->L + iP * Dimension::F;
        kids_complex* R    = this->R + iP * Dimension::FF;
        kids_complex* H    = this->H + iP * Dimension::FF;

        kids_real*    p       = this->p + iP * Dimension::N;
        kids_real*    m       = this->m + iP * Dimension::N;
        int*          occ_nuc = this->occ_nuc + iP;
        kids_complex* rho_ele = this->rho_ele + iP * Dimension::FF;

        switch (representation_type) {
            case RepresentationPolicy::Diabatic: {
                EigenSolve(E, T, V, Dimension::F);
                for (int ik = 0; ik < Dimension::FF; ++ik) H[ik] = V[ik];
                break;
            }
            case RepresentationPolicy::Adiabatic: {
                if (!onthefly) {  // solve adiabatic from diabatic (otherwise E/dE should be provided from Ab Initio

                    for (int i = 0; i < Dimension::FF; ++i) Told[i] = T[i];  // backup old T matrix
                    EigenSolve(E, T, V, Dimension::F);                       // solve new eigen problem

                    if (do_refer) {  ///< refer the sign and order of the previous step
                        // calculate permutation matrix = rountint(T^ * Told)
                        ARRAY_MATMUL_TRANS1(TtTold, T, Told, Dimension::F, Dimension::F, Dimension::F);

                        // ARRAY_SHOW(E, 1, Dimension::F);
                        // ARRAY_SHOW(T, Dimension::F, Dimension::F);
                        // ARRAY_SHOW(TtTold, Dimension::F, Dimension::F);
                        // ARRAY_SHOW(TtTold, Dimension::F, Dimension::F); // @debug

                        if (!basis_switch) {
                            for (int i = 0, ik = 0; i < Dimension::F; ++i) {
                                for (int k = 0; k < Dimension::F; ++k, ++ik) {
                                    TtTold[ik] = (i == k) ? copysign(1.0f, TtTold[ik]) : 0;
                                }
                            }
                        } else {
                            double vset = 0.1 * std::sqrt(1.0e0 / Dimension::F);
                            for (int i = 0; i < Dimension::F; ++i) {
                                double maxnorm = 0;
                                int    csr1 = 0, csr2 = 0, csr12 = 0;
                                for (int k1 = 0, k1k2 = 0; k1 < Dimension::F; ++k1) {
                                    for (int k2 = 0; k2 < Dimension::F; ++k2, ++k1k2) {
                                        // vmax must be larger than sqrt(1/fdim)
                                        if (std::abs(TtTold[k1k2]) > maxnorm) {
                                            maxnorm = std::abs(TtTold[k1k2]);
                                            csr1 = k1, csr2 = k2, csr12 = k1k2;
                                        }
                                    }
                                }
                                double vsign = copysign(1.0f, TtTold[csr12]);
                                for (int k2 = 0, k1k2 = csr1 * Dimension::F;  //
                                     k2 < Dimension::F;                       //
                                     ++k2, ++k1k2) {
                                    TtTold[k1k2] = 0;
                                }
                                for (int k1 = 0, k1k2 = csr2; k1 < Dimension::F; ++k1, k1k2 += Dimension::F) {
                                    TtTold[k1k2] = 0;
                                }
                                TtTold[csr12] = vsign * vset;
                            }
                            for (int i = 0; i < Dimension::FF; ++i) TtTold[i] = round(TtTold[i] / vset);
                        }

                        // ARRAY_SHOW(TtTold, Dimension::F, Dimension::F);

                        // adjust order of eigenvectors & eigenvalues
                        ARRAY_MATMUL(T, T, TtTold, Dimension::F, Dimension::F, Dimension::F);
                        if (basis_switch) {
                            for (int i = 0; i < Dimension::FF; ++i) TtTold[i] = std::abs(TtTold[i]);
                            ARRAY_MATMUL(E, E, TtTold, 1, Dimension::F, Dimension::F);
                        }
                    }

                    if (FORCE_OPT::BATH_FORCE_BILINEAR) {
                        int&    B    = FORCE_OPT::nbath;
                        int&    J    = FORCE_OPT::Nb;
                        int     JFF  = J * Dimension::FF;
                        double* dVb0 = dV;
                        double* dEb0 = dE;
                        for (int b = 0, bb = 0; b < B; ++b, bb += Dimension::Fadd1, dVb0 += JFF, dEb0 += JFF) {
                            ARRAY_MATMUL3_TRANS1(dEb0, T, dVb0, T, Dimension::F, Dimension::F, Dimension::F,
                                                 Dimension::F);
                            for (int j = 0, jik = 0, jbb = bb; j < J; ++j, jbb += Dimension::FF) {
                                double scale = dVb0[jbb] / dVb0[bb];
                                for (int ik = 0; ik < Dimension::FF; ++ik, ++jik) { dEb0[jik] = dEb0[ik] * scale; }
                            }
                        }
                    } else {
                        Eigen::Map<EigMX<double>> Map_T(T, Dimension::F, Dimension::F);
                        Eigen::Map<EigMX<double>> Map_dV(dV, Dimension::NF, Dimension::F);
                        Eigen::Map<EigMX<double>> Map_dE1(dE, Dimension::NF, Dimension::F);
                        Eigen::Map<EigMX<double>> Map_dE2(dE, Dimension::F, Dimension::NF);

                        Map_dE2 = Map_dV.transpose();
                        Map_dE2 = (Map_T.adjoint() * Map_dE2).eval();
                        Map_dE1 = (Map_dE1 * Map_T).eval();
                        Map_dE1 = Map_dE2.transpose().eval();
                    }
                }

                // calc H = E - im * nacv * p / m, and note here nacv_{ij} = dEij / (Ej - Ei)
                for (int i = 0; i < Dimension::N; ++i) ve[i] = p[i] / m[i];
                ARRAY_MATMUL(vedE, ve, dE, 1, Dimension::N, Dimension::FF);

                double Emean = 0.0e0;
                for (int i = 0; i < Dimension::F; ++i) {
                    Emean += E[i];
                    E_copy[i] = E[i];
                }
                Emean /= Dimension::F;

                for (int i = 0, ij = 0; i < Dimension::F; ++i) {
                    for (int j = 0; j < Dimension::F; ++j, ++ij) {  //
                        H[ij] = ((i == j) ? E[i] - Emean : -phys::math::im * vedE[ij] / (E[j] - E[i]));
                    }
                }

                if (phase_correction) {
                    kids_real Ekin = 0;
                    for (int j = 0; j < Dimension::N; ++j) Ekin += 0.5f * p[j] * p[j] / m[j];
                    double Epes = 0.0f;
                    if (Kernel_NADForce::NADForce_type == NADForcePolicy::BO) {
                        Epes = E[*occ_nuc];
                    } else {
                        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
                            Epes += std::real(rho_ele[ii]) * E[i];
                        }
                    }
                    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
                        H[ii] = -2 * Ekin * sqrt(std::max<double>(1.0 + (Epes - E[i]) / Ekin, 0.0f));
                    }
                }
                EigenSolve(L, R, H, Dimension::F);  // R*L*R^ = H
                break;
            }
            case RepresentationPolicy::Force:
            case RepresentationPolicy::Density:
            default:
                break;
        }
    }

    return 0;
}

RepresentationPolicy::_type Kernel_Representation::representation_type;
RepresentationPolicy::_type Kernel_Representation::inp_repr_type;
RepresentationPolicy::_type Kernel_Representation::ele_repr_type;
RepresentationPolicy::_type Kernel_Representation::nuc_repr_type;
RepresentationPolicy::_type Kernel_Representation::tcf_repr_type;
bool                        Kernel_Representation::onthefly;

};  // namespace PROJECT_NS
