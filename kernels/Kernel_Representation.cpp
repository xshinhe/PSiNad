#include "Kernel_Representation.h"

#include "Kernel_Dimension.h"
#include "Kernel_Elec.h"
#include "Kernel_NADForce.h"

namespace PROJECT_NS {

void Kernel_Representation::read_param_impl(Param* PM) {
    std::string rep_string = PM->get<std::string>("representation_flag", LOC(), "Diabatic");
    representation_type    = RepresentationPolicy::_from(rep_string);
    ini_repr_type          = RepresentationPolicy::_from(PM->get<std::string>("ini_repr_flag", LOC(), rep_string));
    ele_repr_type          = RepresentationPolicy::_from(PM->get<std::string>("ele_repr_flag", LOC(), rep_string));
    nuc_repr_type          = RepresentationPolicy::_from(PM->get<std::string>("nuc_repr_flag", LOC(), rep_string));
    tcf_repr_type          = RepresentationPolicy::_from(PM->get<std::string>("tcf_repr_flag", LOC(), rep_string));
    phase_correction       = PM->get<bool>("phase_correction", LOC(), false);
}

void Kernel_Representation::init_data_impl(DataSet* DS) {
    V  = DS->reg<num_real>("model.V", Kernel_Dimension::FF);
    dV = DS->reg<num_real>("model.dV", Kernel_Dimension::NFF);
    // ddV = DS->reg<num_real>("model.ddV", Kernel_Dimension::NNFF);
    E  = DS->reg<num_real>("model.rep.E", Kernel_Dimension::F);
    T  = DS->reg<num_real>("model.rep.T", Kernel_Dimension::FF);
    dE = DS->reg<num_real>("model.rep.dE", Kernel_Dimension::NFF);
    // ddE = DS->reg<num_real>("model.rep.ddE", Kernel_Dimension::NNFF);
    L  = DS->reg<num_real>("model.rep.L", Kernel_Dimension::F);
    R  = DS->reg<num_complex>("model.rep.R", Kernel_Dimension::FF);
    dL = DS->reg<num_complex>("model.rep.dL", Kernel_Dimension::NFF);
    // ddL = DS->reg<num_complex>("model.rep.ddL", Kernel_Dimension::NNFF);
    H  = DS->reg<num_complex>("model.rep.H", Kernel_Dimension::FF);
    dH = DS->reg<num_complex>("model.rep.dH", Kernel_Dimension::NFF);
    // ddH = DS->reg<num_complex>("model.rep.ddH", Kernel_Dimension::NNFF);

    m = DS->reg<num_real>("integrator.m", Kernel_Dimension::N);
    x = DS->reg<num_real>("integrator.x", Kernel_Dimension::N);
    p = DS->reg<num_real>("integrator.p", Kernel_Dimension::N);

    Told     = DS->reg<num_real>("model.rep.Told", Kernel_Dimension::FF);
    Tnew     = DS->reg<num_real>("model.rep.Tnew", Kernel_Dimension::FF);
    Enew     = DS->reg<num_real>("model.rep.Enew", Kernel_Dimension::F);
    TtTold   = DS->reg<num_real>("model.rep.TtTold", Kernel_Dimension::FF);
    Matr_tmp = DS->reg<num_real>("model.rep.Matr", Kernel_Dimension::FF);
}

void Kernel_Representation::init_calc_impl(int stat) {
    do_refer = false;
    exec_kernel(stat);
    do_refer = true;
    _DataSet->set("init.T", T, Kernel_Dimension::FF);
}

int Kernel_Representation::exec_kernel_impl(int stat) {
    if (Kernel_Dimension::F <= 1) return 0;

    switch (representation_type) {
        case RepresentationPolicy::Diabatic: {
            EigenSolve(E, T, V, Kernel_Dimension::F);
            break;
        }
        case RepresentationPolicy::Adiabatic: {
            if (!onthefly) {  // solve adiabatic from diabatic (otherwise E/dE should be provided from Ab Initio Calc)
                // solve E
                for (int i = 0; i < Kernel_Dimension::FF; ++i) Told[i] = T[i];  // backup old T matrix
                EigenSolve(E, T, V, Kernel_Dimension::F);                       // solve new eigen problem

                if (do_refer) {  ///< refer the sign and order of the previous step
                    // calculate permutation matrix = rountint(T^ * Told)
                    ARRAY_MATMUL_TRANS1(TtTold, T, Told, Kernel_Dimension::F, Kernel_Dimension::F, Kernel_Dimension::F);
                    double vset = 0.1 * std::sqrt(1.0e0 / Kernel_Dimension::F);
                    for (int i = 0; i < Kernel_Dimension::F; ++i) {
                        double maxnorm = 0;
                        int csr1 = 0, csr2 = 0, csr12 = 0;
                        for (int k1 = 0, k1k2 = 0; k1 < Kernel_Dimension::F; ++k1) {
                            for (int k2 = 0; k2 < Kernel_Dimension::F; ++k2, ++k1k2) {
                                // vmax must be larger than sqrt(1/fdim)
                                if (std::abs(TtTold[k1k2]) > maxnorm) {
                                    maxnorm = std::abs(TtTold[k1k2]);
                                    csr1 = k1, csr2 = k2, csr12 = k1k2;
                                }
                            }
                        }
                        double vsign = copysign(1.0f, TtTold[csr12]);
                        for (int k2 = 0, k1k2 = csr1 * Kernel_Dimension::F;  //
                             k2 < Kernel_Dimension::F;                       //
                             ++k2, ++k1k2) {
                            TtTold[k1k2] = 0;
                        }
                        for (int k1 = 0, k1k2 = csr2; k1 < Kernel_Dimension::F; ++k1, k1k2 += Kernel_Dimension::F) {
                            TtTold[k1k2] = 0;
                        }
                        TtTold[csr12] = vsign * vset;
                    }
                    for (int i = 0; i < Kernel_Dimension::FF; ++i) TtTold[i] = round(TtTold[i] / vset);

                    // adjust order of eigenvectors
                    ARRAY_MATMUL(Tnew, T, TtTold, Kernel_Dimension::F, Kernel_Dimension::F, Kernel_Dimension::F);
                    for (int i = 0; i < Kernel_Dimension::FF; ++i) T[i] = Tnew[i];

                    // adjust order of eigenvalues
                    for (int i = 0; i < Kernel_Dimension::FF; ++i) TtTold[i] = std::abs(TtTold[i]);
                    ARRAY_MATMUL(Enew, E, TtTold, 1, Kernel_Dimension::F, Kernel_Dimension::F);
                    for (int i = 0; i < Kernel_Dimension::F; ++i) E[i] = Enew[i];
                }

                // solve dE
                num_real *dEi = dE, *dVi = dV;
                for (int i = 0; i < Kernel_Dimension::N; ++i) {
                    // calc dEi = T^*dVi*T
                    ARRAY_MATMUL(Matr_tmp, dVi, T, Kernel_Dimension::F, Kernel_Dimension::F, Kernel_Dimension::F);
                    ARRAY_MATMUL_TRANS1(dEi, T, Matr_tmp, Kernel_Dimension::F, Kernel_Dimension::F,
                                        Kernel_Dimension::F);
                    dVi += Kernel_Dimension::FF, dEi += Kernel_Dimension::FF;
                }
            }

            // calc H = E - im*nacv*np / nm in first step
            for (int i = 0, idx = 0; i < Kernel_Dimension::F; ++i)
                for (int j = 0; j < Kernel_Dimension::F; ++j, ++idx)
                    H[idx] = ((i == j) ? phys::math::iu * E[i] : phys::math::iz);

            // calc H = E - im*nacv*np/nm in second step (nacv = dEjk/(Ek-Ej))
            num_real* dEi = dE;
            for (int i = 0; i < Kernel_Dimension::N; ++i) {
                num_real vi = p[i] / m[i];
                for (int j = 0, jk = 0; j < Kernel_Dimension::F; ++j) {  // @TODO optimize
                    for (int k = 0; k < Kernel_Dimension::F; ++k, ++jk) {
                        if (j == k) continue;
                        H[jk] -= phys::math::im * vi * dEi[jk] / (E[k] - E[j]);
                    }
                }
                dEi += Kernel_Dimension::FF;
            }

            if (phase_correction) {
                num_real Ekin = 0;
                for (int j = 0; j < Kernel_Dimension::N; ++j) Ekin += 0.5f * p[j] * p[j] / m[j];
                double Epes = 0.0f;
                if (Kernel_NADForce::NADForce_type == NADForcePolicy::BO) {
                    Epes = E[*Kernel_Elec::occ_nuc];
                } else {
                    for (int i = 0, ii = 0; i < Kernel_Dimension::F; ++i, ii += Kernel_Dimension::Fadd1) {
                        Epes += std::real(Kernel_Elec::rho_ele[ii]) * E[i];
                    }
                }
                for (int i = 0, ii = 0; i < Kernel_Dimension::F; ++i, ii += Kernel_Dimension::Fadd1) {
                    H[ii] = -2 * Ekin * sqrt(std::max<double>(1.0 + (Epes - E[i]) / Ekin, 0.0f));
                }
            }
            EigenSolve(L, R, H, Kernel_Dimension::F);  // R*L*R^ = H
            break;
        }
        case RepresentationPolicy::Force:
        case RepresentationPolicy::Density:
        default:
            break;
    }
    return 0;
}

RepresentationPolicy::_type Kernel_Representation::representation_type;
RepresentationPolicy::_type Kernel_Representation::ini_repr_type;
RepresentationPolicy::_type Kernel_Representation::ele_repr_type;
RepresentationPolicy::_type Kernel_Representation::nuc_repr_type;
RepresentationPolicy::_type Kernel_Representation::tcf_repr_type;
bool Kernel_Representation::onthefly;

};  // namespace PROJECT_NS
