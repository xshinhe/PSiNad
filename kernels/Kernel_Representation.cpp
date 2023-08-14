#include "Kernel_Representation.h"

#include "Kernel_Declare.h"
#include "Kernel_Elec.h"
#include "Kernel_NADForce.h"

#define ARRAY_SHOW(_A, _n1, _n2)                                                     \
    ({                                                                               \
        std::cout << "Show Array <" << #_A << ">\n";                                 \
        int _idxA = 0;                                                               \
        for (int _i = 0; _i < (_n1); ++_i) {                                         \
            for (int _j = 0; _j < (_n2); ++_j) std::cout << FMT(4) << (_A)[_idxA++]; \
            std::cout << std::endl;                                                  \
        }                                                                            \
    })

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
    V  = DS->reg<num_real>("model.V", Dimension::FF);
    dV = DS->reg<num_real>("model.dV", Dimension::NFF);
    // ddV = DS->reg<num_real>("model.ddV", Dimension::NNFF);
    E        = DS->reg<num_real>("model.rep.E", Dimension::F);
    T        = DS->reg<num_real>("model.rep.T", Dimension::FF);
    dE       = DS->reg<num_real>("model.rep.dE", Dimension::NFF);
    invdiffE = DS->reg<num_real>("model.rep.invdiffE", Dimension::FF);
    vedE     = DS->reg<num_real>("model.rep.vedE", Dimension::FF);
    dEprime  = DS->reg<num_real>("model.rep.dEprime", Dimension::NFF);
    // ddE = DS->reg<num_real>("model.rep.ddE", Dimension::NNFF);
    L  = DS->reg<num_real>("model.rep.L", Dimension::F);
    R  = DS->reg<num_complex>("model.rep.R", Dimension::FF);
    dL = DS->reg<num_complex>("model.rep.dL", Dimension::NFF);
    // ddL = DS->reg<num_complex>("model.rep.ddL", Dimension::NNFF);
    H  = DS->reg<num_complex>("model.rep.H", Dimension::FF);
    dH = DS->reg<num_complex>("model.rep.dH", Dimension::NFF);
    // ddH = DS->reg<num_complex>("model.rep.ddH", Dimension::NNFF);

    m  = DS->reg<num_real>("integrator.m", Dimension::N);
    x  = DS->reg<num_real>("integrator.x", Dimension::N);
    p  = DS->reg<num_real>("integrator.p", Dimension::N);
    ve = DS->reg<num_real>("integrator.ve", Dimension::N);

    Told     = DS->reg<num_real>("model.rep.Told", Dimension::FF);
    Tnew     = DS->reg<num_real>("model.rep.Tnew", Dimension::FF);
    Enew     = DS->reg<num_real>("model.rep.Enew", Dimension::F);
    TtTold   = DS->reg<num_real>("model.rep.TtTold", Dimension::FF);
    Matr_tmp = DS->reg<num_real>("model.rep.Matr", Dimension::FF);
}

void Kernel_Representation::init_calc_impl(int stat) {
    do_refer = false;
    exec_kernel(stat);
    do_refer = true;
    _DataSet->set("init.T", T, Dimension::FF);
}

int Kernel_Representation::exec_kernel_impl(int stat) {
    if (Dimension::F <= 1) return 0;

    switch (representation_type) {
        case RepresentationPolicy::Diabatic: {
            EigenSolve(E, T, V, Dimension::F);
            break;
        }
        case RepresentationPolicy::Adiabatic: {
            if (!onthefly) {  // solve adiabatic from diabatic (otherwise E/dE should be provided from Ab Initio Calc)
                // solve E
                for (int i = 0; i < Dimension::FF; ++i) Told[i] = T[i];  // backup old T matrix
                EigenSolve(E, T, V, Dimension::F);                       // solve new eigen problem

                if (do_refer) {  ///< refer the sign and order of the previous step
                    // calculate permutation matrix = rountint(T^ * Told)
                    ARRAY_MATMUL_TRANS1(TtTold, T, Told, Dimension::F, Dimension::F, Dimension::F);
                    double vset = 0.1 * std::sqrt(1.0e0 / Dimension::F);
                    for (int i = 0; i < Dimension::F; ++i) {
                        double maxnorm = 0;
                        int csr1 = 0, csr2 = 0, csr12 = 0;
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

                    // adjust order of eigenvectors
                    ARRAY_MATMUL(Tnew, T, TtTold, Dimension::F, Dimension::F, Dimension::F);
                    for (int i = 0; i < Dimension::FF; ++i) T[i] = Tnew[i];

                    // adjust order of eigenvalues
                    for (int i = 0; i < Dimension::FF; ++i) TtTold[i] = std::abs(TtTold[i]);
                    ARRAY_MATMUL(Enew, E, TtTold, 1, Dimension::F, Dimension::F);
                    for (int i = 0; i < Dimension::F; ++i) E[i] = Enew[i];
                }

                // Eigen::Map<EigMX<double>> Map_T(T, Dimension::F, Dimension::F);
                // Eigen::Map<EigMX<double>> Map_dV(dV, Dimension::NF, Dimension::F);
                // Eigen::Map<EigMX<double>> Map_dE(dE, Dimension::NF, Dimension::F);
                // Map_dE = ((Map_dV * Map_T).transpose().reshaped(Dimension::NF, Dimension::F) * Map_T)
                //              .reshaped(Dimension::F, Dimension::NF)
                //              .transpose();

                Eigen::Map<EigMX<double>> Map_T(T, Dimension::F, Dimension::F);
                Eigen::Map<EigMX<double>> Map_dV(dV, Dimension::NF, Dimension::F);
                Eigen::Map<EigMX<double>> Map_dE1(dE, Dimension::NF, Dimension::F);
                Eigen::Map<EigMX<double>> Map_dE2(dE, Dimension::F, Dimension::NF);

                Map_dE2 = Map_dV.transpose();
                Map_dE2 = (Map_T.adjoint() * Map_dE2).eval();
                Map_dE1 = (Map_dE1 * Map_T).eval();
                Map_dE1 = Map_dE2.transpose().eval();

                // solve dE (time consuming)
                // num_real *dEi = dE, *dVi = dV;
                // for (int i = 0; i < Dimension::N; ++i) {  // calc dEi = T^*dVi*T
                //     ARRAY_MATMUL3_TRANS1(dEi, T, dVi, T, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
                //     dVi += Dimension::FF, dEi += Dimension::FF;
                // }
            }

            // calc H = E - im*nacv*np / nm in first step
            for (int i = 0, ij = 0; i < Dimension::F; ++i) {
                for (int j = 0; j < Dimension::F; ++j, ++ij) {  //
                    H[ij]        = ((i == j) ? E[i] : phys::math::iz);
                    invdiffE[ij] = (i == j) ? 0.0e0 : 1 / (E[i] - E[j]);
                }
            }

            for (int i = 0; i < Dimension::N; ++i) ve[i] = p[i] / m[i];
            ARRAY_MATMUL(vedE, ve, dE, 1, Dimension::N, Dimension::FF);
            for (int jk = 0; jk < Dimension::FF; ++jk) H[jk] -= phys::math::im * vedE[jk] * (-invdiffE[jk]);

            // calc H = E - im*nacv*np/nm in second step (nacv = dEjk/(Ek-Ej))
            // num_real* dEi = dE;
            // for (int i = 0; i < Dimension::N; ++i) {
            //     num_real vi = p[i] / m[i];
            //     for (int j = 0, jk = 0; j < Dimension::F; ++j) {  // @TODO optimize
            //         for (int k = 0; k < Dimension::F; ++k, ++jk) {
            //             if (j == k) continue;
            //             H[jk] -= phys::math::im * vi * dEi[jk] * (-invdiffE[jk]);  // / (E[k] - E[j]);
            //         }
            //     }
            //     dEi += Dimension::FF;
            // }

            if (phase_correction) {
                num_real Ekin = 0;
                for (int j = 0; j < Dimension::N; ++j) Ekin += 0.5f * p[j] * p[j] / m[j];
                double Epes = 0.0f;
                if (Kernel_NADForce::NADForce_type == NADForcePolicy::BO) {
                    Epes = E[*Kernel_Elec::occ_nuc];
                } else {
                    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
                        Epes += std::real(Kernel_Elec::rho_ele[ii]) * E[i];
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
    return 0;
}

RepresentationPolicy::_type Kernel_Representation::representation_type;
RepresentationPolicy::_type Kernel_Representation::ini_repr_type;
RepresentationPolicy::_type Kernel_Representation::ele_repr_type;
RepresentationPolicy::_type Kernel_Representation::nuc_repr_type;
RepresentationPolicy::_type Kernel_Representation::tcf_repr_type;
bool Kernel_Representation::onthefly;

};  // namespace PROJECT_NS
