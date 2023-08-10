#include "Model_SystemBath.h"

// #include <glog/logging.h>

#include "../core/linalg.h"
#include "../kernels/Kernel_Dimension.h"
#include "../kernels/Kernel_NADForce.h"
#include "../kernels/Kernel_Random.h"
#include "hamiltonian_data.h"

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

void Model_SystemBath::read_param_impl(Param* PM) {
    // size information
    Nb    = _Param->get<int>("Nb", LOC());
    nbath = _Param->get<int>("nbath", LOC());
    //// CHECK_EQ(Nb * nbath, N);  // Nb*nbath must be N

    system_type   = SystemPolicy::_from(_Param->get<std::string>("system_flag", LOC(), "SB"));
    coupling_type = CouplingPolicy::_from(_Param->get<std::string>("coupling_flag", LOC(), "SB"));
}

void Model_SystemBath::init_data_impl(DataSet* DS) {
    /// 1) System
    Hsys = DS->reg<num_real>("model.Hsys", Kernel_Dimension::FF);
    memset(Hsys, 0, Kernel_Dimension::FF * sizeof(num_real));
    L = 1;
    switch (system_type) {
        case SystemPolicy::SB: {
            //// CHECK_EQ(F, 2);
            L             = 2;
            double bias   = _Param->get<double>("bias", LOC(), phys::energy_d, 1.0f);
            double delta  = _Param->get<double>("delta", LOC(), phys::energy_d, 1.0f);
            double HSB[4] = {bias, delta, delta, -bias};
            for (int i = 0; i < Kernel_Dimension::FF; ++i) Hsys[i] = HSB[i];
            break;
        }
        case SystemPolicy::FMO: {
            //// CHECK_EQ(nbath, 7);
            //// CHECK_GE(F, 7);  // F = 7, or F = 8 (include ground state in the last; be careful)
            double tmp_unit = phys::au_2_wn;
            memset(Hsys, 0, Kernel_Dimension::FF * sizeof(num_real));
            for (int i = 0, ik = 0; i < 7; ++i) {
                for (int k = 0; k < 7; ++k, ++ik) { Hsys[i * Kernel_Dimension::F + k] = HFMO_data[ik] / tmp_unit; }
            }
            break;
        }
        case SystemPolicy::SF3a: {
            //// CHECK_EQ(F, 3);
            double tmp_unit = phys::au_2_ev;
            for (int i = 0; i < Kernel_Dimension::FF; ++i) Hsys[i] = HSF3a_data[i] / tmp_unit;
            break;
        }
        case SystemPolicy::SF3b: {
            //// CHECK_EQ(F, 3);
            double tmp_unit = phys::au_2_ev;
            for (int i = 0; i < Kernel_Dimension::FF; ++i) Hsys[i] = HSF3b_data[i] / tmp_unit;
            break;
        }
        case SystemPolicy::SF5a: {
            //// CHECK_EQ(F, 5);
            double tmp_unit = phys::au_2_ev;
            for (int i = 0; i < Kernel_Dimension::FF; ++i) Hsys[i] = HSF5a_data[i] / tmp_unit;
            break;
        }
        case SystemPolicy::SF5b: {
            //// CHECK_EQ(F, 5);
            double tmp_unit = phys::au_2_ev;
            for (int i = 0; i < Kernel_Dimension::FF; ++i) Hsys[i] = HSF5b_data[i] / tmp_unit;
            break;
        }
        case SystemPolicy::FCP: {
            //// CHECK_EQ(nbath, 9);
            //// CHECK_GE(F, 9);  // F = 9, or F = 10 (include ground state; be careful)
            double tmp_unit = phys::au_2_wn;
            for (int i = 0, ik = 0; i < 9; ++i) {
                for (int k = 0; k < 9; ++k, ++ik) { Hsys[i * Kernel_Dimension::F + k] = HFCP_data[ik] / tmp_unit; }
            }
            break;
        }
        case SystemPolicy::AGG: {
            //// CHECK_GE(F, 2);
            double delta = _Param->get<double>("delta", LOC(), phys::energy_d, 1.0f);
            for (int i = 0, ik = 0; i < Kernel_Dimension::F; ++i) {
                for (int k = 0; k < Kernel_Dimension::F; ++k, ++ik)
                    Hsys[ik] = (i - k == 1 || k - i == 1) ? delta : 0.0f;
            }
            break;
        }
        case SystemPolicy::CYC: {
            //// CHECK_GE(F, 2);
            double delta = _Param->get<double>("delta", LOC(), phys::energy_d, 1.0f);
            for (int i = 0, ik = 0; i < Kernel_Dimension::F; ++i) {
                for (int k = 0; k < Kernel_Dimension::F; ++k, ++ik)
                    Hsys[ik] = (i - k == 1 || k - i == 1) ? delta : 0.0f;
            }
            Hsys[0 * Kernel_Dimension::F + (Kernel_Dimension::F - 1)] = delta;
            Hsys[(Kernel_Dimension::F - 1) * Kernel_Dimension::F + 0] = delta;
            break;
        }
        case SystemPolicy::Read: {
            std::string system_readfile = _Param->get<std::string>("system_readfile", LOC(), "system.dat");
            std::ifstream ifs(system_readfile);
            std::string H_unit_str;
            std::string firstline;
            getline(ifs, firstline);
            std::stringstream sstr(firstline);
            sstr >> H_unit_str;  ///< the firstline stores H's unit

            double H_unit = phys::us::conv(phys::au::unit, phys::us::parse(H_unit_str));
            num_real val;
            for (int i = 0; i < Kernel_Dimension::FF; ++i)
                if (ifs >> val) Hsys[i] = val / H_unit;
            ifs.close();
        }
    }

    /// 2) init Bath sub-kernel (declaration & call)
    for (auto pkernel : _kernel_vector) pkernel->init_data(DS);
    omegas  = DS->reg<double>("model.bath.omegas", Nb);
    coeffs  = DS->reg<double>("model.bath.coeffs", Nb);
    x_sigma = DS->reg<double>("model.bath.x_sigma", Nb);
    p_sigma = DS->reg<double>("model.bath.p_sigma", Nb);

    /// 3) bilinear Coupling (saving order: L, nbath, Nb, FF)
    Q   = DS->reg<double>("model.coupling.Q", nbath * Kernel_Dimension::FF);
    CL  = DS->reg<double>("model.coupling.CL", L * Nb);
    QL  = DS->reg<double>("model.coupling.QL", L * nbath * Kernel_Dimension::FF);
    Xnj = DS->reg<double>("model.coupling.Xnj", Kernel_Dimension::NFF);
    switch (coupling_type) {
        case CouplingPolicy::SB: {
            Q[0] = 1.0f, Q[1] = 0.0f, Q[2] = 0.0f, Q[3] = -1.0f;
            break;
        }
        case CouplingPolicy::SE: {
            //// CHECK_LE(nbath, F);
            for (int i = 0, idx = 0; i < nbath; ++i) {
                for (int j = 0; j < Kernel_Dimension::F; ++j) {
                    for (int k = 0; k < Kernel_Dimension::F; ++k, ++idx) Q[idx] = (i == j && i == k) ? 1.0f : 0.0f;
                }
            }
            break;
        }
        default: {
            std::string coupling_readfile = _Param->get<std::string>("coupling_readfile", LOC(), "coupling.dat");
            std::ifstream ifs(coupling_readfile);
            num_real tmp;
            for (int i = 0; i < nbath * Kernel_Dimension::FF; ++i)
                if (ifs >> tmp) Q[i] = tmp;
            ifs.close();
        }
    }

    ARRAY_CLEAR(QL, L * nbath * Kernel_Dimension::FF);
    for (int ibath = 0, idx = 0; ibath < nbath; ++ibath) {
        for (int j = 0; j < Nb; ++j) {
            for (int i = 0, iL = 0; i < Kernel_Dimension::FF; ++i, ++idx) {
                double Qval = Q[ibath * Kernel_Dimension::FF + i];
                if (Qval != 0) {
                    QL[iL * nbath * Kernel_Dimension::FF + ibath * Kernel_Dimension::FF + i] =
                        1;                               // record there is a nonzero value
                    CL[iL * Nb + j] = coeffs[j] * Qval;  // reduce this value to CL
                    iL++;
                    if (iL > L) throw std::runtime_error(" Q shoule be sparsed!");
                }
                Xnj[idx] = coeffs[j] * Q[ibath * Kernel_Dimension::FF + i];  // merge coeffs into Q to obtain Xnj
            }
        }
    }

    // model field
    mass = DS->reg<double>("model.mass", Kernel_Dimension::N);
    for (int j = 0; j < Kernel_Dimension::N; ++j) mass[j] = 1.0f;
    vpes = DS->reg<double>("model.vpes");
    grad = DS->reg<double>("model.grad", Kernel_Dimension::N);
    hess = DS->reg<double>("model.hess", Kernel_Dimension::NN);
    V    = DS->reg<double>("model.V", Kernel_Dimension::FF);
    dV   = DS->reg<double>("model.dV", Kernel_Dimension::NFF);
    // ddV  = DS->reg<double>("model.ddV", Kernel_Dimension::NNFF);

    // init & integrator
    x_init = DS->reg<double>("init.x", Kernel_Dimension::N);
    p_init = DS->reg<double>("init.p", Kernel_Dimension::N);
    x      = DS->reg<double>("integrator.x", Kernel_Dimension::N);
    p      = DS->reg<double>("integrator.p", Kernel_Dimension::N);

    // ARRAY_SHOW(Hsys, Kernel_Dimension::F, Kernel_Dimension::F);
    // ARRAY_SHOW(omegas, 1, Nb);
    // ARRAY_SHOW(coeffs, 1, Nb);
    // ARRAY_SHOW(x_sigma, 1, Nb);
    // ARRAY_SHOW(p_sigma, 1, Nb);
    // exit(-1);
}

void Model_SystemBath::init_calc_impl(int stat) {
    Kernel_Random::rand_gaussian(x_init, Kernel_Dimension::N);
    Kernel_Random::rand_gaussian(p_init, Kernel_Dimension::N);
    for (int ibath = 0, idxR = 0; ibath < nbath; ++ibath) {
        for (int j = 0; j < Nb; ++j, ++idxR) {
            x_init[idxR] = x_init[idxR] * x_sigma[j];
            p_init[idxR] = p_init[idxR] * p_sigma[j];
            // copy "init" field to "integrator"
            x[idxR] = x_init[idxR];
            p[idxR] = p_init[idxR];
        }
    }
}

int Model_SystemBath::exec_kernel_impl(int stat) {
    // note we implement mass = 1

    // calculate nuclear vpes and grad
    double term = 0.0f;
    for (int ibath = 0, idxR = 0; ibath < nbath; ++ibath) {
        for (int j = 0; j < Nb; ++j, ++idxR) {
            term += omegas[j] * omegas[j] * x[idxR] * x[idxR];
            grad[idxR] = omegas[j] * omegas[j] * x[idxR];
        }
    }
    vpes[0] = 0.5 * term;

    // calculate electronic V and dV
    for (int i = 0; i < Kernel_Dimension::FF; ++i) V[i] = Hsys[i];
    switch (coupling_type) {
        case CouplingPolicy::SB: {
            double term = 0;
            for (int ibath = 0, idxR = 0; ibath < nbath; ++ibath) {
                for (int j = 0; j < Nb; ++j, ++idxR) { term += coeffs[j] * x[idxR]; }
            }
            V[0] += term;
            V[3] -= term;
            break;
        }
        case CouplingPolicy::SE: {
            for (int ibath = 0, idxV = 0, idxR = 0; ibath < nbath; ++ibath, idxV += Kernel_Dimension::Fadd1) {
                for (int j = 0; j < Nb; ++j, ++idxR) { V[idxV] += coeffs[j] * x[idxR]; }
            }
            break;
        }
        default: {
            for (int ibath = 0, idxR = 0, idxQ0 = 0; ibath < nbath; ++ibath, idxQ0 += Kernel_Dimension::FF) {
                for (int j = 0; j < Nb; ++j, ++idxR) {
                    double cxj = coeffs[j] * x[idxR];
                    for (int i = 0, idxQ = idxQ0; i < Kernel_Dimension::FF; ++i, ++idxQ) { V[i] += Q[idxQ] * cxj; }
                }
            }
            break;
        }
    }

    if (count_exec == 0) {  // only calculate once time
        ARRAY_CLEAR(dV, Kernel_Dimension::NFF);
        for (int ibath = 0, idxQ0 = 0, idxdV = 0; ibath < nbath; ++ibath, idxQ0 += Kernel_Dimension::FF) {
            for (int j = 0; j < Nb; ++j) {
                for (int i = 0, idxQ = idxQ0; i < Kernel_Dimension::FF; ++i, ++idxQ, ++idxdV) {
                    dV[idxdV] = Q[idxQ] * coeffs[j];
                }
            }
        }
        // ARRAY_SHOW(dV, Kernel_Dimension::N, Kernel_Dimension::FF);
        // exit(-1);
    }

    // if (flag < 2) return 0;

    if (count_exec == 0) {
        ARRAY_CLEAR(hess, Kernel_Dimension::NN);
        for (int j = 0, idx = 0, Nadd1 = Kernel_Dimension::N + 1; j < Kernel_Dimension::N; ++j, idx += Nadd1) {
            hess[idx] = omegas[j % Nb] * omegas[j % Nb];
        }
        // for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    }
    return 0;
}

};  // namespace PROJECT_NS
