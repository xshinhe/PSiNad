#include "Model_SystemBath.h"

// #include <glog/logging.h>

#include "../core/linalg.h"
#include "../kernels/Kernel_Declare.h"
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

    FORCE_OPT::nbath = nbath;
    FORCE_OPT::Nb    = Nb;

    system_type   = SystemPolicy::_from(_Param->get<std::string>("system_flag", LOC(), "SB"));
    coupling_type = CouplingPolicy::_from(_Param->get<std::string>("coupling_flag", LOC(), "SB"));
    nsamp_type    = NSampPolicy::_from(_Param->get<std::string>("nsamp_flag", LOC(), "Wigner"));
}

void Model_SystemBath::init_data_impl(DataSet* DS) {
    /// 1) System
    Hsys = DS->def<kids_real>("model.Hsys", Dimension::FF);
    memset(Hsys, 0, Dimension::FF * sizeof(kids_real));
    L = 1;
    switch (system_type) {
        case SystemPolicy::SB: {
            //// CHECK_EQ(F, 2);
            L             = 2;
            double bias   = _Param->get<double>("bias", LOC(), phys::energy_d, 1.0f);
            double delta  = _Param->get<double>("delta", LOC(), phys::energy_d, 1.0f);
            double HSB[4] = {bias, delta, delta, -bias};
            for (int i = 0; i < Dimension::FF; ++i) Hsys[i] = HSB[i];
            break;
        }
        case SystemPolicy::FMO: {
            //// CHECK_EQ(nbath, 7);
            //// CHECK_GE(F, 7);  // F = 7, or F = 8 (include ground state in the last; be careful)
            double tmp_unit = phys::au_2_wn;
            memset(Hsys, 0, Dimension::FF * sizeof(kids_real));
            for (int i = 0, ik = 0; i < 7; ++i) {
                for (int k = 0; k < 7; ++k, ++ik) { Hsys[i * Dimension::F + k] = HFMO_data[ik] / tmp_unit; }
            }
            break;
        }
        case SystemPolicy::SF3a: {
            //// CHECK_EQ(F, 3);
            double tmp_unit = phys::au_2_ev;
            for (int i = 0; i < Dimension::FF; ++i) Hsys[i] = HSF3a_data[i] / tmp_unit;
            break;
        }
        case SystemPolicy::SF3b: {
            //// CHECK_EQ(F, 3);
            double tmp_unit = phys::au_2_ev;
            for (int i = 0; i < Dimension::FF; ++i) Hsys[i] = HSF3b_data[i] / tmp_unit;
            break;
        }
        case SystemPolicy::SF5a: {
            //// CHECK_EQ(F, 5);
            double tmp_unit = phys::au_2_ev;
            for (int i = 0; i < Dimension::FF; ++i) Hsys[i] = HSF5a_data[i] / tmp_unit;
            break;
        }
        case SystemPolicy::SF5b: {
            //// CHECK_EQ(F, 5);
            double tmp_unit = phys::au_2_ev;
            for (int i = 0; i < Dimension::FF; ++i) Hsys[i] = HSF5b_data[i] / tmp_unit;
            break;
        }
        case SystemPolicy::FCP: {
            //// CHECK_EQ(nbath, 9);
            //// CHECK_GE(F, 9);  // F = 9, or F = 10 (include ground state; be careful)
            double tmp_unit = phys::au_2_wn;
            for (int i = 0, ik = 0; i < 9; ++i) {
                for (int k = 0; k < 9; ++k, ++ik) { Hsys[i * Dimension::F + k] = HFCP_data[ik] / tmp_unit; }
            }
            break;
        }
        case SystemPolicy::AGG: {
            //// CHECK_GE(F, 2);
            double delta = _Param->get<double>("delta", LOC(), phys::energy_d, 1.0f);
            for (int i = 0, ik = 0; i < Dimension::F; ++i) {
                for (int k = 0; k < Dimension::F; ++k, ++ik) Hsys[ik] = (i - k == 1 || k - i == 1) ? delta : 0.0f;
            }
            break;
        }
        case SystemPolicy::CYC: {
            //// CHECK_GE(F, 2);
            double delta = _Param->get<double>("delta", LOC(), phys::energy_d, 1.0f);
            for (int i = 0, ik = 0; i < Dimension::F; ++i) {
                for (int k = 0; k < Dimension::F; ++k, ++ik) Hsys[ik] = (i - k == 1 || k - i == 1) ? delta : 0.0f;
            }
            Hsys[0 * Dimension::F + (Dimension::F - 1)] = delta;
            Hsys[(Dimension::F - 1) * Dimension::F + 0] = delta;
            break;
        }
        case SystemPolicy::Read: {
            std::string   system_readfile = _Param->get<std::string>("system_readfile", LOC(), "system.dat");
            std::ifstream ifs(system_readfile);
            std::string   H_unit_str;
            std::string   firstline;
            getline(ifs, firstline);
            std::stringstream sstr(firstline);
            sstr >> H_unit_str;  ///< the firstline stores H's unit

            double    H_unit = phys::us::conv(phys::au::unit, phys::us::parse(H_unit_str));
            kids_real val;
            for (int i = 0; i < Dimension::FF; ++i)
                if (ifs >> val) Hsys[i] = val / H_unit;
            ifs.close();
        }
    }

    /// 2) init Bath sub-kernel (declaration & call)
    for (auto pkernel : _kernel_vector) pkernel->init_data(DS);
    omegas  = DS->def<kids_real>("model.bath.omegas", Nb);
    coeffs  = DS->def<kids_real>("model.bath.coeffs", Nb);
    x_sigma = DS->def<kids_real>("model.bath.x_sigma", Nb);
    p_sigma = DS->def<kids_real>("model.bath.p_sigma", Nb);

    /// 3) bilinear Coupling (saving order: L, nbath, Nb, FF)
    Q   = DS->def<kids_real>("model.coupling.Q", nbath * Dimension::FF);
    CL  = DS->def<kids_real>("model.coupling.CL", L * Nb);
    QL  = DS->def<kids_real>("model.coupling.QL", L * nbath * Dimension::FF);
    Xnj = DS->def<kids_real>("model.coupling.Xnj", Dimension::NFF);
    switch (coupling_type) {
        case CouplingPolicy::SB: {
            Q[0] = 1.0f, Q[1] = 0.0f, Q[2] = 0.0f, Q[3] = -1.0f;
            break;
        }
        case CouplingPolicy::SE: {
            //// CHECK_LE(nbath, F);
            for (int i = 0, idx = 0; i < nbath; ++i) {
                for (int j = 0; j < Dimension::F; ++j) {
                    for (int k = 0; k < Dimension::F; ++k, ++idx) Q[idx] = (i == j && i == k) ? 1.0f : 0.0f;
                }
            }
            break;
        }
        default: {
            std::string   coupling_readfile = _Param->get<std::string>("coupling_readfile", LOC(), "coupling.dat");
            std::ifstream ifs(coupling_readfile);
            kids_real     tmp;
            for (int i = 0; i < nbath * Dimension::FF; ++i)
                if (ifs >> tmp) Q[i] = tmp;
            ifs.close();
        }
    }

    ARRAY_CLEAR(QL, L * nbath * Dimension::FF);
    for (int ibath = 0, idx = 0; ibath < nbath; ++ibath) {
        for (int j = 0; j < Nb; ++j) {
            for (int i = 0, iL = 0; i < Dimension::FF; ++i, ++idx) {
                double Qval = Q[ibath * Dimension::FF + i];
                if (Qval != 0) {
                    QL[iL * nbath * Dimension::FF + ibath * Dimension::FF + i] = 1;  // record there is a nonzero value
                    CL[iL * Nb + j] = coeffs[j] * Qval;                              // reduce this value to CL
                    iL++;
                    if (iL > L) throw std::runtime_error(" Q shoule be sparsed!");
                }
                Xnj[idx] = coeffs[j] * Q[ibath * Dimension::FF + i];  // merge coeffs into Q to obtain Xnj
            }
        }
    }

    // model field
    mass = DS->def<kids_real>("model.mass", Dimension::N);
    for (int j = 0; j < Dimension::N; ++j) mass[j] = 1.0f;

    vpes = DS->def<kids_real>("model.vpes", Dimension::P);
    grad = DS->def<kids_real>("model.grad", Dimension::PN);
    hess = DS->def<kids_real>("model.hess", Dimension::PNN);
    V    = DS->def<kids_real>("model.V", Dimension::PFF);
    dV   = DS->def<kids_real>("model.dV", Dimension::PNFF);
    // ddV  = DS->def<kids_real>("model.ddV", Dimension::NNFF);

    // init & integrator
    x = DS->def<kids_real>("integrator.x", Dimension::PN);
    p = DS->def<kids_real>("integrator.p", Dimension::PN);

    // ARRAY_SHOW(Hsys, Dimension::F, Dimension::F);
    // ARRAY_SHOW(omegas, 1, Nb);
    // ARRAY_SHOW(coeffs, 1, Nb);
    // ARRAY_SHOW(x_sigma, 1, Nb);
    // ARRAY_SHOW(p_sigma, 1, Nb);
    // exit(-1);
}

void Model_SystemBath::init_calc_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* x = this->x + iP * Dimension::N;
        kids_real* p = this->p + iP * Dimension::N;

        Kernel_Random::rand_gaussian(x, Dimension::N);
        Kernel_Random::rand_gaussian(p, Dimension::N);
        for (int ibath = 0, idxR = 0; ibath < nbath; ++ibath) {
            for (int j = 0; j < Nb; ++j, ++idxR) {
                x[idxR] = x[idxR] * x_sigma[j];
                p[idxR] = p[idxR] * p_sigma[j];
            }
        }

        if (nsamp_type == NSampPolicy::QCT) {
            for (int ibath = 0, idxR = 0; ibath < nbath; ++ibath) {
                for (int j = 0; j < Nb; ++j, ++idxR) {
                    double Eq    = omegas[j] * (0 + 0.5);
                    double Ec    = 0.5 * p[idxR] * p[idxR] + 0.5 * omegas[j] * omegas[j] * x[idxR] * x[idxR];
                    double scale = sqrt(Eq / Ec);
                    x[idxR]      = x[idxR] * scale;
                    p[idxR]      = p[idxR] * scale;
                }
            }
        }
    }

    _DataSet->def("init.x", x, Dimension::PN);
    _DataSet->def("init.p", p, Dimension::PN);
    exec_kernel(stat);
}

int Model_SystemBath::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* x    = this->x + iP * Dimension::N;
        kids_real* vpes = this->vpes + iP;
        kids_real* grad = this->grad + iP * Dimension::N;
        kids_real* hess = this->hess + iP * Dimension::NN;
        kids_real* V    = this->V + iP * Dimension::FF;
        kids_real* dV   = this->dV + iP * Dimension::NFF;

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
        for (int i = 0; i < Dimension::FF; ++i) V[i] = Hsys[i];
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
                for (int ibath = 0, idxV = 0, idxR = 0; ibath < nbath; ++ibath, idxV += Dimension::Fadd1) {
                    for (int j = 0; j < Nb; ++j, ++idxR) { V[idxV] += coeffs[j] * x[idxR]; }
                }
                break;
            }
            default: {
                for (int ibath = 0, idxR = 0, idxQ0 = 0; ibath < nbath; ++ibath, idxQ0 += Dimension::FF) {
                    for (int j = 0; j < Nb; ++j, ++idxR) {
                        double cxj = coeffs[j] * x[idxR];
                        for (int i = 0, idxQ = idxQ0; i < Dimension::FF; ++i, ++idxQ) { V[i] += Q[idxQ] * cxj; }
                    }
                }
                break;
            }
        }

        if (count_exec == 0) {  // only calculate once time
            ARRAY_CLEAR(dV, Dimension::NFF);
            for (int ibath = 0, idxQ0 = 0, idxdV = 0; ibath < nbath; ++ibath, idxQ0 += Dimension::FF) {
                for (int j = 0; j < Nb; ++j) {
                    for (int i = 0, idxQ = idxQ0; i < Dimension::FF; ++i, ++idxQ, ++idxdV) {
                        dV[idxdV] = Q[idxQ] * coeffs[j];
                    }
                }
            }
            // ARRAY_SHOW(dV, Dimension::N, Dimension::FF);
            // exit(-1);
        }

        // if (flag < 2) return 0;

        if (count_exec == 0) {
            ARRAY_CLEAR(hess, Dimension::NN);
            for (int j = 0, idx = 0, Nadd1 = Dimension::N + 1; j < Dimension::N; ++j, idx += Nadd1) {
                hess[idx] = omegas[j % Nb] * omegas[j % Nb];
            }
            // for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
        }
    }
    return 0;
}

};  // namespace PROJECT_NS
