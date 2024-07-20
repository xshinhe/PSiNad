#include "kids/Model_SystemBath.h"

#include "kids/Kernel_NAForce.h"
#include "kids/Kernel_Random.h"
#include "kids/hamiltonian_data.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Model_SystemBath::getName() { return "Model_SystemBath"; }

int Model_SystemBath::getType() const { return utils::hash(FUNCTION_NAME); }

void Model_SystemBath::setInputParam_impl(std::shared_ptr<Param> PM) {
    int N = _param->get_int({"model.N"}, LOC());
    Nb    = _param->get_int({"model.Nb", "model.bath.Nb"}, LOC());
    nbath = _param->get_int({"model.nbath", "model.bath.nbath"}, LOC());
    // assert(N == Nb * nbath);

    FORCE_OPT::Nb    = Nb;
    FORCE_OPT::nbath = nbath;
    Dimension::Nb    = Nb;
    Dimension::nbath = nbath;
    Dimension::shape_Nb.static_build();       // for building "omegas, coeffs, x_sigma, p_sigma"
    Dimension::shape_nbathFF.static_build();  // for building "Q"

    system_type   = SystemPolicy::_from(_param->get_string({"model.system_flag"}, LOC(), "SB"));
    coupling_type = CouplingPolicy::_from(_param->get_string({"model.coupling_flag"}, LOC(), "SB"));

    // @deprecated moved from initializeKernel_impl
    nsamp_type = NSampPolicy::_from(_param->get_string({"model.nsamp_flag"}, LOC(), "Wigner"));
}

void Model_SystemBath::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    /// 1) System
    Hsys = DS->def(DATA::model::Hsys);
    L    = 1;
    switch (system_type) {
        case SystemPolicy::SB: {
            assert(nbath == 1);
            assert(Dimension::F == 2);
            L             = 2;  // sigma_z
            double bias   = _param->get_real({"model.bias"}, LOC(), phys::energy_d, 1.0f);
            double delta  = _param->get_real({"model.delta"}, LOC(), phys::energy_d, 1.0f);
            double HSB[4] = {bias, delta, delta, -bias};
            for (int i = 0; i < Dimension::FF; ++i) Hsys[i] = HSB[i];
            break;
        }
        case SystemPolicy::AGG: {
            assert(Dimension::F >= 2);
            L            = 1;  // |n><n|
            double delta = _param->get_real({"model.delta"}, LOC(), phys::energy_d, 1.0f);
            for (int i = 0, ik = 0; i < Dimension::F; ++i) {
                for (int k = 0; k < Dimension::F; ++k, ++ik) Hsys[ik] = (i - k == 1 || k - i == 1) ? delta : 0.0f;
            }
            bool cyclic = _param->get_bool({"model.agg_cyclic"}, LOC(), false);
            if (cyclic) {
                Hsys[0 * Dimension::F + (Dimension::F - 1)] = delta;
                Hsys[(Dimension::F - 1) * Dimension::F + 0] = delta;
            }
            break;
        }
        case SystemPolicy::SF3a:
        case SystemPolicy::SF3b:
        case SystemPolicy::SF5a:
        case SystemPolicy::SF5b:
        case SystemPolicy::FMO:
        case SystemPolicy::FCP: {
            L                             = 1;  // |n><n|
            auto [Fcheck, Hdata, UnitStr] = Hsys_Dict.at(_param->get_string({"model.system_flag"}, LOC()));
            assert(Dimension::nbath == Fcheck);
            assert(Dimension::F == Fcheck || Dimension::F == Fcheck + 1);  // last is the ground state
            double data_unit = phys::au::as(phys::energy_d, UnitStr);
            for (int i = 0, ik = 0; i < Fcheck; ++i) {
                for (int k = 0; k < Fcheck; ++k, ++ik) { Hsys[i * Dimension::F + k] = Hdata[ik] / data_unit; }
            }
            break;
        }
        case SystemPolicy::Read: {
            std::string   system_file = _param->get_string({"model.system_file"}, LOC(), "system.dat");
            std::ifstream ifs(system_file);
            std::string   data_unit_str;
            std::string   firstline;
            getline(ifs, firstline);
            std::stringstream sstr(firstline);
            sstr >> data_unit_str;  ///< the firstline stores H's unit
            double    data_unit = phys::us::conv(phys::au::unit, phys::us::parse(data_unit_str));
            kids_real val;
            for (int i = 0; i < Dimension::FF; ++i)
                if (ifs >> val) Hsys[i] = val / data_unit;
            ifs.close();
        }
    }
    Dimension::L = L;
    Dimension::shape_LNb.static_build();       // for building "CL"
    Dimension::shape_LnbathFF.static_build();  // for building "QL"

    /// 2) init Bath sub-kernel (declaration & call)
    for (auto pkernel : _child_kernels) pkernel->setInputDataSet(DS);
    omegas  = DS->def(DATA::model::bath::omegas);
    coeffs  = DS->def(DATA::model::bath::coeffs);
    x_sigma = DS->def(DATA::model::bath::x_sigma);
    p_sigma = DS->def(DATA::model::bath::p_sigma);

    /// 3) bilinear Coupling (saving order: L, nbath, Nb, FF)
    Q   = DS->def(DATA::model::coupling::Q);
    CL  = DS->def(DATA::model::coupling::CL);
    QL  = DS->def(DATA::model::coupling::QL);
    Xnj = DS->def(DATA::model::coupling::Xnj);
    switch (coupling_type) {
        case CouplingPolicy::SB: {
            Q[0] = 1.0f, Q[1] = 0.0f, Q[2] = 0.0f, Q[3] = -1.0f;
            break;
        }
        case CouplingPolicy::SE: {
            assert(Dimension::F == Dimension::nbath || Dimension::F == Dimension::nbath + 1);
            for (int i = 0, idx = 0; i < Dimension::nbath; ++i) {
                for (int j = 0; j < Dimension::F; ++j) {
                    for (int k = 0; k < Dimension::F; ++k, ++idx) Q[idx] = (i == j && i == k) ? 1.0f : 0.0f;
                }
            }
            break;
        }
        default: {
            std::string coupling_file =
                _param->get_string({"model.coupling_file", "model.coupling.file"}, LOC(), "coupling.dat");
            std::ifstream ifs(coupling_file);
            kids_real     tmp;
            for (int i = 0; i < Dimension::nbath * Dimension::FF; ++i)
                if (ifs >> tmp) Q[i] = tmp;
            ifs.close();
        }
    }

    for (int ibath = 0, idx = 0; ibath < Dimension::nbath; ++ibath) {
        for (int j = 0; j < Dimension::Nb; ++j) {
            for (int ik = 0, iL = 0; ik < Dimension::FF; ++ik, ++idx) {
                double Qval = Q[ibath * Dimension::FF + ik];
                if (Qval != 0) {
                    QL[iL * Dimension::nbath * Dimension::FF + ibath * Dimension::FF + ik] = 1;
                    CL[iL * Nb + j] = coeffs[j] * Qval;  // reduce this value to CL
                    iL++;
                    if (iL > L) throw kids_error(" Q shoule be more sparsed; please enlarge L!");
                }
                Xnj[idx] = coeffs[j] * Q[ibath * Dimension::FF + ik];  // merge coeffs into Q to obtain Xnj
            }
        }
    }

    // model field
    mass = DS->def(DATA::model::mass);
    for (int j = 0; j < Dimension::N; ++j) mass[j] = 1.0f;
    vpes = DS->def(DATA::model::vpes);
    grad = DS->def(DATA::model::grad);
    hess = DS->def(DATA::model::hess);
    V    = DS->def(DATA::model::V);
    dV   = DS->def(DATA::model::dV);
    // ddV  = DS->def(DATA::model::ddV);
    x = DS->def(DATA::integrator::x);
    p = DS->def(DATA::integrator::p);
}

Status& Model_SystemBath::initializeKernel_impl(Status& stat) {
    executeKernel(stat);
    return stat;  // @todo

    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* x = this->x + iP * Dimension::N;
        kids_real* p = this->p + iP * Dimension::N;

        Kernel_Random::rand_gaussian(x, Dimension::N);
        Kernel_Random::rand_gaussian(p, Dimension::N);
        for (int ibath = 0, idxR = 0; ibath < Dimension::nbath; ++ibath) {
            for (int j = 0; j < Dimension::Nb; ++j, ++idxR) {
                x[idxR] = x[idxR] * x_sigma[j];
                p[idxR] = p[idxR] * p_sigma[j];
            }
        }

        if (nsamp_type == NSampPolicy::QCT) {
            for (int ibath = 0, idxR = 0; ibath < Dimension::nbath; ++ibath) {
                for (int j = 0; j < Dimension::Nb; ++j, ++idxR) {
                    double Eq    = omegas[j] * (0 + 0.5);
                    double Ec    = 0.5 * p[idxR] * p[idxR] + 0.5 * omegas[j] * omegas[j] * x[idxR] * x[idxR];
                    double scale = sqrt(Eq / Ec);
                    x[idxR]      = x[idxR] * scale;
                    p[idxR]      = p[idxR] * scale;
                }
            }
        }
    }

    _dataset->def_real("init.x", x, Dimension::PN);
    _dataset->def_real("init.p", p, Dimension::PN);
    executeKernel(stat);
    return stat;
}

Status& Model_SystemBath::executeKernel_impl(Status& stat) {
    for (int iP = 0; iP < Dimension::P_NOW; ++iP) {
        ///< get one instance from trajectory swarms
        kids_real* x    = this->x + iP * Dimension::N;
        kids_real* vpes = this->vpes + iP;
        kids_real* grad = this->grad + iP * Dimension::N;
        kids_real* hess = this->hess + iP * Dimension::NN;
        kids_real* V    = this->V + iP * Dimension::FF;
        kids_real* dV   = this->dV + iP * Dimension::NFF;

        // calculate nuclear vpes and grad (assume mass = 1)
        if (true) {
            double term = 0.0f;
            for (int ibath = 0, idxR = 0; ibath < Dimension::nbath; ++ibath) {
                for (int j = 0; j < Dimension::Nb; ++j, ++idxR) {
                    term += omegas[j] * omegas[j] * x[idxR] * x[idxR];
                    grad[idxR] = omegas[j] * omegas[j] * x[idxR];
                }
            }
            vpes[0] = 0.5 * term;
        } else {
            vpes[0] = ARRAY_INNER_VMV_TRANS1(x, Kmat, x, Dimension::N, Dimension::N);
            vpes[0] *= 0.5e0;
            ARRAY_MATMUL(grad, Kmat, x, Dimension::N, Dimension::N, 1);
        }

        // calculate electronic V and dV
        for (int i = 0; i < Dimension::FF; ++i) V[i] = Hsys[i];
        switch (coupling_type) {
            case CouplingPolicy::SB: {
                double term = 0;
                for (int ibath = 0, idxR = 0; ibath < Dimension::nbath; ++ibath) {
                    for (int j = 0; j < Dimension::Nb; ++j, ++idxR) { term += coeffs[j] * x[idxR]; }
                }
                V[0] += term;
                V[3] -= term;
                break;
            }
            case CouplingPolicy::SE: {
                for (int ibath = 0, idxV = 0, idxR = 0; ibath < Dimension::nbath; ++ibath, idxV += Dimension::Fadd1) {
                    for (int j = 0; j < Dimension::Nb; ++j, ++idxR) { V[idxV] += coeffs[j] * x[idxR]; }
                }
                break;
            }
            default: {
                for (int ibath = 0, idxR = 0, idxQ0 = 0; ibath < Dimension::nbath; ++ibath, idxQ0 += Dimension::FF) {
                    for (int j = 0; j < Dimension::Nb; ++j, ++idxR) {
                        double cxj = coeffs[j] * x[idxR];
                        for (int i = 0, idxQ = idxQ0; i < Dimension::FF; ++i, ++idxQ) { V[i] += Q[idxQ] * cxj; }
                    }
                }
                break;
            }
        }

        if (count_exec == 0) {  // only calculate once time
            ARRAY_CLEAR(dV, Dimension::NFF);
            for (int ibath = 0, idxQ0 = 0, idxdV = 0; ibath < Dimension::nbath; ++ibath, idxQ0 += Dimension::FF) {
                for (int j = 0; j < Dimension::Nb; ++j) {
                    for (int i = 0, idxQ = idxQ0; i < Dimension::FF; ++i, ++idxQ, ++idxdV) {
                        dV[idxdV] = Q[idxQ] * coeffs[j];
                    }
                }
            }
        }

        // if (flag < 2) return 0;

        if (count_exec == 0) {
            ARRAY_CLEAR(hess, Dimension::NN);
            for (int j = 0, idx = 0, Nadd1 = Dimension::N + 1; j < Dimension::N; ++j, idx += Nadd1) {
                hess[idx] = omegas[j % Dimension::Nb] * omegas[j % Dimension::Nb];
            }
            // for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
        }
    }
    return stat;
}

};  // namespace PROJECT_NS
