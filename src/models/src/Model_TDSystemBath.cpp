#include "kids/Model_TDSystemBath.h"

#include <unistd.h>

#include "kids/Kernel_NAForce.h"
#include "kids/Kernel_Random.h"
#include "kids/hamiltonian_data.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Model_TDSystemBath::getName() { return "Model_TDSystemBath"; }

int Model_TDSystemBath::getType() const { return utils::hash(FUNCTION_NAME); }

void Model_TDSystemBath::setInputParam_impl(std::shared_ptr<Param> PM) {
    int N           = _param->get_int({"model.N"}, LOC());
    Nb              = _param->get_int({"model.Nb", "model.bath.Nb"}, LOC());
    nbath           = _param->get_int({"model.nbath", "model.bath.nbath"}, LOC());
    system_type     = SystemPolicy::_from(_param->get_string({"model.system_flag"}, LOC(), "SB"));
    coupling_type   = CouplingPolicy::_from(_param->get_string({"model.coupling_flag"}, LOC(), "SB"));
    is_et_transform =  //
        _param->get_bool({"model.bath_et_transform", "model.bath.et_transform"}, LOC(), false);

    perx = _param->get_real({"model.perx"}, LOC(), 1.0);
    pery = 1.0e0 - perx;
    freqd = _param->get_real({"model.freqd"}, LOC(), phys::energy_d, 1.0);

    if (nbath <= 0) {
        system_type      = SystemPolicy::Read;
        coupling_type    = CouplingPolicy::Read;
        FORCE_OPT::Nb    = N;
        Dimension::Nb    = N;
        FORCE_OPT::nbath = 1;
        Dimension::nbath = 1;
    } else {
        kids_assert(N == Nb * nbath, "Dimension Error");
        FORCE_OPT::Nb    = Nb;
        Dimension::Nb    = Nb;
        FORCE_OPT::nbath = nbath;
        Dimension::nbath = nbath;
    }
    Dimension::shape_Nb.static_build();       // for building "omegas, coeffs, x_sigma, p_sigma"
    Dimension::shape_nbathFF.static_build();  // for building "Q"
}

void Model_TDSystemBath::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    /// 1) System
    Hsys = DS->def(DATA::model::Hsys);
    L    = 1;
    switch (system_type) {
        case SystemPolicy::SB: {
            kids_assert(nbath == 1, "Dimension Error");
            kids_assert(Dimension::F == 2, "Dimension Error");
            L             = 2;  // sigma_z
            double bias   = _param->get_real({"model.bias"}, LOC(), phys::energy_d, 1.0f);
            double delta  = _param->get_real({"model.delta"}, LOC(), phys::energy_d, 1.0f);
            double HSB[4] = {bias, delta, delta, -bias};
            for (int i = 0; i < Dimension::FF; ++i) Hsys[i] = HSB[i];
            break;
        }
        case SystemPolicy::AGG: {
            kids_assert(Dimension::F >= 2, "Dimension Error");
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
            kids_assert(Dimension::nbath == Fcheck, "Dimension Error");
            kids_assert(Dimension::F == Fcheck || Dimension::F == Fcheck + 1,
                        "Dimension Error");  // last is the ground state
            double data_unit = 1.0e0 / phys::au::as(phys::energy_d, UnitStr);
            for (int i = 0, ik = 0; i < Fcheck; ++i) {
                for (int k = 0; k < Fcheck; ++k, ++ik) { Hsys[i * Dimension::F + k] = Hdata[ik] / data_unit; }
            }
            break;
        }
        case SystemPolicy::Read: {
            try {
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
            } catch (std::runtime_error& e) { throw kids_error("read system.dat fails"); }
        }
    }
    Dimension::L = L;
    Dimension::shape_LNb.static_build();       // for building "CL"
    Dimension::shape_LnbathFF.static_build();  // for building "QL"

    /// 2) init Bath sub-kernel (declaration & call)
    Kmat = DS->def(DATA::model::Kmat);
    for (auto pkernel : _child_kernels) pkernel->setInputDataSet(DS);
    omegas  = DS->def(DATA::model::bath::omegas);
    coeffs  = DS->def(DATA::model::bath::coeffs);
    x_sigma = DS->def(DATA::model::x_sigma);
    p_sigma = DS->def(DATA::model::p_sigma);
    mass    = DS->def(DATA::model::mass);
    for (int j = 0; j < Dimension::N; ++j) mass[j] = 1.0e0;

    /// 3) bilinear Coupling (saving order: L, nbath, Nb, FF)
    Qmat = DS->def(DATA::model::Qmat);
    Q    = DS->def(DATA::model::coupling::Q);
    CL   = DS->def(DATA::model::coupling::CL);
    QL   = DS->def(DATA::model::coupling::QL);
    if (nbath <= 0) {
        try {
            std::string coupling_file =
                _param->get_string({"model.coupling_file", "model.coupling.file"}, LOC(), "coupling.dat");
            std::ifstream ifs(coupling_file);
            kids_real     tmp;
            for (int i = 0; i < Dimension::NFF; ++i)
                if (ifs >> tmp) Qmat[i] = tmp;
            ifs.close();
        } catch (std::runtime_error& e) { throw kids_error("read Qmat from coupling.dat fails"); }
        // keep Q as zero!
    } else {
        // initialize Q, then Qmat (assuming coeffs is prepared)
        switch (coupling_type) {
            case CouplingPolicy::SB: {
                kids_assert(Dimension::F == 2, "Dimension Error");
                Q[0] = 1.0f, Q[1] = 0.0f, Q[2] = 0.0f, Q[3] = -1.0f;
                break;
            }
            case CouplingPolicy::SE: {
                kids_assert(Dimension::F == Dimension::nbath || Dimension::F == Dimension::nbath + 1,
                            "Dimension Error");
                for (int i = 0, idx = 0; i < Dimension::nbath; ++i) {
                    for (int j = 0; j < Dimension::F; ++j) {
                        for (int k = 0; k < Dimension::F; ++k, ++idx) Q[idx] = (i == j && i == k) ? 1.0f : 0.0f;
                    }
                }
                break;
            }
            case CouplingPolicy::Read:
            default: {
                try {
                    std::string coupling_file =
                        _param->get_string({"model.coupling_file", "model.coupling.file"}, LOC(), "coupling.dat");
                    std::ifstream ifs(coupling_file);
                    kids_real     tmp;
                    for (int i = 0; i < Dimension::nbath * Dimension::FF; ++i)
                        if (ifs >> tmp) Q[i] = tmp;
                    ifs.close();
                } catch (std::runtime_error& e) { throw kids_error("read Q from coupling.dat fails"); }
            }
        }
        // check et transfrom
        if (is_et_transform) {  // set Q=0 except for the first
            for (int i = Dimension::FF; i < Dimension::nbath * Dimension::FF; ++i) Q[i] = 0.0e0;
        }
        // decompose Q and calculate Qmat
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
                    Qmat[idx] = coeffs[j] * Q[ibath * Dimension::FF + ik];  // merge coeffs into Q to obtain Xnj
                }
            }
        }
    }
    // model field
    vpes = DS->def(DATA::model::vpes);
    grad = DS->def(DATA::model::grad);
    hess = DS->def(DATA::model::hess);
    V    = DS->def(DATA::model::V);
    dV   = DS->def(DATA::model::dV);
    // ddV  = DS->def(DATA::model::ddV);
    x = DS->def(DATA::integrator::x);
    p = DS->def(DATA::integrator::p);
    t = DS->def(DATA::flowcontrol::t);
}

Status& Model_TDSystemBath::initializeKernel_impl(Status& stat) {
    // executeKernel(stat);
    return stat;
}

Status& Model_TDSystemBath::executeKernel_impl(Status& stat) {
    if (stat.frozen) return stat;

    for (int iP = 0; iP < Dimension::P_NOW; ++iP) {
        ///< get one instance from trajectory swarms
        auto x    = this->x.subspan(iP * Dimension::N, Dimension::N);
        auto vpes = this->vpes.subspan(iP, 1);
        auto grad = this->grad.subspan(iP * Dimension::N, Dimension::N);
        auto hess = this->hess.subspan(iP * Dimension::NN, Dimension::NN);
        auto V    = this->V.subspan(iP * Dimension::FF, Dimension::FF);
        auto dV   = this->dV.subspan(iP * Dimension::NFF, Dimension::NFF);

        if (nbath <= 0) {
            // H = 0.5*p*p + Hsys + 0.5*x*Kmat*x + x*Qmat
            vpes[0] = ARRAY_INNER_VMV_TRANS1(x.data(), Kmat.data(), x.data(), Dimension::N, Dimension::N);
            vpes[0] *= 0.5e0;
            ARRAY_MATMUL(grad.data(), Kmat.data(), x.data(), Dimension::N, Dimension::N, 1);
        } else {
            double term = 0.0f;
            for (int ibath = 0, idxR = 0; ibath < Dimension::nbath; ++ibath) {
                for (int j = 0; j < Dimension::Nb; ++j, ++idxR) {
                    term += omegas[j] * omegas[j] * x[idxR] * x[idxR];
                    grad[idxR] = omegas[j] * omegas[j] * x[idxR];
                }
            }
            vpes[0] = 0.5 * term;
            if (is_et_transform) {
                for (int j = 1; j < Dimension::Nb; ++j) {
                    vpes[0] += coeffs[j] * x[0] * x[j];
                    grad[0] += coeffs[j] * x[j];
                    grad[j] += coeffs[j] * x[0];
                }
            }
        }

        // calculate electronic V and dV
        if (nbath <= 0) {
            // H = 0.5*p*p + Hsys + 0.5*x*Kmat*x + x*Qmat
            ARRAY_MATMUL(V.data(), x.data(), Qmat.data(), 1, Dimension::N, Dimension::FF);
            for (int i = 0; i < Dimension::FF; ++i) V[i] += Hsys[i];

            if (count_exec == 0) {  // only calculate once time
                ARRAY_CLEAR(dV.data(), Dimension::NFF);
                for (int idxdV = 0; idxdV < Dimension::NFF; ++idxdV) dV[idxdV] = Qmat[idxdV];
            }
        } else {
            for (int i = 0; i < Dimension::FF; ++i) V[i] = (i%(Dimension::F+1) == 0)? 
                Hsys[i] * (perx + pery * std::cos(freqd * t[0])) : 
                Hsys[i];

            // std::cout << perx << ", " << pery << std::endl;
            // std::cout << FMT(4) << t[0] << std::endl;
            // for(int i=0; i< Dimension::FF; ++i) std::cout << FMT(4) << V[i];
            // std::cout << std::endl;
            // std::cout << std::endl;

            if (is_et_transform) {  // for optimize
                switch (coupling_type) {
                    case CouplingPolicy::SB: {
                        double term = 0;
                        for (int ibath = 0, idxR = 0; ibath < Dimension::nbath; ++ibath) {
                            term += coeffs[0] * x[ibath * Dimension::Nb];
                        }
                        V[0] += term;
                        V[3] -= term;
                        break;
                    }
                    case CouplingPolicy::SE: {
                        for (int ibath = 0, idxV = 0, idxR = 0; ibath < Dimension::nbath;
                             ++ibath, idxV += Dimension::Fadd1) {
                            V[idxV] += coeffs[0] * x[ibath * Dimension::Nb];
                        }
                        break;
                    }
                    default: {  // TODO OPT
                        for (int ibath = 0, idxR = 0, idxQ0 = 0; ibath < Dimension::nbath;
                             ++ibath, idxQ0 += Dimension::FF) {
                            for (int j = 0; j < Dimension::Nb; ++j, ++idxR) {
                                double cxj = coeffs[j] * x[idxR];
                                for (int i = 0, idxQ = idxQ0; i < Dimension::FF; ++i, ++idxQ) { V[i] += Q[idxQ] * cxj; }
                            }
                        }
                        break;
                    }
                }
            } else {
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
                        for (int ibath = 0, idxV = 0, idxR = 0; ibath < Dimension::nbath;
                             ++ibath, idxV += Dimension::Fadd1) {
                            for (int j = 0; j < Dimension::Nb; ++j, ++idxR) { V[idxV] += coeffs[j] * x[idxR]; }
                        }
                        break;
                    }
                    default: {
                        for (int ibath = 0, idxR = 0, idxQ0 = 0; ibath < Dimension::nbath;
                             ++ibath, idxQ0 += Dimension::FF) {
                            for (int j = 0; j < Dimension::Nb; ++j, ++idxR) {
                                double cxj = coeffs[j] * x[idxR];
                                for (int i = 0, idxQ = idxQ0; i < Dimension::FF; ++i, ++idxQ) { V[i] += Q[idxQ] * cxj; }
                            }
                        }
                        break;
                    }
                }
            }
            if (count_exec == 0) {  // only calculate once time
                ARRAY_CLEAR(dV.data(), Dimension::NFF);
                for (int ibath = 0, idxQ0 = 0, idxdV = 0; ibath < Dimension::nbath; ++ibath, idxQ0 += Dimension::FF) {
                    for (int j = 0; j < Dimension::Nb; ++j) {
                        for (int i = 0, idxQ = idxQ0; i < Dimension::FF; ++i, ++idxQ, ++idxdV) {
                            dV[idxdV] = Q[idxQ] * coeffs[j];
                        }
                    }
                }
            }
        }

        // if (flag < 2) return 0;
        if (count_exec == 0) {
            for (int i = 0; i < Dimension::NN; ++i) hess[i] = Kmat[i];
            // for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
        }
    }
    return stat;
}

};  // namespace PROJECT_NS
