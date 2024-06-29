#include "kids/Sampling_Elec.h"

#include "kids/Kernel_Elec_Utils.h"
#include "kids/Kernel_NAForce.h"
#include "kids/Kernel_Random.h"
#include "kids/Kernel_Representation.h"
#include "kids/debug_utils.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

inline bool isFileExists(const std::string& name) { return std::ifstream{name.c_str()}.good(); }

const std::string Sampling_Elec::getName() { return "Sampling_Elec"; }

int Sampling_Elec::getType() const { return utils::hash(FUNCTION_NAME); }

void Sampling_Elec::setInputParam_impl(std::shared_ptr<Param> PM) {
    // for c & rho_ele
    sampling_type = ElectronicSamplingPolicy::_from(
        _param->get_string({"solver.sampling_ele_flag", "solver.sampling.ele_flag"}, LOC(), "Constraint"));

    occ0 = _param->get_int({"model.occ", "solver.occ"}, LOC(), -1);
    if (occ0 < 0) throw std::runtime_error("occ < 0");
    if (occ0 >= Dimension::F) throw std::runtime_error("occ >= F");

    gamma1 = _param->get_real({"solver.gamma"}, LOC(), elec_utils::gamma_wigner(Dimension::F));
    if (gamma1 < -1.5) gamma1 = elec_utils::gamma_opt(Dimension::F);
    if (gamma1 < -0.5) gamma1 = elec_utils::gamma_wigner(Dimension::F);
    xi1 = (1 + Dimension::F * gamma1);

    // for rho_nuc
    use_cv   = _param->get_bool({"solver.use_cv"}, LOC(), false);
    use_wmm  = _param->get_bool({"solver.use_wmm"}, LOC(), false);
    use_sum  = _param->get_bool({"solver.use_sum"}, LOC(), false);
    use_fssh = _param->get_bool({"solver.use_fssh"}, LOC(), false);
}

void Sampling_Elec::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    c       = DS->def(DATA::integrator::c);
    rho_ele = DS->def(DATA::integrator::rho_ele);
    rho_nuc = DS->def(DATA::integrator::rho_nuc);
    occ_nuc = DS->def(DATA::integrator::occ_nuc);
    w       = DS->def(DATA::integrator::w);
    T       = DS->def(DATA::model::rep::T);
}

Status& Sampling_Elec::initializeKernel_impl(Status& stat) { return stat; }

Status& Sampling_Elec::executeKernel_impl(Status& stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {  // use P insread P_NOW
        kids_complex* c       = this->c + iP * Dimension::F;
        kids_complex* rho_ele = this->rho_ele + iP * Dimension::FF;
        kids_complex* rho_nuc = this->rho_nuc + iP * Dimension::FF;
        kids_complex* w       = this->w + iP;
        kids_real*    T       = this->T + iP * Dimension::FF;
        int*          occ_nuc = this->occ_nuc + iP;

        /////////////////////////////////////////////////////////////////

        int iocc;
        Kernel_Random::rand_catalog(&iocc, 1, true, 0, Dimension::F - 1);
        iocc = ((use_sum) ? iocc : occ0);
        w[0] = (use_sum) ? kids_complex(Dimension::F) : phys::math::iu;

        switch (sampling_type) {
            case ElectronicSamplingPolicy::Focus: {
                elec_utils::c_focus(c, xi1, gamma1, iocc, Dimension::F);
                elec_utils::ker_from_c(rho_ele, c, 1, 0, Dimension::F);
                elec_utils::ker_from_rho(rho_nuc, rho_ele, xi1, gamma1, Dimension::F, use_cv, iocc);
                break;
            }
            case ElectronicSamplingPolicy::GDTWA: {
                elec_utils::c_focus(c, xi1, gamma1, iocc, Dimension::F);  // @useless

                /// GDTWA sampling step 1: discrete random phase
                for (int j = 0; j < Dimension::F; ++j) {
                    if (j == iocc) {
                        rho_ele[j * Dimension::Fadd1] = 1.0e0;
                    } else {
                        double randu;
                        Kernel_Random::rand_uniform(&randu);
                        randu                            = phys::math::halfpi * (int(randu / 0.25f) + 0.5);
                        rho_ele[iocc * Dimension::F + j] = cos(randu) + phys::math::im * sin(randu);
                        rho_ele[j * Dimension::F + iocc] = std::conj(rho_ele[iocc * Dimension::F + j]);
                    }
                }
                for (int i = 0, ij = 0; i < Dimension::F; ++i) {
                    for (int j = 0; j < Dimension::F; ++j, ++ij) {
                        if (i == iocc || j == iocc) continue;
                        rho_ele[ij] = rho_ele[iocc * Dimension::F + j] / rho_ele[iocc * Dimension::F + i];
                    }
                }
                /// GDTWA sampling step 2: set off-diagonal
                double gamma_ou = phys::math::sqrthalf;
                double gamma_uu = 0.0e0;
                for (int i = 0, ij = 0; i < Dimension::F; ++i) {
                    for (int j = 0; j < Dimension::F; ++j, ++ij) {
                        if (i == j) {
                            rho_ele[ij] = (i == iocc) ? phys::math::iu : phys::math::iz;
                        } else if (i == iocc || j == iocc) {
                            rho_ele[ij] *= gamma_ou;
                        } else {
                            rho_ele[ij] *= gamma_uu;
                        }
                    }
                }
                elec_utils::ker_from_rho(rho_nuc, rho_ele, 1, 0, Dimension::F, false, iocc);
                break;
            }
            case ElectronicSamplingPolicy::SQCtri: {
                elec_utils::c_window(c, iocc, ElectronicSamplingPolicy::SQCtri, Dimension::F);
                elec_utils::ker_from_c(rho_ele, c, 1, 0, Dimension::F);
                elec_utils::ker_from_rho(rho_nuc, rho_ele, 1.0, gamma1, Dimension::F, use_cv, iocc);
                break;
            }
            case ElectronicSamplingPolicy::SQCspx: {
                elec_utils::c_sphere(c, Dimension::F);
                for (int i = 0; i < Dimension::F; ++i) c[i] = std::abs(c[i] * c[i]);
                c[iocc] += 1.0e0;
                for (int i = 0; i < Dimension::F; ++i) {
                    kids_real randu;
                    Kernel_Random::rand_uniform(&randu);
                    randu *= phys::math::twopi;
                    c[i] = sqrt(c[i]) * (cos(randu) + phys::math::im * sin(randu));
                }
                elec_utils::ker_from_c(rho_ele, c, 1, 0, Dimension::F);
                elec_utils::ker_from_rho(rho_nuc, rho_ele, 1.0, gamma1, Dimension::F, use_cv, iocc);
                break;
            }
            case ElectronicSamplingPolicy::SQCtest01: {
                elec_utils::c_window(c, iocc, ElectronicSamplingPolicy::SQCtri, Dimension::F);
                elec_utils::ker_from_c(rho_ele, c, 1, 0, Dimension::F);
                double norm = 0.0;
                for (int i = 0; i < Dimension::F; ++i) norm += std::abs(rho_ele[i * Dimension::Fadd1]);
                elec_utils::ker_from_rho(rho_nuc, rho_ele, xi1 / norm, gamma1, Dimension::F, use_cv, iocc);
                break;
            }
            case ElectronicSamplingPolicy::SQCtest02: {
                elec_utils::c_window(c, iocc, ElectronicSamplingPolicy::SQCtri, Dimension::F);
                double norm = 0.0e0;
                for (int i = 0; i < Dimension::F; ++i) norm += std::abs(c[i] * c[i]);
                xi1    = norm;
                gamma1 = (xi1 - 1.0e0) / Dimension::F;
                norm   = sqrt(norm);
                for (int i = 0; i < Dimension::F; ++i) c[i] /= norm;
                elec_utils::ker_from_c(rho_ele, c, 1, 0, Dimension::F);
                norm = 0.0;
                for (int i = 0; i < Dimension::F; ++i) norm += std::abs(rho_ele[i * Dimension::Fadd1]);
                double gmeff = (norm - 1.0) / Dimension::F;
                elec_utils::ker_from_rho(rho_nuc, rho_ele, 1.0, gmeff, Dimension::F, use_cv, iocc);
                break;
            }
            case ElectronicSamplingPolicy::Gaussian: {
                // elec_utils::c_gaussian(c, Dimension::F); /// @debug
                elec_utils::ker_from_c(rho_ele, c, 1, 0, Dimension::F);
                elec_utils::ker_from_rho(rho_nuc, rho_ele, xi1, gamma1, Dimension::F, use_cv, iocc);
                w[0] = kids_complex(Dimension::F);
                break;
            }
            case ElectronicSamplingPolicy::Constraint: {
                elec_utils::c_sphere(c, Dimension::F);
                elec_utils::ker_from_c(rho_ele, c, 1, 0, Dimension::F);
                elec_utils::ker_from_rho(rho_nuc, rho_ele, xi1, gamma1, Dimension::F, use_cv, iocc);
                w[0] = kids_complex(Dimension::F);
                break;
            }
            default: {
                std::string open_file = sampling_file;
                if (!isFileExists(sampling_file)) open_file = utils::concat(sampling_file, stat.icalc, ".ds");
                std::string   stmp, eachline;
                std::ifstream ifs(open_file);
                while (getline(ifs, eachline)) {
                    if (eachline.find("init.c") != eachline.npos) {
                        getline(ifs, eachline);
                        for (int i = 0; i < Dimension::F; ++i) ifs >> c[i];
                    }
                    // if (eachline.find("init.rho_ele") != eachline.npos) {
                    //     getline(ifs, eachline);
                    //     for (int i = 0; i < Dimension::FF; ++i) ifs >> rho_ele[i];
                    // }
                    // if (eachline.find("init.rho_nuc") != eachline.npos) {
                    //     getline(ifs, eachline);
                    //     for (int i = 0; i < Dimension::FF; ++i) ifs >> rho_nuc[i];
                    // }
                }
                elec_utils::ker_from_c(rho_ele, c, 1, 0, Dimension::F);  ///< initial rho_ele
                elec_utils::ker_from_rho(rho_nuc, rho_ele, xi1, gamma1, Dimension::F, use_cv, iocc);
                w[0] = phys::math::iu;
            }
        }

        // BO occupation in adiabatic representation
        Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                         Kernel_Representation::inp_repr_type,  //
                                         Kernel_Representation::nuc_repr_type,  //
                                         SpacePolicy::L);
        occ_nuc[0] = elec_utils::max_choose(rho_nuc);
        if (use_fssh) occ_nuc[0] = elec_utils::pop_choose(rho_nuc);
        Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                         Kernel_Representation::nuc_repr_type,  //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);

        // weight factor in tcf_repr /// moved to Kernel_Elec_Functions
        // Kernel_Representation::transform(rho_ele, T, Dimension::F,              //
        //                                  Kernel_Representation::inp_repr_type,  //
        //                                  RepresentationPolicy::Adiabatic,       //
        //                                  SpacePolicy::L);
        // wz_A[0]        = std::abs(rho_ele[0] - rho_ele[3]);
        // int    max_pop = elec_utils::max_choose(rho_ele);
        // double max_val = std::abs(rho_ele[max_pop * Dimension::Fadd1]);
        // ww_A[0]        = 4.0 - 1.0 / (max_val * max_val);
        // Kernel_Representation::transform(rho_ele, T, Dimension::F,         //
        //                                  RepresentationPolicy::Adiabatic,  //
        //                                  RepresentationPolicy::Diabatic,   //
        //                                  SpacePolicy::L);
        // wz_D[0] = std::abs(rho_ele[0] - rho_ele[3]);
        // max_pop = elec_utils::max_choose(rho_ele);
        // max_val = std::abs(rho_ele[max_pop * Dimension::Fadd1]);
        // ww_D[0] = 4.0 - 1.0 / (max_val * max_val);
        // Kernel_Representation::transform(rho_ele, T, Dimension::F,              //
        //                                  RepresentationPolicy::Diabatic,        //
        //                                  Kernel_Representation::inp_repr_type,  //
        //                                  SpacePolicy::L);
    }
    _dataset->def_complex("init.c", c, Dimension::PF);
    _dataset->def_complex("init.rho_ele", rho_ele, Dimension::PFF);
    _dataset->def_complex("init.rho_nuc", rho_nuc, Dimension::PFF);
    _dataset->def_real("init.T", T, Dimension::PFF);
    return stat;
}

// Status& Sampling_Elec::executeKernel_impl(Status& stat) { return stat; }

};  // namespace PROJECT_NS
