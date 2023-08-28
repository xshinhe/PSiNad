#include "Kernel_NADForce.h"

#include "../core/linalg.h"
#include "Kernel_Declare.h"
#include "Kernel_Elec.h"
#include "Kernel_Representation.h"
#include "Kernel_Update.h"

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

void Kernel_NADForce::read_param_impl(Param* PM) {
    FORCE_OPT::BATH_FORCE_BILINEAR = _Param->get<bool>("BATH_FORCE_BILINEAR", LOC(), false);
};

void Kernel_NADForce::init_data_impl(DataSet* DS) {
    f    = DS->reg<double>("integrator.f", Dimension::PN);
    grad = DS->reg<double>("model.grad", Dimension::PN);
    dV   = DS->reg<double>("model.dV", Dimension::PNFF);
    dE   = DS->reg<double>("model.rep.dE", Dimension::PNFF);
    T    = DS->reg<double>("model.rep.T", Dimension::PFF);

    switch (Kernel_Representation::nuc_repr_type) {
        case RepresentationPolicy::Diabatic:
            Force = dV;
            break;
        case RepresentationPolicy::Adiabatic:
            Force = dE;
            break;
    }
}

void Kernel_NADForce::init_calc_impl(int stat) { exec_kernel(stat); }

int Kernel_NADForce::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        int* occ_nuc         = Kernel_Elec::occ_nuc + iP;
        num_complex* rho_nuc = Kernel_Elec::rho_nuc + iP * Dimension::FF;

        num_real* f     = this->f + iP * Dimension::N;
        num_real* grad  = this->grad + iP * Dimension::N;
        num_real* Force = this->Force + iP * Dimension::NFF;
        num_real* T     = this->T + iP * Dimension::FF;

        /////////////////////////////////////////////////////////////////
        switch (NADForce_type) {
            case NADForcePolicy::BO: {
                for (int j = 0, idxdV0 = 0; j < Dimension::N; ++j, idxdV0 += Dimension::FF)
                    f[j] = Force[j * Dimension::FF + (*occ_nuc) * Dimension::Fadd1];
                break;
            }
            case NADForcePolicy::EHR: {
                Kernel_Representation::transform(rho_nuc, T, Dimension::F, Kernel_Representation::inp_repr_type,
                                                 Kernel_Representation::nuc_repr_type, SpacePolicy::L);
                if (FORCE_OPT::BATH_FORCE_BILINEAR) {  // for both dV & dE (only for FMO-like model)
                    int& B  = FORCE_OPT::nbath;
                    int& J  = FORCE_OPT::Nb;
                    int JFF = J * Dimension::FF;
                    for (int b = 0, bj = 0, b0FF = 0, b0bb = 0; b < B;
                         ++b, b0FF += JFF, b0bb += (JFF + Dimension::Fadd1)) {
                        double* Forceb0 = Force + b0FF;
                        double fb0 = std::real(ARRAY_TRACE2(Kernel_Elec::rho_nuc, Forceb0, Dimension::F, Dimension::F));
                        for (int j = 0, bjbb = b0bb; j < J; ++j, ++bj, bjbb += Dimension::FF) {
                            f[bj] = fb0 * Force[bjbb] / Force[b0bb];
                        }
                    }
                } else {
                    for (int j = 0, jFF = 0; j < Dimension::N; ++j, jFF += Dimension::FF) {
                        double* dVj = Force + jFF;
                        f[j]        = std::real(ARRAY_TRACE2(rho_nuc, dVj, Dimension::F, Dimension::F));
                    }
                }
                Kernel_Representation::transform(rho_nuc, T, Dimension::F, Kernel_Representation::nuc_repr_type,
                                                 Kernel_Representation::inp_repr_type,
                                                 SpacePolicy::L);  // not need
                break;
            }
        }
        for (int j = 0; j < Dimension::N; ++j) f[j] += grad[j];
    }
    return 0;
}

NADForcePolicy::_type Kernel_NADForce::NADForce_type = NADForcePolicy::EHR;

namespace FORCE_OPT {
int nbath                = 1;
int Nb                   = 1;
bool BATH_FORCE_BILINEAR = false;
};  // namespace FORCE_OPT

};  // namespace PROJECT_NS
