#include "Kernel_NADForce.h"

#include "Kernel_Dimension.h"
#include "Kernel_Elec.h"
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
    BATH_FORCE_OPT = _Param->get<bool>("BATH_FORCE_OPT", LOC(), false);
};

void Kernel_NADForce::init_data_impl(DataSet* DS) {
    f    = DS->reg<double>("integrator.f", Kernel_Dimension::N);
    grad = DS->reg<double>("model.grad", Kernel_Dimension::N);
    dV   = DS->reg<double>("model.dV", Kernel_Dimension::NFF);
    dE   = DS->reg<double>("model.rep.dE", Kernel_Dimension::NFF);
    T    = DS->reg<double>("model.rep.T", Kernel_Dimension::FF);

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
    switch (NADForce_type) {
        case NADForcePolicy::BO: {
            for (int j = 0, idxdV0 = 0; j < Kernel_Dimension::N; ++j, idxdV0 += Kernel_Dimension::FF)
                f[j] = Force[j * Kernel_Dimension::FF + (*Kernel_Elec::occ_nuc) * Kernel_Dimension::Fadd1];
            break;
        }
        case NADForcePolicy::EHR: {
            if (Kernel_Representation::ini_repr_type == RepresentationPolicy::Diabatic &&
                Kernel_Representation::nuc_repr_type == RepresentationPolicy::Adiabatic) {
                ARRAY_MATMUL3_TRANS1(Kernel_Elec::rho_nuc, T, Kernel_Elec::rho_nuc, T,  //
                                     Kernel_Dimension::F, Kernel_Dimension::F, Kernel_Dimension::F,
                                     Kernel_Dimension::F);
            }

            if (BATH_FORCE_OPT) {
                int nbath = Kernel_Dimension::F;
                int Nb    = Kernel_Dimension::N / nbath;
                int NbFF  = Nb * Kernel_Dimension::FF;
                for (int ibath = 0, ibj = 0, ib0FF = 0, ib0bb = 0; ibath < nbath;
                     ++ibath, ib0FF += NbFF, ib0bb = (NbFF + Kernel_Dimension::Fadd1)) {
                    double fib0 = 0.0f;
                    for (int i = 0, ii = 0, ib0ii = ib0FF; i < Kernel_Dimension::F;
                         ++i, ii += Kernel_Dimension::Fadd1, ib0ii += Kernel_Dimension::Fadd1) {
                        fib0 += std::real(Kernel_Elec::rho_nuc[ii]) * Force[ib0ii];
                    }
                    for (int j = 0, ibjbb = ib0bb; j < Nb; ++j, ++ibj, ibjbb += Kernel_Dimension::FF) {
                        f[ibj] = fib0 * Force[ibjbb] / Force[ib0bb];
                    }
                }
            } else {
                for (int j = 0, jFF = 0; j < Kernel_Dimension::N; ++j, jFF += Kernel_Dimension::FF) {
                    double* dVj = Force + jFF;
                    f[j] = std::real(ARRAY_TRACE2(Kernel_Elec::rho_nuc, dVj, Kernel_Dimension::F, Kernel_Dimension::F));
                }
            }

            if (Kernel_Representation::ini_repr_type == RepresentationPolicy::Diabatic &&
                Kernel_Representation::nuc_repr_type == RepresentationPolicy::Adiabatic) {
                ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_nuc, T, Kernel_Elec::rho_nuc, T,  //
                                     Kernel_Dimension::F, Kernel_Dimension::F, Kernel_Dimension::F,
                                     Kernel_Dimension::F);
            }
            break;
        }
    }
    // ARRAY_SHOW(f, 1, Kernel_Dimension::N);

    for (int j = 0; j < Kernel_Dimension::N; ++j) f[j] += grad[j];

    // ARRAY_SHOW(f, 1, Kernel_Dimension::N);
    return 0;
}

NADForcePolicy::_type Kernel_NADForce::NADForce_type = NADForcePolicy::EHR;

};  // namespace PROJECT_NS
