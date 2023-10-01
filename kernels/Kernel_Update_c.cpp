#include "Kernel_Update_c.h"

#include "../core/linalg.h"
#include "Kernel_Declare.h"
#include "Kernel_Representation.h"

namespace PROJECT_NS {

void Kernel_Update_c::read_param_impl(Param* PM) {
    dt  = PM->get<double>("dt", LOC(), phys::time_d);  //
    sdt = scale * dt;
}

void Kernel_Update_c::init_data_impl(DataSet* S) {
    E   = S->reg<num_real>("model.rep.E", Dimension::PF);
    T   = S->reg<num_real>("model.rep.T", Dimension::PFF);
    L   = S->reg<num_real>("model.rep.L", Dimension::PF);
    R   = S->reg<num_complex>("model.rep.R", Dimension::PFF);
    U   = S->reg<num_complex>("integrator.U", Dimension::PFF);
    Udt = S->reg<num_complex>("integrator.Udt", Dimension::PFF);

    invexpidiagdt = S->reg<num_complex>("integrator.tmp.invexpidiagdt", Dimension::F);
}

int Kernel_Update_c::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        // local variables for iP-th of swarm
        num_real* E      = this->E + iP * Dimension::F;
        num_real* T      = this->T + iP * Dimension::FF;
        num_real* L      = this->L + iP * Dimension::F;
        num_complex* R   = this->R + iP * Dimension::FF;
        num_complex* U   = this->U + iP * Dimension::FF;
        num_complex* Udt = this->Udt + iP * Dimension::FF;

        switch (Kernel_Representation::ele_repr_type) {
            case RepresentationPolicy::Diabatic: {
                for (int i = 0; i < Dimension::F; ++i) invexpidiagdt[i] = exp(-phys::math::im * E[i] * dt);
                ARRAY_MATMUL3_TRANS2(Udt, T, invexpidiagdt, T, Dimension::F, Dimension::F, 0, Dimension::F);
                ARRAY_MATMUL(U, Udt, U, Dimension::F, Dimension::F, Dimension::F);
                break;
            }
            case RepresentationPolicy::Adiabatic: {
                for (int i = 0; i < Dimension::F; ++i) invexpidiagdt[i] = exp(-phys::math::im * L[i] * dt);
                ARRAY_MATMUL3_TRANS2(Udt, R, invexpidiagdt, R, Dimension::F, Dimension::F, 0, Dimension::F);
                ARRAY_MATMUL(U, Udt, U, Dimension::F, Dimension::F, Dimension::F);
                break;
            }
            default:  // representation_policy::force, representation_policy::density
                      // LOG(FATAL);
                break;
        }
    }
    return 0;
}
};  // namespace PROJECT_NS