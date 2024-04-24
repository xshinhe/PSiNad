#include "kids/Kernel_Update_c.h"

#include "kids/Kernel_Representation.h"
#include "kids/linalg.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

void Kernel_Update_c::setInputDataSet_impl(std::shared_ptr<DataSet>& DS) {
    dt_ptr        = DS->def(DATA::iter::dt);
    E             = DS->def(DATA::model::rep::E);
    T             = DS->def(DATA::model::rep::T);
    L             = DS->def(DATA::model::rep::L);
    R             = DS->def(DATA::model::rep::R);
    U             = DS->def(DATA::integrator::U);
    Udt           = DS->def(DATA::integrator::Udt);
    succ_ptr      = DS->def(DATA::iter::succ);
    frez_ptr      = DS->def(DATA::iter::frez);
    invexpidiagdt = DS->def(DATA::integrator::tmp::invexpidiagdt);
}

Status& Kernel_Update_c::executeKernel_impl(Status& stat) {
    if (frez_ptr[0]) return 0;

    for (int iP = 0; iP < Dimension::P; ++iP) {
        // local variables for iP-th of swarm
        kids_real*    E   = this->E + iP * Dimension::F;
        kids_real*    T   = this->T + iP * Dimension::FF;
        kids_real*    L   = this->L + iP * Dimension::F;
        kids_complex* R   = this->R + iP * Dimension::FF;
        kids_complex* U   = this->U + iP * Dimension::FF;
        kids_complex* Udt = this->Udt + iP * Dimension::FF;

        switch (Kernel_Representation::ele_repr_type) {
            case RepresentationPolicy::Diabatic: {
                for (int i = 0; i < Dimension::F; ++i)
                    invexpidiagdt[i] = exp(-phys::math::im * E[i] * scale * dt_ptr[0]);
                ARRAY_MATMUL3_TRANS2(Udt, T, invexpidiagdt, T, Dimension::F, Dimension::F, 0, Dimension::F);
                ARRAY_MATMUL(U, Udt, U, Dimension::F, Dimension::F, Dimension::F);
                break;
            }
            case RepresentationPolicy::Adiabatic: {
                for (int i = 0; i < Dimension::F; ++i)
                    invexpidiagdt[i] = exp(-phys::math::im * L[i] * scale * dt_ptr[0]);
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