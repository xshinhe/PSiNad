#include "Kernel_Update_c.h"

#include "../core/linalg.h"
#include "Kernel_Declare.h"
#include "Kernel_Representation.h"

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

void Kernel_Update_c::init_data_impl(DataSet* DS) {
    dt_ptr = DS->def<kids_real>("iter.dt");

    E        = DS->def<kids_real>("model.rep.E", Dimension::PF);
    T        = DS->def<kids_real>("model.rep.T", Dimension::PFF);
    L        = DS->def<kids_real>("model.rep.L", Dimension::PF);
    R        = DS->def<kids_complex>("model.rep.R", Dimension::PFF);
    U        = DS->def<kids_complex>("integrator.U", Dimension::PFF);
    Udt      = DS->def<kids_complex>("integrator.Udt", Dimension::PFF);
    succ_ptr = DS->def<bool>("iter.succ");
    frez_ptr = DS->def<bool>("iter.frez");

    invexpidiagdt = DS->def<kids_complex>("integrator.tmp.invexpidiagdt", Dimension::F);
}

int Kernel_Update_c::exec_kernel_impl(int stat) {
    if (!succ_ptr[0] || frez_ptr[0]) return 0;

    for (int iP = 0; iP < Dimension::P; ++iP) {
        // local variables for iP-th of swarm
        kids_real* E      = this->E + iP * Dimension::F;
        kids_real* T      = this->T + iP * Dimension::FF;
        kids_real* L      = this->L + iP * Dimension::F;
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