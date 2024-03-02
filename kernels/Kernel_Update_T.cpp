#include "Kernel_Update_T.h"

#include "Kernel_Declare.h"
#include "Kernel_Random.h"

namespace kids {

void Kernel_Update_T::read_param_impl(Param* PM) { gammal = PM->get<double>("gammal", LOC(), 0.1); }

void Kernel_Update_T::init_data_impl(DataSet* DS) {
    dt_ptr = DS->def<kids_real>("iter.dt");
    m      = DS->def<kids_real>("integrator.m", Dimension::PN);
    p      = DS->def<kids_real>("integrator.p", Dimension::PN);

    // if Langevin dynamics, set optimal c1 & c2p
    c1  = DS->def<kids_real>("integrator.c1", Dimension::PN);
    c2p = DS->def<kids_real>("integrator.c2p", Dimension::PN);
    for (int i = 0; i < Dimension::PN; ++i) {
        c1[i]  = exp(-gammal * scale * dt_ptr[0]);
        c2p[i] = sqrt(1.0 - c1[i] * c1[i]);
    }

    // if for NHC; registeration for auxiliary variables
    // ...
}

int Kernel_Update_T::exec_kernel_impl(int stat) {
    for (int i = 0; i < Dimension::PN; ++i) {
        Kernel_Random::rand_gaussian(&randu);
        p[i] = c1[i] * p[i] + c2p[i] * sqrt(m[i] / beta) * randu;
    }
    return 0;
}

};  // namespace kids
