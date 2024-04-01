#include "Kernel_Update_T.h"

#include "../core/vars_list.h"
#include "Kernel_Declare.h"
#include "Kernel_Random.h"

namespace PROJECT_NS {

void Kernel_Update_T::read_param_impl(Param* PM) { gammal = PM->get<double>("gammal", LOC(), 0.1); }

void Kernel_Update_T::init_data_impl(DataSet* DS) {
    dt_ptr = DS->def(DATA::iter::dt);
    m      = DS->def(DATA::integrator::m);
    p      = DS->def(DATA::integrator::p);

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

};  // namespace PROJECT_NS
