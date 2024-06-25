#include "kids/Kernel_Update_T.h"

#include "kids/Kernel_Random.h"
#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_Update_T::getName() { return "Kernel_Update_T"; }

int Kernel_Update_T::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Update_T::setInputParam_impl(std::shared_ptr<Param> PM) { gammal = PM->get_double("gammal", LOC(), 0.1); }

void Kernel_Update_T::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    dt_ptr = DS->def(DATA::iter::dt);
    m      = DS->def(DATA::integrator::m);
    p      = DS->def(DATA::integrator::p);

    // if Langevin dynamics, set optimal c1 & c2p
    c1  = DS->def_real("integrator.c1", Dimension::PN);
    c2p = DS->def_real("integrator.c2p", Dimension::PN);
    for (int i = 0; i < Dimension::PN; ++i) {
        c1[i]  = exp(-gammal * scale * dt_ptr[0]);
        c2p[i] = sqrt(1.0 - c1[i] * c1[i]);
    }

    // if for NHC; registeration for auxiliary variables
    // ...
}

Status& Kernel_Update_T::executeKernel_impl(Status& stat) {
    for (int i = 0; i < Dimension::PN; ++i) {
        Kernel_Random::rand_gaussian(&randu);
        p[i] = c1[i] * p[i] + c2p[i] * sqrt(m[i] / beta) * randu;
    }
    return stat;
}

};  // namespace PROJECT_NS
