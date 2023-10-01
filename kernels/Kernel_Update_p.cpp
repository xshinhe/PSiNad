#include "Kernel_Update_p.h"

#include "Kernel_Declare.h"

namespace PROJECT_NS {

void Kernel_Update_p::read_param_impl(Param* PM) {
    dt  = PM->get<double>("dt", LOC(), phys::time_d);
    sdt = scale * dt;
}

void Kernel_Update_p::init_data_impl(DataSet* DS) {
    f = DS->reg<num_real>("integrator.f", Dimension::PN);
    p = DS->reg<num_real>("integrator.p", Dimension::PN);
}

int Kernel_Update_p::exec_kernel_impl(int stat) {
    for (int i = 0; i < Dimension::PN; ++i) p[i] -= f[i] * sdt;
    return 0;
}

};  // namespace PROJECT_NS