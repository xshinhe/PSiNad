#include "Kernel_Iter.h"

namespace PROJECT_NS {

void Kernel_Iter::init_data_impl(DataSet* DS) {
    istep_ptr = DS->reg<int>("timer.istep");
    nstep_ptr = DS->reg<int>("timer.nstep");
}

int Kernel_Iter::exec_kernel_impl(int stat) {
    int& istep_ref = *istep_ptr;
    int& nstep_ref = *nstep_ptr;
    while (istep_ref < nstep_ref) {
        for (auto& pkernel : _kernel_vector) { pkernel->exec_kernel(stat); }
    }
    return 0;
}
};  // namespace PROJECT_NS
