#include "kids/Kernel_Iter.h"

#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_Iter::getName() { return "Kernel_Iter"; }

int Kernel_Iter::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Iter::setInputParam_impl(std::shared_ptr<Param> PM) {
    t0    = PM->get_double("t0", LOC(), phys::time_d, 0.0f);
    tend  = PM->get_double("tend", LOC(), phys::time_d, 1.0f);
    dt    = PM->get_double("dt", LOC(), phys::time_d, 0.1f);
    sstep = PM->get_int("sstep", LOC(), 1);
    nstep = sstep * (int((tend - t0) / (sstep * dt)) + 1);  // @bug?
    nsamp = nstep / sstep + 1;
}

void Kernel_Iter::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    t_ptr                         = DS->def(DATA::iter::t);
    dt_ptr                        = DS->def(DATA::iter::dt);
    istep_ptr                     = DS->def(DATA::iter::istep);
    isamp_ptr                     = DS->def(DATA::iter::isamp);
    at_samplingstep_initially_ptr = DS->def(DATA::iter::at_samplingstep_initially);
    at_samplingstep_finally_ptr   = DS->def(DATA::iter::at_samplingstep_finally);
    succ_ptr                      = DS->def(DATA::iter::succ);
    // initializarion
    DS->def_int("iter.sstep", &sstep);
    DS->def_int("iter.nstep", &nstep);
    DS->def_int("iter.nsamp", &nsamp);
}

Status& Kernel_Iter::initializeKernel_impl(Status& stat) {
    t_ptr[0]     = t0;
    dt_ptr[0]    = dt;
    istep_ptr[0] = 0;
    isamp_ptr[0] = 0;
    succ_ptr[0]  = true;
    return stat;
}

Status& Kernel_Iter::executeKernel_impl(Status& stat) {
    while (istep_ptr[0] < nstep) {
        at_samplingstep_finally_ptr[0]   = ((istep_ptr[0] + 1) % sstep == 0);
        at_samplingstep_initially_ptr[0] = (istep_ptr[0] % sstep == 0);

        for (auto& pkernel : _child_kernels) { pkernel->executeKernel(stat); }

        t_ptr[0] += dt_ptr[0];

        istep_ptr[0]++;
        isamp_ptr[0] = istep_ptr[0] / sstep;
    }
    at_samplingstep_initially_ptr[0] = true;  // only record!
    dt_ptr[0]                        = 0;     // no-dynamics!
    return stat;
}
};  // namespace PROJECT_NS
