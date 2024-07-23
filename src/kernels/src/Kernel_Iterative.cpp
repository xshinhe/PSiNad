#include "kids/Kernel_Iterative.h"

#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_Iterative::getName() { return "Kernel_Iterative"; }

int Kernel_Iterative::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Iterative::setInputParam_impl(std::shared_ptr<Param> PM) {
    t0    = _param->get_real({"model.t0", "solver.t0"}, LOC(), phys::time_d, 0.0f);
    tend  = _param->get_real({"model.tend", "solver.tend"}, LOC(), phys::time_d, 1.0f);
    dt0   = _param->get_real({"model.dt", "solver.dt"}, LOC(), phys::time_d, 0.1f);
    sstep = _param->get_int({"solver.sstep"}, LOC(), 1);

    // set time grids
    nstep = sstep * (int((tend - t0) / (sstep * dt0)));  // @bug? (try new algo for nstep)
    nsamp = nstep / sstep + 1;

    // Dimension::sstep = sstep;
    // Dimension::nstep = nstep;
    // Dimension::nsamp = nsamp;
}

void Kernel_Iterative::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    t            = DS->def(DATA::flowcontrol::t);
    dt           = DS->def(DATA::flowcontrol::dt);
    istep        = DS->def(DATA::flowcontrol::istep);
    isamp        = DS->def(DATA::flowcontrol::isamp);
    at_condition = DS->def(DATA::flowcontrol::at_condition);
    DS->def_int("flowcontrol.sstep", &sstep);
    DS->def_int("flowcontrol.nstep", &nstep);
    DS->def_int("flowcontrol.nsamp", &nsamp);
}

Status& Kernel_Iterative::initializeKernel_impl(Status& stat) {
    t[0]     = t0;
    dt[0]    = dt0;
    istep[0] = 0;
    isamp[0] = 0;
    return stat;
}

Status& Kernel_Iterative::executeKernel_impl(Status& stat) {
    stat.first_step = true;
    while (istep[0] <= nstep) {
        if (istep[0] == nstep) dt[0] = 0;  // set dt=0 to remove dynamics! only record in last step
        at_condition[0] = (istep[0] % sstep == 0);
        isamp[0]        = istep[0] / sstep;
        for (auto& pkernel : _child_kernels) { pkernel->executeKernel(stat); }
        t[0] += dt[0];
        istep[0]++;
        stat.first_step = false;
    }
    return stat;
}
};  // namespace PROJECT_NS
