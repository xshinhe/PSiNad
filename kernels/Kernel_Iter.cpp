#include "Kernel_Iter.h"

namespace PROJECT_NS {

void Kernel_Iter::read_param_impl(Param* PM) {
    t0    = PM->get<double>("t0", LOC(), phys::time_d, 0.0f);
    tend  = PM->get<double>("tend", LOC(), phys::time_d, 1.0f);
    dt    = PM->get<double>("dt", LOC(), phys::time_d, 0.1f);
    sstep = PM->get<int>("sstep", LOC(), 1);
    nstep = sstep * (int((tend - t0) / (sstep * dt)) + 1);  // @bug?
    nsamp = nstep / sstep + 1;
}

void Kernel_Iter::init_data_impl(DataSet* DS) {
    t_ptr                         = DS->def<kids_real>("iter.t");
    dt_ptr                        = DS->def<kids_real>("iter.dt");
    istep_ptr                     = DS->def<int>("iter.istep");
    isamp_ptr                     = DS->def<int>("iter.isamp");
    at_samplingstep_initially_ptr = DS->def<bool>("iter.at_samplingstep_initially");
    at_samplingstep_finally_ptr   = DS->def<bool>("iter.at_samplingstep_finally");
    // initializarion
    DS->def<int>("iter.sstep", &sstep);
    DS->def<int>("iter.nstep", &nstep);
    DS->def<int>("iter.nsamp", &nsamp);
}

void Kernel_Iter::init_calc_impl(int stat) {
    t_ptr[0]     = t0;
    dt_ptr[0]    = dt;
    istep_ptr[0] = 0;
    isamp_ptr[0] = 0;
}

int Kernel_Iter::exec_kernel_impl(int stat) {
    while (istep_ptr[0] < nstep) {
        at_samplingstep_finally_ptr[0]   = ((istep_ptr[0] + 1) % sstep == 0);
        at_samplingstep_initially_ptr[0] = (istep_ptr[0] % sstep == 0);

        for (auto& pkernel : _kernel_vector) { pkernel->exec_kernel(stat); }

        t_ptr[0] += dt_ptr[0];

        istep_ptr[0]++;
        isamp_ptr[0] = istep_ptr[0] / sstep;
    }
    at_samplingstep_initially_ptr[0] = true;  // only record!
    dt_ptr[0]                        = 0;     // no-dynamics!
    return 0;
}
};  // namespace PROJECT_NS
