#include "Kernel_Iter.h"

namespace kids {

void Kernel_Iter::read_param_impl(Param* PM) {
    t0    = PM->get<double>("t0", LOC(), phys::time_d, 0.0f);
    tend  = PM->get<double>("tend", LOC(), phys::time_d, 1.0f);
    dt    = PM->get<double>("dt", LOC(), phys::time_d, 0.1f);
    sstep = PM->get<int>("sstep", LOC(), 1);
    nstep = sstep * (int((tend - t0) / (sstep * dt)) + 1);  // @bug?
    nsamp = nstep / sstep + 1;
}

void Kernel_Iter::init_data_impl(DataSet* DS) {
    t_ptr       = DS->def<kids_real>("iter.t");
    dt_ptr      = DS->def<kids_real>("iter.dt");
    tsec_ptr    = DS->def<kids_real>("iter.tsec");
    tend_ptr    = DS->def<kids_real>("iter.tend");
    sstep_ptr   = DS->def<int>("iter.sstep");
    istep_ptr   = DS->def<int>("iter.istep");
    nstep_ptr   = DS->def<int>("iter.nstep");
    isamp_ptr   = DS->def<int>("iter.isamp");
    nsamp_ptr   = DS->def<int>("iter.nsamp");
    do_recd_ptr = DS->def<int>("iter.do_recd");
    do_prec_ptr = DS->def<int>("iter.do_prec");
    // initial
    dt_ptr[0]    = dt;
    sstep_ptr[0] = sstep;
    nstep_ptr[0] = nstep;
    nsamp_ptr[0] = nsamp;
    tsec_ptr[0]  = sstep * dt;
}

void Kernel_Iter::init_calc_impl(int stat) {
    istep_ptr[0] = 0;
    isamp_ptr[0] = 0;
    t_ptr[0]     = t0;
    dt_ptr[0]    = dt;
}

int Kernel_Iter::exec_kernel_impl(int stat) {
    while (istep_ptr[0] < nstep_ptr[0]) {
        do_prec_ptr[0] = ((istep_ptr[0] + 1) % sstep == 0) ? 1 : 0;
        do_recd_ptr[0] = (istep_ptr[0] % sstep == 0) ? 1 : 0;
        for (auto& pkernel : _kernel_vector) { pkernel->exec_kernel(stat); }
        t_ptr[0] += dt_ptr[0];
        istep_ptr[0]++;
        isamp_ptr[0] = istep_ptr[0] / sstep;
    }
    do_recd_ptr[0] = (istep_ptr[0] % sstep == 0) ? 1 : 0;  // if record last step?
    dt_ptr[0]      = 0;
    return 0;
}
};  // namespace kids
