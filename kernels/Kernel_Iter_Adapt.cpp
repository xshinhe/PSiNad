#include "Kernel_Iter_Adapt.h"

#include "Kernel_Declare.h"

namespace PROJECT_NS {

void Kernel_Iter_Adapt::read_param_impl(Param* PM) {
    t0      = PM->get<double>("t0", LOC(), phys::time_d, 0.0f);
    tend    = PM->get<double>("tend", LOC(), phys::time_d, 1.0f);
    dt      = PM->get<double>("dt", LOC(), phys::time_d, 0.1f);
    sstep   = PM->get<int>("sstep", LOC(), 1);
    msize   = PM->get<int>("msize", LOC(), 128);
    nbackup = PM->get<int>("nbackup", LOC(), 1);
    nstep   = sstep * (int((tend - t0) / (sstep * dt)) + 1);  // @bug?
    nsamp   = nstep / sstep + 1;
}

void Kernel_Iter_Adapt::init_data_impl(DataSet* DS) {
    t_ptr      = DS->def<kids_real>("iter.t");
    dt_ptr     = DS->def<kids_real>("iter.dt");
    istep_ptr  = DS->def<int>("iter.istep");
    isamp_ptr  = DS->def<int>("iter.isamp");
    tsize_ptr  = DS->def<int>("iter.tsize");
    dtsize_ptr = DS->def<int>("iter.dtsize");

    succ_ptr    = DS->def<bool>("iter.succ");
    do_recd_ptr = DS->def<bool>("iter.do_recd");
    do_prec_ptr = DS->def<bool>("iter.do_prec");

    // initializarion
    DS->def<int>("iter.sstep", &sstep);
    DS->def<int>("iter.nstep", &nstep);
    DS->def<int>("iter.nsamp", &nsamp);
    DS->def<int>("iter.msize", &msize);
}

void Kernel_Iter_Adapt::init_calc_impl(int stat) {
    t_ptr[0]      = t0;
    dt_ptr[0]     = dt;
    isamp_ptr[0]  = 0;
    istep_ptr[0]  = 0;
    tsize_ptr[0]  = 0;
    dtsize_ptr[0] = msize;
}

int Kernel_Iter_Adapt::exec_kernel_impl(int stat) {
    while (istep_ptr[0] < nstep) {
        t_ptr[0]  = t0 + dt * (tsize_ptr[0] / ((double) msize));
        dt_ptr[0] = dt * (dtsize_ptr[0] / ((double) msize));

        do_recd_ptr[0] = tsize_ptr[0] % (sstep * msize) == 0;
        do_prec_ptr[0] = (tsize_ptr[0] + dtsize_ptr[0]) % (sstep * msize) == 0;

        // backups
        for (auto& fname : backup_fields) {
            for (int bto = nbackup, bfrom = bto - 1; bto > 1; --bto, --bfrom) {
                _DataSet->_def(utils::concat("backup.", bto, ".", fname), utils::concat("backup.", bfrom, ".", fname));
            }
            _DataSet->_def(utils::concat("backup.", 1, ".", fname),  //
                           utils::concat("integrator.", fname));
        }

        // each loop
        for (auto& pkernel : _kernel_vector) { pkernel->exec_kernel(stat); }

        // check success & update timer
        // double* V = _DataSet->def<double>("model.V", Dimension::PFF);
        // if (fabs(V[0] - V[4]) * dt_ptr[0] < 0.01) succ_ptr[0] = 1;

        if (istep_ptr[0] % 2 == 0) {
            if ((tsize_ptr[0] + dtsize_ptr[0]) % msize == 0) istep_ptr[0]++;
            tsize_ptr[0] += dtsize_ptr[0];

            int extend_dtsize = 2 * dtsize_ptr[0];
            int remain_dtsize = msize - (tsize_ptr[0] % msize);
            dtsize_ptr[0]     = std::min({msize, extend_dtsize, remain_dtsize});

        } else {
            if (dtsize_ptr[0] % 2 == 0) {
                dtsize_ptr[0] /= 2;
                tsize_ptr[0] += 0;
                // recover backups
                for (auto& fname : backup_fields) {
                    _DataSet->_def(utils::concat("integrator.", fname),     //
                                   utils::concat("backup.", 1, ".", fname)  //
                    );
                    for (int bto = nbackup, bfrom = bto - 1; bto > 1; --bto, --bfrom) {
                        _DataSet->_def(utils::concat("backup.", bfrom, ".", fname),
                                       utils::concat("backup.", bto, ".", fname));
                    }
                }

            } else {
                if ((tsize_ptr[0] + dtsize_ptr[0]) % msize == 0) istep_ptr[0]++;
                tsize_ptr[0] += dtsize_ptr[0];
                dtsize_ptr[0] = dtsize_ptr[0];

                std::cout << "current dt_dynamic / dt = " << dtsize_ptr[0] / ((double) msize) << "\n";
                std::cout << "exceed minial dt!\n";
            }
        }

        std::cout << "dt ===== " << dt_ptr[0] << "\n";

        isamp_ptr[0] = istep_ptr[0] / sstep;
    }
    do_recd_ptr[0] = true;
    dt_ptr[0]      = 0;
    return 0;
}

};  // namespace PROJECT_NS