#include "Kernel_Iter_Adapt.h"

#include "Kernel_Declare.h"

namespace kids {

void Kernel_Iter_Adapt::read_param_impl(Param* PM) {
    t0        = PM->get<double>("t0", LOC(), phys::time_d, 0.0f);
    tend      = PM->get<double>("tend", LOC(), phys::time_d, 1.0f);
    dt        = PM->get<double>("dt", LOC(), phys::time_d, 0.1f);
    int sstep = PM->get<int>("sstep", LOC(), 10);
    tsec      = PM->get<double>("tsec", LOC(), phys::time_d, sstep * dt);
}

void Kernel_Iter_Adapt::init_data_impl(DataSet* DS) {
    t_ptr       = DS->def<double>("iter.t");
    dt_ptr      = DS->def<double>("iter.dt");
    tend_ptr    = DS->def<double>("iter.tend");
    tsec_ptr    = DS->def<double>("iter.tsec");
    succ_ptr    = DS->def<int>("iter.succ");
    nsamp_ptr   = DS->def<int>("iter.nsamp");
    isamp_ptr   = DS->def<int>("iter.isamp");
    do_recd_ptr = DS->def<int>("iter.do_recd");
    do_prec_ptr = DS->def<int>("iter.do_prec");

    *nsamp_ptr  = round(tend / tsec) + 2;  // @debug to be remove 1
    (*tend_ptr) = tend;
    (*tsec_ptr) = tsec;
}

void Kernel_Iter_Adapt::init_calc_impl(int stat) {
    (*t_ptr)     = t0;
    (*dt_ptr)    = dt;
    (*isamp_ptr) = 0;
}

int Kernel_Iter_Adapt::exec_kernel_impl(int stat) {
    while (isamp_ptr[0] < nsamp_ptr[0]) {
        // update indicators
        double t_over_tsec        = (t_ptr[0] - t0) / tsec;
        double dt_over_tsec       = dt_ptr[0] / tsec;
        double t_add_dt_over_tsec = t_over_tsec + dt_over_tsec;

        do_recd_ptr[0] = fabs(t_over_tsec - round(t_over_tsec)) < 0.5f * dt_over_tsec ? 1 : 0;
        do_prec_ptr[0] = fabs(t_add_dt_over_tsec - round(t_add_dt_over_tsec)) < 0.5f * dt_over_tsec ? 1 : 0;

        // backups
        // for (auto& fname : backup_fields) {
        //     for (int bto = nbackup, bfrom = bto - 1; bto > 1; --bto, --bfrom) {
        //         _DataSet->copy_field(utils::concat("backup.", bto, ".", fname),
        //                              utils::concat("backup.", bfrom, ".", fname));
        //     }
        //     _DataSet->copy_field(utils::concat("backup.", 1, ".", fname),  //
        //                          utils::concat("integrator.", fname));
        // }

        // each loop
        for (auto& pkernel : _kernel_vector) { pkernel->exec_kernel(stat); }

        // check success & update timer
        // double* V = _DataSet->def<double>("model.V", Dimension::PFF);
        // if (fabs(V[0] - V[4]) * dt_ptr[0] < 0.01) succ_ptr[0] = 1;

        if (succ_ptr[0] == 1 || true) {
            t_ptr[0] += dt_ptr[0];
            double extend_dt = 2 * dt_ptr[0];
            double remain_dt = dt + int((t_ptr[0] - t0) / dt) * dt - (t_ptr[0] - t0);
            if (remain_dt < dt / 1024) remain_dt = dt;

            // suggested dt for next step
            dt_ptr[0] = std::min({extend_dt, dt, remain_dt});
        } else {
            t_ptr[0] += 0.0e0;
            dt_ptr[0] = std::max({dt_ptr[0] / 2, dt / 1024});

            // recover backups
            // for (auto& fname : backup_fields) {
            //     _DataSet->copy_field(utils::concat("integrator.", fname),      //
            //                          utils::concat("backup.", 1, ".", fname),  //
            //     );
            //     for (int bto = 2, bfrom = bto - 1; bto > 1; --bto, --bfrom) {
            //         _DataSet->copy_field(utils::concat("backup.", bfrom, ".", fname),
            //                              utils::concat("backup.", bto, ".", fname));
            //     }
            // }

            std::cout << dt_ptr[0] << "\n";
            // exit(0);
        }

        std::cout << "dt=" << dt_ptr[0] << "\n";

        if (do_prec_ptr[0] == 1) isamp_ptr[0]++;
    }

    double t_over_tsec        = (t_ptr[0] - t0) / tsec;
    double dt_over_tsec       = dt_ptr[0] / tsec;
    double t_add_dt_over_tsec = t_over_tsec + dt_over_tsec;
    do_recd_ptr[0]            = fabs(t_over_tsec - round(t_over_tsec)) < 0.5f * dt_over_tsec ? 1 : 0;
    dt_ptr[0]                 = 0;

    return 0;
}

};  // namespace kids