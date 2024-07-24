#include "kids/Kernel_Iterative_Adapt.h"

#include <algorithm>

#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

#define FMTF(X)                                                      \
    " " << std::setiosflags(std::ios::fixed) /*scientific notation*/ \
        << std::setprecision(X)              /*precision*/           \
        << std::setw(X + 4)                  /*precision*/

namespace PROJECT_NS {

const std::string Kernel_Iterative_Adapt::getName() { return "Kernel_Iterative_Adapt"; }

int Kernel_Iterative_Adapt::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Iterative_Adapt::setInputParam_impl(std::shared_ptr<Param> PM) {
    t0        = _param->get_real({"model.t0", "solver.t0"}, LOC(), phys::time_d, 0.0f);
    tend      = _param->get_real({"model.tend", "solver.tend"}, LOC(), phys::time_d, 1.0f);
    dt0       = _param->get_real({"model.dt", "solver.dt"}, LOC(), phys::time_d, 0.1f);
    time_unit = _param->get_real({"model.time_unit", "solver.time_unit"}, LOC(), phys::time_d, 1.0f);
    sstep     = _param->get_int({"solver.sstep"}, LOC(), 1);
    msize     = _param->get_int({"solver.msize"}, LOC(), 128);
    nbackup   = _param->get_int({"solver.nbackup"}, LOC(), 1);

    nstep = sstep * (int((tend - t0) / (sstep * dt0)));  // @bug? (try new algo for nstep)
    nsamp = nstep / sstep + 1;
}

void Kernel_Iterative_Adapt::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    t            = DS->def(DATA::flowcontrol::t);
    dt           = DS->def(DATA::flowcontrol::dt);
    istep        = DS->def(DATA::flowcontrol::istep);
    isamp        = DS->def(DATA::flowcontrol::isamp);
    tsize        = DS->def(DATA::flowcontrol::tsize);
    dtsize       = DS->def(DATA::flowcontrol::dtsize);
    at_condition = DS->def(DATA::flowcontrol::at_condition);

    succ_ptr         = DS->def(DATA::flowcontrol::succ);
    last_attempt_ptr = DS->def(DATA::flowcontrol::last_attempt);
    frez_ptr         = DS->def(DATA::flowcontrol::frez);
    fail_type_ptr    = DS->def(DATA::flowcontrol::fail_type);

    // initializarion
    DS->def_int("flowcontrol.sstep", &sstep);
    DS->def_int("flowcontrol.nstep", &nstep);
    DS->def_int("flowcontrol.nsamp", &nsamp);
    DS->def_int("flowcontrol.msize", &msize);
}

Status& Kernel_Iterative_Adapt::initializeKernel_impl(Status& stat) {
    t[0]      = t0;
    dt[0]     = dt0;
    isamp[0]  = 0;
    istep[0]  = 0;
    tsize[0]  = 0;
    dtsize[0] = msize;

    stat.succ         = true;
    stat.last_attempt = false;
    stat.frozen       = false;
    stat.fail_type    = 0;
    return stat;
}

Status& Kernel_Iterative_Adapt::executeKernel_impl(Status& stat) {
    int last_tried_dtsize = msize;

    std::cout << "S|: "                       //
              << std::setw(10) << "Progress"  //
              << std::setw(10) << "Time"      //
              << std::setw(10) << "tsize"     //
              << std::setw(10) << "dtsize"    //
              << std::setw(10) << "try"       //
              << "\n";

    while (istep[0] <= nstep) {
        int  tsize_before_loop     = tsize[0];              ///< current time-point tick
        int  tsize_after_loop      = tsize[0] + dtsize[0];  ///< next time-point tick after loop
        bool at_fullstep_initially = tsize_before_loop % (msize) == 0;
        bool at_fullstep_finally   = tsize_after_loop % (msize) == 0;
        at_condition[0]            = tsize_before_loop % (sstep * msize) == 0;
        t[0]                       = t0 + dt0 * (tsize[0] / ((double) msize));
        dt[0]                      = dt0 * (dtsize[0] / ((double) msize));

        if (istep[0] == nstep) {
            at_condition[0] = true;
            dt[0]           = 0;  // frozen dynamics
        }

        // backup before loop
        for (auto& fname : backup_fields) {
            for (int bto = nbackup, bfrom = bto - 1; bto > 1; --bto, --bfrom) {
                _dataset->_def(utils::concat("backup.", bto, ".", fname), utils::concat("backup.", bfrom, ".", fname));
            }
            _dataset->_def(utils::concat("backup.", 1, ".", fname),  //
                           utils::concat("integrator.", fname));
        }

        // for each loop, we always previously set succ as true
        stat.succ = true;
        for (auto& pkernel : _child_kernels) { pkernel->executeKernel(stat); }

        char statc = '?';
        if (stat.succ && !stat.last_attempt) statc = 'T';
        if (stat.succ && stat.last_attempt) statc = 'R';
        if (!stat.succ && stat.last_attempt) statc = 'X';
        if (!stat.succ && !stat.last_attempt && dtsize[0] > 1) statc = 'F';
        if (!stat.succ && !stat.last_attempt && dtsize[0] == 1) statc = 'L';
        if (stat.frozen) statc = 'Z';

        if (std::ifstream{"X_STAT"}.good() && !stat.frozen) statc = 'X';
        if (std::ifstream{utils::concat("X_STAT", stat.icalc)}.good() && !stat.frozen) statc = 'X';

        switch (statc) {
            case 'X': {
                // save breakdown information
                std::ofstream ofs{utils::concat(directory, "/fail", stat.icalc, "-", istep[0], ".ds")};
                _dataset->dump(ofs);
                ofs.close();

                stat.frozen = true;
                // not break here!
            }
            case 'T':
            case 'R':
            case 'Z': {
                stat.last_attempt = false;

                if (at_fullstep_finally) istep[0]++;
                tsize[0] += dtsize[0];

                int extend_dtsize = (at_fullstep_finally) ? 2 * last_tried_dtsize : 2 * dtsize[0];
                int remain_dtsize = msize - (tsize[0] % msize);
                int new_dtsize    = std::min({msize, extend_dtsize, remain_dtsize});
                last_tried_dtsize = dtsize[0];
                dtsize[0]         = new_dtsize;
                break;
            }
            case 'L':
            case 'F': {
                stat.last_attempt = (dtsize[0] == 1);
                last_tried_dtsize = dtsize[0];
                dtsize[0]         = (stat.last_attempt) ? dtsize[0] : dtsize[0] / 2;
                tsize[0] += 0;

                for (auto& fname : backup_fields) {
                    _dataset->_def(utils::concat("integrator.", fname),     //
                                   utils::concat("backup.", 1, ".", fname)  //
                    );
                    for (int bto = nbackup, bfrom = bto - 1; bto > 1; --bto, --bfrom) {
                        _dataset->_def(utils::concat("backup.", bfrom, ".", fname),
                                       utils::concat("backup.", bto, ".", fname));
                    }
                }
                break;
            }
        }
        std::cout << statc << "|:"  //
                  << std::resetiosflags(std::ios::scientific) << std::setiosflags(std::ios::fixed)
                  << std::setprecision(2) << std::setw(10) << 100 * t[0] / tend << "%"  //
                  << std::resetiosflags(std::ios::fixed) << std::setiosflags(std::ios::scientific)
                  << std::setprecision(2) << std::setw(10) << t[0] / time_unit  //
                  << std::setw(10) << tsize_before_loop                         //
                  << std::setw(10) << last_tried_dtsize                         //
                  << std::setw(10) << dtsize[0] << "\n";
        isamp[0] = istep[0] / sstep;
    }
    return stat;
}

};  // namespace PROJECT_NS
