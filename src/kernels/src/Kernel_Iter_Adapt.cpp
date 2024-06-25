#include "kids/Kernel_Iter_Adapt.h"

#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

#define FMTF(X)                                                      \
    " " << std::setiosflags(std::ios::fixed) /*scientific notation*/ \
        << std::setprecision(X)              /*precision*/           \
        << std::setw(X + 4)                  /*precision*/

namespace PROJECT_NS {

const std::string Kernel_Iter_Adapt::getName() { return "Kernel_Iter_Adapt"; }

int Kernel_Iter_Adapt::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Iter_Adapt::setInputParam_impl(std::shared_ptr<Param> PM) {
    t0   = PM->get_double("t0", LOC(), phys::time_d, 0.0f);
    tend = PM->get_double("tend", LOC(), phys::time_d, 1.0f);
    dt   = PM->get_double("dt", LOC(), phys::time_d, 0.1f);

    time_unit = PM->get_double("time_unit", LOC(), phys::time_d, 1.0f);

    sstep   = PM->get_int("sstep", LOC(), 1);
    msize   = PM->get_int("msize", LOC(), 128);
    nbackup = PM->get_int("nbackup", LOC(), 1);
    nstep   = sstep * (int((tend - t0) / (sstep * dt)) + 1);  // @bug?
    nsamp   = nstep / sstep + 1;
}

void Kernel_Iter_Adapt::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    t_ptr      = DS->def(DATA::iter::t);
    dt_ptr     = DS->def(DATA::iter::dt);
    istep_ptr  = DS->def(DATA::iter::istep);
    isamp_ptr  = DS->def(DATA::iter::isamp);
    tsize_ptr  = DS->def(DATA::iter::tsize);
    dtsize_ptr = DS->def(DATA::iter::dtsize);

    at_samplingstep_initially_ptr = DS->def(DATA::iter::at_samplingstep_initially);
    at_samplingstep_finally_ptr   = DS->def(DATA::iter::at_samplingstep_finally);

    succ_ptr         = DS->def(DATA::iter::succ);
    last_attempt_ptr = DS->def(DATA::iter::last_attempt);
    frez_ptr         = DS->def(DATA::iter::frez);
    fail_type_ptr    = DS->def(DATA::iter::fail_type);

    // initializarion
    DS->def_int("iter.sstep", &sstep);
    DS->def_int("iter.nstep", &nstep);
    DS->def_int("iter.nsamp", &nsamp);
    DS->def_int("iter.msize", &msize);
}

Status& Kernel_Iter_Adapt::initializeKernel_impl(Status& stat) {
    t_ptr[0]      = t0;
    dt_ptr[0]     = dt;
    isamp_ptr[0]  = 0;
    istep_ptr[0]  = 0;
    tsize_ptr[0]  = 0;
    dtsize_ptr[0] = msize;

    succ_ptr[0]         = true;
    last_attempt_ptr[0] = false;
    frez_ptr[0]         = false;
    fail_type_ptr[0]    = 0;
    return stat;
}

Status& Kernel_Iter_Adapt::executeKernel_impl(Status& stat) {
    int last_tried_dtsize = msize;

    std::cout << "S|: "                       //
              << std::setw(10) << "Progress"  //
              << std::setw(10) << "Time"      //
              << std::setw(10) << "tsize"     //
              << std::setw(10) << "dtsize"    //
              << std::setw(10) << "try"       //
              << "\n";

    while (istep_ptr[0] < nstep) {
        int tsize_before_loop = tsize_ptr[0];                  ///< current time-point tick
        int tsize_after_loop  = tsize_ptr[0] + dtsize_ptr[0];  ///< next time-point tick after loop

        bool at_fullstep_initially       = tsize_before_loop % (msize) == 0;
        bool at_fullstep_finally         = tsize_after_loop % (msize) == 0;
        at_samplingstep_initially_ptr[0] = tsize_before_loop % (sstep * msize) == 0;
        at_samplingstep_finally_ptr[0]   = tsize_after_loop % (sstep * msize) == 0;
        t_ptr[0]                         = t0 + dt * (tsize_ptr[0] / ((double) msize));
        dt_ptr[0]                        = dt * (dtsize_ptr[0] / ((double) msize));

        // backup before loop
        for (auto& fname : backup_fields) {
            for (int bto = nbackup, bfrom = bto - 1; bto > 1; --bto, --bfrom) {
                _dataset->_def(utils::concat("backup.", bto, ".", fname), utils::concat("backup.", bfrom, ".", fname));
            }
            _dataset->_def(utils::concat("backup.", 1, ".", fname),  //
                           utils::concat("integrator.", fname));
        }

        // for each loop, we always previously set succ_ptr[0] = true
        succ_ptr[0] = true;
        for (auto& pkernel : _child_kernels) { pkernel->executeKernel(stat); }

        char statc = '?';
        if (succ_ptr[0] && !last_attempt_ptr[0]) statc = 'T';
        if (succ_ptr[0] && last_attempt_ptr[0]) statc = 'R';
        if (!succ_ptr[0] && last_attempt_ptr[0]) statc = 'X';
        if (!succ_ptr[0] && !last_attempt_ptr[0] && dtsize_ptr[0] > 1) statc = 'F';
        if (!succ_ptr[0] && !last_attempt_ptr[0] && dtsize_ptr[0] == 1) statc = 'L';
        if (frez_ptr[0]) statc = 'Z';

        if (std::ifstream{"X_STAT"}.good() && !frez_ptr[0]) statc = 'X';
        if (std::ifstream{utils::concat("X_STAT", stat.icalc)}.good() && !frez_ptr[0]) statc = 'X';

        switch (statc) {
            case 'X': {
                // save breakdown information
                std::string   directory = _param->get_string("directory", LOC());
                std::ofstream ofs{utils::concat(directory, "/fail", stat.icalc, "-", istep_ptr[0], ".ds")};
                _dataset->dump(ofs);
                ofs.close();

                frez_ptr[0] = true;
                // not break here!
            }
            case 'T':
            case 'R':
            case 'Z': {
                last_attempt_ptr[0] = false;
                if (at_fullstep_finally) istep_ptr[0]++;
                tsize_ptr[0] += dtsize_ptr[0];

                int extend_dtsize = (at_fullstep_finally) ? 2 * last_tried_dtsize : 2 * dtsize_ptr[0];
                int remain_dtsize = msize - (tsize_ptr[0] % msize);
                int new_dtsize    = std::min({msize, extend_dtsize, remain_dtsize});
                last_tried_dtsize = dtsize_ptr[0];
                dtsize_ptr[0]     = new_dtsize;
                break;
            }
            case 'L':
            case 'F': {
                last_attempt_ptr[0] = (dtsize_ptr[0] == 1);
                last_tried_dtsize   = dtsize_ptr[0];
                dtsize_ptr[0]       = (last_attempt_ptr[0]) ? dtsize_ptr[0] : dtsize_ptr[0] / 2;
                tsize_ptr[0] += 0;

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
                  << std::setprecision(2) << std::setw(10) << 100 * t_ptr[0] / tend << "%"  //
                  << std::resetiosflags(std::ios::fixed) << std::setiosflags(std::ios::scientific)
                  << std::setprecision(2) << std::setw(10) << t_ptr[0] / time_unit  //
                  << std::setw(10) << tsize_before_loop                             //
                  << std::setw(10) << last_tried_dtsize                             //
                  << std::setw(10) << dtsize_ptr[0] << "\n";
        isamp_ptr[0] = istep_ptr[0] / sstep;
    }
    at_samplingstep_initially_ptr[0] = true;
    dt_ptr[0]                        = 0;
    return stat;
}

};  // namespace PROJECT_NS
