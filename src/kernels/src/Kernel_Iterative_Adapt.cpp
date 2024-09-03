#include "kids/Kernel_Iterative_Adapt.h"

#include <unistd.h>

#include <algorithm>
#include <chrono>

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
    t0            = _param->get_real({"model.t0", "solver.t0"}, LOC(), phys::time_d, 0.0f);
    tend          = _param->get_real({"model.tend", "solver.tend"}, LOC(), phys::time_d, 1.0f);
    dt0           = _param->get_real({"model.dt", "solver.dt"}, LOC(), phys::time_d, 0.1f);
    time_unit     = _param->get_real({"model.time_unit", "solver.time_unit"}, LOC(), phys::time_d, 1.0f);
    sstep         = _param->get_int({"solver.sstep"}, LOC(), 1);
    msize         = _param->get_int({"solver.msize"}, LOC(), 128);
    nbackup       = _param->get_int({"solver.nbackup"}, LOC(), 1);
    use_exchange  = _param->get_bool({"solver.use_exchange", "use_exchange"}, LOC(), false);
    exchange_root = _param->get_int({"solver.exchange_root", "exchange_root"}, LOC(), -1);
    exchange_num  = _param->get_int({"solver.exchange_num", "exchange_num"}, LOC(), 100);
    exchange_time = _param->get_real({"solver.exchange_time", "exchange_time"}, LOC(), 600.0);  // in second
    nstep         = sstep * (int((tend - t0) / (sstep * dt0)));  // @bug? (try new algo for nstep)
    nsamp         = nstep / sstep + 1;
}

void Kernel_Iterative_Adapt::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    t                 = DS->def(DATA::flowcontrol::t);
    dt                = DS->def(DATA::flowcontrol::dt);
    istep             = DS->def(DATA::flowcontrol::istep);
    isamp             = DS->def(DATA::flowcontrol::isamp);
    tsize             = DS->def(DATA::flowcontrol::tsize);
    dtsize            = DS->def(DATA::flowcontrol::dtsize);
    last_tried_dtsize = DS->def(DATA::flowcontrol::last_tried_dtsize);
    at_condition      = DS->def(DATA::flowcontrol::at_condition);

    // initializarion
    DS->def(VARIABLE<kids_int>("flowcontrol.sstep", &Dimension::shape_1, "@"))[0] = sstep;
    DS->def(VARIABLE<kids_int>("flowcontrol.nstep", &Dimension::shape_1, "@"))[0] = nstep;
    DS->def(VARIABLE<kids_int>("flowcontrol.nsamp", &Dimension::shape_1, "@"))[0] = nsamp;
    DS->def(VARIABLE<kids_int>("flowcontrol.msize", &Dimension::shape_1, "@"))[0] = msize;
}

Status& Kernel_Iterative_Adapt::initializeKernel_impl(Status& stat) {
    if (_param->get_bool({"restart"}, LOC(), false)) {  //
        std::string loadfile = _param->get_string({"load"}, LOC(), "NULL");
        if (loadfile == "NULL" || loadfile == "" || loadfile == "null") loadfile = "restart.ds";

        if (std::ifstream{"X_STAT"}.good()) remove("X_STAT");
        if (std::ifstream{utils::concat("X_STAT", stat.icalc)}.good()) {
            std::string rmfile = utils::concat("X_STAT", stat.icalc);
            remove(rmfile.c_str());
        }
        istep[0]             = _dataset->def_int("restart.istep", 1)[0];
        tsize[0]             = _dataset->def_int("restart.tsize", 1)[0];
        dtsize[0]            = _dataset->def_int("restart.dtsize", 1)[0];
        last_tried_dtsize[0] = _dataset->def_int("restart.last_tried_dtsize", 1)[0];

        stat.succ         = true;
        stat.last_attempt = false;
        stat.frozen       = false;
        stat.fail_type    = 0;
        return stat;
    }
    if (use_exchange) {
        ex_begin          = std::chrono::steady_clock::now();
        exchange_fulltime = exchange_time;
    }

    t[0]                 = t0;
    dt[0]                = dt0;
    isamp[0]             = 0;
    istep[0]             = 0;
    tsize[0]             = 0;
    dtsize[0]            = msize;
    last_tried_dtsize[0] = msize;

    stat.succ         = true;
    stat.last_attempt = false;
    stat.frozen       = false;
    stat.fail_type    = 0;

    return stat;
}

Status& Kernel_Iterative_Adapt::exchange(Status& stat) {
    auto ex_end     = std::chrono::steady_clock::now();
    auto total_time = static_cast<std::chrono::duration<double>>(ex_end - ex_begin).count();
    // std::cout << LOC() << "total time = " << total_time << " <?> ex time = " << exchange_fulltime << "\n";

    std::string signal_file   = "EXCHANGE";
    std::string nosignal_file = "NOEXCHANGE";
    std::string allowed_file  = "EXCHANGE_ALLOWED";
    std::string flag_file     = utils::concat("EXCHANGE_THREAD_", stat.icalc);
    if (stat.icalc == exchange_root) {
        remove("EXCHANGE_COUNT");      // for root
        remove(allowed_file.c_str());  // for root
    }
    remove(flag_file.c_str());  // for all

    bool do_exchange = false;
    if (std::ifstream{signal_file}.good()) do_exchange = true;
    if (total_time > exchange_fulltime && stat.icalc == exchange_root) do_exchange = true;
    if (!do_exchange) return stat;
    if (std::ifstream{nosignal_file}.good()) return stat;

    // std::cout << LOC() << "TRY EXCHANGE!\n";
    exchange_fulltime += exchange_time;

    bool                             raised    = false;
    bool                             locked    = false;
    bool                             engaged   = true;
    bool                             exchanged = false;
    int                              maxcyc = 99, icyc = 0;
    int                              idx1, idx2;
    std::vector<std::pair<int, int>> list1;
    std::vector<std::pair<int, int>> list2;
    while (engaged || locked) {
        icyc++;
        // for all: lock statuc is associated with allowed_file
        locked = std::ifstream{allowed_file}.good();
        raised = std::ifstream{flag_file}.good();

        if (raised && !locked && icyc > maxcyc &&  //
            stat.icalc != exchange_root && !std::ifstream{signal_file}.good()) {
            remove(flag_file.c_str());
            std::cout << "EXIT AFTER ROOT IS FREED!\n";
            break;
        }

        // for all: if locked, check if it in allowed_file
        // only in locked block, the engaged status can be changed
        if (raised && locked) {
            std::ifstream ifs(allowed_file);
            std::string   eachline;
            engaged = false;
            while (getline(ifs, eachline)) {
                std::stringstream ss(eachline);
                ss >> idx1 >> idx2;
                if (idx1 == stat.icalc) {
                    engaged = true;
                    break;
                }
            }
            if (engaged == false) {  // for NO-ROOT: if not paritcipated in exchange, just free them
                remove(flag_file.c_str());
                if (stat.icalc != exchange_root) break;
            } else {  // for all: if locked and engaged. then do exchange.
                if (idx1 != idx2 && !exchanged) {
                    // dump idx2
                    std::string dump_file = utils::concat("exchange-", idx1, ".ds");
                    std::string load_file = utils::concat("exchange-", idx2, ".ds");
                    _dataset->def_int("restart.istep", istep.data(), 1);
                    _dataset->def_int("restart.tsize", tsize.data(), 1);
                    _dataset->def_int("restart.dtsize", dtsize.data(), 1);
                    _dataset->def_int("restart.last_tried_dtsize", last_tried_dtsize.data(), 1);
                    std::ofstream ofs{dump_file};
                    _dataset->dump(ofs);
                    ofs.close();

                    // wait for idx2 to be prepared (including root)
                    while (!std::ifstream{load_file}.good()) { usleep(5000); }

                    // load idx2
                    std::ifstream ifs{load_file};
                    _dataset->load(ifs);
                    ifs.close();
                    istep[0]             = _dataset->def_int("restart.istep", 1)[0];
                    tsize[0]             = _dataset->def_int("restart.tsize", 1)[0];
                    dtsize[0]            = _dataset->def_int("restart.dtsize", 1)[0];
                    last_tried_dtsize[0] = _dataset->def_int("restart.last_tried_dtsize", 1)[0];
                    stat.succ            = true;
                    stat.last_attempt    = false;
                    stat.frozen          = false;
                    stat.fail_type       = 0;
                    std::cout << "IDX " << stat.icalc << " has loaded for dataset " << load_file << "\n";
                    remove(load_file.c_str());  // avoid reload by other times
                }
                exchanged = true;  // including self-exchange
                // delete file and wait for leaving
                remove(flag_file.c_str());
                if (stat.icalc != exchange_root) {
                    while (std::ifstream{allowed_file}.good()) usleep(5000);
                    break;
                }
            }
        }
        if (raised && exchanged && !locked) break;

        // to generate allow file by root
        // 1) for ROOT: signal for exchange
        if (stat.icalc == exchange_root && !std::ifstream{signal_file}.good()) {
            std::ofstream ofs{signal_file};
            ofs << "PROPOSED BY ROOT = " << stat.icalc;
            ofs.close();
        }
        // 2) for all: raise for paritcipation
        if (!exchanged && std::ifstream{signal_file}.good() && !std::ifstream{flag_file}.good()) {
            std::ofstream ofs{flag_file};
            ofs << stat.icalc << " " << istep[0] << "\n";
            ofs.close();
        }
        // 3) for ROOT: generate allowed_file for exchange
        if (stat.icalc == exchange_root && !std::ifstream{allowed_file}.good()) {
            system("cat EXCHANGE_THREAD_* > EXCHANGE_COUNT");
            list1.clear();
            std::ifstream ifs("EXCHANGE_COUNT");
            std::string   eachline;
            while (getline(ifs, eachline)) {
                std::stringstream ss(eachline);
                int               idx1, ival;
                ss >> idx1 >> ival;
                list1.push_back(std::pair<int, int>(idx1, ival));
            }
            if (list1.size() >= exchange_num) {
                // ordering icalc number
                std::sort(list1.begin(), list1.end(),                                     //
                          [](const std::pair<int, int>& a, const std::pair<int, int>& b)  //
                          { return a.first < b.first; });
                // ordering istep
                list2 = list1;
                std::sort(list2.begin(), list2.end(),                                     //
                          [](const std::pair<int, int>& a, const std::pair<int, int>& b)  //
                          { return a.second < b.second; });

                std::ofstream ofs(allowed_file);
                for (int i = 0; i < list1.size(); ++i) { ofs << list1[i].first << " " << list2[i].first << "\n"; }
                ofs.close();

                std::cout << "SUMMARY OF EXCHANGE\n";
                for (int i = 0; i < list1.size(); ++i) {
                    std::cout << "ID " << list1[i].first << " (with istep " << list1[i].second << ")"  //
                              << "-> ID " << list2[i].first << " (with istep " << list2[i].second << ")\n";
                }
            } else {
                if (icyc > maxcyc) {
                    remove(flag_file.c_str());
                    remove(signal_file.c_str());
                    std::cout << "ROOT: GIVE UP GATHERING FOR THE NOW RAISED NO. = " << list1.size()  //
                              << " < " << exchange_num << " REQUIRED\n";
                    break;
                }
            }
        }
        // 4) for ROOT: discontruct allowed_file for exchange
        if (stat.icalc == exchange_root && std::ifstream{allowed_file}.good()) {
            bool all_flag_is_deleted = true;
            for (int i = 0; i < list1.size(); ++i) {
                if (std::ifstream{utils::concat("EXCHANGE_THREAD_", i)}.good()) {
                    all_flag_is_deleted = false;
                    break;
                }
            }
            // only all flag file in allowed_file are removed, then allowed_file is removed
            if (all_flag_is_deleted) {
                std::cout << "END OF EXCHANGE\n";
                remove(signal_file.c_str());
                remove(allowed_file.c_str());
                engaged = false;
                break;
            }
        }
        usleep(100000);
    }
    remove(flag_file.c_str());
    if (stat.icalc == exchange_root) {
        remove(signal_file.c_str());
        remove(allowed_file.c_str());
    }
    if (icyc >= maxcyc && stat.icalc == exchange_root) {
        // if you want recover exchange during simulation. just delete the NOCHANGED file.
        system(utils::concat("touch ", nosignal_file).c_str());
        std::cout << "EXCAHNGE FAILED! SHUT OFF EXCHANGE! ";
        std::cout << "PlEASE CHECK YOUR SETUP!\n";
        // use_exchange = false; // keep use_exchange = true, use NOCHANGED to block temporary exchange
    }

    // ex_end            = std::chrono::steady_clock::now();
    // total_time        = static_cast<std::chrono::duration<double>>(ex_end - ex_begin).count();
    // exchange_fulltime = total_time + exchange_time;

    return stat;
}

Status& Kernel_Iterative_Adapt::executeKernel_impl(Status& stat) {
    std::cout << "S|: "                       // Status of Adaptive Integrator
              << std::setw(10) << "Progress"  // Progress before one step
              << std::setw(10) << "Time"      // Time before one step
              << std::setw(10) << "tsize"     // tsize before one step
              << std::setw(10) << "dtsize"    // stepsize before this step
              << std::setw(10) << "try"       // stepsize after this step
              << std::endl;
    if (_param->get_bool({"restart"}, LOC(), false)) {
        stat.first_step = false;
    } else {
        stat.first_step = true;
    }
    isamp[0] = istep[0] / sstep;
    while (istep[0] <= nstep) {
        stat.istep = istep[0];  // @TODO CAUTION!!!
        if (use_exchange) exchange(stat);

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
            stat.frozen     = true;
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
        // stat.fail_type = 0; // we may use the value in last last step, so don't reset it
        for (auto& pkernel : _child_kernels) { pkernel->executeKernel(stat); }

        char statc = '?';
        if (stat.succ && !stat.last_attempt) statc = 'T';                     // SUCC=True
        if (stat.succ && stat.last_attempt) statc = 'R';                      // Recovered
        if (!stat.succ && stat.last_attempt) statc = 'X';                     // END HERE
        if (!stat.succ && !stat.last_attempt && dtsize[0] > 1) statc = 'F';   // FAILED
        if (!stat.succ && !stat.last_attempt && dtsize[0] == 1) statc = 'L';  // FAILED and TODO TRY LAST
        if (stat.frozen) statc = 'Z';                                         // FROZEN THEN

        if (std::ifstream{"X_STAT"}.good() && !stat.frozen) statc = 'X';
        if (std::ifstream{utils::concat("X_STAT", stat.icalc)}.good() && !stat.frozen) statc = 'X';

        switch (statc) {
            case 'X': {
                for (auto& fname : backup_fields) {
                    _dataset->_def(utils::concat("integrator.", fname),     //
                                   utils::concat("backup.", 1, ".", fname)  //
                    );
                    for (int bto = nbackup, bfrom = bto - 1; bto > 1; --bto, --bfrom) {
                        _dataset->_def(utils::concat("backup.", bfrom, ".", fname),
                                       utils::concat("backup.", bto, ".", fname));
                    }
                }

                _dataset->def_int("restart.istep", istep.data(), 1);
                _dataset->def_int("restart.tsize", tsize.data(), 1);
                _dataset->def_int("restart.dtsize", dtsize.data(), 1);
                _dataset->def_int("restart.last_tried_dtsize", last_tried_dtsize.data(), 1);

                // save breakdown information
                std::ofstream ofs{utils::concat(directory, "/frozen", stat.icalc, "-", tsize[0], ".ds")};
                _dataset->dump(ofs);
                ofs.close();
                stat.frozen = true;
            }
            case 'T':
            case 'R':
            case 'Z': {
                stat.last_attempt = false;
                stat.first_step   = false;

                if (at_fullstep_finally) istep[0]++;
                tsize[0] += dtsize[0];

                int extend_dtsize    = (at_fullstep_finally) ? 2 * last_tried_dtsize[0] : 2 * dtsize[0];
                int remain_dtsize    = msize - (tsize[0] % msize);
                int new_dtsize       = std::min({msize, extend_dtsize, remain_dtsize});
                last_tried_dtsize[0] = dtsize[0];
                dtsize[0]            = new_dtsize;
                break;
            }
            case 'L':
            case 'F': {
                stat.last_attempt    = (dtsize[0] == 1);
                last_tried_dtsize[0] = dtsize[0];
                dtsize[0]            = (stat.last_attempt) ? dtsize[0] : dtsize[0] / 2;
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
                  << std::setw(10) << last_tried_dtsize[0]                      //
                  << std::setw(10) << dtsize[0] << std::endl; // flush into log
        isamp[0] = istep[0] / sstep;
    }
    return stat;
}

};  // namespace PROJECT_NS
