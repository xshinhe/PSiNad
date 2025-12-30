#include "psnd/Kernel_Iterative_Adapt.h"

#include <unistd.h>

#include <algorithm>
#include <chrono>

#include "psnd/hash_fnv1a.h"
#include "psnd/macro_utils.h"
#include "psnd/vars_list.h"

// @blame FMT is in psnd/fmt.h 不要重复造轮子 @hexin
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

    // add by hclu251026
    dump_in_string = _param->get_bool({"solver.dump_in_string"}, LOC(), false);
}

void Kernel_Iterative_Adapt::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    t                 = DS->def(DATA::control::t);
    dt                = DS->def(DATA::control::dt);
    istep             = DS->def(DATA::control::istep);
    isamp             = DS->def(DATA::control::isamp);  // 通过istep实时计算
    tsize             = DS->def(DATA::control::tsize);
    dtsize            = DS->def(DATA::control::dtsize);
    last_tried_dtsize = DS->def(DATA::control::last_tried_dtsize);
    at_condition      = DS->def(DATA::control::at_condition);

    // initializarion
    DS->def(VARIABLE<psnd_int>("control.sstep", &Dimension::shape_1, "@"))[0] = sstep;
    DS->def(VARIABLE<psnd_int>("control.nstep", &Dimension::shape_1, "@"))[0] = nstep;
    DS->def(VARIABLE<psnd_int>("control.nsamp", &Dimension::shape_1, "@"))[0] = nsamp;
    DS->def(VARIABLE<psnd_int>("control.msize", &Dimension::shape_1, "@"))[0] = msize;

    // backup for dt, used for recover
    DS->def(VARIABLE<psnd_real>("control.dt_backup", &Dimension::shape_1, "@"))[0] = dt0;
}

Status& Kernel_Iterative_Adapt::initializeKernel_impl(Status& stat) {
    if (_param->get_string({"load", "solver.load"}, LOC(), "").find(":continue") != std::string::npos) {  //
        if (_dataset_load == nullptr) throw psnd_error(utils::concat(LOC(), ": DataSet Load error"));
        if (std::ifstream{"X_STAT"}.good()) remove("X_STAT");
        if (std::ifstream{utils::concat("X_STAT", stat.icalc)}.good()) {
            std::string rmfile = utils::concat("X_STAT", stat.icalc);
            remove(rmfile.c_str());
        }

        // exactly copy from _dataset_load to _dataset
        // load istep即可 isamp通过istep实时计算

        // 方法1：直接检查recover节点是否存在
        if (_dataset_load->haskey("recover")) {
            std::cout << "[Kernel_Iterative_Adapt] Found recover node" << std::endl;
            istep[0]             = _dataset_load->def_int("recover.istep", 1)[0];
            tsize[0]             = _dataset_load->def_int("recover.tsize", 1)[0];
            dtsize[0]            = _dataset_load->def_int("recover.dtsize", 1)[0];
            last_tried_dtsize[0] = _dataset_load->def_int("recover.last_tried_dtsize", 1)[0];
        } else {
            std::cout << "[Kernel_Iterative_Adapt] No recover node found, the recover node is loaded from control "
                         "node. Be careful, when the trajectory finish sucessfully, control.dt will be set to ZERO."
                      << std::endl;
            istep[0]             = _dataset_load->def_int("control.istep", 1)[0];
            tsize[0]             = _dataset_load->def_int("control.tsize", 1)[0];
            dtsize[0]            = _dataset_load->def_int("control.dtsize", 1)[0];
            last_tried_dtsize[0] = _dataset_load->def_int("control.last_tried_dtsize", 1)[0];
        }


        stat.succ         = true;
        stat.last_attempt = false;
        stat.frozen       = false;
        stat.fail_type    = 0;
        return stat;
    }
    if (_param->get_string({"load", "solver.load"}, LOC(), "").find(":restart") != std::string::npos) {  //
        if (_dataset_load == nullptr) throw psnd_error(utils::concat(LOC(), ": DataSet Load error"));
        if (std::ifstream{"X_STAT"}.good()) remove("X_STAT");
        if (std::ifstream{utils::concat("X_STAT", stat.icalc)}.good()) {
            std::string rmfile = utils::concat("X_STAT", stat.icalc);
            remove(rmfile.c_str());
        }

        // exactly copy from _dataset_load to _dataset
        t[0]  = _dataset_load->def_real("control.t", 1)[0];
        dt[0] = dt0;

        isamp[0]             = 0;
        istep[0]             = 0;
        tsize[0]             = _dataset_load->def_int("control.tsize", 1)[0];
        dtsize[0]            = msize;
        last_tried_dtsize[0] = msize;
        stat.succ            = true;
        stat.last_attempt    = false;
        stat.frozen          = false;
        stat.fail_type       = 0;
        return stat;
    }
    if (use_exchange) {
        ex_begin          = std::chrono::steady_clock::now();
        exchange_fulltime = exchange_time;
    }

    // 普通初始化
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
                    _dataset->def_int("recover.istep", istep.data(), 1);
                    _dataset->def_int("recover.tsize", tsize.data(), 1);
                    _dataset->def_int("recover.dtsize", dtsize.data(), 1);
                    _dataset->def_int("recover.last_tried_dtsize", last_tried_dtsize.data(), 1);
                    std::ofstream ofs{dump_file};
                    _dataset->dump(ofs);
                    ofs.close();

                    // wait for idx2 to be prepared (including root)
                    while (!std::ifstream{load_file}.good()) { usleep(5000); }

                    // load idx2
                    std::ifstream ifs{load_file};
                    _dataset->load(ifs);
                    ifs.close();
                    istep[0]             = _dataset->def_int("recover.istep", 1)[0];
                    tsize[0]             = _dataset->def_int("recover.tsize", 1)[0];
                    dtsize[0]            = _dataset->def_int("recover.dtsize", 1)[0];
                    last_tried_dtsize[0] = _dataset->def_int("recover.last_tried_dtsize", 1)[0];
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
    if (_param->get_bool({"verbose"}, LOC(), false)) {
        std::cout << "S|: "                       // Status of Adaptive Integrator
                  << std::setw(10) << "Progress"  // Progress before one step
                  << std::setw(10) << "Time"      // Time before one step
                  << std::setw(10) << "tsize"     // tsize before one step
                  << std::setw(10) << "dtsize"    // stepsize before this step
                  << std::setw(10) << "try"       // stepsize after this step
                  << std::endl;
    }
    if (_param->get_string({"load", "solver.load"}, LOC(), "").find(":restart") != std::string::npos) {
        stat.first_step = false;
    } else {
        stat.first_step = true;
    }

    int count_fail_type1 = 0;

    // tsize 是一直累加的 一直到结束
    // dtsize 是每一步的步长大小，初始为msize，每一步根据情况调整

    isamp[0]          = istep[0] / sstep;
    int backup_dtsize = 0;
    while (istep[0] <= nstep) {
        stat.istep = istep[0];  // @TODO CAUTION!!!
        if (use_exchange) exchange(stat);

        int  tsize_before_loop     = tsize[0];              ///< current time-point tick
        int  tsize_after_loop      = tsize[0] + dtsize[0];  ///< next time-point tick after loop
        bool at_fullstep_initially = tsize_before_loop % (msize) == 0;
        bool at_fullstep_finally   = tsize_after_loop % (msize) == 0;  // 判断这一步跑完之后是不是恰好到一个整的格点
        at_condition[0]            = tsize_before_loop % (sstep * msize) == 0;
        t[0]                       = t0 + dt0 * (tsize[0] / ((double) msize));
        dt[0]                      = dt0 * (dtsize[0] / ((double) msize));

        // 当轨线结束时，dtsize和dt都被设为0， stat.frozen被设为true
        // 当stat.frozen被设为true时，里面的每个子kernel都会跳过计算
        // 跳出while循环之后，会把最后的dtsize和dt恢复回来
        if (istep[0] == nstep) {
            at_condition[0] = true;
            backup_dtsize   = dtsize[0];
            dtsize[0]       = 0;  // frozen dynamics
            dt[0]           = 0; //dt0 * (dtsize[0] / ((double) msize));
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
        if (!stat.succ && !stat.last_attempt && dtsize[0] > 1) statc = 'F';   // FAILED
        if (!stat.succ && !stat.last_attempt && dtsize[0] == 1) statc = 'L';  // FAILED and TODO TRY LAST
        if (!stat.succ && stat.last_attempt) statc = 'X';                     // END HERE
        if (!stat.succ && stat.fail_type == 1) {
            if (count_fail_type1 >= 2) {
                statc = 'X';
            } else {
                count_fail_type1++;  // don't change statc
            }
        }
        if (stat.fail_type != 1) count_fail_type1 = 0;  // reset to zero
        if (stat.frozen) statc = 'Z';                   // FROZEN THEN

        if (std::ifstream{"X_STAT"}.good() && !stat.frozen) statc = 'X';
        if (std::ifstream{utils::concat("X_STAT_", stat.icalc)}.good() && !stat.frozen) statc = 'X';

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

                _dataset->def_int("recover.istep", istep.data(), 1);
                _dataset->def_int("recover.tsize", tsize.data(), 1);
                _dataset->def_int("recover.dtsize", dtsize.data(), 1);
                _dataset->def_int("recover.last_tried_dtsize", last_tried_dtsize.data(), 1);

                // save breakdown information
                std::ofstream ofs{utils::concat(directory, "/frozen", stat.icalc, "-", tsize[0], ".ds")};
                _dataset->dump(ofs);
                ofs.close();
                std::ofstream ofs2{utils::concat(directory, "/frozen_last.ds")};
                _dataset->dump(ofs2);
                ofs2.close();
                stat.frozen = true;

                return stat;
            }
            case 'T':
            case 'R':
            case 'Z': {
                stat.last_attempt = false;
                stat.first_step   = false;

                if (at_fullstep_finally) istep[0]++;
                tsize[0] += dtsize[0];

                int extend_dtsize = (at_fullstep_finally) ? 2 * last_tried_dtsize[0] : 2 * dtsize[0];

                // in some case, don't enlarge the time size
                // if (statc == 'T' && stat.fail_type == -1) extend_dtsize    = (at_fullstep_finally) ?
                // last_tried_dtsize[0] : dtsize[0];

                int remain_dtsize    = msize - (tsize[0] % msize);
                int new_dtsize       = std::min({msize, extend_dtsize, remain_dtsize});
                last_tried_dtsize[0] = dtsize[0];
                dtsize[0]            = new_dtsize;
                break;
            }
            case 'L':
            case 'F': {
                if (!stat.succ) {
                    if (stat.fail_type == 1) {
                        count_fail_type1++;
                    } else {
                        count_fail_type1 = 0;
                    }
                }
                // save last tried
                stat.last_attempt    = (dtsize[0] == 1);
                last_tried_dtsize[0] = dtsize[0];

                // suggest new dt (don't minimize dt for fail_type==1)
                dtsize[0] = (stat.last_attempt || stat.fail_type == 1) ? dtsize[0] : dtsize[0] / 2;
                // tsize[0] += 0;

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
        if (_param->get_bool({"verbose"}, LOC(), false)) {
            std::cout << statc << "|:"  //
                      << std::resetiosflags(std::ios::scientific) << std::setiosflags(std::ios::fixed)
                      << std::setprecision(2) << std::setw(10) << 100 * t[0] / tend << "%"  //
                      << std::resetiosflags(std::ios::fixed) << std::setiosflags(std::ios::scientific)
                      << std::setprecision(2) << std::setw(10) << t[0] / time_unit  //
                      << std::setw(10) << tsize_before_loop                         //
                      << std::setw(10) << last_tried_dtsize[0]                      //
                      << std::setw(10) << dtsize[0] << std::endl;                   // flush into log
        }
        isamp[0] = istep[0] / sstep;


        // Dump dataset every successful step (or every N steps)
        // 为了每一步都能有一个checkpoint, 方便重跑

        // 每一步都记录一下recover字段，方便dump
        _dataset->def_int("recover.istep", istep.data(), 1);
        _dataset->def_int("recover.tsize", tsize.data(), 1);
        _dataset->def_int("recover.dtsize", dtsize.data(), 1);
        _dataset->def_int("recover.last_tried_dtsize", last_tried_dtsize.data(), 1);

        bool should_dump    = false;
        int  dump_frequency = _param->get_int({"solver.dump_frequency", "dump_frequency"}, LOC(), 0);

        if (dump_frequency > 0 && istep[0] % dump_frequency == 0) { should_dump = true; }

        std::string dump_filename = "dump";

        if (should_dump) {
            if (dump_in_string){
                try {
                    std::ofstream ofs{utils::concat(directory, "/", dump_filename, "-calc", stat.icalc, ".ds")};
                    _dataset->dump(ofs);
                    ofs.close();
                    std::cout << "Dumped step " << istep[0] << " to " << dump_filename << "-calc" << stat.icalc << ".ds"
                            << std::endl;

                } catch (std::runtime_error& e) {
                    std::cerr << "Warning: Failed to dump at step " << istep[0] << ": " << e.what() << std::endl;
                }
            } 

            try{
                _dataset->save_binary_filter(utils::concat(directory, "/", dump_filename, "-calc", stat.icalc, ".bin.ds"));
            } catch (std::runtime_error& e) { throw psnd_error(utils::concat(directory, "/", dump_filename, "-calc", stat.icalc, ".bin.ds")); }

        }
    }
    // 由于最后一步运行完了dt被设为0，所以需要恢复最后的dt
    dtsize[0] = backup_dtsize;  // frozen dynamics
    dt[0]     = dt0 * (dtsize[0] / ((double) msize));
    return stat;
}

};  // namespace PROJECT_NS
