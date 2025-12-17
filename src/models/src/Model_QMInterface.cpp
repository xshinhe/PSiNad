#include "psnd/Model_QMInterface.h"

#include <unistd.h>

#include <cstdlib>
#include <sstream>

#include "ghc/filesystem.hpp"
#include "psnd/Kernel_Representation.h"
#include "psnd/chem.h"
#include "psnd/debug_utils.h"
#include "psnd/hash_fnv1a.h"
#include "psnd/linalg.h"
#include "psnd/macro_utils.h"
#include "psnd/vars_list.h"

static std::string toLower(const std::string& input) {
    std::string result = input;  // 创建一个副本
    std::transform(result.begin(), result.end(), result.begin(), [](unsigned char c) { return std::tolower(c); });
    return result;
}


inline int removeFile(const std::string& filename) { return remove(filename.c_str()); }

inline void clearFile(const std::string& filename) { std::ofstream clear(filename, std::ios::trunc); }

inline void closeOFS(std::ofstream& ofs) {
    if (ofs.is_open()) ofs.close();
}

inline bool isFileExists(const std::string& name) { return std::ifstream{name.c_str()}.good(); }

namespace PROJECT_NS {

const std::string Model_QMInterface::getName() { return "Model_QMInterface"; }

int Model_QMInterface::getType() const { return utils::hash(FUNCTION_NAME); }

void Model_QMInterface::setInputParam_impl(std::shared_ptr<Param> PM) {
    Kernel_Representation::onthefly = true;

    qm_string       = _param->get_string({"model.qm_flag"}, LOC(), "MNDO");
    qm_config_in    = _param->get_string({"model.qm_config"}, LOC(), "QM.in");
    save_every_calc = _param->get_bool({"model.qm_save_every_calc"}, LOC(), true);
    save_every_step = _param->get_bool({"model.qm_save_every_step"}, LOC(), false);
    sstep_dataset   = _param->get_int({"model.sstep_dataset"}, LOC(), 0);
    use_state_detection = _param->get_bool({"model.use_state_detection"}, LOC(), false); // 动力学减小计算量

    nac_threshold = _param->get_real({"model.nac_threshold"}, LOC(), 1.0);

    char* p = getenv("PSND_PYTHON");
    if (p != nullptr) pypsnd_path = p;
    if (pypsnd_path == "" || !isFileExists(utils::concat(pypsnd_path, "/", "QM.py")))
        throw psnd_error("please correctly setup env: PSND_PYTHON");

    if (!isFileExists(qm_config_in))
        throw psnd_error("QM config not found, please set model.qm_config = <your config file>");

    // read temperature
    double temperature = _param->get_real({"model.temperature"}, LOC(), phys::temperature_d, 1.0f);
    beta               = 1.0f / (phys::au::k * temperature);  // don't ignore k_Boltzman
    // read task
    init_nuclinp = _param->get_string({"model.init_nuclinp"}, LOC(), "#hess");
    time_unit    = _param->get_real({"model.time_unit", "solver.time_unit"}, LOC(), phys::time_d, 1.0f);
}

void Model_QMInterface::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    x = DS->def(DATA::integrator::x);
    p = DS->def(DATA::integrator::p);

    x0      = DS->def(DATA::model::x0);
    p0      = DS->def(DATA::model::p0);
    w       = DS->def(DATA::model::w);
    x_sigma = DS->def(DATA::model::x_sigma);
    p_sigma = DS->def(DATA::model::p_sigma);

    // model field
    atoms            = DS->def(DATA::model::atoms);
    mass             = DS->def(DATA::model::mass);
    vpes             = DS->def(DATA::model::vpes);
    grad             = DS->def(DATA::model::grad);
    // hess             = DS->def(DATA::model::hess);
    // Tmod             = DS->def(DATA::model::Tmod);
    // f_r              = DS->def(DATA::model::f_r);
    // f_p              = DS->def(DATA::model::f_p);
    // f_rp             = DS->def(DATA::model::f_rp);
    osc_strength = DS->def_real("model.osc_strength", Dimension::F, "Oscillator strengths for transitions");
    V                = DS->def(DATA::model::V);
    dV               = DS->def(DATA::model::dV);
    eig              = DS->def(DATA::model::rep::eig);
    T                = DS->def(DATA::model::rep::T);
    dE               = DS->def(DATA::model::rep::dE);
    nac              = DS->def(DATA::model::rep::nac);
    nac_prev         = DS->def(DATA::model::rep::nac_prev);
    succ_ptr         = DS->def(DATA::control::succ);
    last_attempt_ptr = DS->def(DATA::control::last_attempt);
    frez_ptr         = DS->def(DATA::control::frez);
    fail_type_ptr    = DS->def(DATA::control::fail_type);
    dt_ptr           = DS->def(DATA::control::dt);
    t_ptr            = DS->def(DATA::control::t);
    istep_ptr        = DS->def(DATA::control::istep);

    ARRAY_EYE(T.data(), Dimension::F);

    
    double        dtmp;
    int           itmp;
    std::string   stmp;
    std::ifstream ifs(qm_config_in);
    getline(ifs, stmp, '\n');
    getline(ifs, stmp, '\n');
    getline(ifs, stmp, '\n');
    if (std::stringstream{stmp} >> itmp) natom = itmp;
    std::cout << "number of atom: " << stmp << std::endl;
    psnd_assert(natom * 3 == Dimension::N, "Dimension Error");
    getline(ifs, stmp, '\n');
    for (int iatom = 0, idx = 0; iatom < natom; ++iatom) {
        if (ifs >> stmp) atoms[iatom] = chem::getElemIndex(stmp);
        for (int a = 0; a < 3; ++a) {
            if (ifs >> dtmp) x0[idx] = dtmp / phys::au_2_ang;
            mass[idx] = chem::getElemMass(atoms[iatom]) / phys::au_2_amu;  // 计算原子质量  提前
            idx++;
        }
    }
    std::cout << "x0: " << x0[0] << " " << x0[1] << " " << x0[2] << "\n";

    config_content = "";
    while (getline(ifs, stmp, '\n')) config_content += stmp + "\n";
    ifs.close();

    // 初始化读取hessian 提前
    std::string read_hess = _param->get_string({"model.read_hess", "solver.read_hess"}, LOC(), "NULL");
    if (read_hess != "NULL") {  // used for sampling
        if (!isFileExists(read_hess)) throw psnd_error("cannot open hess as .ds file");
        hess             = DS->def(DATA::model::hess);
        Tmod             = DS->def(DATA::model::Tmod);
        ARRAY_EYE(Tmod.data(), Dimension::N);
        f_r              = DS->def(DATA::model::f_r);
        f_p              = DS->def(DATA::model::f_p);
        f_rp             = DS->def(DATA::model::f_rp);
        std::ifstream ifs(read_hess);  // prepare for NMA
        bool          read_w = false;
        bool          read_H = false;
        bool          read_T = false;
        std::string   eachline;
        while (getline(ifs, eachline)) {
            if (eachline.find("model.vpes") != eachline.npos) {
                getline(ifs, eachline);
                for (int i = 0; i < 1; ++i) ifs >> ener_refered;
            }
            if (eachline.find("model.x0") != eachline.npos && false) {  // read from config_content
                getline(ifs, eachline);
                for (int i = 0; i < Dimension::N; ++i) ifs >> x0[i];
            }
            if (eachline.find("model.p0") != eachline.npos) {
                getline(ifs, eachline);
                for (int i = 0; i < Dimension::N; ++i) ifs >> p0[i];
            }
            if (eachline.find("model.hess") != eachline.npos) {
                read_H = true;
                getline(ifs, eachline);
                for (int i = 0; i < Dimension::NN; ++i) ifs >> hess[i];
            }
            if (eachline.find("model.w") != eachline.npos) {
                read_w = true;
                getline(ifs, eachline);
                for (int i = 0; i < Dimension::N; ++i) ifs >> w[i];
            }
            if (eachline.find("model.Tmod") != eachline.npos) {
                read_T = true;
                getline(ifs, eachline);
                for (int i = 0; i < Dimension::NN; ++i) ifs >> Tmod[i];
            }
        }
        ifs.close();
        if (read_H && !read_w) { EigenSolve(w.data(), Tmod.data(), hess.data(), Dimension::N); }
        if (!read_H && !read_w && !read_T) throw psnd_error("cannot read hess from ds");
    } else {
        for (int j = 0; j < Dimension::N; ++j) x_sigma[j] = 0.0e0, p_sigma[j] = 0.0e0;
    }
}

Status& Model_QMInterface::initializeKernel_impl(Status& stat) {
    try_level = 0;
    return stat;
}

Status& Model_QMInterface::executeKernel_impl(Status& stat) {
    if (stat.frozen) return stat;

    // prepare path for calculation
    // @bug: ghc::filesystem is not reliable in MPI ENV
    std::string path_str;
    if (save_every_calc) {
        path_str = utils::concat(directory, "/QM-", stat.icalc);
    } else {
        path_str = directory + "/QM";
    }
    ghc::filesystem::path path(path_str);
    if (!ghc::filesystem::is_directory(path)) { ghc::filesystem::create_directory(path); }

    // detect if there is a STOP file
    if (isFileExists(utils::concat(path_str, "/STOP"))) {
        stat.succ      = false;  // force stop
        stat.fail_type = 1;      // force stop
        // DON'T FROZEN HERE!!!
        return stat;
    }

    // prepare input run for calculation
    std::string tmp_input;
    if (save_every_step) {
        tmp_input = utils::concat(path_str, "/QM.run.", istep_ptr[0]);
    } else {
        tmp_input = utils::concat(path_str, "/QM.run");
    }

    // convert AU to Angstrom
    for (int i = 0; i < Dimension::N; ++i) x[i] *= phys::au_2_ang;


    // check if x has NaN
    for (int i = 0; i < Dimension::N; ++i) {
        if (std::isnan(x[i])) { throw psnd_error("x has NaN, please check your input file or QM code"); }
    }


    // write input
    std::ofstream ofs(tmp_input);
    ofs << "[GEOM]" << "\n";
    ofs << "xyz = \"\"\"\n";
    ofs << natom << "\n";
    ofs << "Autogenerated at t = " << t_ptr[0] / time_unit << "\n";
    for (int iatom = 0, idx = 0; iatom < natom; ++iatom) {
        ofs << chem::getElemLabel(atoms[iatom]);  //
        for (int a = 0; a < 3; ++a) {
            ofs << FMT(8) << x[idx];
            idx++;
        }
        ofs << "\n";
    }
    ofs << "\n";
    ofs << config_content;
    ofs.close();

    // determine the level of calculation (larger level, more loose convergence criteria)
    if (stat.last_attempt && stat.fail_type == 1) {
        try_level++;
    } else {
        try_level = 0;
    }

    // call python executation
    //     // first lower the qm_String
    std::string  qm_string_lower = toLower(qm_string);
    auto         occ_nuc         = _dataset->def(DATA::integrator::occ_nuc);
    const double deltaE_thres    = nac_threshold / phys::au_2_ev;  // default in 1.0 eV in atomic unit
    // 计算上一步中其他态跟占据态 occ 之间的能量差，能量差小于nac_threshold eV才计算这两个态的耦合，注意eig是原子单位，需要转换为eV 
    std::vector<std::pair<int, int>> occ_pairs;
    for (int i = 0; i < Dimension::F; ++i) {
        int j = occ_nuc[0];
        if (i != j) {
            double deltaE = eig[i] - eig[j];
            if (std::abs(deltaE) < deltaE_thres) {  // 0.3 eV
                std::cout << "[QMInterface] Coupling between state " << i << " and occupied state " << j
                          << " is considered, ΔE = " << deltaE * 27.2114 << " eV\n";
                occ_pairs.emplace_back(i, j);
            }
        }
    }
    // 把occ pairs 转化为字符串的形式 如(0,1), (1,2) -> "01 12"
    // @ambiguous: Does '112' represent (11,2) or (1,12)?
    // @hexin
    std::string occ_pairs_str;
    for (const auto& pair : occ_pairs) { occ_pairs_str += utils::concat(pair.first, ",", pair.second, " "); }

    std::string qm_call_str;
    if (use_state_detection){
        qm_call_str = utils::concat("python ", pypsnd_path, "/QM.py -t ", try_level,  //
                                            " -d ", path_str, " -i ", tmp_input, " -qm ", qm_string_lower, " -occ ",
                                            occ_nuc[0], " -ncouple ", occ_pairs_str);
    } else {
        qm_call_str = utils::concat("python ", pypsnd_path, "/QM.py -t ", try_level,  //
                                            " -d ", path_str, " -i ", tmp_input, " -qm ", qm_string_lower);
    }

    int         s           = system(qm_call_str.c_str());

    // checkout the result
    if (s == 0 && isFileExists(utils::concat(path_str, "/interface.ds"))) {
        std::ifstream ifs;
        int           stat_number = 1;
        std::string   eachline;

        // all quantities are needed in AU
        ifs.open(utils::concat(path_str, "/interface.ds"));
        while (getline(ifs, eachline)) {
            if (eachline.find("interface.stat") != eachline.npos) {
                getline(ifs, eachline);
                for (int i = 0; i < 1; ++i) ifs >> stat_number;
            }
            if (eachline.find("interface.eig") != eachline.npos) {
                getline(ifs, eachline);
                for (int i = 0; i < Dimension::F; ++i) ifs >> eig[i];
            }
            if (eachline.find("interface.dE") != eachline.npos) {
                getline(ifs, eachline);
                for (int j = 0, jFF = 0; j < Dimension::N; ++j, jFF += Dimension::FF) {
                    for (int i = 0, jii = jFF; i < Dimension::F; ++i, jii += Dimension::Fadd1) ifs >> dE[jii];
                }
            }
            if (eachline.find("interface.nac") != eachline.npos) {
                getline(ifs, eachline);
                for (int jik = 0; jik < Dimension::NFF; ++jik) ifs >> nac[jik];
            } else {
                for (int jik = 0; jik < Dimension::NFF; ++jik) nac[jik] = 0;
            }
            if (eachline.find("interface.strength") != eachline.npos) {
                getline(ifs, eachline);
                for (int i = 0; i < Dimension::F; ++i) ifs >> osc_strength[i];
            }
        }
        std::string command;
        if (sstep_dataset > 0 && istep_ptr[0] % sstep_dataset == 0 && stat_number == 0) {
            command = utils::concat("cp ", path_str, "/interface.ds ",  //
                                    path_str, "/interface-", istep_ptr[0], ".ds ");
            system(command.c_str());
            if (!save_every_step) {  // also save structure
                command = utils::concat("cp ", tmp_input, " ", tmp_input, ".", istep_ptr[0]);
                system(command.c_str());
            }
        }
        if (stat_number == 0) {
            command = utils::concat("mv ", path_str, "/interface.ds ", path_str, "/interface-old.ds ");
            system(command.c_str());
        }

        if (stat_number != 0) {
            stat.succ      = false;
            stat.fail_type = 1;
        } else {
            // in last_attempt, even though succeed, we don't reset fail_type
            if (stat.fail_type == 1 && !stat.last_attempt) stat.fail_type = 0;
        }
    } else {
        if (s != 0) std::cout << "psnd external shell status bug\n";
        if (!isFileExists(utils::concat(path_str, "/interface.ds"))) std::cout << "interface.ds is not generated\n";
        stat.succ      = false;
        stat.fail_type = 1;
    }

    if (stat.succ) {
        if (!stat.first_step) track_nac_sign();  // @note track_nac_sign is important
        for (int i = 0, idx = 0; i < Dimension::N; ++i) {
            for (int j = 0; j < Dimension::F; ++j) {
                for (int k = 0; k < Dimension::F; ++k, ++idx) {
                    if (j == k) continue;
                    dE[idx] = nac[idx] * (eig[k] - eig[j]);
                }
            }
        }
    }
    for (int i = 0; i < Dimension::N; ++i) x[i] /= phys::au_2_ang;
    return stat;
}

int Model_QMInterface::track_nac_sign() {
    for (int i = 0; i < Dimension::F; ++i) {
        for (int j = 0; j < Dimension::F; ++j) {  // check if NAC(:,i,j) should flip its sign
            if (i == j) continue;

            const double norm_eps = 10e-14;
            double       norm_old = 0.0f;
            double       norm_new = 0.0f;
            double       cosangle = 0.0f;
            int          IJ       = i * Dimension::F + j;
            for (int k = 0, idx = IJ; k < Dimension::N; ++k, idx += Dimension::FF) {
                norm_old += nac_prev[idx] * nac_prev[idx];
                norm_new += nac[idx] * nac[idx];
                cosangle += nac_prev[idx] * nac[idx];
            }
            norm_old = sqrt(norm_old);
            norm_new = sqrt(norm_new);
            if (norm_old < norm_eps || norm_new < norm_eps) {
                cosangle = 1.0f;
            } else {
                cosangle = cosangle / (norm_old * norm_new);
            }

            if (norm_new > 10e12 || norm_old > 10e12) {
                for (int k = 0; k < Dimension::N; ++k) {
                    nac[k * Dimension::FF + i * Dimension::F + j] =
                        copysign(nac[k * Dimension::FF + i * Dimension::F + j],
                                 nac_prev[k * Dimension::FF + i * Dimension::F + j]);
                }
            } else if (cosangle < 0) {  // in this case we flip the sign of NAC(:,i,j)
                for (int k = 0; k < Dimension::N; ++k) { nac[k * Dimension::FF + i * Dimension::F + j] *= -1; }
            }
        }
    }
    for (int i = 0; i < Dimension::NFF; ++i) nac_prev[i] = nac[i];  // save a copy
    return 0;
}
};  // namespace PROJECT_NS
