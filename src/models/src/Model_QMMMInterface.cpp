#include "psnd/Model_QMMMInterface.h"

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

const std::string Model_QMMMInterface::getName() { return "Model_QMMMInterface"; }

int Model_QMMMInterface::getType() const { return utils::hash(FUNCTION_NAME); }

void Model_QMMMInterface::setInputParam_impl(std::shared_ptr<Param> PM) {
    Kernel_Representation::onthefly = true;

    qm_string       = _param->get_string({"model.qmmm_flag"}, LOC(), "MNDO");
    mm_string       = _param->get_string({"model.mm_flag"}, LOC(), "amber");
    qmmm_config_in  = _param->get_string({"model.qmmm_config"}, LOC(), "QMMM.in");
    qmmm_layer_info = _param->get_string({"model.qmmm_layer_info"}, LOC(), "layer_real.xyz");
    qmmm_top_file   = _param->get_string({"model.qmmm_top_file"}, LOC(), "real.top,model-H.top");
    crd_input       = _param->get_string({"model.qmmm_crd_input"}, LOC(), "real.crd");
    save_every_calc = _param->get_bool({"model.qmmm_save_every_calc"}, LOC(), true);
    save_every_step = _param->get_bool({"model.qmmm_save_every_step"}, LOC(), false);
    sstep_dataset   = _param->get_int({"model.sstep_dataset"}, LOC(), 0);

    char* p = getenv("PSNDQMMM_PYTHON");
    if (p != nullptr) psndqmmm_path = p;
    if (psndqmmm_path == "" || !isFileExists(utils::concat(psndqmmm_path, "/", "psndqmmm.py")))
        throw psnd_error("please correctly setup env: PSNDQMMM_PYTHON");

    if (!isFileExists(qmmm_config_in))
        throw psnd_error("QMMM config not found, please set model.qm_config = <your config file>");

    // read temperature
    double temperature = _param->get_real({"model.temperature"}, LOC(), phys::temperature_d, 1.0);
    beta               = 1.0 / (phys::au::k * temperature);  // don't ignore k_Boltzman
    // read task
    init_nuclinp = _param->get_string({"model.init_nuclinp"}, LOC(), "#hess");
    time_unit    = _param->get_real({"model.time_unit", "solver.time_unit"}, LOC(), phys::time_d, 1.0);
}

void Model_QMMMInterface::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    x = DS->def(DATA::integrator::x);
    p = DS->def(DATA::integrator::p);

    x0      = DS->def(DATA::model::x0);
    p0      = DS->def(DATA::model::p0);
    w       = DS->def(DATA::model::w);
    x_sigma = DS->def(DATA::model::x_sigma);
    p_sigma = DS->def(DATA::model::p_sigma);

    // model field
    atoms      = DS->def(DATA::model::atoms);
    layer_type = DS->def(DATA::model::layer_type);
    mass       = DS->def(DATA::model::mass);
    vpes       = DS->def(DATA::model::vpes);
    grad       = DS->def(DATA::model::grad);
    // hess             = DS->def(DATA::model::hess);
    // Tmod             = DS->def(DATA::model::Tmod);
    f_r              = DS->def(DATA::model::f_r);
    f_p              = DS->def(DATA::model::f_p);
    f_rp             = DS->def(DATA::model::f_rp);
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

    ARRAY_EYE(T.data(), Dimension::F);  // 为什么要在这儿初始化T矩阵？
    // ARRAY_EYE(Tmod.data(), Dimension::N);

    if (!isFileExists(qmmm_layer_info)) { throw psnd_error("qmm_layer_info is needed!"); }

    double        dtmp;
    int           itmp;
    std::string   stmp;
    std::ifstream ifs(qmmm_layer_info);
    int           iatom = 0;
    while (getline(ifs, stmp, '\n')) {
        std::stringstream ss{stmp};
        std::string       sym, type;
        double            charg;
        ss >> sym >> charg >> x0[3 * iatom] >> x0[3 * iatom + 1] >> x0[3 * iatom + 2] >> type;
        atoms[iatom] = chem::getElemIndex(sym);
        if (type == "H") layer_type[iatom] = 0;
        if (type == "M") layer_type[iatom] = 1;
        if (type == "L") layer_type[iatom] = 2;
        for (int a = 0; a < 3; ++a) {
            x0[3 * iatom + a] /= phys::au_2_ang;
            mass[3 * iatom + a] = chem::getElemMass(atoms[iatom]) / phys::au_2_amu;
        }
        iatom++;
    }
    // PRINT_ARRAY(mass, 1, 12);
    ifs.close();
    natom = iatom;
    if (Dimension::N != 3 * natom) throw psnd_error("mismatch real_layers size with N");
    for (int j = 0; j < Dimension::N; ++j) x_sigma[j] = 0.0e0, p_sigma[j] = 0.0e0;
}

Status& Model_QMMMInterface::initializeKernel_impl(Status& stat) {
    try_level = 0;
    return stat;
}

Status& Model_QMMMInterface::executeKernel_impl(Status& stat) {
    if (stat.frozen) return stat;

    // prepare path for calculation
    // @bug: ghc::filesystem is not reliable in MPI ENV
    std::string path_str;
    if (save_every_calc) {
        path_str = utils::concat(directory, "/QMMM-", stat.icalc);
    } else {
        path_str = utils::concat(directory, "/QMMM");
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
    if (save_every_step) {
        crd_input = utils::concat(path_str, "/real", istep_ptr[0], ".crd");
    } else {
        crd_input = utils::concat(path_str, "/real.crd");
    }
    std::cout << "Current working directory: " << ghc::filesystem::current_path() << std::endl;
    std::cout << "crd_input: " << crd_input << std::endl;
    std::cout << "QMMM -d path_str: " << path_str << std::endl;

    // convert AU to Angstrom
    for (int i = 0; i < Dimension::N; ++i) x[i] *= phys::au_2_ang;
    // PRINT_ARRAY(x, 1, 12);
    // std::cout << std::endl;

    // write input
    std::ofstream ofs(crd_input);
    ofs << "\n" << natom << "\n";

    for (int iatom = 0, idx = 0; iatom < natom; ++iatom) {
        // ofs << chem::getElemLabel(atoms[iatom]);  //
        for (int a = 0; a < 3; ++a) {
            ofs << std::fixed << std::setprecision(7) << std::setw(12) << x[idx] << " ";
            idx++;
        }
        if (iatom % 2 == 1) ofs << "\n";
    }
    ofs << "\n";
    ofs.close();

    // determine the level of calculation (larger level, more loose convergence criteria)
    if (stat.last_attempt && stat.fail_type == 1) {
        try_level++;
    } else {
        try_level = 0;
    }

    // call python executation
    // std::string qm_call_str =
    //     utils::concat("python ", psndqmmm_path, "/psndqmmm.py -t ", try_level,  //
    //                   " -d ", path_str, " -i ", qmmm_config_in, " -c ", crd_input, " > ", path_str, "/log");
    std::string qm_string_lower = toLower(qm_string);
    std::string qm_call_str =
        utils::concat("python ", psndqmmm_path, "/psndqmmm.py -t ", try_level, " -mm ", mm_string, " -qm ",
                      qm_string_lower, " -d ", path_str, " -i ", qmmm_config_in, " -c ", crd_input, " -l ",
                      qmmm_layer_info, " -top ", qmmm_top_file, " > ", path_str, "/log 2>&1");
    int s = system(qm_call_str.c_str());

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
            }
        }
        std::string command;
        if (sstep_dataset > 0 && istep_ptr[0] % sstep_dataset == 0 && stat_number == 0) {
            command = utils::concat("cp ", path_str, "/interface.ds ",  //
                                    path_str, "/interface-", istep_ptr[0], ".ds ");
            system(command.c_str());
            if (!save_every_step) {  // also save structure
                command = utils::concat("cp ", crd_input, " ", crd_input, ".", istep_ptr[0]);
                system(command.c_str());
                command = utils::concat("cp ", path_str, "/log ", path_str, "/log.", istep_ptr[0]);
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

int Model_QMMMInterface::track_nac_sign() {
    for (int i = 0; i < Dimension::F; ++i) {
        for (int j = 0; j < Dimension::F; ++j) {  // check if NAC(:,i,j) should flip its sign
            if (i == j) continue;

            const double norm_eps = 10e-14;
            double       norm_old = 0.0;
            double       norm_new = 0.0;
            double       cosangle = 0.0;
            int          IJ       = i * Dimension::F + j;
            for (int k = 0, idx = IJ; k < Dimension::N; ++k, idx += Dimension::FF) {
                norm_old += nac_prev[idx] * nac_prev[idx];
                norm_new += nac[idx] * nac[idx];
                cosangle += nac_prev[idx] * nac[idx];
            }
            norm_old = sqrt(norm_old);
            norm_new = sqrt(norm_new);
            if (norm_old < norm_eps || norm_new < norm_eps) {
                cosangle = 1.0;
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
