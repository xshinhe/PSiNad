#include "kids/Model_QMInterface.h"

#include <cstdlib>
#include <sstream>

#include "ghc/filesystem.hpp"
#include "kids/Kernel_Representation.h"
#include "kids/chem.h"
#include "kids/debug_utils.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

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

    std::string qm_string = _param->get_string({"model.qm_flag"}, LOC(), "MNDO");
    qm_config_in          = _param->get_string({"model.qm_config"}, LOC(), "QM.in");
    save_every_calc       = _param->get_bool({"model.qm_save_every_calc"}, LOC(), true);
    save_every_step       = _param->get_bool({"model.qm_save_every_step"}, LOC(), false);
    qm_type               = QMPolicy::_from(qm_string);

    char* p = getenv("KIDS_PYTHON");
    if (p != nullptr) pykids_path = p;
    if (pykids_path == "" || !isFileExists(utils::concat(pykids_path, "/", "QM.py")))
        throw kids_error("please correctly setup env: KIDS_PYTHON");

    if (!isFileExists(qm_config_in))
        throw kids_error("QM config not found, please set model.qm_config = <your config file>");

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
    hess             = DS->def(DATA::model::hess);
    Tmod             = DS->def(DATA::model::Tmod);
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
    succ_ptr         = DS->def(DATA::flowcontrol::succ);
    last_attempt_ptr = DS->def(DATA::flowcontrol::last_attempt);
    frez_ptr         = DS->def(DATA::flowcontrol::frez);
    fail_type_ptr    = DS->def(DATA::flowcontrol::fail_type);
    dt_ptr           = DS->def(DATA::flowcontrol::dt);
    t_ptr            = DS->def(DATA::flowcontrol::t);
    istep_ptr        = DS->def(DATA::flowcontrol::istep);

    ARRAY_EYE(T, Dimension::F);
    ARRAY_EYE(Tmod, Dimension::N);

    double        dtmp;
    int           itmp;
    std::string   stmp;
    std::ifstream ifs(qm_config_in);
    getline(ifs, stmp, '\n');
    if (std::stringstream{stmp} >> itmp) natom = itmp;
    kids_assert(natom * 3 == Dimension::N, "Dimension Error");
    getline(ifs, stmp, '\n');
    for (int iatom = 0, idx = 0; iatom < natom; ++iatom) {
        if (ifs >> stmp) atoms[iatom] = chem::getElemIndex(stmp);
        for (int a = 0; a < 3; ++a) {
            if (ifs >> dtmp) x0[idx] = dtmp / phys::au_2_ang;
            mass[idx] = chem::getElemMass(atoms[iatom]) / phys::au_2_amu;
            idx++;
        }
    }
    config_content = "";
    while (getline(ifs, stmp, '\n')) config_content += stmp + "\n";
    ifs.close();

    if (isFileExists("HESS.in")) {  // used for sampling
        std::ifstream ifs("HESS.in");
        for (int i = 0; i < Dimension::F; ++i) {
            if (ifs >> dtmp) ener_refered = dtmp;
        }
        for (int i = 0, ik = 0; i < Dimension::N; ++i) {
            for (int k = 0; k < Dimension::N; ++k, ++ik) {
                if (ifs >> dtmp) hess[ik] = dtmp;
            }
        }
        for (int i = 0; i < Dimension::N; ++i) {
            if (ifs >> dtmp) w[i] = dtmp;
        }
        for (int i = 0, ik = 0; i < Dimension::N; ++i) {
            for (int k = 0; k < Dimension::N; ++k, ++ik) {
                if (ifs >> dtmp) Tmod[ik] = dtmp;
            }
        }
        // // from hessian & temperature to prepare initial sampling
        // for (int j = 0; j < Dimension::N; ++j) {
        //     if (j < 6 || w[j] < 1.0e-10) {  // cutoff low frequency
        //         x_sigma[j] = 0.0f;
        //         p_sigma[j] = 0.0f;
        //     } else {  // NOTE: it's for normal-mode!
        //         double Qoverbeta = 0.5f * w[j] / std::tanh(0.5f * beta * w[j]);
        //         if (classical_bath) Qoverbeta = 1.0e0 / beta;
        //         x_sigma[j] = std::sqrt(Qoverbeta / (w[j] * w[j]));
        //         p_sigma[j] = std::sqrt(Qoverbeta);
        //     }
        // }
        ifs.close();
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

    for (int i = 0; i < Dimension::N; ++i) x[i] *= phys::au_2_ang;

    std::string path_str;
    if (save_every_calc) {
        path_str = utils::concat(directory, "/QM-", stat.icalc);
    } else {
        path_str = directory + "/QM";
    }

    ghc::filesystem::path path(path_str);
    if (!ghc::filesystem::is_directory(path)) {
        std::cout << LOC() << "\n";
        ghc::filesystem::create_directory(path);
    }

    if (isFileExists(utils::concat(path_str, "/STOP"))) {
        stat.frozen = true;  // force stop
        return stat;
    }

    std::string tmp_input;
    if (save_every_step) {
        tmp_input = utils::concat(path_str, "/QM.run.", istep_ptr[0]);
    } else {
        tmp_input = utils::concat(path_str, "/QM.run");
    }

    std::ofstream ofs(tmp_input);
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

    if (stat.last_attempt && stat.fail_type == 1) {
        try_level++;
    } else {
        try_level = 0;
    }

    std::string qm_call_str = utils::concat("python ", pykids_path, "/QM.py -t ", try_level,  //
                                            " -d ", path_str, " -i ", tmp_input);
    int         s           = system(qm_call_str.c_str());
    if (s != 0) {
        stat.succ = false;
        std::cout << LOC() << "\n";
        exit(0);
        return stat;
    }

    if (!isFileExists(utils::concat(path_str, "/stat.dat")) ||
        !isFileExists(utils::concat(path_str, "/energy.dat")) ||    //
        !isFileExists(utils::concat(path_str, "/gradient.dat")) ||  //
        !isFileExists(utils::concat(path_str, "/nacv.dat"))) {
        stat.succ = false;
        throw kids_error("DEBUG TEST");
        return stat;
    }

    std::ifstream ifs;
    ifs.open(utils::concat(path_str, "/stat.dat"));
    int         stat_number;
    std::string error_msg;
    ifs >> stat_number >> error_msg;
    ifs.close();
    if (stat_number == 0) {
        stat.succ = true;
        if (stat.last_attempt && stat.fail_type == 1) {
            std::cout << "survive in last try mndo\n";
        } else if (stat.last_attempt && stat.fail_type == 2) {
            std::cout << "mndo pass first, see next\n";
        }
        if (stat.fail_type == 1) stat.fail_type = 0;
    } else {
        stat.succ      = false;
        stat.fail_type = 1;  // failture from QM
        std::cout << "fail in calling MNDO! " << error_msg << "\n";
    }

    ifs.open(utils::concat(path_str, "/energy.dat"));
    for (int i = 0; i < Dimension::F; ++i) ifs >> eig[i];
    ifs.close();

    ifs.open(utils::concat(path_str, "/gradient.dat"));
    for (int j = 0, jFF = 0; j < Dimension::N; ++j, jFF += Dimension::FF) {
        for (int i = 0, jii = jFF; i < Dimension::F; ++i, jii += Dimension::Fadd1) ifs >> dE[jii];
    }
    ifs.close();

    ifs.open(utils::concat(path_str, "/nacv.dat"));
    for (int jik = 0; jik < Dimension::NFF; ++jik) ifs >> nac[jik];
    ifs.close();

    if (!stat.first_step) track_nac_sign();  // @note track_nac_sign is important
    for (int i = 0, idx = 0; i < Dimension::N; ++i) {
        for (int j = 0; j < Dimension::F; ++j) {
            for (int k = 0; k < Dimension::F; ++k, ++idx) {
                if (j == k) continue;
                dE[idx] = nac[idx] * (eig[k] - eig[j]);
            }
        }
    }
    // exit(0);

    for (int i = 0; i < Dimension::N; ++i) x[i] /= phys::au_2_ang;
    // output unit conversion is performed in python. everthing is au unit now.
    // for (int i = 0; i < Dimension::F; ++i) eig[i] /= phys::au_2_kcal_1mea;  ///< convert kcalpmol to Hartree
    // for (int i = 0; i < Dimension::NFF; ++i)
    //     dE[i] /= (phys::au_2_kcal_1mea / phys::au_2_ang);  ///< convert to Hartree/Bohr
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
