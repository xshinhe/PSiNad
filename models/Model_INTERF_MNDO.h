#ifndef MODEL_INTERF_MNDO_H
#define MODEL_INTERF_MNDO_H


#include <unistd.h>

#include <cmath>

#include "../core/Kernel.h"

struct MNDO99KW {
    std::string key;
    int val;
};

using MNDO99KW_map = std::map<std::string, int>;

class Model_Interf_MNDO final : public Nad_ForceField {
   public:
    inline virtual const std::string name() { return "Model_Interf_MNDO"; }

    Model_Interf_MNDO(){};

    virtual ~Model_Interf_MNDO(){};


   private:
    num_real temp, beta;
    bool diff_nac;
    std::string init_nuclinp, savename;
    // std::vector<std::string> mndo99_keyword, mndo99_comment, mndo99_data, mndo99_addition;
    std::vector<std::string> mndo99_data;
    std::string mndo99_keyword, mndo99_comment, mndo99_addition;

    std::vector<MNDO99KW> keyword;  // keyword wrapper

    std::string savefile_traj, savefile_ener, savefile_grad, savefile_nac;
    std::string savefile_RP, savefile_eac;

    num_real* Hsys;
    num_real* x_sigma;
    num_real* p_sigma;

    // integrator
    num_real *x, *p, *m, *w;

    // model
    num_real* mass;
    num_real *vpes, *grad, *hess;
    num_real *V, *dV, *ddV, *nacv;

    int natom;
    int read_flag;
    int ncigrd;
    int iroot;

    DEFINE_POINTER_PROTECTED(int, atoms);
    DEFINE_POINTER_PROTECTED(num_real, mod_Hess);
    DEFINE_POINTER_PROTECTED(num_real, mod_Tmat);
    DEFINE_POINTER_PROTECTED(num_real, nr_samp);
    DEFINE_POINTER_PROTECTED(num_real, np_samp);
    DEFINE_POINTER_PROTECTED(num_real, nac_prev);
    DEFINE_POINTER_PROTECTED(num_real, nac);
    DEFINE_POINTER_PROTECTED(num_real, ener);
    DEFINE_POINTER_PROTECTED(num_real, ener2);
    DEFINE_POINTER_PROTECTED(num_real, grad);


    void read_param_impl(Param* PM);

    void init_data_impl(DataSet* DS);

    void init_calc_impl(int stat);

    int exec_kernel_impl(int stat = -1);

    int parse_mndo99(const std::string& mndo99inp);
    std::string new_keyword(const MNDO99KW_map& newkeyword);
    int new_task(num_real* R, const std::string& file, const std::string& task_flag, const int& rdim);
    int track_nac_sign();
    int parse_standard(const std::string& log);
    int parse_hessian(const std::string& log);
    int calc_normalmode(num_real* R, const int& rdim);
    int calc_samp();
    int calc_scan();
};


Model_Interf_MNDO::read_param_impl(Param* PM) {
    type = ForceFieldOnTheFly;  // one-the-fly adiabatic forcefield

    tag = name() + "_" + tag;

    // reset unit convertor
    Param_Reset(iou.leng, phys::au_2_ang);
    Param_Reset(iou.mass, phys::au_2_amu);
    Param_Reset(iou.ener, phys::au_2_kcal_1mea);
    Param_Reset(iou.time, phys::au_2_fs);

    // parse mndo99 input
    std::string mndo99inp = PM->get<std::string>("mndo99inp", LOC(), "null");

    if (mndo99inp != "null") {
        natom = parse_mndo99(mndo99inp);
    } else {
        throw std::runtime_error("lack of file input for mndo");
    }

    CHECK_GT(F, 0) << "fdim should > 0";
    CHECK_EQ(N, 3 * natom) << "rdim is mismatch with natom";
    CHECK_LE(F, ncigrd) << "fdim should equal to (at least less than) ncigrd";
}

Model_Interf_MNDO::init_data_impl(DataSet* DS) {
    ALLOCATE_PTR_TO_VECTOR(atoms, natom);

    DS->reg<int>("model.atoms", Dimension::N);
    for (int i = 0, idx = 0, idxR = 0; i < natom; ++i) {
        atoms[i] = stoi(mndo99_data[idx++]);
        for (int j = 0; j < 3; ++j, ++idxR) {
            mod_M[idxR]  = ELEMENTS_MASS[atoms[i]] / iou.mass;   // convert amu to au
            mod_R0[idxR] = stod(mndo99_data[idx++]) / iou.leng;  // convert angstrom into au
            idx++;                                               // skip fixflag
        }
    }

    // parse temperature
    temp = Param_GetV(temp, PM, 300.0f) / iou.temp;
    beta = 1 / temp;

    // try to read/calculate hessian
    ALLOCATE_PTR_TO_VECTOR(nr_samp, N);
    ALLOCATE_PTR_TO_VECTOR(np_samp, N);
    ALLOCATE_PTR_TO_VECTOR(mod_Hess, NN);
    ALLOCATE_PTR_TO_VECTOR(mod_Tmat, NN);
    for (int i = 0; i < NN; ++i) mod_Hess[i] = 0.0f, mod_Tmat[i] = 0.0f;
    init_nuclinp = Param_GetT(std::string, PM, "init_nuclinp", "#hess");
    Param_GetV(read_flag, PM, 0);  // for init_nuclinp
    ALLOCATE_PTR_TO_VECTOR(nac_prev, NFF);
    ALLOCATE_PTR_TO_VECTOR(nac, NFF);
    ALLOCATE_PTR_TO_VECTOR(ener, F);
    ALLOCATE_PTR_TO_VECTOR(ener2, F);
    ALLOCATE_PTR_TO_VECTOR(grad, NF);

    if (init_nuclinp == "#calc") {
        calc_normalmode(mod_R0, N);
        exit(0);
    }
    if (init_nuclinp == "#scan") {
        calc_scan();
        exit(0);
    }
    if (init_nuclinp == "#samp") {
        calc_samp();
        exit(0);
    }
    if (init_nuclinp == "#hess") {
        std::ifstream ifs(".hess");  // read from .hess (saves in au)
        num_real tmp;
        if (ifs) {
            for (int i = 0, idx = 0; i < N; ++i)  // hessian
                for (int j = 0; j < N; ++j, ++idx) {
                    if (ifs >> tmp) mod_Hess[idx] = tmp;         // Hessian without mass-weighted (in au)
                    mod_Hess[idx] /= sqrt(mod_M[i] * mod_M[j]);  // make it mass-weighted Hessian
                }
            for (int i = 0; i < N; ++i)  // read frequency (in au)
                if (ifs >> tmp) mod_W[i] = tmp;
            for (int i = 0, idx = 0; i < N; ++i)  // read normal-mode vectors (in au)
                for (int j = 0; j < N; ++j, ++idx) {
                    if (ifs >> tmp) mod_Tmat[idx] = tmp;
                }
        } else {
            LOG(FATAL) << "loss of .hess in currect directory";
        }
        ifs.close();

        // from hessian & temperature to prepare initial sampling
        for (int j = 0; j < N; ++j) {
            if (j < 6 || mod_W[j] < 1.0e-10) {  // cutoff low frequency
                mod_sigmaR[j] = 0.0f;
                mod_sigmaP[j] = 0.0f;
            } else {  // NOTE: it's for normal-mode!
                double Qoverbeta = 0.5f * mod_W[j] / std::tanh(0.5f * beta * mod_W[j]);
                mod_sigmaR[j]    = std::sqrt(Qoverbeta / (mod_W[j] * mod_W[j]));
                mod_sigmaP[j]    = std::sqrt(Qoverbeta);
            }
        }
    }
}

/**
 * @brief ForceField_init for mndo99
 * @param
 *     nr:  nuclear configuration
 *     np:  nuclear momentum
 *     nm:  nuclear mass
 *   erho:  electronic density
 *   eeac:  electronic amplititude
 *   eocc:  electronic occupation
 *   rdim:  nuclear freedoms (=3*natom)
 *   fdim:  electronic freedoms
 *  itraj:  index of trajectory
 * @return int
 * @bug none
 */
int Model_Interf_MNDO::init_calc_impl(int stat) {
    if (itraj < 0) {
        eocc = mod_occ;
        for (int i = 0; i < N; ++i) nm[i] = mod_M[i];
        return 0;
    }

    if (init_nuclinp == "#hess") {  // assuming .hess is given

        // sampling normal-mode
        rand_gaussian(nr_samp, N);
        rand_gaussian(np_samp, N);
        for (int i = 0; i < 6; ++i) nr_samp[i] = 0.0f, np_samp[i] = 0.0f;
        for (int i = 0; i < N; ++i) {
            nm[i] = mod_M[i];
            nr_samp[i] *= mod_sigmaR[i];
            np_samp[i] *= mod_sigmaP[i];
        }

        // transfrom normal-mode to cartesian coordinates
        ARRAY_MATMUL(nr, mod_Tmat, nr_samp, N, N, 1);
        ARRAY_MATMUL(np, mod_Tmat, np_samp, N, N, 1);
        for (int i = 0; i < N; ++i) {
            nr[i] = nr[i] / std::sqrt(nm[i]) + mod_R0[i];
            np[i] = np[i] * std::sqrt(nm[i]) + mod_P0[i];
        }
    } else if (init_nuclinp == "#fix") {  // for initial md
        for (int i = 0; i < N; ++i) nr[i] = mod_R0[i], np[i] = 0.0f, nm[i] = mod_M[i];
    } else {
        num_real tmp;
        std::string stmp;

        LOG(INFO) << "opening " << utils::concat(init_nuclinp, itraj);
        std::ifstream ifs(utils::concat(init_nuclinp, itraj));
        if (!ifs.is_open()) LOG(FATAL) << "init_nuclinp cannot open";

        getline(ifs, stmp);  // skip atoms number
        getline(ifs, stmp);  // skip comments
        for (int iatom = 0, idx1 = 0, idx2 = 0; iatom < natom; ++iatom) {
            ifs >> stmp;                           // skip atomic flag
            for (int i = 0; i < 3; ++i, ++idx1) {  // read in [angstrom]
                if (ifs >> tmp) nr[idx1] = tmp / phys::au_2_ang;
            }
            for (int i = 0; i < 3; ++i, ++idx2) {  // read in [angstrom/ps]
                if (ifs >> tmp) np[idx2] = tmp / phys::au_2_angoverps * mod_M[idx2];
            }
        }
        for (int i = 0; i < N; ++i) nm[i] = mod_M[i];

        // ARRAY_SHOW(nr, N / 3, 3);
        // LOG(FATAL);
    }

    Nad_ForceField::ForceField_init_elec(erho, eeac, eocc, F, itraj);
    diff_nac = false;
    if (itraj >= 0) {  // open file
        savefile_traj = utils::concat("traj-", itraj, ".xyz");
        savefile_ener = utils::concat("ener-", itraj, ".dat");
        savefile_grad = utils::concat("grad-", itraj, ".dat");
        savefile_nac  = utils::concat("nac-", itraj, ".dat");
        savefile_RP   = utils::concat("RV-", itraj, ".dat");
        savefile_eac  = utils::concat("eac-", itraj, ".dat");
        if (FLAGS_r == "") {  // if not restart, clear old files
            utils::clearFile(savefile_traj);
            utils::clearFile(savefile_ener);
            utils::clearFile(savefile_grad);
            utils::clearFile(savefile_nac);
            utils::clearFile(savefile_RP);
            utils::clearFile(savefile_eac);
        }
    }
    return 0;
}

int Model_Interf_MNDO::exec_kernel_impl(int stat) {
    // convert atomic unit
    for (int i = 0; i < N; ++i) R[i] *= phys::au_2_ang;  ///< convert Bohr to Angstrom
    for (int i = 0; i < F; ++i) E[i] = 0;                // clear previous data
    for (int i = 0; i < NFF; ++i) dE[i] = 0;             // clear previous data

    // prepare a mndo99 input file
    std::string inpfile = utils::concat(".mndo99inp.", mpi_rank);
    std::string outfile = utils::concat(".mndo99out.", mpi_rank);
    new_task(R, inpfile, "def", rdim);  ///< generate a default task
    std::string cmd_exe = utils::concat("mndo99 < ", inpfile, " > ", outfile);
    int stat            = system(cmd_exe.c_str());
    if (stat != 0) return stat;

    parse_standard(outfile);  // parse in MNDO's units

    /**
     * @note track_nac_sign is important!
     */
    track_nac_sign();

    /**
     * @brief Construct Force Matrix
     * @details dE[:,i,j] = nacv[:,i,j] * (E[j] - E[i]), Force Matrix is useful to
     * evolve nuclear trajectories
     */
    for (int i = 0; i < F; ++i) E[i] = ener[i];
    for (int i = 0, idx = 0; i < N; ++i) {
        for (int j = 0; j < F; ++j) {
            for (int k = 0; k < F; ++k, ++idx) {
                if (j == k) {
                    dE[idx] = grad[i * F + j];
                } else {
                    dE[idx] = nac[idx] * (ener[k] - ener[j]);
                }
            }
        }
    }

    // convert units
    for (int i = 0; i < N; ++i) R[i] /= phys::au_2_ang;                              ///< convert Angstrom to Bohr
    for (int i = 0; i < F; ++i) E[i] /= phys::au_2_kcal_1mea;                        ///< convert kcalpmol to Hartree
    for (int i = 0; i < NFF; ++i) dE[i] /= (phys::au_2_kcal_1mea / phys::au_2_ang);  ///< convert to Hartree/Bohr

    if (!ARRAY_ISFINITE(R, N) || !ARRAY_ISFINITE(E, F) || !ARRAY_ISFINITE(dE, NFF)) {
        LOG(WARNING) << "Problem for calling forcefield";
        stat = -1;
    }
    return stat;
}


namespace MNDO99_PARSE {
enum _enum { KEYWORD, COMMENT, DATA, ADDITION };
};  // namespace MNDO99_PARSE

/**
 * @brief this function parse mndo99 input (only support cartesian format)
 * @return int: total atom number parsed
 * @todo:   support zmat format parser
 */
int Model_Interf_MNDO::parse_mndo99(const std::string& mndo99inp) {
    int parse_type = MNDO99_PARSE::KEYWORD;

    std::stringstream mndo99_keyword_sstr;
    std::stringstream mndo99_comment_sstr;
    std::stringstream mndo99_addition_sstr;

    int count_atom = 0;
    std::ifstream ifs(mndo99inp);
    std::string eachline;
    while (getline(ifs, eachline)) {
        eachline.erase(eachline.find_last_not_of(" ") + 1);  // remove last blank
        switch (parse_type) {
            case MNDO99_PARSE::KEYWORD: {
                mndo99_keyword_sstr << eachline << std::endl;
                std::string data, key, val;
                std::istringstream input(eachline);
                while (input >> data) {
                    if (data == "+") continue;
                    std::istringstream pair(data);  // "key=val" pattern
                    getline(pair, key, '=');
                    getline(pair, val, '=');
                    std::for_each(key.begin(), key.end(), [](char& c) { c = ::tolower(c); });

                    if (key == "ncigrd") ncigrd = std::stoi(val);
                    if (key == "iroot") iroot = std::stoi(val);

                    keyword.push_back({key, std::stoi(val)});
                }
                if (eachline.back() != '+') parse_type = MNDO99_PARSE::COMMENT;
                break;
            }
            case MNDO99_PARSE::COMMENT: {
                mndo99_comment_sstr << eachline << std::endl;
                getline(ifs, eachline);
                mndo99_comment_sstr << eachline << std::endl;
                parse_type = MNDO99_PARSE::DATA;
                break;
            }
            case MNDO99_PARSE::DATA: {
                int ndata = 0;
                std::string data, databuf[16];
                std::istringstream input(eachline);
                while (input >> data) databuf[ndata++] = data;
                /** @comment
                    for cartesian format: ndata should be 7
                    for zmat format: ndata should be 10
                */
                if (std::stoi(databuf[0]) != 0) {  // normal data
                    count_atom++;
                    for (int i = 0; i < ndata; ++i) mndo99_data.push_back(databuf[i]);
                } else {
                    parse_type = MNDO99_PARSE::ADDITION;
                    mndo99_addition_sstr << eachline << std::endl;
                }
                break;
            }
            case MNDO99_PARSE::ADDITION: {
                mndo99_addition_sstr << eachline << std::endl;
                break;
            }
        }
    }
    ifs.close();
    mndo99_keyword  = mndo99_keyword_sstr.str();
    mndo99_comment  = mndo99_comment_sstr.str();
    mndo99_addition = mndo99_addition_sstr.str();

    CHECK_GT(count_atom, 0) << "Atom number should larger than 0";
    int data_size = mndo99_data.size();
    if (7 * count_atom == data_size) LOG(INFO) << "read from cartesian format";
    if (10 * count_atom == data_size) LOG(FATAL) << "cannot read from zmat format";
    return count_atom;
}


/**
 * @brief this function generates a set of keywords
 * @return std::string: describe keywords
 * @bug none
 */
std::string Model_Interf_MNDO::new_keyword(const MNDO99KW_map& newkeyword) {
    std::stringstream sstr;
    int i                  = 0;
    const int ntermperline = 8;

    for (auto iter = keyword.begin(); iter != keyword.end(); ++iter, ++i) {
        MNDO99KW& kw = *iter;
        if (i != 0 && i % ntermperline == 0) sstr << "+" << std::endl;
        if (newkeyword.find(kw.key) != newkeyword.end()) {
            sstr << kw.key << "=" << newkeyword.at(kw.key) << " ";
        } else {
            sstr << kw.key << "=" << kw.val << " ";
        }
    }
    sstr << std::endl;
    return sstr.str();
}


namespace MNDO99_TASK {
enum _enum { SP, FORCE, NAD, HESS, DEF };
const std::map<std::string, _enum> _dict = {{"sp", SP}, {"force", FORCE}, {"nad", NAD}, {"hess", HESS}, {"def", DEF}};
};  // namespace MNDO99_TASK

/**
 * @brief this function generates a input file
 * @return int
 * @bug none
 */
int Model_Interf_MNDO::new_task(num_real* R, const std::string& file, const std::string& task_flag, const int& rdim) {
    // parse task_flag
    int task_type = MNDO99_TASK::_dict.at(task_flag);

    std::string revised_keyword, revised_addition;
    switch (task_type) {
        case MNDO99_TASK::SP:
            revised_keyword = new_keyword({{"jop", -1}, {"icross", 0}});  ///< @bug: cannot find energies
            break;
        case MNDO99_TASK::FORCE:
            /** @comment
                for icross != 0. additional lines needed in mndo99_addition part
                be careful.
            */
            revised_keyword = new_keyword({{"jop", -1}, {"icross", 1}});
            break;
        case MNDO99_TASK::NAD:
            /** @comment
                for icross != 0. additional lines needed in mndo99_addition part
                be careful.
            */
            revised_keyword = new_keyword({{"jop", -1}, {"icross", 7}, {"mprint", 1}});
            break;
        case MNDO99_TASK::HESS:
            revised_keyword = new_keyword({{"jop", 2}, {"icross", 0}, {"kprint", 1}});
            break;
        default:
            revised_keyword = mndo99_keyword;
            break;
    }
    const int fixflag = 0;
    std::ofstream ofs(file);
    ofs << revised_keyword;
    ofs << mndo99_comment;
    for (int i = 0, idx = 0; i < natom; ++i) {
        ofs << FMT(4) << atoms[i]                       // Atom index
            << FMT(8) << R[idx++] << FMT(4) << fixflag  // Xk, flag_Xk
            << FMT(8) << R[idx++] << FMT(4) << fixflag  // Yk, flag_Yk
            << FMT(8) << R[idx++] << FMT(4) << fixflag  // Zk, flag_Zk
            << std::endl;
    }
    ofs << mndo99_addition;
    ofs.close();
    return 0;
}

/**
 * @brief this part track the sign of NAC between current step with the last step,
 *  in order to make sure nac changes continuously with time.
 * @return int
 * @bug: attention when nac cross zero line
 */
int Model_Interf_MNDO::track_nac_sign() {
    if (diff_nac) {  ///< skip starting time
        for (int i = 0; i < F; ++i) {
            for (int j = 0; j < F; ++j) {  // check if NAC(:,i,j) should flip its sign
                if (i == j) continue;

                const double norm_eps = 10e-14;
                double norm_old       = 0.0f;
                double norm_new       = 0.0f;
                double cosangle       = 0.0f;
                int IJ                = i * F + j;
                for (int k = 0, idx = IJ; k < N; ++k, idx += FF) {
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
                    for (int k = 0; k < N; ++k) {
                        nac[k * FF + i * F + j] = copysign(nac[k * FF + i * F + j], nac_prev[k * FF + i * F + j]);
                    }
                } else if (cosangle < 0) {  // in this case we flip the sign of NAC(:,i,j)
                    for (int k = 0; k < N; ++k) { nac[k * FF + i * F + j] *= -1; }
                }
            }
        }
    }
    diff_nac = true;                                     ///< diff_nac is set to true after starting
    for (int i = 0; i < NFF; ++i) nac_prev[i] = nac[i];  // save a copy

    return 0;
}


/**
 * @brief parse energy/gradients/nac from output
 * @param
 *    log:  log file to parse
 * @return status
 * @bug none
 */
int Model_Interf_MNDO::parse_standard(const std::string& log) {
    int stat = -1;
    int istate, jstate;

    std::ifstream ifs(log);
    std::string stmp, eachline;
    while (getline(ifs, eachline)) {
        /**
         * @brief find energy surfaces (for icross=0, this section is missed)
         */
        if (eachline.find("SUMMARY OF MULTIPLE CI CALCULATIONS") != eachline.npos) {
            for (int i = 0; i < 4; ++i) getline(ifs, eachline);
            for (int i = 0; i < F; ++i) {  ///< E in [kcalpmol]
                ifs >> stmp >> stmp >> ener[i] >> stmp >> stmp;
            }
            stat = 0;
        }
        /**
         * @brief find gradients
         */
        else if (eachline.find("CI CALCULATION FOR STATE:") != eachline.npos) {
            std::istringstream sstr(eachline);
            sstr >> stmp >> stmp >> stmp >> stmp >> istate;
            istate--;
            if (istate < F) {  // skip additional state
                while (getline(ifs, eachline)) {
                    // if (eachline.find("Mult. 1") != eachline.npos) {  ///< ener2 saved in [ev]
                    //     std::istringstream sstr(eachline);
                    //     sstr >> stmp >> jstate >> stmp >> stmp >> stmp >> stmp >> stmp >> stmp;
                    //     jstate--;
                    //     sstr >> ener[jstate];
                    // }
                    if (eachline.find("GRADIENTS (KCAL/(MOL*ANGSTROM))") != eachline.npos) {
                        for (int i = 0; i < 8; ++i) ifs >> stmp;
                        for (int i = 0, idx = 0; i < natom; ++i) {       ///< grad in kcalpmol/angstrom
                            ifs >> stmp >> stmp >> stmp >> stmp >> stmp  // skips
                                >> grad[(idx++) * F + istate]            // grad(Xk,i,i)
                                >> grad[(idx++) * F + istate]            // grad(Yk,i,i)
                                >> grad[(idx++) * F + istate];           // grad(Zk,i,i)
                        }
                        break;
                    }
                }
            }
            stat = 1;
        }
        /**
         * @brief find non-adiabatic coupling(NAC) terms
         *  note we use the complete formalism (let `MPRINT=1`) for NAC in MNDO99
         */
        else if (eachline.find("CI CALCULATION FOR INTERSTATE "
                               "COUPLING OF STATES:") != eachline.npos) {
            std::istringstream sstr(eachline);
            sstr >> stmp >> stmp >> stmp >> stmp >> stmp >> stmp >> stmp >> istate >> jstate;
            istate--, jstate--;
            if (istate < F && jstate < F) {  // skip additional NAC
                int IJ = istate * F + jstate;
                int JI = jstate * F + istate;
                while (getline(ifs, eachline)) {
                    if (eachline.find("COMPLETE EXPRESSION.") != eachline.npos) {  // let `MPRINT=1`
                        for (int i = 0, idx = 0; i < natom; ++i) {                 ///< found nacv in 1/angstrom
                            ifs >> stmp                                            // skips
                                >> nac[(idx++) * FF + IJ]                          // nac(Xk, i, j)
                                >> nac[(idx++) * FF + IJ]                          // nac(Yk, i, j)
                                >> nac[(idx++) * FF + IJ];                         // nac(Zk, i, j)
                        }
                        for (int idx = 0; idx < N; ++idx) {  // copy to the other half-side nacv
                            nac[idx * FF + JI] = -nac[idx * FF + IJ];
                        }
                        break;
                    }
                }
            }
            stat = 2;
        }
    }
    ifs.close();
    return stat;
}

/**
 * @brief this function parse frequency calculation (JOP=2)
 * @param
 *    log:  log file to parse
 * @return int: if it is on equilibrium postion
 *      0: yes
 *     -1: no
 * @bug none
 */
int Model_Interf_MNDO::parse_hessian(const std::string& log) {
    std::ifstream ifs(log);

    // initialization of arrays
    for (int i = 0; i < N; ++i) mod_W[i] = 0.0f;
    for (int i = 0; i < NN; ++i) mod_Hess[i] = 0.0f, mod_Tmat[i] = 0.0f;
    std::string stmp, eachline;
    int v1 = N / 10, v2 = N % 10, eqstat = 0, itmp;
    double dtmp = 0.0f;
    while (getline(ifs, eachline)) {
        if (eachline.find("GRADIENT NORM =") != eachline.npos) {
            double norm;
            std::string stmp;
            std::istringstream sstr(eachline);
            sstr >> stmp >> stmp >> stmp >> norm;
            eqstat = (norm < 10.0f) ? 0 : -1;
        }
        if (eachline.find("FORCE CONSTANT MATRIX AFTER SYMMETRIZATION.") != eachline.npos) {
            for (int i = 0; i < v1; ++i) {
                getline(ifs, eachline);                    // skip a blankline
                for (int k = 0; k < 10; ++k) ifs >> itmp;  // skip header
                for (int j = 0; j < N; ++j) {
                    ifs >> itmp;                                                       // skip index
                    for (int k = 0; k < 10; ++k) ifs >> mod_Hess[N * j + i * 10 + k];  // read
                }
            }
            getline(ifs, eachline);                    // skip a blankline
            for (int k = 0; k < v2; ++k) ifs >> itmp;  // skip header
            for (int j = 0; j < N; ++j) {
                ifs >> itmp;                                                        // skip index
                for (int k = 0; k < v2; ++k) ifs >> mod_Hess[N * j + v1 * 10 + k];  // read
            }
        }
        if (eachline.find("EIGENVECTORS OF THE MASS-WEIGHTED") != eachline.npos) {
            for (int i = 0; i < v1; ++i) {
                getline(ifs, eachline);                                         // skip a blankline
                for (int k = 0; k < 10; ++k) ifs >> itmp;                       // skip header
                for (int k = 0; k < 10; ++k) {                                  // read frequency
                    if (ifs >> dtmp) mod_W[i * 10 + k] = dtmp / phys::au_2_wn;  // [convert wn to au]
                }
                for (int j = 0; j < N; ++j) {  // read normal-mode
                    ifs >> itmp;               // skip index
                    for (int k = 0; k < 10; ++k) {
                        if (ifs >> dtmp) mod_Tmat[N * j + i * 10 + k] = dtmp;
                    }
                }
            }
            getline(ifs, eachline);                                          // skip a blankline
            for (int k = 0; k < v2; ++k) ifs >> itmp;                        // skip header
            for (int k = 0; k < v2; ++k) {                                   // read frequency
                if (ifs >> dtmp) mod_W[v1 * 10 + k] = dtmp / phys::au_2_wn;  // [convert wn to au]
            }
            for (int j = 0; j < N; ++j) {  // read normal-mode
                ifs >> itmp;               // skip index
                for (int k = 0; k < v2; ++k) {
                    if (ifs >> dtmp) mod_Tmat[N * j + v1 * 10 + k] = dtmp;
                }
            }
        }
    }
    ifs.close();

    std::ofstream ofs(".hess");
    for (int i = 0, idx = 0; i < N; ++i) {  // hessian (in au)
        for (int j = 0; j < N; ++j, ++idx) ofs << FMT(8) << mod_Hess[idx];
        ofs << std::endl;
    }
    for (int i = 0, idx = 0; i < N; ++i) {  // frequancey (in au)
        ofs << FMT(8) << mod_W[i];
    }
    ofs << std::endl;
    for (int i = 0, idx = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j, ++idx) ofs << FMT(8) << mod_Tmat[idx];
        ofs << std::endl;
    }
    ofs.close();
    return eqstat;
}

/**
 * @brief this function generates normalmode trajectories
 * @param
 *      R: equilibrium configuration
 * @return int
 * @bug none
 */
int Model_Interf_MNDO::calc_normalmode(num_real* R, const int& rdim) {
    // update mndo99hess input/output files
    if (!utils::isFileExists(".mndo99hess.out")) {
        new_task(R, ".mndo99hess.in", "hess", rdim);
        int stat = system("mndo99 < .mndo99hess.in > .mndo99hess.out");
        if (stat != 0) return stat;
    }
    int stat = parse_hessian(".mndo99hess.out");
    CHECK_EQ(stat, 0) << "parse_hessian fails";

    /**
     * @brief generate normal-mode trajectories
     * @details
     *      r = Xeq + M^(-1/2) * T * Q
     *      p = M^(1/2) * T * P
     *      here, r, p are cartesian postion & momentum, while Q, P are normal-mode postion & momentum
     */
    int ntime = 50;                   // report 50 points in [0,2*pi] period
    for (int iw = 0; iw < N; ++iw) {  // iw-the mode
        std::string save = "normalmode";
        std::ofstream ofs(save + std::to_string(iw) + ".xyz");
        for (int itime = 0; itime < ntime; ++itime) {  // itime
            for (int i = 0; i < N; ++i) {
                nr_samp[i] = mod_R0[i] + mod_Tmat[N * i + iw] / std::sqrt(mod_M[i]) *
                                             std::sin(itime * phys::math::twopi / ntime) / std::sqrt(2.0f * mod_W[iw]);
            }
            ofs << natom << std::endl;
            ofs << "time: " << (itime * phys::math::twopi) / (mod_W[iw] * ntime) << std::endl;
            for (int i = 0, idx = 0; i < natom; ++i) {
                ofs << FMT(8) << ELEMENTS_LABEL[atoms[i]]   //
                    << FMT(8) << nr_samp[idx++] * iou.leng  //
                    << FMT(8) << nr_samp[idx++] * iou.leng  //
                    << FMT(8) << nr_samp[idx++] * iou.leng << std::endl;
            }
        }
        ofs.close();
    }
    return 0;
}

/**
 * @brief this function generates initialization configuration
 * @param
 *      R: equilibrium configuration
 * @return int
 * @bug none
 */
int Model_Interf_MNDO::calc_samp() {
    int stat;

    if (!utils::isFileExists(".mndo99hess.out")) {
        LOG(WARNING);
        new_task(mod_R0, ".mndo99hess.in", "hess", N);
        stat = system("mndo99 < .mndo99hess.in > .mndo99hess.out");
        if (stat != 0) return stat;
    }
    stat = parse_hessian(".mndo99hess.out");
    CHECK_EQ(stat, 0) << "parse_hessian fails";

    // ARRAY_SHOW(mod_Hess, N, N);
    // LOG(FATAL);

    /**
     * mass-weighted Hessian to obtain normalmode
     */
    for (int i = 0, idx = 0; i < N; ++i)
        for (int j = 0; j < N; ++j, ++idx) mod_Hess[idx] /= sqrt(mod_M[i] * mod_M[j]);


    // ARRAY_SHOW(mod_W, N / 3, 3);

    double Qeff = 6.0;
    for (int j = 0; j < N; ++j) {
        if (j < 6 || mod_W[j] < 1.0e-3) {  // cutoff low frequency
            mod_sigmaR[j] = 0.0f;
            mod_sigmaP[j] = 0.0f;
        } else {                             // NOTE: it's for normal-mode!
            double Qoverbeta = 1.0f / beta;  // 0.5f * mod_W[j] / std::tanh(0.5f * beta * mod_W[j]);
            Qeff += Qoverbeta * beta;
            mod_sigmaR[j] = std::sqrt(Qoverbeta / (mod_W[j] * mod_W[j]));
            mod_sigmaP[j] = std::sqrt(Qoverbeta);
        }
    }
    Qeff /= N;
    Qeff = 0.5f * beta * mod_W[N - 1] / std::tanh(0.5f * beta * mod_W[N - 1]);

    LOG(WARNING) << "Qeff = " << Qeff;

    // sampling normal-mode
    double* tmpwork = new double[N];
    num_real* E     = new num_real[F];
    num_real* dE    = new num_real[NFF];
    num_real* ddE;

    ForceField_epes(E, dE, ddE, mod_R0, 1, N, F, 0, 0);
    double ref_ener = E[0];
    LOG(WARNING) << ref_ener;

    for (int isamp = 0; isamp < 500; ++isamp) {
        rand_gaussian(nr_samp, N);
        rand_gaussian(np_samp, N);
        for (int i = 0; i < 6; ++i) nr_samp[i] = 0.0f, np_samp[i] = 0.0f;
        for (int i = 0; i < N; ++i) {
            nr_samp[i] *= mod_sigmaR[i];
            np_samp[i] *= mod_sigmaP[i];
        }
        // transfrom normal-mode to cartesian coordinates
        ARRAY_MATMUL(tmpwork, mod_Tmat, nr_samp, N, N, 1);
        for (int i = 0; i < N; ++i) nr_samp[i] = tmpwork[i];
        ARRAY_MATMUL(tmpwork, mod_Tmat, np_samp, N, N, 1);
        for (int i = 0; i < N; ++i) np_samp[i] = tmpwork[i];
        for (int i = 0; i < N; ++i) {
            nr_samp[i] = nr_samp[i] / std::sqrt(mod_M[i]) + mod_R0[i];
            np_samp[i] = np_samp[i] * std::sqrt(mod_M[i]);
        }

        // fluctuation of Kinetic energy
        double Ekin = 0.0;
        for (int i = 0; i < N; ++i) Ekin += 0.5f * np_samp[i] * np_samp[i] / mod_M[i];

        // fluctuation of potential energy
        ForceField_epes(E, dE, ddE, nr_samp, 1, N, F, 0, isamp);
        double Epot = (E[0] - ref_ener);  // convert to au
        LOG(WARNING) << Epot;
        LOG(WARNING) << "::: " << beta * Epot / Qeff << " " << beta * Ekin / Qeff;

        if (beta * Epot > 10 || beta * Ekin > 10) {
            LOG(WARNING) << "failed sample once";
            std::ofstream ofs(utils::concat("sampfail", isamp));
            ofs << FMT(8) << natom << std::endl;  // number of atoms
            ofs << FMT(8) << isamp << FMT(8) << beta << FMT(8) << Epot << FMT(8) << Ekin << std::endl;
            for (int i = 0, idx1 = 0, idx2 = 0; i < natom; ++i) {
                // output configuration ([angstrom])
                ofs << FMT(8) << ELEMENTS_LABEL[atoms[i]];
                for (int k = 0; k < 3; ++k, ++idx1) ofs << FMT(8) << nr_samp[idx1] * phys::au_2_ang;
                // output velocity ([angstrom/ps])
                for (int k = 0; k < 3; ++k, ++idx2) ofs << FMT(8) << np_samp[idx2] / mod_M[idx2] * phys::au_2_angoverps;
                ofs << std::endl;
            }
            ofs.close();
            isamp--;
            continue;  // for resampling
        };

        LOG(WARNING) << "test: " << beta * Epot << " " << beta * Ekin;

        std::ofstream ofs(utils::concat("samp", isamp));
        ofs << FMT(8) << natom << std::endl;  // number of atoms
        ofs << FMT(8) << isamp << std::endl;
        for (int i = 0, idx1 = 0, idx2 = 0; i < natom; ++i) {
            // output configuration ([angstrom])
            ofs << FMT(8) << ELEMENTS_LABEL[atoms[i]];
            for (int k = 0; k < 3; ++k, ++idx1) ofs << FMT(8) << nr_samp[idx1] * phys::au_2_ang;
            // output velocity ([angstrom/ps])
            for (int k = 0; k < 3; ++k, ++idx2) ofs << FMT(8) << np_samp[idx2] / mod_M[idx2] * phys::au_2_angoverps;
            ofs << std::endl;
        }
        ofs.close();
    }

    delete[] tmpwork, delete[] E, delete[] dE;
    return 0;
}

int Model_Interf_MNDO::calc_scan() {
    int istep = 0, readn;
    num_real tmp;
    std::string eachline;
    num_real* E  = new num_real[F];
    num_real* dE = new num_real[NFF];
    num_real* ddE;

    savefile_traj = utils::concat("traj-", 0, ".xyz");
    savefile_ener = utils::concat("ener-", 0, ".dat");
    savefile_grad = utils::concat("grad-", 0, ".dat");
    savefile_nac  = utils::concat("nac-", 0, ".dat");
    utils::clearFile(savefile_traj);
    utils::clearFile(savefile_ener);
    utils::clearFile(savefile_grad);
    utils::clearFile(savefile_nac);

    std::ifstream ifs("scan.xyz");
    if (!ifs.is_open()) LOG(FATAL) << "scan.xyz cannot open";

    while (getline(ifs, eachline)) {
        std::istringstream isstr(eachline);
        isstr >> readn;
        CHECK_EQ(natom, readn);
        getline(ifs, eachline);  // skip comments
        for (int iatom = 0, idx1 = 0; iatom < natom; ++iatom) {
            ifs >> eachline;                       // skip atomic flag
            for (int i = 0; i < 3; ++i, ++idx1) {  // read in [angstrom]
                if (ifs >> tmp) nr_samp[idx1] = tmp / phys::au_2_ang;
            }
        }
        getline(ifs, eachline);  // a line!!!
        // ARRAY_SHOW(mod_R0, N / 3, 3);
        ForceField_epes(E, dE, ddE, nr_samp, 1, N, F, 0, istep);
        istep++;
    }
    ifs.close();
    delete[] E, delete[] dE;
    return 0;
}


#endif  // MODEL_INTERF_MNDO_H