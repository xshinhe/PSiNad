#include "Model_Interf_MNDO.h"

#include "../core/Element.h"
#include "../core/linalg.h"
#include "../kernels/Kernel_Declare.h"
#include "../kernels/Kernel_Random.h"
#include "../kernels/Kernel_Representation.h"
// #include "../mpi_utils.h"

#define ARRAY_SHOW(_A, _n1, _n2)                                                            \
    ({                                                                                      \
        std::cout << "Show Array <" << #_A << ">\n";                                        \
        int _idxA = 0;                                                                      \
        for (int _i = 0; _i < (_n1); ++_i) {                                                \
            for (int _j = 0; _j < (_n2); ++_j) std::cout << FMT(4) << (_A)[_idxA++] << ","; \
            std::cout << std::endl;                                                         \
        }                                                                                   \
    })

inline int removeFile(const std::string& filename) { return remove(filename.c_str()); }

inline void clearFile(const std::string& filename) { std::ofstream clear(filename, std::ios::trunc); }

inline void closeOFS(std::ofstream& ofs) {
    if (ofs.is_open()) ofs.close();
}

inline bool isFileExists(const std::string& name) { return std::ifstream{name.c_str()}.good(); }

namespace PROJECT_NS {

void Model_Interf_MNDO::read_param_impl(Param* PM) {
    Kernel_Representation::onthefly = true;

    // parse mndo input
    exec_file = PM->get<std::string>("exec_file", LOC(), "mndo");
    directory = PM->get<std::string>("directory", LOC());

    std::string mndoinp = PM->get<std::string>("mndoinp", LOC(), "null");
    natom               = parse_mndo(mndoinp);

    assert(Dimension::N == 3 * natom);
    assert(Dimension::F <= ncigrd);
    assert(Dimension::F <= nciref);
    assert(directory != "");
}

void Model_Interf_MNDO::init_data_impl(DataSet* DS) {
    x = DS->def<double>("integrator.x", Dimension::N);
    p = DS->def<double>("integrator.p", Dimension::N);

    x0      = DS->def<double>("model.x0", Dimension::N);
    p0      = DS->def<double>("model.p0", Dimension::N);
    w       = DS->def<double>("model.w", Dimension::N);
    x_sigma = DS->def<double>("model.x_sigma", Dimension::N);
    p_sigma = DS->def<double>("model.p_sigma", Dimension::N);

    // model field
    atoms = DS->def<int>("model.atoms", Dimension::N);
    mass  = DS->def<double>("model.mass", Dimension::N);
    vpes  = DS->def<double>("model.vpes");
    grad  = DS->def<double>("model.grad", Dimension::N);
    hess  = DS->def<double>("model.hess", Dimension::NN);
    Tmod  = DS->def<double>("model.Tmod", Dimension::NN);
    f_r   = DS->def<double>("model.f_r", nciref);
    f_p   = DS->def<double>("model.f_p", nciref);
    f_rp  = DS->def<double>("model.f_rp", nciref);

    V  = DS->def<double>("model.V", Dimension::FF);
    dV = DS->def<double>("model.dV", Dimension::NFF);
    // ddV  = DS->def<double>("model.ddV", Dimension::NNFF);
    E  = DS->def<double>("model.rep.E", Dimension::F);
    T  = DS->def<double>("model.rep.T", Dimension::FF);  // diabatic-to-adiabatic matrix. set to I as default
    dE = DS->def<double>("model.rep.dE", Dimension::NFF);
    // ddE  = DS->def<double>("model.rep.ddE", Dimension::NNFF);
    nac      = DS->def<double>("model.rep.nac", Dimension::NFF);
    nac_prev = DS->def<double>("model.rep.nac_prev", Dimension::NFF);

    for (int i = 0, ik = 0; i < Dimension::F; ++i) {
        for (int k = 0; k < Dimension::F; ++k, ++ik) T[ik] = (i == k) ? 1.0e0 : 0.0e0;
    }

    // read z index
    for (int i = 0, idx = 0, idxR = 0; i < natom; ++i) {
        atoms[i] = stoi(mndo_data[idx++]);
        for (int j = 0; j < 3; ++j, ++idxR) {
            mass[idxR] = ELEMENTS_MASS[atoms[i]] / phys::au_2_amu;  // convert amu to au
            x0[idxR]   = stod(mndo_data[idx++]) / phys::au_2_ang;   // convert angstrom into au
            idx++;                                                  // skip fixflag
        }
    }

    // read temperature
    double temperature = _Param->get<double>("temperature", LOC(), phys::temperature_d, 1.0f);
    beta               = 1.0f / (phys::au::k * temperature);  // don't ignore k_Boltzman

    // read task
    init_nuclinp = _Param->get<std::string>("init_nuclinp", LOC(), "#hess");

    if (init_nuclinp == "#normalmode") {
        calc_normalmode();
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

    // read hessian from hessian calculation (jop=2 kprint=1)
    if (init_nuclinp == "#hess") {
        std::string hess_log = _Param->get<std::string>("hess_log", LOC(), "hess.out");
        if (!isFileExists(hess_log))
            throw std::runtime_error(utils::concat("hess_log file [", hess_log, "] is missed!"));

        parse_hessian(hess_log);

        // from hessian & temperature to prepare initial sampling
        for (int j = 0; j < Dimension::N; ++j) {
            if (j < 6 || w[j] < 1.0e-10) {  // cutoff low frequency
                x_sigma[j] = 0.0f;
                p_sigma[j] = 0.0f;
            } else {  // NOTE: it's for normal-mode!
                double Qoverbeta = 0.5f * w[j] / std::tanh(0.5f * beta * w[j]);
                x_sigma[j]       = std::sqrt(Qoverbeta / (w[j] * w[j]));
                p_sigma[j]       = std::sqrt(Qoverbeta);
            }
        }
    }
    if (init_nuclinp == "#hess2") {
        std::string hess_mol = _Param->get<std::string>("hess_mol", LOC(), "hess_molden.dat");
        if (!isFileExists(hess_mol))
            throw std::runtime_error(utils::concat("hess_mol file [", hess_mol, "] is missed!"));

        parse_hessian2(hess_mol);
    }
}

/**
 * @brief ForceField_init for mndo
 * @param
 *     nr:  nuclear configuration
 *     np:  nuclear momentum
 *     m:   nuclear mass
 *   erho:  electronic density
 *   eeac:  electronic amplititude
 *   eocc:  electronic occupation
 *   fdim:  electronic freedoms
 *  itraj:  index of trajectory
 * @return int
 * @bug none
 */
void Model_Interf_MNDO::init_calc_impl(int stat) {
    if (init_nuclinp == "#hess" || init_nuclinp == "#hess2") {  // assuming .hess is given
        // sampling normal-mode
        Kernel_Random::rand_gaussian(x, Dimension::N);
        Kernel_Random::rand_gaussian(p, Dimension::N);
        for (int i = 0; i < 6; ++i) x[i] = 0.0f, p[i] = 0.0f;
        for (int i = 0; i < Dimension::N; ++i) {
            x[i] *= x_sigma[i];
            p[i] *= p_sigma[i];
        }
        // transfrom normal-mode to cartesian coordinates
        ARRAY_MATMUL(x, Tmod, x, Dimension::N, Dimension::N, 1);  // .eval()!
        ARRAY_MATMUL(p, Tmod, p, Dimension::N, Dimension::N, 1);  // .eval()!
        for (int i = 0; i < Dimension::N; ++i) {
            x[i] = x[i] / std::sqrt(mass[i]) + x0[i];
            p[i] = p[i] * std::sqrt(mass[i]) + p0[i];
        }
    } else if (init_nuclinp == "#fix") {  // for initial md
        for (int i = 0; i < Dimension::N; ++i) x[i] = x0[i], p[i] = 0.0f;
    } else {  // init_nuclinp as dataset from which we read x and p
        std::string stmp, eachline;
        std::ifstream ifs(utils::concat(init_nuclinp, stat, ".ds"));
        while (getline(ifs, eachline)) {
            if (eachline.find("init.x") != eachline.npos) {
                getline(ifs, eachline);
                for (int i = 0; i < Dimension::N; ++i) ifs >> x[i];
            }
            if (eachline.find("init.p") != eachline.npos) {
                getline(ifs, eachline);
                for (int i = 0; i < Dimension::N; ++i) ifs >> p[i];
            }
        }

        // ARRAY_SHOW(x, 1, Dimension::N);
        // ARRAY_SHOW(p, 1, Dimension::N);
    }

    _DataSet->def("init.x", x, Dimension::N);
    _DataSet->def("init.p", p, Dimension::N);

    std::string hdlr_str = _Param->get<std::string>("handler", LOC());

    refer        = false;
    task_control = (hdlr_str == "sampling") ? "samp" : "nad";
    removeFile(utils::concat(directory, "/imomap.dat"));
    exec_kernel_impl(stat);
    refer        = true;
    task_control = "nad";
}

int Model_Interf_MNDO::exec_kernel_impl(int stat_in) {
    // convert atomic unit
    for (int i = 0; i < Dimension::N; ++i) x[i] *= phys::au_2_ang;  ///< convert Bohr to Angstrom
    ARRAY_CLEAR(E, Dimension::F);
    ARRAY_CLEAR(dE, Dimension::NFF);

    // prepare a mndo input file
    std::string inpfile = utils::concat(".mndoinp.", stat_in);
    std::string outfile = utils::concat(".mndoout.", stat_in);

    new_task(utils::concat(directory, "/", inpfile), task_control);

    std::string cmd_exe = utils::concat("cd ", directory, " && ", exec_file, " < ", inpfile, " > ", outfile);
    int stat            = system(cmd_exe.c_str());

    parse_standard(utils::concat(directory, "/", outfile), stat_in);  // parse in MNDO's units

    // ARRAY_SHOW(x, 1, Dimension::N);
    // ARRAY_SHOW(E, 1, Dimension::F);
    // ARRAY_SHOW(dE, 1, Dimension::F);
    // ARRAY_SHOW(nac, Dimension::N, Dimension::FF);

    track_nac_sign();  // @note track_nac_sign is important
    for (int i = 0, idx = 0; i < Dimension::N; ++i) {
        for (int j = 0; j < Dimension::F; ++j) {
            for (int k = 0; k < Dimension::F; ++k, ++idx) {
                if (j == k) continue;
                dE[idx] = nac[idx] * (E[k] - E[j]);
            }
        }
    }

    // convert units
    for (int i = 0; i < Dimension::N; ++i) x[i] /= phys::au_2_ang;        ///< convert Angstrom to Bohr
    for (int i = 0; i < Dimension::F; ++i) E[i] /= phys::au_2_kcal_1mea;  ///< convert kcalpmol to Hartree
    for (int i = 0; i < Dimension::NFF; ++i)
        dE[i] /= (phys::au_2_kcal_1mea / phys::au_2_ang);  ///< convert to Hartree/Bohr

    // @TODO
    // if (!ARRAY_ISFINITE(R, Dimension::N) || !ARRAY_ISFINITE(E, F) || !ARRAY_ISFINITE(dE, NFF)) {
    //     std::cout << "Problem for calling forcefield";
    //     stat = -1;
    // }
    return stat;
}

/**
 * @brief this function parse mndo input (only support cartesian format)
 * @return int: total atom number parsed
 * @todo:   support zmat format parser
 */
int Model_Interf_MNDO::parse_mndo(const std::string& mndoinp) {
    enum Stage { KEYWORD, COMMENT, DATA, ADDITION };
    Stage istage = KEYWORD;

    std::stringstream mndo_keyword_sstr;
    std::stringstream mndo_comment_sstr;
    std::stringstream mndo_addition_sstr;

    int count_atom = 0;
    std::ifstream ifs(mndoinp);
    std::string eachline;
    while (getline(ifs, eachline)) {
        eachline.erase(eachline.find_last_not_of(" ") + 1);  // remove last blank
        switch (istage) {
            case KEYWORD: {
                mndo_keyword_sstr << eachline << std::endl;
                std::string data, key, val;
                std::istringstream input(eachline);
                while (input >> data) {
                    if (data == "+") continue;
                    std::istringstream pair(data);  // "key=val" pattern
                    getline(pair, key, '=');
                    getline(pair, val, '=');
                    std::for_each(key.begin(), key.end(), [](char& c) { c = ::tolower(c); });

                    if (key == "ncigrd") ncigrd = std::stoi(val);
                    if (key == "nciref") nciref = std::stoi(val);
                    if (key == "iroot") iroot = std::stoi(val);

                    keyword.push_back({key, std::stoi(val)});
                }
                if (eachline.back() != '+') istage = COMMENT;
                break;
            }
            case COMMENT: {
                mndo_comment_sstr << eachline << std::endl;
                getline(ifs, eachline);
                mndo_comment_sstr << eachline << std::endl;
                istage = DATA;
                break;
            }
            case DATA: {
                int ndata = 0;
                std::string data, databuf[16];
                std::istringstream input(eachline);
                while (input >> data) databuf[ndata++] = data;
                if (std::stoi(databuf[0]) != 0) {
                    count_atom++;
                    for (int i = 0; i < ndata; ++i) mndo_data.push_back(databuf[i]);
                } else {
                    istage = ADDITION;
                    mndo_addition_sstr << eachline << std::endl;
                }
                break;
            }
            case ADDITION: {
                mndo_addition_sstr << eachline << std::endl;
                break;
            }
        }
    }
    ifs.close();
    mndo_keyword  = mndo_keyword_sstr.str();
    mndo_comment  = mndo_comment_sstr.str();
    mndo_addition = mndo_addition_sstr.str();
    assert(count_atom > 0);
    int data_size = mndo_data.size();
    if (7 * count_atom == data_size) std::cout << "read from cartesian format";
    if (10 * count_atom == data_size) throw std::runtime_error("cannot read from zmat format");
    return count_atom;
}


/**
 * generate a new set of keywords
 */
std::string Model_Interf_MNDO::new_keyword(const MNDOKW_map& newkeyword) {
    std::stringstream sstr;
    int i                  = 0;
    const int ntermperline = 8;
    for (auto iter = keyword.begin(); iter != keyword.end(); ++iter, ++i) {
        MNDOKW& kw = *iter;
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


/**
 * generate input file for a new task based on the template
 */
int Model_Interf_MNDO::new_task(const std::string& file, const std::string& task_flag) {
    std::string revised_keyword, revised_addition;
    if (task_flag == "sp") {                                          // single point calculation
        revised_keyword = new_keyword({{"jop", -1}, {"icross", 0}});  ///< @bug: cannot find energies
    } else if (task_flag == "samp") {                                 // single point calculation
        revised_keyword =
            new_keyword({{"jop", -1}, {"icross", 7}, {"mprint", 1}, {"iuvcd", 2}, {"imomap", 0}, {"mapthr", 90}});
    } else if (task_flag == "force") {  // force calculation
        revised_keyword = new_keyword({{"jop", -1}, {"icross", 1}});
        // revised_addition = ...; // additional lines
    } else if (task_flag == "nad") {  // non-adiabatic coupling calculation
        revised_keyword = new_keyword({{"jop", -1}, {"icross", 7}, {"mprint", 1}, {"imomap", 3}, {"mapthr", 95}});
        // revised_addition = ...; // additional lines
    } else if (task_flag == "hess") {  // hessian calculation
        revised_keyword = new_keyword({{"jop", 2}, {"icross", 0}, {"kprint", 1}});
    } else {  // default
        revised_keyword = mndo_keyword;
    }

    revised_addition  = mndo_addition;
    const int fixflag = 0;
    std::ofstream ofs(file);
    ofs << revised_keyword;
    ofs << mndo_comment;
    for (int i = 0, idx = 0; i < natom; ++i) {
        ofs << FMT(4) << atoms[i]                       // Atom index
            << FMT(8) << x[idx++] << FMT(4) << fixflag  // Xk, flag_Xk
            << FMT(8) << x[idx++] << FMT(4) << fixflag  // Yk, flag_Yk
            << FMT(8) << x[idx++] << FMT(4) << fixflag  // Zk, flag_Zk
            << std::endl;
    }
    ofs << revised_addition;
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
    if (refer) {  ///< skip starting time
        for (int i = 0; i < Dimension::F; ++i) {
            for (int j = 0; j < Dimension::F; ++j) {  // check if NAC(:,i,j) should flip its sign
                if (i == j) continue;

                const double norm_eps = 10e-14;
                double norm_old       = 0.0f;
                double norm_new       = 0.0f;
                double cosangle       = 0.0f;
                int IJ                = i * Dimension::F + j;
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
    }
    for (int i = 0; i < Dimension::NFF; ++i) nac_prev[i] = nac[i];  // save a copy
    return 0;
}


/**
 * @brief parse energy/gradients/nac from output
 * @param
 *    log:  log file to parse
 * @return status
 * @bug none
 */
int Model_Interf_MNDO::parse_standard(const std::string& log, int stat_in) {
    int stat         = -1;
    int istate_force = 0;
    int istate, jstate;

    std::ifstream ifs(log);
    std::string stmp, eachline;
    while (getline(ifs, eachline, '\n')) {
        /**
         * @brief find energy surfaces (for icross=0, this section is missed)
         */
        if (eachline.find("Properties of transitions   1 -> #") != eachline.npos) {
            for (int i = 0; i < 2; ++i) getline(ifs, eachline);  // blankline + headline
            for (int i = 1; i < nciref; ++i) {                   ///< E in [kcalpmol]
                ifs >> stmp >> stmp >> stmp >> stmp >> stmp >> stmp >> f_r[i] >> f_p[i] >> f_rp[i] >> stmp;
            }
            stat = 0;
        }
        /**
         * @brief find energy surfaces (for icross=0, this section is missed)
         */
        else if (eachline.find("SUMMARY OF MULTIPLE CI CALCULATIONS") != eachline.npos) {
            for (int i = 0; i < 4; ++i) getline(ifs, eachline);
            for (int i = 0; i < Dimension::F; ++i) {  ///< E in [kcalpmol]
                ifs >> stmp >> stmp >> E[i] >> stmp >> stmp;
            }
            stat = 0;
        }
        /**
         * @brief find gradients
         */
        else if (eachline.find("CI CALCULATION FOR STATE:") != eachline.npos) {
            std::istringstream sstr(eachline);
            sstr >> stmp >> stmp >> stmp >> stmp >> istate_force;
            istate_force--;
        }
        //
        else if (eachline.find("GRADIENTS (KCAL/(MOL*ANGSTROM))") != eachline.npos) {
            for (int i = 0; i < 8; ++i) ifs >> stmp;
            for (int i = 0, idx = 0; i < natom; ++i) {                                 ///< grad in kcalpmol/angstrom
                ifs >> stmp >> stmp >> stmp >> stmp >> stmp                            // skips
                    >> dE[(idx++) * Dimension::FF + istate_force * Dimension::Fadd1]   // grad(Xk,i,i)
                    >> dE[(idx++) * Dimension::FF + istate_force * Dimension::Fadd1]   // grad(Yk,i,i)
                    >> dE[(idx++) * Dimension::FF + istate_force * Dimension::Fadd1];  // grad(Zk,i,i)
            }
        }
        /**
         * @brief find non-adiabatic coupling(NAC) terms
         *  note we use the complete formalism (let `MPRINT=1`) in MNDO
         */
        else if (eachline.find("CI CALCULATION FOR INTERSTATE "
                               "COUPLING OF STATES:") != eachline.npos) {
            std::istringstream sstr(eachline);
            sstr >> stmp >> stmp >> stmp >> stmp >> stmp >> stmp >> stmp >> istate >> jstate;
            istate--, jstate--;
            if (istate < Dimension::F && jstate < Dimension::F && istate != jstate) {  // skip additional NAC
                int IJ = istate * Dimension::F + jstate;
                int JI = jstate * Dimension::F + istate;
                while (getline(ifs, eachline)) {
                    if (eachline.find("COMPLETE EXPRESSION.") != eachline.npos) {  // let `MPRINT=1`
                        for (int i = 0, idx = 0; i < natom; ++i) {                 ///< found nacv in 1/angstrom
                            ifs >> stmp                                            // skips
                                >> nac[(idx++) * Dimension::FF + IJ]               // nac(Xk, i, j)
                                >> nac[(idx++) * Dimension::FF + IJ]               // nac(Yk, i, j)
                                >> nac[(idx++) * Dimension::FF + IJ];              // nac(Zk, i, j)
                        }
                        for (int idx = 0; idx < Dimension::N; ++idx) {  // copy to the other half-side nacv
                            nac[idx * Dimension::FF + JI] = -nac[idx * Dimension::FF + IJ];
                        }
                        break;
                    }
                }
            }
            stat = 2;
        }
    }
    ifs.close();
    if (stat != 2) {
        bool* succ_ptr = _DataSet->def<bool>("iter.succ");
        succ_ptr[0]    = false;
        std::cerr << "fail in calling MNDO!\n";

        int* istep_ptr      = _DataSet->def<int>("iter.istep");
        std::string cmd_exe = utils::concat("cp ", directory, "/.mndoinp.", stat_in, "  ", directory, "/.mndoinp.",
                                            stat_in, ".fail.", istep_ptr[0]);
        system(cmd_exe.c_str());
        cmd_exe = utils::concat("cp ", directory, "/.mndoout.", stat_in, "  ", directory, "/.mndoout.", stat_in,
                                ".fail.", istep_ptr[0]);
    } else {
        bool* succ_ptr = _DataSet->def<bool>("iter.succ");
        succ_ptr[0]    = true;
    }
    return stat;
}

/**
 * parse frequency calculation log
 * where the frequency should be calculated by specifying (JOP=2), (ICROSS=0) and (KPRINT=1)
 */
int Model_Interf_MNDO::parse_hessian(const std::string& log) {
    std::ifstream ifs(log);

    // initialization of arrays
    ARRAY_CLEAR(w, Dimension::N);
    ARRAY_CLEAR(hess, Dimension::NN);
    ARRAY_CLEAR(Tmod, Dimension::NN);

    std::string stmp, eachline;
    int v1 = Dimension::N / 10, v2 = Dimension::N % 10, eqstat = 0, itmp;
    double dtmp = 0.0f;
    while (getline(ifs, eachline)) {
        if (eachline.find("GRADIENT NORM =") != eachline.npos) {
            double norm;
            std::string stmp;
            std::istringstream sstr(eachline);
            sstr >> stmp >> stmp >> stmp >> norm;
            eqstat = (norm < 10.0f) ? 0 : -1;
            if (eqstat != 0) std::cerr << "Warning Hessian is not used under equilibrium!\n";
        }
        if (eachline.find("FORCE CONSTANT MATRIX AFTER SYMMETRIZATION.") != eachline.npos) {
            for (int i = 0; i < v1; ++i) {
                getline(ifs, eachline);                    // skip a blankline
                for (int k = 0; k < 10; ++k) ifs >> itmp;  // skip header
                for (int j = 0; j < Dimension::N; ++j) {
                    ifs >> itmp;                                                              // skip index
                    for (int k = 0; k < 10; ++k) ifs >> hess[Dimension::N * j + i * 10 + k];  // read
                }
            }
            if (v2 == 0) continue;
            getline(ifs, eachline);                    // skip a blankline
            for (int k = 0; k < v2; ++k) ifs >> itmp;  // skip header
            for (int j = 0; j < Dimension::N; ++j) {
                ifs >> itmp;                                                               // skip index
                for (int k = 0; k < v2; ++k) ifs >> hess[Dimension::N * j + v1 * 10 + k];  // read
            }
        }
        if (eachline.find("EIGENVECTORS OF THE MASS-WEIGHTED") != eachline.npos) {
            for (int i = 0; i < v1; ++i) {
                getline(ifs, eachline);                                     // skip a blankline
                for (int k = 0; k < 10; ++k) ifs >> itmp;                   // skip header
                for (int k = 0; k < 10; ++k) {                              // read frequency
                    if (ifs >> dtmp) w[i * 10 + k] = dtmp / phys::au_2_wn;  // [convert wn to au]
                }
                for (int j = 0; j < Dimension::N; ++j) {  // read normal-mode
                    ifs >> itmp;                          // skip index
                    for (int k = 0; k < 10; ++k) {
                        if (ifs >> dtmp) Tmod[Dimension::N * j + i * 10 + k] = dtmp;
                    }
                }
            }
            if (v2 == 0) continue;
            getline(ifs, eachline);                                      // skip a blankline
            for (int k = 0; k < v2; ++k) ifs >> itmp;                    // skip header
            for (int k = 0; k < v2; ++k) {                               // read frequency
                if (ifs >> dtmp) w[v1 * 10 + k] = dtmp / phys::au_2_wn;  // [convert wn to au]
            }
            for (int j = 0; j < Dimension::N; ++j) {  // read normal-mode
                ifs >> itmp;                          // skip index
                for (int k = 0; k < v2; ++k) {
                    if (ifs >> dtmp) Tmod[Dimension::N * j + v1 * 10 + k] = dtmp;
                }
            }
        }
    }
    ifs.close();
    return 0;
}

/**
 * parse frequency calculation (JOP=2) and (KPRINT=1)
 */
int Model_Interf_MNDO::parse_hessian2(const std::string& molden_file) {
    std::ifstream ifs(molden_file);

    // initialization of arrays
    ARRAY_CLEAR(w, Dimension::N);
    ARRAY_CLEAR(hess, Dimension::NN);
    ARRAY_CLEAR(Tmod, Dimension::NN);

    std::string stmp, eachline;
    int itmp;
    double dtmp = 0.0f;
    while (getline(ifs, eachline)) {
        if (eachline.find("[FREQ]") != eachline.npos) {
            for (int j = 6; j < Dimension::N; ++j) {
                getline(ifs, eachline);
                std::stringstream sstr{eachline};
                if (sstr >> dtmp) w[j] = dtmp / phys::au_2_wn;
            }
        }
        if (eachline.find("[FR-NORM-COORD]") != eachline.npos) {
            for (int i = 6; i < Dimension::N; ++i) {
                ifs >> stmp >> itmp;
                for (int j = 0; j < Dimension::N; ++j) {
                    if (ifs >> dtmp) Tmod[j * Dimension::N + i] = dtmp * sqrt(mass[j] * phys::au_2_amu);
                }
            }
        }
    }
    ifs.close();

    // ARRAY_SHOW(w, 1, Dimension::N);
    // ARRAY_SHOW(mass, 1, Dimension::N);
    // ARRAY_SHOW(Tmod, Dimension::N, Dimension::N);
    // std::cout << FMT(8) << phys::au_2_amu << "\n";
    return 0;
}


/**
 * @brief this function generates normalmode trajectories
 * @param
 *      R: equilibrium configuration
 * @return int
 * @bug none
 */
int Model_Interf_MNDO::calc_normalmode() {
    // update mndohess input/output files
    if (isFileExists(".mndohess.out")) {
        new_task(".mndohess.in", "hess");
        int stat = system("mndo < .mndohess.in > .mndohess.out");
        if (stat != 0) return stat;
    }
    int stat = parse_hessian(".mndohess.out");
    assert(stat == 0);

    /**
     * @brief generate normal-mode trajectories
     * @details
     *      r = Xeq + M^(-1/2) * T * Q
     *      p = M^(1/2) * T * P
     *      here, r, p are cartesian postion & momentum, while Q, P are normal-mode postion & momentum
     */
    int ntime = 50;                              // report 50 points in [0,2*pi] period
    for (int iw = 0; iw < Dimension::N; ++iw) {  // iw-the mode
        std::string save = "normalmode";
        std::ofstream ofs(save + std::to_string(iw) + ".xyz");
        for (int itime = 0; itime < ntime; ++itime) {  // itime
            for (int i = 0; i < Dimension::N; ++i) {
                x[i] = x0[i] + Tmod[Dimension::N * i + iw] / std::sqrt(mass[i]) *
                                   std::sin(itime * phys::math::twopi / ntime) / std::sqrt(2.0f * w[iw]);
            }
            ofs << natom << std::endl;
            ofs << "time: " << (itime * phys::math::twopi) / (w[iw] * ntime) << std::endl;
            for (int i = 0, idx = 0; i < natom; ++i) {
                ofs << FMT(8) << ELEMENTS_LABEL[atoms[i]]   //
                    << FMT(8) << x[idx++] * phys::au_2_ang  //
                    << FMT(8) << x[idx++] * phys::au_2_ang  //
                    << FMT(8) << x[idx++] * phys::au_2_ang << std::endl;
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

    if (!isFileExists(".mndohess.out")) {
        new_task(".mndohess.in", "hess");
        stat = system("mndo < .mndohess.in > .mndohess.out");
        if (stat != 0) return stat;
    }
    stat = parse_hessian(".mndohess.out");
    assert(stat == 0);

    /**
     * mass-weighted Hessian to obtain normalmode
     */
    for (int i = 0, idx = 0; i < Dimension::N; ++i)
        for (int j = 0; j < Dimension::N; ++j, ++idx) hess[idx] /= sqrt(mass[i] * mass[j]);

    double Qeff = 6.0;
    for (int j = 0; j < Dimension::N; ++j) {
        if (j < 6 || w[j] < 1.0e-3) {  // cutoff low frequency
            x_sigma[j] = 0.0f;
            p_sigma[j] = 0.0f;
        } else {                             // NOTE: it's for normal-mode!
            double Qoverbeta = 1.0f / beta;  // 0.5f * w[j] / std::tanh(0.5f * beta * w[j]);
            Qeff += Qoverbeta * beta;
            x_sigma[j] = std::sqrt(Qoverbeta / (w[j] * w[j]));
            p_sigma[j] = std::sqrt(Qoverbeta);
        }
    }
    Qeff /= Dimension::N;
    Qeff = 0.5f * beta * w[Dimension::N - 1] / std::tanh(0.5f * beta * w[Dimension::N - 1]);

    std::cout << "Qeff = " << Qeff;

    for (int i = 0; i < Dimension::N; ++i) x[i] = x0[i];
    exec_kernel_impl();
    double ref_ener = E[0];
    std::cout << ref_ener;

    for (int isamp = 0; isamp < 500; ++isamp) {
        Kernel_Random::rand_gaussian(x, Dimension::N);
        Kernel_Random::rand_gaussian(p, Dimension::N);
        for (int i = 0; i < 6; ++i) x[i] = 0.0f, p[i] = 0.0f;
        for (int i = 0; i < Dimension::N; ++i) {
            x[i] *= x_sigma[i];
            p[i] *= p_sigma[i];
        }
        // transfrom normal-mode to cartesian coordinates
        ARRAY_MATMUL(x, Tmod, x, Dimension::N, Dimension::N, 1);
        ARRAY_MATMUL(p, Tmod, p, Dimension::N, Dimension::N, 1);
        for (int i = 0; i < Dimension::N; ++i) {
            x[i] = x[i] / std::sqrt(mass[i]) + x0[i];
            p[i] = p[i] * std::sqrt(mass[i]);
        }

        // fluctuation of Kinetic energy
        double Ekin = 0.0;
        for (int i = 0; i < Dimension::N; ++i) Ekin += 0.5f * p[i] * p[i] / mass[i];

        // fluctuation of potential energy
        exec_kernel_impl();
        double Epot = (E[0] - ref_ener);  // convert to au
        std::cout << Epot;
        std::cout << "::: " << beta * Epot / Qeff << " " << beta * Ekin / Qeff;

        if (beta * Epot > 10 || beta * Ekin > 10) {
            std::cout << "failed sample once";
            std::ofstream ofs(utils::concat("sampfail", isamp));
            ofs << FMT(8) << natom << std::endl;  // number of atoms
            ofs << FMT(8) << isamp << FMT(8) << beta << FMT(8) << Epot << FMT(8) << Ekin << std::endl;
            for (int i = 0, idx1 = 0, idx2 = 0; i < natom; ++i) {
                // output configuration ([angstrom])
                ofs << FMT(8) << ELEMENTS_LABEL[atoms[i]];
                for (int k = 0; k < 3; ++k, ++idx1) ofs << FMT(8) << x[idx1] * phys::au_2_ang;
                // output velocity ([angstrom/ps])
                for (int k = 0; k < 3; ++k, ++idx2) ofs << FMT(8) << p[idx2] / mass[idx2] * phys::au_2_angoverps;
                ofs << std::endl;
            }
            ofs.close();
            isamp--;
            continue;  // for resampling
        };

        std::cout << "test: " << beta * Epot << " " << beta * Ekin;

        std::ofstream ofs(utils::concat("samp", isamp));
        ofs << FMT(8) << natom << std::endl;  // number of atoms
        ofs << FMT(8) << isamp << std::endl;
        for (int i = 0, idx1 = 0, idx2 = 0; i < natom; ++i) {
            // output configuration ([angstrom])
            ofs << FMT(8) << ELEMENTS_LABEL[atoms[i]];
            for (int k = 0; k < 3; ++k, ++idx1) ofs << FMT(8) << x[idx1] * phys::au_2_ang;
            // output velocity ([angstrom/ps])
            for (int k = 0; k < 3; ++k, ++idx2) ofs << FMT(8) << p[idx2] / mass[idx2] * phys::au_2_angoverps;
            ofs << std::endl;
        }
        ofs.close();
    }
    return 0;
}

int Model_Interf_MNDO::calc_scan() {
    int istep = 0, readn;
    kids_real tmp;
    std::string eachline;

    // savefile_traj = utils::concat("traj-", 0, ".xyz");
    // savefile_ener = utils::concat("ener-", 0, ".dat");
    // savefile_grad = utils::concat("grad-", 0, ".dat");
    // savefile_nac  = utils::concat("nac-", 0, ".dat");
    // clearFile(savefile_traj);
    // clearFile(savefile_ener);
    // clearFile(savefile_grad);
    // clearFile(savefile_nac);

    std::ifstream ifs("scan.xyz");
    if (!ifs.is_open()) throw std::runtime_error("scan.xyz cannot open");

    while (getline(ifs, eachline)) {
        std::istringstream isstr(eachline);
        isstr >> readn;
        assert(natom == readn);
        getline(ifs, eachline);  // skip comments
        for (int iatom = 0, idx1 = 0; iatom < natom; ++iatom) {
            ifs >> eachline;                       // skip atomic flag
            for (int i = 0; i < 3; ++i, ++idx1) {  // read in [angstrom]
                if (ifs >> tmp) x[idx1] = tmp / phys::au_2_ang;
            }
        }
        getline(ifs, eachline);  // a line!!!
        exec_kernel_impl();
        istep++;
    }
    ifs.close();
    return 0;
}
};  // namespace PROJECT_NS