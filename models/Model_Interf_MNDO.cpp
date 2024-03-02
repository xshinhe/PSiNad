#include "Model_Interf_MNDO.h"

#include "../core/Element.h"
#include "../core/linalg.h"
#include "../kernels/Kernel_Declare.h"
#include "../kernels/Kernel_Random.h"
#include "../kernels/Kernel_Representation.h"
// #include "../mpi_utils.h"

#define ARRAY_SHOW(_A, _n1, _n2)                                                     \
    ({                                                                               \
        std::cout << "Show Array <" << #_A << ">\n";                                 \
        int _idxA = 0;                                                               \
        for (int _i = 0; _i < (_n1); ++_i) {                                         \
            for (int _j = 0; _j < (_n2); ++_j) std::cout << FMT(4) << (_A)[_idxA++]; \
            std::cout << std::endl;                                                  \
        }                                                                            \
    })

inline int removeFile(std::string& filename) { return remove(filename.c_str()); }

inline void clearFile(std::string& filename) { std::ofstream clear(filename, std::ios::trunc); }

inline void closeOFS(std::ofstream& ofs) {
    if (ofs.is_open()) ofs.close();
}

inline bool isFileExists(const std::string& name) { return std::ifstream{name.c_str()}.good(); }

namespace kids {

void Model_Interf_MNDO::read_param_impl(Param* PM) {
    Kernel_Representation::onthefly = true;

    // parse mndo99 input
    std::string mndo99inp = PM->get<std::string>("mndo99inp", LOC(), "null");
    natom                 = parse_mndo99(mndo99inp);

    assert(Dimension::N == 3 * natom);
    assert(Dimension::F <= ncigrd);
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

    V  = DS->def<double>("model.V", Dimension::FF);
    dV = DS->def<double>("model.dV", Dimension::NFF);
    // ddV  = DS->def<double>("model.ddV", Dimension::NNFF);
    E  = DS->def<double>("model.rep.E", Dimension::F);
    dE = DS->def<double>("model.rep.dE", Dimension::NFF);
    // ddE  = DS->def<double>("model.rep.ddE", Dimension::NNFF);
    nac      = DS->def<double>("model.rep.nac", Dimension::NFF);
    nac_prev = DS->def<double>("model.rep.nac_prev", Dimension::NFF);

    // read z index
    for (int i = 0, idx = 0, idxR = 0; i < natom; ++i) {
        atoms[i] = stoi(mndo99_data[idx++]);
        for (int j = 0; j < 3; ++j, ++idxR) {
            mass[idxR] = ELEMENTS_MASS[atoms[i]] / phys::au_2_amu;   // convert amu to au
            x0[idxR]   = stod(mndo99_data[idx++]) / phys::au_2_ang;  // convert angstrom into au
            idx++;                                                   // skip fixflag
        }
    }

    // read temperature
    double temperature = _Param->get<double>("temperature", LOC(), phys::temperature_d, 1.0f);
    beta               = 1.0f / (phys::au::k * temperature);  // don't ignore k_Boltzman

    // read task
    init_nuclinp = _Param->get<std::string>("init_nuclinp", LOC(), "#hess");
    read_flag    = _Param->get<int>("read_flag", LOC(), 0);

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

    if (init_nuclinp == "#hess") {
        std::ifstream ifs(".hess");  // read from .hess (saves in au)
        kids_real tmp;
        if (ifs) {
            for (int i = 0, idx = 0; i < Dimension::N; ++i)  // hessian
                for (int j = 0; j < Dimension::N; ++j, ++idx) {
                    if (ifs >> tmp) hess[idx] = tmp;       // Hessian without mass-weighted (in au)
                    hess[idx] /= sqrt(mass[i] * mass[j]);  // make it mass-weighted Hessian
                }
            for (int i = 0; i < Dimension::N; ++i)  // read frequency (in au)
                if (ifs >> tmp) w[i] = tmp;
            for (int i = 0, idx = 0; i < Dimension::N; ++i)  // read normal-mode vectors (in au)
                for (int j = 0; j < Dimension::N; ++j, ++idx) {
                    if (ifs >> tmp) Tmod[idx] = tmp;
                }
        } else {
            throw std::runtime_error("loss of .hess in currect directory");
        }
        ifs.close();

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
}

/**
 * @brief ForceField_init for mndo99
 * @param
 *     nr:  nuclear configuration
 *     np:  nuclear momentum
 *     m:  nuclear mass
 *   erho:  electronic density
 *   eeac:  electronic amplititude
 *   eocc:  electronic occupation
 *   fdim:  electronic freedoms
 *  itraj:  index of trajectory
 * @return int
 * @bug none
 */
void Model_Interf_MNDO::init_calc_impl(int stat) {
    if (init_nuclinp == "#hess") {  // assuming .hess is given
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

        ARRAY_SHOW(x, 1, Dimension::N);
        ARRAY_SHOW(p, 1, Dimension::N);

    } else if (init_nuclinp == "#fix") {  // for initial md
        for (int i = 0; i < Dimension::N; ++i) x[i] = x0[i], p[i] = 0.0f;
    } else {
        kids_real tmp;
        std::string stmp;
        std::ifstream ifs(utils::concat(init_nuclinp, stat));
        if (!ifs.is_open()) throw std::runtime_error("init_nuclinp cannot open");
        getline(ifs, stmp);  // skip atoms number
        getline(ifs, stmp);  // skip comments
        for (int iatom = 0, idx1 = 0, idx2 = 0; iatom < natom; ++iatom) {
            ifs >> stmp;                           // skip atomic flag
            for (int i = 0; i < 3; ++i, ++idx1) {  // read in [angstrom]
                if (ifs >> tmp) x[idx1] = tmp / phys::au_2_ang;
            }
            for (int i = 0; i < 3; ++i, ++idx2) {  // read in [angstrom/ps]
                if (ifs >> tmp) p[idx2] = tmp / phys::au_2_angoverps * mass[idx2];
            }
        }
    }
    _DataSet->set("init.x", x, Dimension::N);
    _DataSet->set("init.p", p, Dimension::N);

    refer = false;
    exec_kernel_impl(stat);
    refer = true;
}

int Model_Interf_MNDO::exec_kernel_impl(int stat_in) {
    // convert atomic unit
    for (int i = 0; i < Dimension::N; ++i) x[i] *= phys::au_2_ang;  ///< convert Bohr to Angstrom
    ARRAY_CLEAR(E, Dimension::F);
    ARRAY_CLEAR(dE, Dimension::NFF);

    // prepare a mndo99 input file
    std::string inpfile = utils::concat(".mndo99inp.", stat_in);
    std::string outfile = utils::concat(".mndo99out.", stat_in);

    new_task(inpfile, "def");  ///< generate a default task

    std::string cmd_exe = utils::concat("mndo99 < ", inpfile, " > ", outfile);
    int stat            = system(cmd_exe.c_str());
    if (stat != 0) return stat;

    parse_standard(outfile);  // parse in MNDO's units

    ARRAY_SHOW(E, 1, Dimension::F);
    // ARRAY_SHOW(dE, 1, Dimension::F);
    ARRAY_SHOW(nac, Dimension::N, Dimension::FF);

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

    // if (!ARRAY_ISFINITE(R, Dimension::N) || !ARRAY_ISFINITE(E, F) || !ARRAY_ISFINITE(dE, NFF)) {
    //     std::cout << "Problem for calling forcefield";
    //     stat = -1;
    // }
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

    assert(count_atom > 0);
    int data_size = mndo99_data.size();
    if (7 * count_atom == data_size) std::cout << "read from cartesian format";
    if (10 * count_atom == data_size) throw std::runtime_error("cannot read from zmat format");
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
int Model_Interf_MNDO::new_task(const std::string& file, const std::string& task_flag) {
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
            << FMT(8) << x[idx++] << FMT(4) << fixflag  // Xk, flag_Xk
            << FMT(8) << x[idx++] << FMT(4) << fixflag  // Yk, flag_Yk
            << FMT(8) << x[idx++] << FMT(4) << fixflag  // Zk, flag_Zk
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
            sstr >> stmp >> stmp >> stmp >> stmp >> istate;
            istate--;
            if (istate < Dimension::F) {  // skip additional state
                while (getline(ifs, eachline)) {
                    // if (eachline.find("Mult. 1") != eachline.npos) {
                    //     std::istringstream sstr(eachline);
                    //     sstr >> stmp >> jstate >> stmp >> stmp >> stmp >> stmp >> stmp >> stmp;
                    //     jstate--;
                    //     sstr >> E[jstate]; ///< saved in [ev]
                    // }
                    if (eachline.find("GRADIENTS (KCAL/(MOL*ANGSTROM))") != eachline.npos) {
                        for (int i = 0; i < 8; ++i) ifs >> stmp;
                        for (int i = 0, idx = 0; i < natom; ++i) {       ///< grad in kcalpmol/angstrom
                            ifs >> stmp >> stmp >> stmp >> stmp >> stmp  // skips
                                >> dE[(idx++) * Dimension::FF + istate * Dimension::Fadd1]   // grad(Xk,i,i)
                                >> dE[(idx++) * Dimension::FF + istate * Dimension::Fadd1]   // grad(Yk,i,i)
                                >> dE[(idx++) * Dimension::FF + istate * Dimension::Fadd1];  // grad(Zk,i,i)
                        }
                        break;
                    }
                }
            }
            stat = 1;
        }
        /**
         * @brief find non-adiabatic coupling(NAC) terms
         *  note we use the complete formalism (let `MPRINT=1`) in MNDO99
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
    using namespace Dimension;

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
        }
        if (eachline.find("FORCE CONSTANT MATRIX AFTER SYMMETRIZATIODimension::N.") != eachline.npos) {
            for (int i = 0; i < v1; ++i) {
                getline(ifs, eachline);                    // skip a blankline
                for (int k = 0; k < 10; ++k) ifs >> itmp;  // skip header
                for (int j = 0; j < Dimension::N; ++j) {
                    ifs >> itmp;                                                              // skip index
                    for (int k = 0; k < 10; ++k) ifs >> hess[Dimension::N * j + i * 10 + k];  // read
                }
            }
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

    std::ofstream ofs(".hess");
    for (int i = 0, idx = 0; i < Dimension::N; ++i) {  // hessian (in au)
        for (int j = 0; j < Dimension::N; ++j, ++idx) ofs << FMT(8) << hess[idx];
        ofs << std::endl;
    }
    for (int i = 0, idx = 0; i < Dimension::N; ++i) {  // frequancey (in au)
        ofs << FMT(8) << w[i];
    }
    ofs << std::endl;
    for (int i = 0, idx = 0; i < Dimension::N; ++i) {
        for (int j = 0; j < Dimension::N; ++j, ++idx) ofs << FMT(8) << Tmod[idx];
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
int Model_Interf_MNDO::calc_normalmode() {
    // update mndo99hess input/output files
    if (isFileExists(".mndo99hess.out")) {
        new_task(".mndo99hess.in", "hess");
        int stat = system("mndo99 < .mndo99hess.in > .mndo99hess.out");
        if (stat != 0) return stat;
    }
    int stat = parse_hessian(".mndo99hess.out");
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

    if (!isFileExists(".mndo99hess.out")) {
        new_task(".mndo99hess.in", "hess");
        stat = system("mndo99 < .mndo99hess.in > .mndo99hess.out");
        if (stat != 0) return stat;
    }
    stat = parse_hessian(".mndo99hess.out");
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
};  // namespace kids