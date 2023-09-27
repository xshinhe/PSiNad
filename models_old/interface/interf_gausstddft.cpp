#include "interf_gausstddft.h"

#include <unistd.h>

#include <algorithm>

#include "../../utils/definitions.h"
#include "../../utils/elements.h"

using namespace ARRAY_EG;


GAUSS16_ForceField::GAUSS16_ForceField(const Param& iparm) : Nad_ForceField(iparm) {
    type = ForceFieldOnTheFly;  // one-the-fly adiabatic forcefield!

    // reset unit convertor

    Param_Reset(iou.leng, phys::au_2_ang);
    Param_Reset(iou.mass, phys::au_2_amu);
    Param_Reset(iou.ener, phys::au_2_kcal_1mea);

    tag = name() + "_" + tag;

    // parse g16inp & determine system size
    F                  = Param_GetV(F, parm, 0);
    std::string g16inp = Param_GetT(std::string, parm, "g16inp", "null");
    if (g16inp == "null") LOG(FATAL) << "Cannot find <g16inp>" << g16inp << std::endl;
    natom = parse_g16(g16inp);

    CHECK_GT(natom, 0);
    CHECK_EQ(2 + 4 * natom,
             g16_data.size());  // xyz-cart format, for zmat: 10*natom =
                                // g16_data.size()
    CHECK_GT(F, 0);
    CHECK_EQ(N, 3 * natom);

    ALLOCATE_PTR_TO_VECTOR(atoms, natom);
    for (int i = 0, idxd = 2, idxR = 0; i < natom; ++i) {
        atoms[i] = Elements::ELEMENTS_ZNUM(g16_data[idxd++].c_str());
        for (int j = 0; j < 3; ++j) {
            mod_R0[idxR++] = stod(g16_data[idxd++]) / iou.leng;  // convert angstrom into au
        }
    }
    for (int i = 0, idx = 0; i < natom; ++i) {
        for (int j = 0; j < 3; ++j) mod_M[idx++] = ELEMENTS_MASS[atoms[i]] / iou.mass;  // covert amu into au
    }

    ALLOCATE_PTR_TO_VECTOR(mod_Hess, NN);
    ALLOCATE_PTR_TO_VECTOR(mod_Tmat, NN);
    ARRAY_CLEAR(mod_Hess, NN);
    ARRAY_CLEAR(mod_Tmat, NN);
    init_nuclinp = Param_GetT(std::string, parm, "init_nuclinp", "#hess");
    if (init_nuclinp == "#calc") {
        calc_hess(mod_R0, N);
        exit(0);  // prepare hessian file and exit normally
    }
    if (init_nuclinp == "#hess") {  // try to read hessian (saves in au)
        std::ifstream ifs(".hess");
        num_real tmp;
        if (ifs) {
            for (int i = 0, idx = 0; i < N; ++i)
                for (int j = 0; j < N; ++j, ++idx) {
                    if (ifs >> tmp) mod_Hess[idx] = tmp;         // in au (not mass-weighted!)
                    mod_Hess[idx] /= sqrt(mod_M[i] * mod_M[j]);  // make it mass-weighted (hessian
                                                                 // is indirectly used!!!)
                }
            for (int i = 0; i < N; ++i)
                if (ifs >> tmp) mod_W[i] = tmp;  // in au
            for (int i = 0, idx = 0; i < N; ++i)
                for (int j = 0; j < N; ++j, ++idx) {
                    if (ifs >> tmp) mod_Tmat[idx] = tmp;  // in au
                }
        } else {
            LOG(FATAL) << "You should prepare .hess file in currect directory" << std::endl;
        }
        ifs.close();
        // ARRAY_SHOW(mod_W, 1, N);

        // from hessian & temperature to prepare initial sampling
        temp = Param_GetV(temp, parm, 300.0f) / iou.temp;
        beta = 1 / temp;
        ALLOCATE_PTR_TO_VECTOR(nr_samp, N);
        ALLOCATE_PTR_TO_VECTOR(np_samp, N);
        for (int j = 0; j < N; ++j) {
            if (j < 6 || mod_W[j] < 1.0e-10) {  // cutoff low frequency
                mod_sigmaR[j] = 0.0f;
                mod_sigmaP[j] = 0.0f;
            } else {                             // M=1 for normal-mode
                double Qoverbeta = 1.0f / beta;  // classical
                // double Qoverbeta = 0.5f * mod_W[j] / std::tanh(0.5f * beta * mod_W[j]);
                mod_sigmaR[j] = std::sqrt(Qoverbeta / (mod_W[j] * mod_W[j]));
                mod_sigmaP[j] = std::sqrt(Qoverbeta);
            }
        }
    }

    Param_GetV(read_flag, parm, 0);  // for init_nuclinp

    // FINAL REPORT
    CheckForceField();
    LOG(INFO) << "gs keyword" << std::endl << keyword;
}

int GAUSS16_ForceField::ForceField_init(double* nr, double* np, double* nm, num_complex* erho, num_complex* eeac,
                                        int& eocc, const int& rdim, const int& fdim, const int& icycle) {
    savename = GAUSS16_ForceField::name() + "-" + std::to_string(icycle);
    std::ofstream ofsxyz(savename + ".xyz");
    ofsxyz << "";
    ofsxyz.close();
    std::ofstream ofsele(savename + ".dat");
    ofsele << "";
    ofsele.close();

    if (init_nuclinp == "#hess") {  // assuming .hess is given
        rand_gaussian(nr_samp, N);
        for (int i = 0; i < 6; ++i) nr_samp[i] = 0.0f;
        rand_gaussian(np_samp, N);
        for (int i = 0; i < 6; ++i) np_samp[i] = 0.0f;
        for (int i = 0; i < N; ++i) {
            nm[i] = mod_M[i];
            nr_samp[i] *= mod_sigmaR[i];
            np_samp[i] *= mod_sigmaP[i];
        }
        ARRAY_MATMUL(nr, mod_Tmat, nr_samp, N, N,
                     1);  // transfrom normal-mode to cart, nr can not used
                          // in nr_samp
        ARRAY_MATMUL(np, mod_Tmat, np_samp, N, N, 1);
        for (int i = 0; i < N; ++i) {
            nr[i] = nr[i] / std::sqrt(nm[i]) + mod_R0[i];
            np[i] = np[i] * std::sqrt(nm[i]) + mod_P0[i];
        }
    } else {  // read nr, np, directly from init_nuclinp
        num_real tmp;
        if (read_flag == 0) {
            std::ifstream ifs(init_nuclinp);
            for (int ic = 0; ic < icycle; ++ic)
                for (int i = 0; i < 2 * N; ++i) ifs >> tmp;
            for (int i = 0; i < N; ++i)
                if (ifs >> tmp) nr[i] = tmp;
            for (int i = 0; i < N; ++i)
                if (ifs >> tmp) np[i] = tmp;
            ifs.close();
        } else if (read_flag == 1) {  // rst file
            std::ifstream ifs(init_nuclinp);
            ifs >> tmp >> tmp >> tmp >> tmp;
            for (int i = 0; i < F; ++i) ifs >> tmp >> tmp;
            for (int i = 0; i < N; ++i) ifs >> nr[i] >> np[i];
            ifs.close();
        } else if (read_flag == 2) {  // break file
            std::ifstream ifs(init_nuclinp);
            ifs >> tmp >> tmp >> tmp >> tmp;
            for (int i = 0; i < N + 2 * F; ++i) ifs >> tmp >> tmp;
            for (int i = 0; i < N; ++i) ifs >> nr[i] >> np[i];
            ifs.close();
        }
        for (int i = 0; i < N; ++i) nm[i] = mod_M[i];
    }
    Nad_ForceField::ForceField_init_elec(erho, eeac, eocc, F, icycle);
    return 0;
}


int GAUSS16_ForceField::ForceField_epes(num_real* E, num_real* dE, num_real* ddE,
                                        num_real* R,  // input in au
                                        const int& flag, const int& rdim, const int& fdim) {
    if (flag > 1) LOG(FATAL) << "Not support ddE" << std::endl;
    // convert au into angstrom
    for (int i = 0; i < N; ++i) R[i] *= phys::au_2_ang;

    // update g16inp file
    std::ofstream ofs(std::string(".g16.") + std::to_string(mpi_rank) + ".com");
    for (auto line : g16_keyword) {
        if (line.size() > 6 && line.substr(1, 3) == "chk") {
            ofs << "%chk=check" + std::to_string(mpi_rank) + ".chk" << std::endl;
        } else {
            ofs << line << std::endl;
        }
    }
    for (auto line : g16_comment) ofs << line << std::endl;
    ofs << g16_data[0] << " " << g16_data[1] << std::endl;
    for (int i = 0, idx = 0; i < natom; ++i) {  // output position
        ofs << std::setprecision(12) << std::setiosflags(std::ios::scientific) << std::setiosflags(std::ios::left)
            << ELEMENTS_LABEL[atoms[i]] << " " << R[idx++] << " " << R[idx++] << " " << R[idx++] << " " << std::endl;
    }
    for (auto line : g16_addition) {
        if (line.size() > 6 && line.substr(1, 3) == "chk") {
            ofs << "%chk=check" + std::to_string(mpi_rank) + ".chk" << std::endl;
        } else {
            ofs << line << std::endl;
        }
    }
    ofs.close();

    // begin g16 calculation
    std::string cmd_exe = std::string("g16 .g16.") + std::to_string(mpi_rank) + ".com";
    int stat            = system(cmd_exe.c_str());
    if (stat != 0) return stat;

    for (int i = 0; i < F; ++i) E[i] = 0;
    for (int i = 0; i < NFF; ++i) dE[i] = 0;

    // parse the calculation
    std::ifstream ifs(std::string(".g16.") + std::to_string(mpi_rank) + ".log");
    std::string stmp, eachline;
    stat        = -1;
    int energyi = 0, forcei = 0, naci = 0;
    if (ifs)
        while (getline(ifs, eachline)) {
            // find energy
            if (eachline.find("SCF Done:") != eachline.npos && energyi == 0) {  // for ground energy
                // LOG(WARNING) << "FIND SCF:" << eachline;
                std::istringstream sstr(eachline);
                sstr >> stmp >> stmp >> stmp >> stmp >> E[energyi++];
            } else if (eachline.find("Total Energy, E(TD-HF/TD-DFT)") != eachline.npos) {  // excited state
                std::istringstream sstr(eachline);
                // LOG(WARNING) << "FIND TOT:" << eachline;
                sstr >> stmp >> stmp >> stmp >> stmp >> E[energyi++];  // au
            } else if (eachline.find("Forces (Hartrees/Bohr)") != eachline.npos) {
                // LOG(WARNING) << "FIND GAD:" << eachline;
                getline(ifs, eachline);  // skip
                getline(ifs, eachline);  // skip
                std::istringstream sstr(eachline);
                for (int i = 0, idx = 0; i < natom; ++i) {
                    ifs >> stmp >> stmp >> dE[(idx++) * FF + forcei * (F + 1)] >> dE[(idx++) * FF + forcei * (F + 1)] >>
                        dE[(idx++) * FF + forcei * (F + 1)];  // !< in hartree/angstrom
                }
                for (int i = 0; i < 5; ++i) getline(ifs, eachline);  // skip some lines
                forcei++;
            } else if (eachline.find("Nonadiabatic Coup. (Bohr^-1)") != eachline.npos) {
                // LOG(WARNING) << "FIND NAD:" << eachline;
                int istate = 0, jstate = naci;
                while ((istate + jstate) >= F - 1) jstate -= (F - 1 - istate), istate++;
                jstate += istate + 1;

                getline(ifs, eachline);  // skip
                getline(ifs, eachline);  // skip
                for (int i = 0, idx = 0; i < natom; ++i) {
                    ifs >> stmp >> stmp >> dE[(idx++) * FF + istate * F + jstate] >>
                        dE[(idx++) * FF + istate * F + jstate] >>
                        dE[(idx++) * FF + istate * F + jstate];  // !< in hartree/angstrom
                }
                for (int i = 0; i < N; ++i) dE[i * FF + jstate * F + istate] = -dE[i * FF + istate * F + jstate];
                naci++;
                stat = 0;
            }
        }

    std::ofstream ofsxyz(savename + ".xyz", std::ofstream::app);
    ofsxyz << natom << std::endl << "frame" << std::endl;
    for (int i = 0, idx = 0; i < natom; ++i) {  // output position
        ofsxyz << std::setprecision(12) << std::setiosflags(std::ios::scientific) << std::setiosflags(std::ios::left)
               << ELEMENTS_LABEL[atoms[i]] << " " << R[idx++] << " " << R[idx++] << " " << R[idx++] << " " << std::endl;
    }
    ofsxyz.close();
    std::ofstream ofsele(savename + ".dat", std::ofstream::app);
    ofsele << N << " " << F << std::endl;
    for (int i = 0; i < F; ++i)
        ofsele << std::setprecision(12) << std::setiosflags(std::ios::scientific) << std::setiosflags(std::ios::left)
               << E[i] << " ";
    ofsele << std::endl;
    for (int i = 0, idx = 0; i < N; ++i) {
        for (int j = 0; j < FF; ++j, ++idx)
            ofsele << std::setprecision(12) << std::setiosflags(std::ios::scientific)
                   << std::setiosflags(std::ios::left) << dE[idx] << " ";
        ofsele << std::endl;
    }
    ofsele.close();

    // dE = nacv * (E - E)
    for (int i = 0, idx = 0; i < N; ++i) {
        for (int j = 0; j < F; ++j) {
            for (int k = 0; k < F; ++k, ++idx) {
                if (j == k) {
                    dE[idx] = -dE[idx];  // force to grad
                    continue;
                }
                dE[idx] = dE[idx] * (E[k] - E[j]);
            }
        }
    }

    for (int i = 0; i < N; ++i) R[i] /= phys::au_2_ang;
    // ARRAY_SHOW(E, F, 1);
    // ARRAY_SHOW(dE, N, FF);
    // LOG(FATAL) << "test:" << std::endl;
    if (!ARRAY_ISFINITE(R, N) || !ARRAY_ISFINITE(E, F) || !ARRAY_ISFINITE(dE, NFF)) {
        LOG(WARNING) << "Problem for calling forcefield" << std::endl;
        stat = -1;  // CHECK if broken
    }
    // LOG(WARNING) << "stat" << stat;
    return stat;
}



int GAUSS16_ForceField::parse_g16(const std::string& g16inp) {  // save g16 input in vector
    enum { GAUSS16_SEC_KEYWORD, GAUSS16_SEC_COMMENT, GAUSS16_SEC_DATA, GAUSS16_SEC_ADDITION };
    int section = GAUSS16_SEC_KEYWORD, count_atom = 0;
    std::ifstream ifs(g16inp);
    if (ifs) {
        std::string eachline;
        while (getline(ifs, eachline)) {
            eachline.erase(eachline.find_last_not_of(" ") + 1);  // remove last blank
            // std::cout << eachline << std::endl;
            if (section == GAUSS16_SEC_KEYWORD) {
                if (eachline[0] != '%' && eachline[0] != '#') {
                    section = GAUSS16_SEC_COMMENT;
                    g16_comment.push_back(eachline);  // blank line
                    continue;
                }
                g16_keyword.push_back(eachline);  // save of keyword-lines
                std::string data, key, val;
                std::istringstream input(eachline.erase(0, 1));
                while (input >> data) {
                    std::istringstream pair(data);  // "key=val" pattern
                    getline(pair, key, '=');
                    getline(pair, val, '=');
                    std::for_each(key.begin(), key.end(), [](char& c) { c = ::tolower(c); });
                    keyword[key] = val;
                }
            } else if (section == GAUSS16_SEC_COMMENT) {
                g16_comment.push_back(eachline);  // comment
                getline(ifs, eachline);
                g16_comment.push_back(eachline);  // second blank
                section = GAUSS16_SEC_DATA;
            } else if (section == GAUSS16_SEC_DATA) {  //! note!: only support cart format
                if (eachline.size() == 0) {            // blankline, note file must in unix-encode (see doc2unix)
                    section = GAUSS16_SEC_ADDITION;
                    g16_addition.push_back(eachline);
                    continue;
                }
                // push data
                int ndata = 0;
                std::string data, databuf[16];
                std::istringstream input(eachline);
                while (input >> data) databuf[ndata++] = data;
                for (int i = 0; i < ndata; ++i) g16_data.push_back(databuf[i]);
                if (ndata == 4) count_atom++;
            } else {  // GAUSS16_SEC_ADDITION
                g16_addition.push_back(eachline);
            }
        }
    }
    ifs.close();

    // for (auto line : g16_keyword) LOG(INFO) << line;
    // LOG(WARNING);
    // for (auto line : g16_comment) LOG(INFO) << line;
    // LOG(WARNING);
    // for (auto line : g16_data) LOG(INFO) << line;
    // LOG(WARNING);
    // for (auto line : g16_addition) LOG(INFO) << line;

    return count_atom;
}

int GAUSS16_ForceField::calc_hess(num_real* R, const int& rdim) {
    // update g16hess input/output files
    if (!utils::isFileExists(".g16hess.log")) {
        std::ofstream ofs(".g16hess.com");
        for (auto line : g16_keyword) ofs << line << std::endl;  // how to revise keyword @TODO
        for (auto line : g16_comment) ofs << line << std::endl;
        for (int i = 0, idx = 0; i < natom; ++i) {  // output position
            ofs << std::setprecision(12) << std::setiosflags(std::ios::scientific) << std::setiosflags(std::ios::left)
                << atoms[i] << " " << R[idx++] * iou.leng << " 0 " << R[idx++] * iou.leng << " 0 "
                << R[idx++] * iou.leng << " 0 " << std::endl;
        }
        for (auto line : g16_addition) ofs << line << std::endl;
        ofs.close();
        // begin g16 calculation
        int stat = system("g16 .g16hess.com");
        if (stat != 0) return stat;
    }

    // initialization of arrays
    for (int i = 0; i < NN; ++i) mod_Hess[i] = 0.0f, mod_Tmat[i] = 0.0f;
    for (int i = 0; i < N; ++i) mod_W[i] = 0.0f;
    std::ifstream ifs(".g16hess.log");
    std::string stmp, eachline;
    int v1 = N / 10, v2 = N % 10, stat = -1, itmp;
    if (ifs)
        while (getline(ifs, eachline)) {
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
                    getline(ifs, eachline);                                 // skip a blankline
                    for (int k = 0; k < 10; ++k) ifs >> itmp;               // skip header
                    for (int k = 0; k < 10; ++k) ifs >> mod_W[i * 10 + k];  // read frequency
                    for (int j = 0; j < N; ++j) {
                        ifs >> itmp;  // skip index
                        for (int k = 0; k < 10; ++k) ifs >> mod_Tmat[N * j + i * 10 + k];
                    }
                }
                getline(ifs, eachline);                                  // skip a blankline
                for (int k = 0; k < v2; ++k) ifs >> itmp;                // skip header
                for (int k = 0; k < v2; ++k) ifs >> mod_W[v1 * 10 + k];  // read frequency
                for (int j = 0; j < N; ++j) {
                    ifs >> itmp;                                                        // skip index
                    for (int k = 0; k < v2; ++k) ifs >> mod_Tmat[N * j + v1 * 10 + k];  // read
                }
                stat = 0;
            }
        }
    ifs.close();
    for (int i = 0; i < N; ++i) mod_W[i] /= phys::au_2_wn;  // convert to au!

    std::ofstream ofs(".hess");
    for (int i = 0, idx = 0; i < N; ++i) {  // hessian (in au)
        for (int j = 0; j < N; ++j, ++idx)
            ofs << std::setprecision(12) << std::setiosflags(std::ios::scientific) << mod_Hess[idx] << " ";
        ofs << std::endl;
    }
    for (int i = 0, idx = 0; i < N; ++i) {  // frequancey (in au)
        ofs << std::setprecision(12) << std::setiosflags(std::ios::scientific) << mod_W[i] << " ";
    }
    ofs << std::endl;
    for (int i = 0, idx = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j, ++idx)
            ofs << std::setprecision(12) << std::setiosflags(std::ios::scientific) << mod_Tmat[idx] << " ";
        ofs << std::endl;
    }
    ofs.close();

    int ntime = 50;  // report 50 points

    ALLOCATE_PTR_TO_VECTOR(workr, N);
    num_real* nr = workr;
    for (int iw = 0; iw < N; ++iw) {  // iw-the mode
        std::string save = "normalmode";
        std::ofstream ofs(save + std::to_string(iw) + ".xyz");
        for (int itime = 0; itime < ntime; ++itime) {  // itime
            for (int i = 0; i < N; ++i) {
                nr[i] = mod_R0[i] + mod_Tmat[N * i + iw] / std::sqrt(mod_M[i]) *
                                        std::sin(itime * phys::math::twopi / ntime) / std::sqrt(2.0f * mod_W[iw]);
            }
            ofs << natom << std::endl;
            std::string comment_line = g16_comment[1];
            comment_line.erase(0, comment_line.find_first_not_of(" "));
            comment_line.erase(comment_line.find_last_not_of(" ") + 1);
            ofs << "comment: " << comment_line << ", time: " << (itime * phys::math::twopi) / (mod_W[iw] * ntime)
                << std::endl;
            for (int i = 0, idx = 0; i < natom; ++i) {
                ofs << std::setprecision(12) << std::setiosflags(std::ios::scientific) << ELEMENTS_LABEL[atoms[i]]
                    << " " << nr[idx++] * iou.leng << " " << nr[idx++] * iou.leng << " " << nr[idx++] * iou.leng
                    << std::endl;
            }
        }
        ofs.close();
    }
    return 0;
}
