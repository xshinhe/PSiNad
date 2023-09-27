#include "systembath.h"

#include "../../utils/definitions.h"
#include "../bath/bath.h"
#include "hamiltonian_data.h"

using namespace ARRAY_EG;

namespace PROJECT_NS {
void SystemBath_ForceField::read_param_impl(Param* P) {
    // read lambda flag
    lambda_option.flag = P->get<std::string>("lambda_flag", LOC(), "lambda");
    lambda_option.type = lambda_policy::_dict.at(lambda_option.flag);

    // read coupling flag
    coupling_option.flag = P->get<std::string>("coupling_flag", LOC(), "#se");
    coupling_option.type = coupling_policy::_dict.at(coupling_option.flag);

    // read spectrum flag
    spectrum_option.flag = P->get<std::string>("spectrum_flag", LOC(), "read");
    spectrum_option.type = spectrum_policy::_dict.at(spectrum_option.flag);

    // read Hsys flag
    Hsys_option.flag = P->get<std::string>("Hsys_flag", LOC(), "read");
    Hsys_option.type = Hsys_policy::_dict.at(Hsys_option.flag);


    Nb    = P->get<int>("Nb", LOC());
    nbath = P->get<int>("nbath", LOC());
    CHECK_EQ(Nb * nbath, N);  // Nb*nbath must be N


    bias  = P->get<double>("bias", LOC(), phys::energy_d, 1.0f);
    delta = P->get<double>("delta", LOC(), phys::energy_d, 1.0f);

    // read omegac
    omegac = P->get<double>("omegac", LOC(), phys::energy_d, 1.0f);
    CHECK_GT(omegac, 0);


    switch (lambda_option.type) {
        case lambda_policy::lambUsed: {
            double strength = P->get<double>("strength", LOC(), phys::energy_d, 1.0f);
            lambda          = strength;
            break;
        }
        case lambda_policy::alphaUsed: {
            double strength = P->get<double>("strength", LOC(), 1.0f);
            lambda          = 0.5f * omegac * strength;
            break;
        }
        case lambda_policy::etaUsed: {
            double strength = P->get<double>("strength", LOC(), phys::energy_d, 1.0f);
            lambda          = 0.5f * strength;
            break;
        }
        case lambda_policy::ergUsed: {
            double strength = P->get<double>("strength", LOC(), phys::energy_d, 1.0f);
            lambda          = 0.25f * strength;
            break;
        }
        default:
            LOG(FATAL) << "Unknown st_flag";
    }
    CHECK_GT(lambda, 0);

    // read temperature
    double temp = P->get<double>("temp", LOC(), phys::temperature_d, 1.0f);
    beta        = 1.0f / (phys::au::k * temp);  // don't ignore k_Boltzman
    CHECK_GT(beta, 0);



    // allocate Q tensor (the coupling ternsor) & initialization

    switch (coup_type) {
        case BathPolicy::SpinBosonCoup: {
            CHECK_EQ(nbath, 1);
            CHECK_EQ(F, 2);                                       // nbath must be 1, F must be 2
            Q[0] = 1.0f, Q[1] = 0.0f, Q[2] = 0.0f, Q[3] = -1.0f;  // sigmaz
            break;
        }
        case BathPolicy::SiteExcitonCoup: {
            // always be careful with nbath and F
            CHECK_LE(nbath, F);  // nbath should equal or less than F. (a ground state can be placed in last)
            for (int i = 0, idx = 0; i < nbath; ++i) {
                for (int j = 0; j < F; ++j) {
                    for (int k = 0; k < F; ++k, ++idx) Q[idx] = (i == j && i == k) ? 1.0f : 0.0f;
                }
            }
            break;
        }
        case BathPolicy::GeneralCoup:
        default: {
            std::ifstream ifs(coup_flag);
            num_real tmp;
            for (int i = 0; i < nbath * F * F; ++i)
                if (ifs >> tmp) Q[i] = tmp;
            ifs.close();
        }
    }
}

void SystemBath_ForceField::init_data_impl(State* S) {
    Hsys = S->reg<num_real>("model.Hsys", F * F);

    switch (Hsys_option.type) {
        case Hsys_policy::SpinBosonHsys: {
            CHECK_EQ(F, 2);
            Param_Reset(L, 2);
            double HSB[4] = {bias, delta, delta, -bias};
            for (int i = 0; i < FF; ++i) Hsys[i] = HSB[i];
            break;
        }
        case Hsys_policy::FMOHsys: {
            CHECK_EQ(nbath, 7);
            CHECK_GE(F, 7);  // F = 7, or F = 8 (include ground state in the last; be careful)
            Param_Reset(L, 1);
            double tmp_unit = phys::au_2_wn;
            memset(Hsys, 0, FF * sizeof(num_real));
            for (int i = 0, ik = 0; i < 7; ++i) {
                for (int k = 0; k < 7; ++k, ++ik) { Hsys[i * F + k] = HFMO_data[ik] / tmp_unit; }
            }
            break;
        }
        case Hsys_policy::SF3aHsys: {
            CHECK_EQ(F, 3);
            Param_Reset(L, 1);
            double tmp_unit = phys::au_2_ev;
            for (int i = 0; i < FF; ++i) Hsys[i] = HSF3a_data[i] / tmp_unit;
            break;
        }
        case Hsys_policy::SF3bHsys: {
            CHECK_EQ(F, 3);
            Param_Reset(L, 1);
            double tmp_unit = phys::au_2_ev;
            for (int i = 0; i < FF; ++i) Hsys[i] = HSF3b_data[i] / tmp_unit;
            break;
        }
        case Hsys_policy::SF5aHsys: {
            CHECK_EQ(F, 5);
            Param_Reset(L, 1);
            double tmp_unit = phys::au_2_ev;
            for (int i = 0; i < FF; ++i) Hsys[i] = HSF5a_data[i] / tmp_unit;
            break;
        }
        case Hsys_policy::SF5bHsys: {
            CHECK_EQ(F, 5);
            Param_Reset(L, 1);
            double tmp_unit = phys::au_2_ev;
            for (int i = 0; i < FF; ++i) Hsys[i] = HSF5b_data[i] / tmp_unit;
            break;
        }
        case Hsys_policy::FCPHsys: {
            CHECK_EQ(nbath, 9);
            CHECK_GE(F, 9);  // F = 9, or F = 10 (include ground state; be careful)
            Param_Reset(L, 1);
            double tmp_unit = phys::au_2_wn;
            memset(Hsys, 0, FF * sizeof(num_real));
            for (int i = 0, ik = 0; i < 9; ++i) {
                for (int k = 0; k < 9; ++k, ++ik) { Hsys[i * F + k] = HFCP_data[ik] / tmp_unit; }
            }
            break;
        }
        case Hsys_policy::AGGHsys: {
            CHECK_GE(F, 2);
            Param_Reset(L, 1);
            double delta = Param_GetQ(phys::energy_d, "delta", parm, 1.0f);
            for (int i = 0, ik = 0; i < F; ++i) {
                for (int k = 0; k < F; ++k, ++ik) Hsys[ik] = (i - k == 1 || k - i == 1) ? delta : 0.0f;
            }
            break;
        }
        case Hsys_policy::CYCHsys: {
            CHECK_GE(F, 2);
            Param_Reset(L, 1);
            double delta = Param_GetQ(phys::energy_d, "delta", parm, 1.0f);
            for (int i = 0, ik = 0; i < F; ++i) {
                for (int k = 0; k < F; ++k, ++ik) Hsys[ik] = (i - k == 1 || k - i == 1) ? delta : 0.0f;
            }
            Hsys[0 * F + (F - 1)] = delta;
            Hsys[(F - 1) * F + 0] = delta;
            break;
        }
        case Hsys_policy::GivenHsys:
        default: {
            std::ifstream ifs("Hsys.dat");

            std::string H_unit_str;
            std::string firstline;
            getline(ifs, firstline);
            std::stringstream sstr(firstline);
            sstr >> H_unit_str;  ///< the firstline stores H's unit

            // parse unit
            double H_unit = phys::us::conv(phys::au::unit, phys::us::parse(H_unit_str));

            num_real val;
            for (int i = 0; i < FF; ++i)
                if (ifs >> val) Hsys[i] = val / H_unit;
            ifs.close();
        }
    }

    //
    // 2): initialize bath spectrum (save order: L > nbath > Nb > FF)
    ALLOCATE_PTR_TO_VECTOR(omegas, Nb);
    ALLOCATE_PTR_TO_VECTOR(coeffs, Nb);
    ALLOCATE_PTR_TO_VECTOR(Q, nbath * FF);
    ALLOCATE_PTR_TO_VECTOR(CL, L * Nb);
    ALLOCATE_PTR_TO_VECTOR(QL, L * nbath * FF);
    ALLOCATE_PTR_TO_VECTOR(Xnj, NFF);

    // get discretized spectrum
    mybath    = new Bath(parm, nbath, F);
    coup_type = mybath->coup_type;
    mybath->Discretization(omegas, coeffs, Nb);


    // add reorgainization energy to thermostated site
    bool addl = Param_GetT(int, parm, "addl", false);
    if (mybath->coup_type == BathPolicy::SiteExcitonCoup && addl) {
        for (int i = 0, idxQ = 0; i < nbath; ++i) {
            for (int j = 0, idxH = 0; j < F; ++j) {
                for (int k = 0; k < F; ++k, ++idxH, ++idxQ) {
                    Hsys[idxH] += mybath->lambda * Q[idxQ] * Q[idxQ];  // correct as long as Q is correct
                }
            }
        }
    }

    // copy mass and frequency for trajectory-based methods
    for (int i = 0; i < N; ++i) { mod_M[i] = 1.0f, mod_W[i] = omegas[i % Nb]; }

    // copy Q & QL for functional-based methods (HEOM & SSE)
    for (int i = 0; i < nbath * FF; ++i) Q[i] = mybath->Q[i];

    ARRAY_CLEAR(QL, L * nbath * FF);
    for (int ibath = 0, idx = 0; ibath < nbath; ++ibath) {
        for (int j = 0; j < Nb; ++j) {
            for (int i = 0, iL = 0; i < FF; ++i, ++idx) {
                double Qval = Q[ibath * FF + i];
                if (Qval != 0) {
                    QL[iL * nbath * FF + ibath * FF + i] = 1;                 // record there is a nonzero value
                    CL[iL * Nb + j]                      = coeffs[j] * Qval;  // reduce this value to CL
                    iL++;
                    if (iL > L) LOG(FATAL);  // Q shoule be sparsed!
                }
                Xnj[idx] = coeffs[j] * Q[ibath * FF + i];  // merge coeffs into Q to obtain Xnj
            }
        }
    }

    // utils to output some correlation functions // @move
    if (Hsys_flag == "#rub") {
        // isotopic effect in rubene
        double qcoeff = Param_GetV(qcoeff, parm, 1.0f);
        for (int i = 0; i < Nb; ++i) coeffs[i] *= qcoeff, omegas[i] *= qcoeff;
        for (int i = 0; i < N; ++i) mod_W[i] *= qcoeff;

        ALLOCATE_PTR_TO_VECTOR(workr, 5 * FF);  // enlarge space if necessary

        num_real* tmpE    = workr;  // only F-size need
        num_real* tmpT    = workr + FF;
        num_real* tmpJ    = workr + 2 * FF;
        num_real* tmpRho  = workr + 3 * FF;
        num_real* tmpJRho = workr + 4 * FF;

        /* @deprecated
        double lambda_re = 0.0f;
        for (int i = 0; i < Nb; ++i) { lambda_re += 0.5f * coeffs[i] * coeffs[i] / (omegas[i] * omegas[i]); }
        */
        double betaeff = mybath->beta * exp(-mybath->beta * mybath->lambda / 3);
        EigenSolve(tmpE, tmpT, Hsys, F);
        for (int i = 0; i < F; ++i) tmpE[i] = exp(-betaeff * tmpE[i]);
        ARRAY_MATMUL3_TRANS2(tmpRho, tmpT, tmpE, tmpT, F, F, 0, F);
        double tr = 0.0f;
        for (int i = 0; i < F; ++i) tr += tmpRho[i * (F + 1)];
        for (int i = 0; i < FF; ++i) tmpRho[i] /= tr;

        std::ofstream ofs("rhomat");
        for (int i = 0, idx = 0; i < F; ++i) {
            for (int j = 0; j < F; ++j, ++idx) { ofs << FMT(8) << tmpRho[idx]; }
            ofs << std::endl;
        }
        ofs.close();

        double delta = Hsys[1];
        for (int i = 0, idx = 0; i < F; ++i) {
            for (int j = 0; j < F; ++j, ++idx) {
                if (i + 1 == j || (i == F - 1 && j == 0)) {
                    tmpJ[idx] = 1.0f * delta;
                } else if (i - 1 == j || (i == 0 && j == F - 1)) {
                    tmpJ[idx] = -1.0f * delta;
                } else {
                    tmpJ[idx] = 0.0f;
                }
            }
        }
        ARRAY_MATMUL(tmpJRho, tmpJ, tmpRho, F, F, F);

        ofs.open("tcf0");
        for (int i = 0, idx = 0; i < F; ++i) {
            for (int j = 0; j < F; ++j, ++idx) { ofs << FMT(8) << -tmpJRho[idx]; }
            ofs << std::endl;
        }
        ofs.close();

        ofs.open("tcft");
        for (int i = 0, idx = 0; i < F; ++i) {
            for (int j = 0; j < F; ++j, ++idx) { ofs << FMT(8) << tmpJ[idx]; }
            ofs << std::endl;
        }
        ofs.close();
    } else if (Hsys_flag == "#fcp-sp") {
        std::ofstream ofs;
        ofs.open("tcf0");
        for (int i = 0, ik = 0; i < F; ++i) {
            for (int k = 0; k < F; ++k, ++ik) {
                if (i == F - 1 && k != F - 1 || i != F - 1 && k == F - 1)
                    ofs << FMT(8) << 1;
                else
                    ofs << FMT(8) << 0;
            }
            ofs << std::endl;
        }
        ofs.close();

        ofs.open("tcft");
        for (int i = 0, ik = 0; i < F; ++i) {
            for (int k = 0; k < F; ++k, ++ik) {
                if (i == F - 1 && k != F - 1 || i != F - 1 && k == F - 1)
                    ofs << FMT(8) << 1;
                else
                    ofs << FMT(8) << 0;
            }
            ofs << std::endl;
        }
        ofs.close();
    }

    BO_ForceField::ForceField_init_default_build(mybath->beta, N);
    CheckForceField();
    tag = name() + "_" + tag;
}

};  // namespace PROJECT_NS


SystemBath_ForceField::SystemBath_ForceField(const Param& iparm, const int& child)
    : Nad_ForceField(iparm){};  // for child to call

SystemBath_ForceField::SystemBath_ForceField(const Param& iparm) : Nad_ForceField(iparm) {
    Nb    = Param_GetT(int, parm, "Nb");     // Num. for Discretization of one bath
    nbath = Param_GetT(int, parm, "nbath");  // Num. for total bath
    CHECK_EQ(Nb * nbath, N);                 // Nb*nbath must be N

    // 1) initial system Hamiltonian
    ALLOCATE_PTR_TO_VECTOR(Hsys, FF);
    Hsys_flag = Param_GetT(std::string, parm, "Hsys_flag", "read");
    Hsys_type = Hsys_policy::_dict_Hsys.at(Hsys_flag);

    switch (Hsys_type) {
        case Hsys_policy::SpinBosonHsys: {
            CHECK_EQ(F, 2);
            Param_Reset(L, 2);
            double bias   = Param_GetQ(phys::energy_d, parm, "bias", 1.0f);
            double delta  = Param_GetQ(phys::energy_d, parm, "delta", 1.0f);
            double HSB[4] = {bias, delta, delta, -bias};
            for (int i = 0; i < FF; ++i) Hsys[i] = HSB[i];
            break;
        }
        case Hsys_policy::FMOHsys: {
            CHECK_EQ(nbath, 7);
            CHECK_GE(F, 7);  // F = 7, or F = 8 (include ground state in the last; be careful)
            Param_Reset(L, 1);
            double tmp_unit = phys::au_2_wn;
            memset(Hsys, 0, FF * sizeof(num_real));
            for (int i = 0, ik = 0; i < 7; ++i) {
                for (int k = 0; k < 7; ++k, ++ik) { Hsys[i * F + k] = HFMO_data[ik] / tmp_unit; }
            }
            break;
        }
        case Hsys_policy::SF3aHsys: {
            CHECK_EQ(F, 3);
            Param_Reset(L, 1);
            double tmp_unit = phys::au_2_ev;
            for (int i = 0; i < FF; ++i) Hsys[i] = HSF3a_data[i] / tmp_unit;
            break;
        }
        case Hsys_policy::SF3bHsys: {
            CHECK_EQ(F, 3);
            Param_Reset(L, 1);
            double tmp_unit = phys::au_2_ev;
            for (int i = 0; i < FF; ++i) Hsys[i] = HSF3b_data[i] / tmp_unit;
            break;
        }
        case Hsys_policy::SF5aHsys: {
            CHECK_EQ(F, 5);
            Param_Reset(L, 1);
            double tmp_unit = phys::au_2_ev;
            for (int i = 0; i < FF; ++i) Hsys[i] = HSF5a_data[i] / tmp_unit;
            break;
        }
        case Hsys_policy::SF5bHsys: {
            CHECK_EQ(F, 5);
            Param_Reset(L, 1);
            double tmp_unit = phys::au_2_ev;
            for (int i = 0; i < FF; ++i) Hsys[i] = HSF5b_data[i] / tmp_unit;
            break;
        }
        case Hsys_policy::FCPHsys: {
            CHECK_EQ(nbath, 9);
            CHECK_GE(F, 9);  // F = 9, or F = 10 (include ground state; be careful)
            Param_Reset(L, 1);
            double tmp_unit = phys::au_2_wn;
            memset(Hsys, 0, FF * sizeof(num_real));
            for (int i = 0, ik = 0; i < 9; ++i) {
                for (int k = 0; k < 9; ++k, ++ik) { Hsys[i * F + k] = HFCP_data[ik] / tmp_unit; }
            }
            break;
        }
        case Hsys_policy::AGGHsys: {
            CHECK_GE(F, 2);
            Param_Reset(L, 1);
            double delta = Param_GetQ(phys::energy_d, "delta", parm, 1.0f);
            for (int i = 0, ik = 0; i < F; ++i) {
                for (int k = 0; k < F; ++k, ++ik) Hsys[ik] = (i - k == 1 || k - i == 1) ? delta : 0.0f;
            }
            break;
        }
        case Hsys_policy::CYCHsys: {
            CHECK_GE(F, 2);
            Param_Reset(L, 1);
            double delta = Param_GetQ(phys::energy_d, "delta", parm, 1.0f);
            for (int i = 0, ik = 0; i < F; ++i) {
                for (int k = 0; k < F; ++k, ++ik) Hsys[ik] = (i - k == 1 || k - i == 1) ? delta : 0.0f;
            }
            Hsys[0 * F + (F - 1)] = delta;
            Hsys[(F - 1) * F + 0] = delta;
            break;
        }
        case Hsys_policy::GivenHsys:
        default: {
            std::ifstream ifs("Hsys.dat");

            std::string H_unit_str;
            std::string firstline;
            getline(ifs, firstline);
            std::stringstream sstr(firstline);
            sstr >> H_unit_str;  ///< the firstline stores H's unit

            // parse unit
            double H_unit = phys::us::conv(phys::au::unit, phys::us::parse(H_unit_str));

            num_real val;
            for (int i = 0; i < FF; ++i)
                if (ifs >> val) Hsys[i] = val / H_unit;
            ifs.close();
        }
    }


    // 2): initialize bath spectrum (save order: L > nbath > Nb > FF)
    ALLOCATE_PTR_TO_VECTOR(omegas, Nb);
    ALLOCATE_PTR_TO_VECTOR(coeffs, Nb);
    ALLOCATE_PTR_TO_VECTOR(Q, nbath * FF);
    ALLOCATE_PTR_TO_VECTOR(CL, L * Nb);
    ALLOCATE_PTR_TO_VECTOR(QL, L * nbath * FF);
    ALLOCATE_PTR_TO_VECTOR(Xnj, NFF);

    // get discretized spectrum
    mybath    = new Bath(parm, nbath, F);
    coup_type = mybath->coup_type;
    mybath->Discretization(omegas, coeffs, Nb);


    // add reorgainization energy to thermostated site
    bool addl = Param_GetT(int, parm, "addl", false);
    if (mybath->coup_type == BathPolicy::SiteExcitonCoup && addl) {
        for (int i = 0, idxQ = 0; i < nbath; ++i) {
            for (int j = 0, idxH = 0; j < F; ++j) {
                for (int k = 0; k < F; ++k, ++idxH, ++idxQ) {
                    Hsys[idxH] += mybath->lambda * Q[idxQ] * Q[idxQ];  // correct as long as Q is correct
                }
            }
        }
    }

    // copy mass and frequency for trajectory-based methods
    for (int i = 0; i < N; ++i) { mod_M[i] = 1.0f, mod_W[i] = omegas[i % Nb]; }

    // copy Q & QL for functional-based methods (HEOM & SSE)
    for (int i = 0; i < nbath * FF; ++i) Q[i] = mybath->Q[i];

    ARRAY_CLEAR(QL, L * nbath * FF);
    for (int ibath = 0, idx = 0; ibath < nbath; ++ibath) {
        for (int j = 0; j < Nb; ++j) {
            for (int i = 0, iL = 0; i < FF; ++i, ++idx) {
                double Qval = Q[ibath * FF + i];
                if (Qval != 0) {
                    QL[iL * nbath * FF + ibath * FF + i] = 1;                 // record there is a nonzero value
                    CL[iL * Nb + j]                      = coeffs[j] * Qval;  // reduce this value to CL
                    iL++;
                    if (iL > L) LOG(FATAL);  // Q shoule be sparsed!
                }
                Xnj[idx] = coeffs[j] * Q[ibath * FF + i];  // merge coeffs into Q to obtain Xnj
            }
        }
    }

    // utils to output some correlation functions // @move
    if (Hsys_flag == "#rub") {
        // isotopic effect in rubene
        double qcoeff = Param_GetV(qcoeff, parm, 1.0f);
        for (int i = 0; i < Nb; ++i) coeffs[i] *= qcoeff, omegas[i] *= qcoeff;
        for (int i = 0; i < N; ++i) mod_W[i] *= qcoeff;

        ALLOCATE_PTR_TO_VECTOR(workr, 5 * FF);  // enlarge space if necessary

        num_real* tmpE    = workr;  // only F-size need
        num_real* tmpT    = workr + FF;
        num_real* tmpJ    = workr + 2 * FF;
        num_real* tmpRho  = workr + 3 * FF;
        num_real* tmpJRho = workr + 4 * FF;

        /* @deprecated
        double lambda_re = 0.0f;
        for (int i = 0; i < Nb; ++i) { lambda_re += 0.5f * coeffs[i] * coeffs[i] / (omegas[i] * omegas[i]); }
        */
        double betaeff = mybath->beta * exp(-mybath->beta * mybath->lambda / 3);
        EigenSolve(tmpE, tmpT, Hsys, F);
        for (int i = 0; i < F; ++i) tmpE[i] = exp(-betaeff * tmpE[i]);
        ARRAY_MATMUL3_TRANS2(tmpRho, tmpT, tmpE, tmpT, F, F, 0, F);
        double tr = 0.0f;
        for (int i = 0; i < F; ++i) tr += tmpRho[i * (F + 1)];
        for (int i = 0; i < FF; ++i) tmpRho[i] /= tr;

        std::ofstream ofs("rhomat");
        for (int i = 0, idx = 0; i < F; ++i) {
            for (int j = 0; j < F; ++j, ++idx) { ofs << FMT(8) << tmpRho[idx]; }
            ofs << std::endl;
        }
        ofs.close();

        double delta = Hsys[1];
        for (int i = 0, idx = 0; i < F; ++i) {
            for (int j = 0; j < F; ++j, ++idx) {
                if (i + 1 == j || (i == F - 1 && j == 0)) {
                    tmpJ[idx] = 1.0f * delta;
                } else if (i - 1 == j || (i == 0 && j == F - 1)) {
                    tmpJ[idx] = -1.0f * delta;
                } else {
                    tmpJ[idx] = 0.0f;
                }
            }
        }
        ARRAY_MATMUL(tmpJRho, tmpJ, tmpRho, F, F, F);

        ofs.open("tcf0");
        for (int i = 0, idx = 0; i < F; ++i) {
            for (int j = 0; j < F; ++j, ++idx) { ofs << FMT(8) << -tmpJRho[idx]; }
            ofs << std::endl;
        }
        ofs.close();

        ofs.open("tcft");
        for (int i = 0, idx = 0; i < F; ++i) {
            for (int j = 0; j < F; ++j, ++idx) { ofs << FMT(8) << tmpJ[idx]; }
            ofs << std::endl;
        }
        ofs.close();
    } else if (Hsys_flag == "#fcp-sp") {
        std::ofstream ofs;
        ofs.open("tcf0");
        for (int i = 0, ik = 0; i < F; ++i) {
            for (int k = 0; k < F; ++k, ++ik) {
                if (i == F - 1 && k != F - 1 || i != F - 1 && k == F - 1)
                    ofs << FMT(8) << 1;
                else
                    ofs << FMT(8) << 0;
            }
            ofs << std::endl;
        }
        ofs.close();

        ofs.open("tcft");
        for (int i = 0, ik = 0; i < F; ++i) {
            for (int k = 0; k < F; ++k, ++ik) {
                if (i == F - 1 && k != F - 1 || i != F - 1 && k == F - 1)
                    ofs << FMT(8) << 1;
                else
                    ofs << FMT(8) << 0;
            }
            ofs << std::endl;
        }
        ofs.close();
    }

    BO_ForceField::ForceField_init_default_build(mybath->beta, N);
    CheckForceField();
    tag = name() + "_" + tag;
};

SystemBath_ForceField::~SystemBath_ForceField() {
    // if (mybath != nullptr) { delete mybath; }
}

int SystemBath_ForceField::ForceField_init(double* nr, double* np, double* nm, num_complex* erho, num_complex* eeac,
                                           int& eocc, const int& rdim, const int& fdim, const int& icycle) {
    Nad_ForceField::ForceField_init_default(nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle);
    return 0;
}

int SystemBath_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                           const int& rdim) {
    plFunction();

    double Vnucl = 0.0f;
    for (int j = 0; j < N; ++j) {  // @profiling: self-time from 20% => 13%
        double mwRj = mod_W[j] * R[j];
        Vnucl += mod_M[j] * mwRj * mwRj;  // @note: don't add kenitic energy, that is not include in PIMD
        dV[j] = mod_W[j] * mwRj;
    }
    V[0] = 0.5f * Vnucl;

    if (flag < 2) return 0;
    ARRAY_CLEAR(ddV, NN);
    for (int j = 0, idx = 0, add = N + 1; j < N; ++j, idx += add) ddV[idx] = mod_M[j] * mod_W[j] * mod_W[j];
    return 0;
}

int SystemBath_ForceField::ForceField_epes(double* V, double* dV, double* ddV, double* R, const int& flag,
                                           const int& rdim, const int& fdim, const int& itraj, const int& isamp) {
    switch (coup_type) {
        case BathPolicy::SpinBosonCoup:
            ForceField_epes_SpinBoson(V, dV, ddV, R, flag, rdim, fdim, itraj, isamp);
            break;
        case BathPolicy::SiteExcitonCoup:
            ForceField_epes_SiteExciton(V, dV, ddV, R, flag, rdim, fdim, itraj, isamp);
            break;
        case BathPolicy::GeneralCoup:
        default:
            ForceField_epes_General(V, dV, ddV, R, flag, rdim, fdim, itraj, isamp);
    }
    return 0;
}

int SystemBath_ForceField::ForceField_epes_SpinBoson(double* V, double* dV, double* ddV, double* R, const int& flag,
                                                     const int& rdim, const int& fdim, const int& itraj,
                                                     const int& isamp) {
    // V
    for (int i = 0; i < FF; ++i) V[i] = Hsys[i];
    for (int j = 0; j < N; ++j) {
        V[0] += coeffs[j] * R[j];
        V[3] -= coeffs[j] * R[j];
    }
    if (flag < 1) return 0;

    // dV
    if (isamp < 1) {  // independent with R, so only calculate in the first time
        int idxdV = 0;
        for (int j = 0; j < N; ++j) {
            dV[idxdV++] = coeffs[j];
            dV[idxdV++] = 0.0f;
            dV[idxdV++] = 0.0f;
            dV[idxdV++] = -coeffs[j];
        }
    }
    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int SystemBath_ForceField::ForceField_epes_SiteExciton(double* V, double* dV, double* ddV, double* R, const int& flag,
                                                       const int& rdim, const int& fdim, const int& itraj,
                                                       const int& isamp) {
    // V
    for (int i = 0; i < FF; ++i) V[i] = Hsys[i];
    int idxV = 0, idxR = 0, Fadd1 = F + 1;
    for (int ibath = 0; ibath < nbath; ++ibath) {
        for (int j = 0; j < Nb; ++j) { V[idxV] += coeffs[j] * R[idxR++]; }
        idxV += Fadd1;
    }
    if (flag < 1) return 0;

    // dV
    if (isamp < 1) {  // it might run multi-trajectories at same time, thus first call flag is not reliable.
        ARRAY_CLEAR(dV, NFF);

        int idxdV = 0;
        for (int ibath = 0; ibath < nbath; ++ibath) {
            for (int j = 0; j < Nb; ++j) {
                for (int i = 0; i < F; ++i)
                    for (int k = 0; k < F; ++k) { dV[idxdV++] = (i == k && i == ibath) ? coeffs[j] : 0.0f; }
            }
        }
    }

    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int SystemBath_ForceField::ForceField_epes_General(double* V, double* dV, double* ddV, double* R, const int& flag,
                                                   const int& rdim, const int& fdim, const int& itraj,
                                                   const int& isamp) {
    // V
    for (int i = 0; i < FF; ++i) V[i] = Hsys[i];
    int idxR = 0;
    for (int ibath = 0; ibath < nbath; ++ibath) {
        for (int j = 0; j < Nb; ++j) {
            int idxQ = ibath * FF;
            for (int i = 0; i < FF; ++i) { V[i] += Q[idxQ++] * coeffs[j] * R[idxR]; }
            idxR++;
        }
    }
    if (flag < 1) return 0;

    // dV
    int idxdV = 0;
    if (isamp < 1) {
        for (int ibath = 0; ibath < nbath; ++ibath) {
            for (int j = 0; j < Nb; ++j) {
                int idxQ = ibath * FF;
                for (int i = 0; i < FF; ++i) { dV[idxdV++] = Q[idxQ++] * coeffs[j]; }
            }
        }
    }

    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int SystemBath_ForceField::reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim,
                                        const int& fdim) {
    num_real* pdHj = dH;
    int Fadd1      = F + 1;
    switch (coup_type) {
        case BathPolicy::SpinBosonCoup: {
            double val = 0.0f;
            // only diagonal
            for (int i = 0, idxd = 0; i < F; ++i, idxd += Fadd1) val += REAL_OF(rho[idxd]) * pdHj[idxd];
            if (coeffs[0] != 0) {
                for (int j = 0; j < N; ++j) fx[j] = val * coeffs[j] / coeffs[0];
            }
            break;
        }
        case BathPolicy::SiteExcitonCoup: {
            // Nad_ForceField::reduce_force(fx, rho, dH, rdim, fdim);
            for (int ibath = 0, idxR = 0; ibath < nbath; ++ibath) {
                double val = 0.0f;
                pdHj       = dH + idxR * FF;
                // only diagonal
                for (int i = 0, idxd = 0; i < F; ++i, idxd += Fadd1) val += REAL_OF(rho[idxd]) * pdHj[idxd];
                if (coeffs[0] != 0) {
                    for (int j = 0; j < Nb; ++j, ++idxR) fx[idxR] = val * coeffs[j] / coeffs[0];
                }
            }
            break;
        }
        case BathPolicy::GeneralCoup: {  // time consuming without optimization!
            Nad_ForceField::reduce_force(fx, rho, dH, rdim, fdim);
            break;
        }
    }
    return 0;
}

int SystemBath_ForceField::get_nbath() { return nbath; }
int SystemBath_ForceField::get_Nb() { return Nb; }
