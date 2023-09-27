#include "solver_sse.h"

#include "../utils/definitions.h"
#include "../utils/nad_utils.h"

using namespace ARRAY_EG;

SSE_Solver::SSE_Solver(Param iparm, Model* pM) : Solver(iparm, pM) {
    plFunction();

    pForceField = dynamic_cast<SystemBath_ForceField*>(pM);  // casted into child class Nad_ForceField
    LOG(WARNING);

    F     = pForceField->get_F();
    N     = pForceField->get_N();
    nbath = pForceField->get_nbath();
    Nb    = pForceField->get_Nb();  // N = Nb * nbath
    L     = pForceField->L;
    LOG(WARNING) << "L = " << L;
    // CHECK_EQ(N, nbath * Nb);

    FF = F * F, NN = N * N, NF = N * F;
    NFF = N * F * F, NNF = N * N * F;
    NNFF = N * N * F * F, FFFF = F * F * F * F;
    NbFF = Nb * FF;

    std::string ini_flag = Param_GetT(std::string, parm, "ini_flag", "#occ");
    ini_type             = elec_init::_dict.at(ini_flag);

    std::string tcf_flag = Param_GetT(std::string, parm, "tcf_flag", "#rho");
    tcf_type             = nad_tcf::_dict.at(tcf_flag);

    nspec = pForceField->nspec();  // get no. of specified dividings
    // get numerical simulation details
    Param_GetV(tend, iparm, 1.0f);
    Param_GetV(dt, iparm, -1.0f);
    Param_GetV(sstep, iparm, 1);
    Param_GetV(ntraj, iparm, 1);
    if (pForceField->Suggest_tend() > 0.0f) tend = pForceField->Suggest_tend();
    if (pForceField->Suggest_dt() > 0.0f) dt = pForceField->Suggest_dt();
    nstep = sstep * (int(tend / (sstep * dt)) + 1);
    nsamp = nstep / sstep + 1;
    CHECK_GT(tend, 0);
    CHECK_GT(dt, 0);
    CHECK_GT(sstep, 0);
    CHECK_GT(ntraj, 0);
    CHECK_GT(nstep, 0);
    CHECK_GT(nsamp, 0);
    LOG(WARNING);


    try {  // @bugs: size must be determined before
        ALLOCATE_PTR_TO_VECTOR(rand1, N);
        ALLOCATE_PTR_TO_VECTOR(rand2, N);
        ALLOCATE_PTR_TO_VECTOR(Eele, F);
        ALLOCATE_PTR_TO_VECTOR(DE, FF);
        ALLOCATE_PTR_TO_VECTOR(w, Nb);
        ALLOCATE_PTR_TO_VECTOR(tanhqwb, Nb);
        ALLOCATE_PTR_TO_VECTOR(xicoeff1, Nb);
        ALLOCATE_PTR_TO_VECTOR(xicoeff2, Nb);
        ALLOCATE_PTR_TO_VECTOR(alpha_pref, Nb);
        ALLOCATE_PTR_TO_VECTOR(CL, L * Nb);
        ALLOCATE_PTR_TO_VECTOR(DEpW, Nb * FF);
        ALLOCATE_PTR_TO_VECTOR(DEmW, Nb * FF);
        ALLOCATE_PTR_TO_VECTOR(coeff_DEpW, Nb * FF);
        ALLOCATE_PTR_TO_VECTOR(coeff_DEmW, Nb * FF);
        ALLOCATE_PTR_TO_VECTOR(Hele, FF);
        ALLOCATE_PTR_TO_VECTOR(Tele, FF);
        ALLOCATE_PTR_TO_VECTOR(T0, FF);
        ALLOCATE_PTR_TO_VECTOR(Xzero, NFF);
        ALLOCATE_PTR_TO_VECTOR(Xtran, NFF);
        ALLOCATE_PTR_TO_VECTOR(Qzero, L * nbath * FF);
        ALLOCATE_PTR_TO_VECTOR(Qtran, L * nbath * FF);

        ALLOCATE_PTR_TO_VECTOR(BL, L * nbath);
        ALLOCATE_PTR_TO_VECTOR(tmpbarr1, nbath * Nb);
        ALLOCATE_PTR_TO_VECTOR(tmpbarr2, nbath * Nb);
        ALLOCATE_PTR_TO_VECTOR(eac, F);
        ALLOCATE_PTR_TO_VECTOR(eac_adia, F);
        ALLOCATE_PTR_TO_VECTOR(eac0, F);
        ALLOCATE_PTR_TO_VECTOR(rho0, FF);
        ALLOCATE_PTR_TO_VECTOR(rhot, FF);
        ALLOCATE_PTR_TO_VECTOR(fact1, F);
        ALLOCATE_PTR_TO_VECTOR(fact2, F);
        ALLOCATE_PTR_TO_VECTOR(fact3, F);
        ALLOCATE_PTR_TO_VECTOR(fact4, F);
        ALLOCATE_PTR_TO_VECTOR(eactmp, F);
        ALLOCATE_PTR_TO_VECTOR(eactmp1, F);
        ALLOCATE_PTR_TO_VECTOR(eactmp2, F);
        ALLOCATE_PTR_TO_VECTOR(eacadd1, F);
        ALLOCATE_PTR_TO_VECTOR(eacadd2, F);
        ALLOCATE_PTR_TO_VECTOR(eactran, F);
        ALLOCATE_PTR_TO_VECTOR(crand1, N);
        ALLOCATE_PTR_TO_VECTOR(crand2, N);
        ALLOCATE_PTR_TO_VECTOR(invexpiEt, F);
        ALLOCATE_PTR_TO_VECTOR(U, FF);
        ALLOCATE_PTR_TO_VECTOR(expiwt, Nb);
        ALLOCATE_PTR_TO_VECTOR(expiDEt, FF);
        ALLOCATE_PTR_TO_VECTOR(invexpiwt, Nb);
        ALLOCATE_PTR_TO_VECTOR(invexpiDEt, FF);
        ALLOCATE_PTR_TO_VECTOR(expiwdht, Nb);
        ALLOCATE_PTR_TO_VECTOR(invexpiwdht, Nb);
        ALLOCATE_PTR_TO_VECTOR(expiwt_now, Nb);
        ALLOCATE_PTR_TO_VECTOR(invexpiwt_now, Nb);
        ALLOCATE_PTR_TO_VECTOR(Xtmpj, FF);
        ALLOCATE_PTR_TO_VECTOR(Xtmp2j, FF);
        ALLOCATE_PTR_TO_VECTOR(xi, N);
        ALLOCATE_PTR_TO_VECTOR(time_ker, Nb * FF);
        ALLOCATE_PTR_TO_VECTOR(mem_arr, (2 * nstep + 1) * FF);

        ALLOCATE_PTR_TO_VECTOR(workr, FF);
        ALLOCATE_PTR_TO_VECTOR(workc, FF);
    } catch (std::bad_alloc& e) { LOG(FATAL) << e.what(); }

    vscale = 1.0f;
    if (dt > 0.1f) {
        LOG(WARNING) << "Large dt will cause unconvergece for RK-4 expansion, so do scaling on dt=0.1";
        vscale = 0.01f / dt;   // remembering [M] ~ [L-1] ~ [T-1]
        dt     = vscale * dt;  // = 0.01
        tend   = vscale * tend;
    }
    for (int i = 0; i < FF; ++i) Hele[i] = pForceField->Hsys[i] / vscale;
    EigenSolve(Eele, Tele, Hele, F);  // solve and save pure electronic probelm
    for (int i = 0, idx = 0, idx2 = 0; i < F; ++i) {
        for (int j = 0; j < F; ++j, ++idx) DE[idx] = Eele[i] - Eele[j];
    }

    switch (ini_type) {
        case elec_init::occ:
        case elec_init::eac:
        case elec_init::rho:
            ARRAY_EYE(T0, F);
            break;
        case elec_init::eig:
        case elec_init::d2a:
            for (int i = 0; i < FF; ++i) T0[i] = Tele[i];
            break;
        case elec_init::rot: {
            CHECK_EQ(F, 2);
            double rotangle = Param_GetV(rotangle, parm, 0.0f) * phys::math::twopi;
            T0[0] = cos(rotangle), T0[1] = -sin(rotangle);
            T0[2] = sin(rotangle), T0[3] = cos(rotangle);
            break;
        }
        case elec_init::def: {
            try {
                std::ifstream ifs(ini_flag);
                double tmp;
                for (int i = 0; i < FF; ++i)
                    if (ifs >> tmp) T0[i] = tmp;
                ifs.close();
            } catch (std::runtime_error& e) { LOG(FATAL) << e.what(); }
            break;
        }
        default:
            LOG(FATAL);
    }
    LOG(WARNING);

    beta = pForceField->mybath->beta * vscale;  // beta directly read in bath, not in systembath !
    for (int j = 0, iEW = 0; j < Nb; ++j) {
        double wj = pForceField->omegas[j] / vscale;
        double nj = 1.0f / (exp(beta * wj) - 1.0f);
        // approach 1
        w[j]           = wj;
        expiwdht[j]    = exp(0.5e0 * phys::math::im * wj * dt);
        invexpiwdht[j] = exp(-0.5e0 * phys::math::im * wj * dt);

        tanhqwb[j]    = tanh(0.25f * wj * beta);
        xicoeff1[j]   = sqrt(0.5f * (nj + 1)) / sqrt(2.0f * wj);
        xicoeff2[j]   = sqrt(0.5f * nj) / sqrt(2.0f * wj);
        alpha_pref[j] = 1.0f / (2.0f * wj);

        for (int i = 0; i < FF; ++i, ++iEW) {  // DE \pm w
            DEpW[iEW]       = (DE[i] + wj);
            DEmW[iEW]       = (DE[i] - wj);
            coeff_DEpW[iEW] = 0.5e0 * alpha_pref[j] * (tanhqwb[j] + 1) / (DE[i] + wj);  // @TODO
            coeff_DEmW[iEW] = 0.5e0 * alpha_pref[j] * (tanhqwb[j] - 1) / (DE[i] - wj);
        }
    }

    for (int i = 0; i < L * Nb; ++i) CL[i] = pForceField->CL[i] / (vscale * sqrt(vscale));
    for (int i = 0; i < L * nbath * FF; ++i) Qzero[i] = pForceField->QL[i];
    for (int j = 0; j < L * nbath; ++j) {
        num_real* Q0j = Qzero + j * FF;
        num_real* Qxj = Qtran + j * FF;
        ARRAY_MATMUL3_TRANS1(Qxj, Tele, Q0j, Tele, F, F, F, F);
    }

    // X = CLB * QLMFF
    for (int ibath = 0, jj = 0, idxQ = 0, idxX = 0; ibath < nbath; ++ibath) {
        // @deprecated later
        for (int j = 0; j < Nb; ++j, ++jj) {
            for (int k = 0; k < FF; ++k, ++idxX) { Xzero[idxX] = pForceField->Xnj[idxX] / (vscale * sqrt(vscale)); }
            // Xtranj = Tele^ * Xzeroj * Tele
            num_real* X0j = Xzero + jj * FF;
            num_real* Xxj = Xtran + jj * FF;
            ARRAY_MATMUL3_TRANS1(Xxj, Tele, X0j, Tele, F, F, F, F);  // note X0j is diagonal
        }
    }

    for (int i = 0; i < Nb; ++i)
        if (w[i] * dt > 0.2) LOG(WARNING) << "Too large timestep for " << i << "-th <w*dt>=" << w[i] * dt;

    LOG(WARNING);

    save = SSE_Solver::name() + "_" + pForceField->tag;
};

SSE_Solver::~SSE_Solver(){};

int SSE_Solver::init(const int& itraj) {
    plFunction();

    pForceField->ForceField_init_elec(rho0, eac0, occ0, F, itraj);
    CHECK_GE(occ0, -1);
    if (occ0 >= 0 && occ0 < F) {
        for (int i = 0; i < F; ++i) eac[i] = (i == occ0) ? phys::math::iu : phys::math::iz;
    } else if (occ0 == -1) {
        for (int i = 0; i < F; ++i) eac[i] = eac0[i];
    } else {
        LOG(FATAL);
    }

    switch (ini_type) {
        case elec_init::occ:
            break;
        default:
            LOG(WARNING);
            ARRAY_MATMUL_TRANS1(workc, T0, eac, F, F, 1);
            for (int i = 0; i < F; ++i) eac[i] = workc[i];
            break;
    }
    rand_gaussian(rand1, N);
    rand_gaussian(rand2, N);
    for (int j = 0; j < N; ++j)
        crand1[j] = rand1[j] + phys::math::im * rand2[j], crand2[j] = rand1[j] - phys::math::im * rand2[j];
    return 0;
}

int SSE_Solver::correlation(const int& isamp, NAD_TCFer& tcfer) {
    if (isamp == 0) rho_eac(rho0, eac, F);
    rho_eac(rhot, eac, F);
    for (int idx = 0, i0 = 0; i0 < FF; ++i0)
        for (int it = 0; it < FF; ++it)
            if (tcfer.ifrecord(i0, it)) tcfer.val[idx++] = rho0[i0] * rhot[it];
    return 0;
}

int SSE_Solver::sampler(const int& isamp, NAD_TCFer& tcfer) {  // default sampler for population dynamics
    correlation(isamp, tcfer);
    ispec = pForceField->ForceField_spec(nullptr, nullptr, nullptr, N, F);
    tcfer.Count(isamp, ispec);
    return 0;
}

// EOM on adiabatic represnetation
int SSE_Solver::time_kernel(num_complex* Ker, const num_real& t) {
    plFunction();

    // time kernel is Nb*FF array
    if (false) {  // slow
        plScope("A");

#define EXP_LF(w, t) ((exp(-phys::math::im * (w) * (t)) - phys::math::iu) / (w))
        for (int j = 0, idx = 0; j < Nb; ++j) {
            for (int i = 0; i < FF; ++i, ++idx) {
                Ker[idx] = (0.5e0 * phys::math::im) * alpha_pref[j] *
                           ((tanhqwb[j] - 1) * EXP_LF(DEmW[idx], t) + (tanhqwb[j] + 1) * EXP_LF(DEpW[idx], t));
            }
        }
#undef EXP_LF
    }

    {
        plScope("D");
        for (int j = 0; j < Nb; ++j)
            expiwt[j] = exp(phys::math::im * w[j] * t), invexpiwt[j] = phys::math::iu / expiwt[j];
        for (int i = 0; i < FF; ++i)
            expiDEt[i] = exp(phys::math::im * DE[i] * t), invexpiDEt[i] = phys::math::iu / expiDEt[i];
        for (int j = 0, idx = 0; j < Nb; ++j) {
            for (int i = 0; i < FF; ++i, ++idx) {
                Ker[idx] = phys::math::im * ((invexpiDEt[i] * expiwt[j] - phys::math::iu) * coeff_DEmW[idx] +
                                             (invexpiDEt[i] * invexpiwt[j] - phys::math::iu) * coeff_DEpW[idx]);
            }
        }
    }
    return 0;
}

int SSE_Solver::call_memory_array_adia(const num_real& dt) {
    for (int istep_dummy = 0; istep_dummy <= 2 * nstep; ++istep_dummy) {
        double t_now = istep_dummy * 0.5e0 * dt;

        num_complex* mem_now = mem_arr + istep_dummy * FF;
        time_kernel(time_ker, t_now);

        for (int i = 0; i < FF; ++i) mem_now[i] = phys::math::iz;
        for (int ibath = 0, jj = 0; ibath < nbath; ++ibath) {
            for (int j = 0; j < Nb; ++j, ++jj) {
                num_real* Xxj          = Xtran + jj * FF;
                num_complex* time_kerj = time_ker + j * FF;

                for (int i = 0; i < FF; ++i) Xtmpj[i] = Xxj[i] * time_kerj[i];
                ARRAY_MATMUL(Xtmp2j, Xxj, Xtmpj, F, F, F);

                for (int i = 0; i < FF; ++i) mem_now[i] += Xtmp2j[i];
            }
        }
        // num_complex* expitw_now = expitw + istep_dummy * Nb;
        // num_complex* invexpitw_now = invexpitw + istep_dummy * Nb;
        // for(int j=0; j<Nb; ++j){
        //      expitw_now[j] = exp(phys::math::im*w[j]*t_now), invexpitw_now[j] = phys::math::iu/expitw_now[j];
        // }
    }
    std::ofstream ofs;
    ofs.open("ker.check");
    for (int istep_dummy = 0, idx = 0; istep_dummy <= 2 * nstep; ++istep_dummy) {
        ofs << FMT(8) << istep_dummy * 0.5e0 * dt;
        for (int i = 0; i < FF; ++i, ++idx) {
            ofs << FMT(8) << REAL_OF(mem_arr[idx]) << FMT(8) << IMAG_OF(mem_arr[idx]);
        }
        ofs << std::endl;
    }
    ofs.close();
    return 0;
}

// EOM on adiabatic represnetation
int SSE_Solver::action_on_wavafunction_adia(num_complex* eacnew, num_complex* eac, const num_real& t) {
    plFunction();

    for (int i = 0; i < Nb; ++i) expiwt[i] = exp(phys::math::im * w[i] * t), invexpiwt[i] = phys::math::iu / expiwt[i];
    for (int i = 0; i < FF; ++i)
        expiDEt[i] = exp(phys::math::im * DE[i] * t), invexpiDEt[i] = phys::math::iu / expiDEt[i];

    for (int i = 0; i < F; ++i) eacnew[i] = Eele[i] * eac[i];

    for (int ibath = 0, jj = 0; ibath < nbath; ++ibath) {
        for (int j = 0, iEW = 0; j < Nb; ++j, ++jj) {  // Hamiltonian from stochastic part
            xi[jj] = xicoeff1[j] * crand1[jj] * expiwt[j] + xicoeff2[j] * crand2[jj] * invexpiwt[j];

            num_real* Xxj = Xtran + jj * FF;

            ARRAY_MATMUL(eacadd1, Xxj, eac, F, F, 1);

            double coswjt = REAL_OF(expiwt[j]), sinwjt = IMAG_OF(expiwt[j]);
            for (int i = 0; i < FF; ++i, ++iEW) {  // note Y = Xtmpj
                Xtmpj[i] = Xxj[i] * alpha_pref[j] * (-0.5e0 * phys::math::im) *
                           ((tanhqwb[j] - 1) * (phys::math::iu - invexpiDEt[i] * expiwt[j]) * coeff_DEmW[iEW] +
                            (tanhqwb[j] + 1) * (phys::math::iu - invexpiDEt[i] * invexpiwt[j]) * coeff_DEpW[iEW]);
            }
            ARRAY_MATMUL(eactmp1, Xtmpj, eac, F, F, 1);
            ARRAY_MATMUL(eacadd2, Xxj, eactmp1, F, F, 1);
            for (int i = 0; i < F; ++i) eacnew[i] += xi[jj] * eacadd1[i] - phys::math::im * eacadd2[i];
        }
    }
    for (int i = 0; i < F; ++i) eacnew[i] *= -phys::math::im;  // update eac
    // LOG(FATAL);
    return 0;
}

// EOM on adiabatic represnetation
int SSE_Solver::action_on_wavafunction_adia_mem(num_complex* eacnew, num_complex* eacold, const int& hstep,
                                                const bool& hupdate) {
    plFunction();

    MapMXc Map_expiwt(expiwt_now, Nb, 1);
    MapMXc Map_invexpiwt(invexpiwt_now, Nb, 1);
    MapMXc Map_expiwdht(expiwdht, Nb, 1);
    if (hupdate) {  // faster
        plScope("a");
        Map_expiwt.array() *= Map_expiwdht.array();
        Map_invexpiwt.array() = Map_expiwt.array().inverse();

        // for (int j = 0; j < Nb; ++j) expiwt_now[j] *= expiwdht[j], invexpiwt_now[j] = phys::math::iu / expiwt[j];
        // //@slow
    }
    for (int i = 0; i < F; ++i) eacnew[i] = Eele[i] * eacold[i];

    {  // @faster, that we first compute Xeff (Xtmpj) on eac
        plScope("b-rev");
        ARRAY_CLEAR(Xtmpj, FF);

        // // @decopistion: X(NFF) = CL(LB) * Q(LM,FF)
        // // xicoeff(B,1) * crand(M,B) * U(B,1) * [CL(LB) * Q(LM,FF)]
        MapMXr Map_Qtran(Qtran, L * nbath, FF);
        MapMXr Map_CL(CL, L, Nb);
        MapMXc Map_BL1(BL, L, nbath);
        MapMXc Map_BL2(BL, 1, L * nbath);
        MapMXr Map_xicoeff1(xicoeff1, Nb, 1);
        MapMXr Map_xicoeff2(xicoeff2, Nb, 1);
        MapMXc Map_Xtmpj(Xtmpj, 1, FF);
        MapMXc Map_crand1(crand1, nbath, Nb);
        MapMXc Map_crand2(crand2, nbath, Nb);
        {
            plScope("tough2");
            Map_BL1 = Map_CL *
                      ((Map_expiwt.array() * Map_xicoeff1.array()).matrix().asDiagonal() * Map_crand1.transpose() +
                       (Map_invexpiwt.array() * Map_xicoeff2.array()).matrix().asDiagonal() * Map_crand2.transpose());
            Map_Xtmpj = Map_BL2 * Map_Qtran;
        }

        // if(false){ // @though slightly faster in nbath=1 case, but might slow for large nbath case!!!
        //     for (int ibath = 0, jj=0, idx=0; ibath < nbath; ++ibath, idx+=Nb) {
        //         for (int j = 0; j < Nb; ++j, ++jj) {
        //             xi[jj] = xicoeff1[j] * crand1[jj] * expiwt_now[j] + xicoeff2[j] * crand2[jj] * invexpiwt_now[j];
        //         }
        //     }
        //     plScope("tough");
        //     ARRAY_MATMUL(Xtmpj, xi, Xtran, 1, N, FF);
        // }

        ARRAY_MATMUL(eacadd1, Xtmpj, eacold, F, F, 1);
        for (int i = 0; i < F; ++i) eacnew[i] += eacadd1[i];
    }
    {
        plScope("c");
        num_complex* mem_now = mem_arr + hstep * FF;
        ARRAY_MATMUL(eacadd2, mem_now, eacold, F, F, 1);  //
        for (int i = 0; i < F; ++i) eacnew[i] -= phys::math::im * eacadd2[i];
    }
    for (int i = 0; i < F; ++i) eacnew[i] *= -phys::math::im;  // update eac

    return 0;
}


int SSE_Solver::action_on_wavafunction(num_complex* eacnew, num_complex* eac, const num_real& t) {
    plFunction();

    for (int i = 0; i < Nb; ++i) expiwt[i] = exp(phys::math::im * w[i] * t), invexpiwt[i] = phys::math::iu / expiwt[i];
    for (int i = 0; i < FF; ++i)
        expiDEt[i] = exp(phys::math::im * DE[i] * t), invexpiDEt[i] = phys::math::iu / expiDEt[i];
    ARRAY_MATMUL(eacnew, Hele, eac, F, F, 1);          // const part Hele * eac
    ARRAY_MATMUL_TRANS1(eactran, Tele, eac, F, F, 1);  // Tele * eac @bugs: ==> Tele^ eac

    for (int ibath = 0, jj = 0; ibath < nbath; ++ibath) {
        for (int j = 0, iEW = 0; j < Nb; ++j, ++jj) {  // Hamiltonian from stochastic part
            xi[jj] = xicoeff1[j] * crand1[jj] * expiwt[j] + xicoeff2[j] * crand2[jj] * invexpiwt[j];

            num_real* X0j = Xzero + jj * FF;
            num_real* Xxj = Xtran + jj * FF;

            // eacadd1  = X0j * eac, note X0j is diagonal
            for (int i = 0, idx = 0, Fadd1 = F + 1; i < F; ++i, idx += Fadd1) eacadd1[i] = X0j[idx] * eac[i];

            double coswjt = REAL_OF(expiwt[j]), sinwjt = IMAG_OF(expiwt[j]);
            for (int i = 0; i < FF; ++i, ++iEW) {  // note Y = Xtmpj
                Xtmpj[i] = Xxj[i] * alpha_pref[j] * (-0.5e0 * phys::math::im) *
                           ((tanhqwb[j] - 1) * (phys::math::iu - invexpiDEt[i] * expiwt[j]) * coeff_DEmW[iEW] +
                            (tanhqwb[j] + 1) * (phys::math::iu - invexpiDEt[i] * invexpiwt[j]) * coeff_DEpW[iEW]);
                // Xtmpj[i] = Xxj[i] * alpha_pref[j] * phys::math::im * invexpiDEt[i] * coeff_DEmW[iEW] *
                // coeff_DEpW[iEW] *
                //            (tanhqwb[j] * (-DE[i] * expiDEt[i] + DE[i] * coswjt + phys::math::im * w[j] * sinwjt) -
                //             (-w[j] * expiDEt[i] + w[j] * coswjt + phys::math::im * DE[i] * sinwjt));
            }
            ARRAY_MATMUL(eactmp1, Xtmpj, eactran, F, F, 1);  // Y*T^*eac
            ARRAY_MATMUL(eactmp2, Tele, eactmp1, F, F, 1);   // T*Y*T^ * eac
            for (int i = 0, idx = 0, Fadd1 = F + 1; i < F; ++i, idx += Fadd1) eacadd2[i] = X0j[idx] * eactmp2[i];
            // add to eacnew
            for (int i = 0; i < F; ++i) eacnew[i] += xi[jj] * eacadd1[i] - phys::math::im * eacadd2[i];
        }
    }
    for (int i = 0; i < F; ++i) eacnew[i] *= -phys::math::im;  // update eac
    // LOG(WARNING);
    // LOG(FATAL);
    return 0;
}


int SSE_Solver::sse_adia(NAD_TCFer& tcfer, const int& N, const int& F)  // gamma-modified mean-field trajectory
{
    plFunction();

    num_real h1 = dt / 6.0f, h2 = dt / 3.0f, h3 = dt / 3.0f, h4 = dt / 6.0f;
    init(itraj);  // initialization procedure

    ARRAY_MATMUL_TRANS1(eac_adia, Tele, eac, F, F, 1);

    int succ = sampler(0, tcfer);

    for (int j = 0; j < Nb; ++j) expiwt_now[j] = phys::math::iu, invexpiwt_now[j] = phys::math::iu;

    for (int istep_dummy = 0; istep_dummy < nstep; ++istep_dummy) {
        int hstep = 2 * istep_dummy;

        if (succ == 0) istep = istep_dummy;  // otherwise froze istep
        {                                    // stochastic electronic propagation
            if (succ == 0) succ = action_on_wavafunction_adia_mem(fact1, eac_adia, hstep, false);  // istep * dt
            for (int i = 0; i < F; ++i) eactmp[i] = eac_adia[i] + 0.5f * dt * fact1[i];
            if (succ == 0) succ = action_on_wavafunction_adia_mem(fact2, eactmp, hstep + 1, true);  // (istep+0.5) * dt
            for (int i = 0; i < F; ++i) eactmp[i] = eac_adia[i] + 0.5f * dt * fact2[i];
            if (succ == 0) succ = action_on_wavafunction_adia_mem(fact3, eactmp, hstep + 1, false);  // (istep+0.5) * dt
            for (int i = 0; i < F; ++i) eactmp[i] = eac_adia[i] + 1.0f * dt * fact3[i];
            if (succ == 0) succ = action_on_wavafunction_adia_mem(fact4, eactmp, hstep + 2, true);  // (istep+1) * dt

            for (int i = 0; i < F; ++i) eac_adia[i] += (fact1[i] * h1 + fact2[i] * h2 + fact3[i] * h3 + fact4[i] * h4);
        }
        if ((istep_dummy + 1) % sstep == 0) {
            isamp = (istep_dummy + 1) / sstep;

            // convert into adiabatic represnetation
            ARRAY_MATMUL(eac, Tele, eac_adia, F, F, 1);

            if (!ARRAY_ISFINITE(eac, F)) break;
            if (succ == 0) succ = sampler(isamp, tcfer);
        }
    }
    return 0;
}

int SSE_Solver::sse(NAD_TCFer& tcfer, const int& N, const int& F)  // gamma-modified mean-field trajectory
{
    plFunction();

    num_real h1 = dt / 6.0f, h2 = dt / 3.0f, h3 = dt / 3.0f, h4 = dt / 6.0f;
    init(itraj);  // initialization procedure

    int succ = sampler(0, tcfer);
    for (int istep_dummy = 0; istep_dummy < nstep; ++istep_dummy) {
        if (succ == 0) istep = istep_dummy;  // otherwise froze istep
        {                                    // stochastic electronic propagation
            if (succ == 0) succ = action_on_wavafunction(fact1, eac, istep * dt);
            for (int i = 0; i < F; ++i) eactmp[i] = eac[i] + 0.5f * dt * fact1[i];
            if (succ == 0) succ = action_on_wavafunction(fact2, eactmp, (istep + 0.5f) * dt);
            for (int i = 0; i < F; ++i) eactmp[i] = eac[i] + 0.5f * dt * fact2[i];
            if (succ == 0) succ = action_on_wavafunction(fact3, eactmp, (istep + 0.5f) * dt);
            for (int i = 0; i < F; ++i) eactmp[i] = eac[i] + 1.0f * dt * fact3[i];
            if (succ == 0) succ = action_on_wavafunction(fact4, eactmp, (istep + 1.0f) * dt);
            for (int i = 0; i < F; ++i) eac[i] += (fact1[i] * h1 + fact2[i] * h2 + fact3[i] * h3 + fact4[i] * h4);
        }
        if ((istep_dummy + 1) % sstep == 0) {
            isamp = (istep_dummy + 1) / sstep;
            if (!ARRAY_ISFINITE(eac, F)) break;
            if (succ == 0) succ = sampler(isamp, tcfer);
        }
    }
    return 0;
}

int SSE_Solver::run_parallel() {  // basic run_impl() for nadtraj

    int nsave = (ntraj > mpi_nprocs) ? ntraj / mpi_nprocs : 1;
    if (nsave > FLAGS_nsave_mpi) nsave = FLAGS_nsave_mpi;
    double sampunit = dt * sstep * iou.time / vscale;

    init(-1);  // get occ0
    NAD_TCFer coll    = NAD_TCFer(tcf_type, nspec, nsamp, N, F, occ0);
    NAD_TCFer collsum = NAD_TCFer(tcf_type, nspec, nsamp, N, F, occ0);
    NAD_TCFer collmpi = NAD_TCFer(tcf_type, nspec, nsamp, N, F, occ0);

    call_memory_array_adia(dt);

    for (int isave = 0; isave < nsave; ++isave) {
        int eachstart = (isave * ntraj) / nsave, eachend = ((isave + 1) * ntraj) / nsave, istart, iend;
        mpi_range(eachstart, eachend, mpi_nprocs, mpi_rank, istart, iend);
        CHECK_EQ(ntraj % (nsave * mpi_nprocs), 0);
        LOG(INFO) << "During [" << eachstart << ", " << eachend << "), "
                  << "mpi-" << mpi_rank << " cycle in [" << istart << ", " << iend << ")";

        MPI_Barrier(MPI_COMM_WORLD);
        for (int icycle = istart; icycle < iend; ++icycle) {
            itraj = icycle;
            coll.Clear();
            try {
                sse_adia(coll, N, F);
            } catch (std::runtime_error& e) {  // if some error, output currect results
                LOG(WARNING) << "runtime_error cause breakdown for traj=" << itraj;
            }
            collsum.Amount(coll);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        collmpi.MPIAmount(collsum);

        // MPI_Reduce(Collsum, Collmpi, ncoll, MPI::DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
        // MPI_Reduce(Statsum, Statmpi, nsamp, MPI::INT, MPI_SUM, 0, MPI_COMM_WORLD);

        if (mpi_rank == 0) collmpi.report(save + "-cache" + std::to_string(isave), sampunit);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    collmpi.MPIAmount(collsum);
    if (mpi_rank == 0) collmpi.report(save, sampunit);
    collsum.report(save + "-mpi" + std::to_string(mpi_rank), sampunit);

    return 0;
}
