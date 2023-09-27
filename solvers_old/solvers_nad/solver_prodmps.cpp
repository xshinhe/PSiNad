#include "solver_prodmps.h"

#include "../models/nad_forcefield/manysite_models.h"
#include "../solvers_nad/nadtraj.h"
#include "../utils/nad_utils.h"

using namespace ARRAY_EG;

ProductMPS_Solver::ProductMPS_Solver(Param iparm, Model* pM) : NadTraj_Solver(iparm, pM) {
    pForceField = dynamic_cast<ManySite_ForceField*>(pM);

    Param_GetV(scale, parm, 1.0f);

    std::string flag = Param_GetT(std::string, parm, "flag", "#twa");
    type             = ProductMPSPolicy::_dict.at(flag);

    // get dimensions
    F   = pForceField->get_F();
    FF  = F * F;
    M   = pForceField->get_M();
    MF  = M * F;
    MFF = M * FF;

    CHECK_GE(F, 2);
    CHECK_GE(M, 1);

    eom_type = elec_eom::rho;
    tcf_type = nad_tcf::redtcf;

    std::ofstream ofs;
    ofs.open("tcf0");  // sigmax
    ofs << FMT(8) << 0 << FMT(8) << 1 << std::endl << FMT(8) << 1 << FMT(8) << 0 << std::endl;
    ofs.close();

    ofs.open("tcft");  // sigmax
    ofs << FMT(8) << 0 << FMT(8) << 1 << std::endl << FMT(8) << 1 << FMT(8) << 0 << std::endl;
    ofs.close();

    tcf_weight = 0.5e0;

    ALLOCATE_PTR_TO_VECTOR(rhos, MFF);
    ALLOCATE_PTR_TO_VECTOR(rhos0, MFF);
    ALLOCATE_PTR_TO_VECTOR(rhost, MFF);
    ALLOCATE_PTR_TO_VECTOR(Hs, MFF);

    std::string suffix = std::string("mstraj") + flag;
    if (type == ProductMPSPolicy::MMF) suffix += std::to_string(scale);

    save = name() + suffix + "_" + pForceField->tag;
}

ProductMPS_Solver::~ProductMPS_Solver(){};

int ProductMPS_Solver::init(const int& itraj) {
    int randi;
    for (int is = 0; is < M; ++is) {
        num_complex* rhoi = rhos + is * FF;
        // NadTraj_Solver::init(itraj);
        // rho_eac(rhoi, eac, F);

        switch (type) {
            case ProductMPSPolicy::TWA: {
                num_real x = 1;              // {only 1}
                rand_catalog(&randi);        // {0,1}
                num_real y = 1 - 2 * randi;  // {-1,1}
                rand_catalog(&randi);        // {0,1}
                num_real z = 1 - 2 * randi;  // {-1,1}
                rhoi[0]    = z * phys::math::iu;
                rhoi[1]    = x * phys::math::iu + y * phys::math::im;
                rhoi[2]    = x * phys::math::iu - y * phys::math::im;
                rhoi[3]    = -z * phys::math::iu;
                break;
            }
            case ProductMPSPolicy::MMF: {
                rand_catalog(&randi);
                if (randi == 0) {
                    rhoi[0] = phys::math::iu, rhoi[3] = phys::math::iz;
                } else {
                    rhoi[0] = phys::math::iz, rhoi[3] = phys::math::iu;
                }
                num_real randu;
                rand_uniform(&randu);
                randu *= phys::math::twopi;

                rhoi[1] = scale * phys::math::sqrthalf * (cos(randu) + phys::math::im * sin(randu));
                rhoi[2] = CONJ_OF(rhoi[1]);

                break;
            }
            case ProductMPSPolicy::CMM: {
                num_real gamma0 = GAMMA_WIGNER(F);
                num_real xi0    = 1.0f;
                num_real totact = (1 + F * gamma0) / xi0;

                samp_mvc_sphere(mvc, 2.0f * totact, F);
                eac_mvc(eac0, mvc, F);
                rho_eac(rhoi, eac0, F);
                for (int i = 0; i < F; ++i) rhoi[i * (F + 1)] -= gamma0;

                // LOG(WARNING) << rhoi[0] + rhoi[3];
                // ARRAY_SHOW(rhoi, F, F);
                break;
            }
        }
    }
    return 0;
}

int ProductMPS_Solver::traj(NAD_TCFer& tcfer, const int& N, const int& F) {
    init(itraj);
    int succ = sampler(0, tcfer);
    for (int istep_dummy = 0; istep_dummy < nstep; ++istep_dummy) {
        if (succ == 0) istep = istep_dummy;  // otherwise froze istep
        pForceField->ForceField_heff(Hs, rhos, M, F);

        // ARRAY_SHOW(Hs, M, FF);
        // LOG(FATAL);

        for (int is = 0; is < M; ++is) {
            num_complex* Heff = Hs + is * FF;
            num_complex* rhoi = rhos + is * FF;

            EigenSolve(L, S, Heff, F);
            for (int i = 0; i < F; ++i) workc[i] = cos(L[i] * dt) - phys::math::im * sin(L[i] * dt);
            ARRAY_MATMUL3_TRANS2(U, S, workc, S, F, F, 0, F);

            update_rho(rhoi, U, N, F, workr, workc);  // N=1 fixed with nomeaning
        }

        if ((istep_dummy + 1) % sstep == 0) {
            plScope("sample");
            isamp = (istep_dummy + 1) / sstep;
            if (check_break(succ) != 0) break;
            if (succ == 0) succ = sampler(isamp, tcfer);
        }
    }
    final(itraj);
    return 0;
}

int ProductMPS_Solver::correlation(const int& isamp, NAD_TCFer& tcfer) {
    num_complex paulix[4] = {0.0e0, 1.0e0, 1.0e0, 0.0e0};

    switch (type) {
        case ProductMPSPolicy::TWA: {
            if (isamp == 0) {
                for (int is = 0, idx = 0; is < M; ++is) {
                    rhos0[idx++] = phys::math::iz;
                    rhos0[idx++] = phys::math::iu;
                    rhos0[idx++] = phys::math::iu;
                    rhos0[idx++] = phys::math::iz;
                }
                // for (int i = 0; i < MFF; ++i) rhos0[i] = rhos[i];
            }
            for (int i = 0; i < MFF; ++i) rhost[i] = rhos[i];
            if (tcfer.tcf_reduced) {
                tcfer.val[0] = phys::math::iz;
                for (int is = 0; is < M; ++is) {
                    num_complex* rho0 = rhos0 + is * FF;
                    num_complex* rhot = rhost + is * FF;
                    tcfer.val[0] += (1.0e0 / M) * 0.5e0 * ARRAY_TRACE2(rho0, rhot, F, F);
                }
            }
            break;
        }
        case ProductMPSPolicy::MMF:
        case ProductMPSPolicy::CMM: {
            if (isamp == 0) {
                for (int i = 0; i < MFF; ++i) rhos0[i] = rhos[i];
            }
            for (int i = 0; i < MFF; ++i) rhost[i] = rhos[i];
            if (tcfer.tcf_reduced) {
                tcfer.val[0] = phys::math::iz;
                for (int is = 0; is < M; ++is) {
                    num_complex* rho0 = rhos0 + is * FF;
                    num_complex* rhot = rhost + is * FF;
                    tcfer.val[0] +=
                        (1.0e0 / M) * num_real(F) * ARRAY_TRACE2(rho0, paulix, F, F) * ARRAY_TRACE2(rhot, paulix, F, F);
                }
            }
            break;
        }
    }

    return 0;
}
