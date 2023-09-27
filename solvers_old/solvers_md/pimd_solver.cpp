#include "pimd_solver.h"

using namespace ARRAY_EG;

PIMDTraj_Solver::PIMDTraj_Solver(Param iparm, Model* pM) : Traj_Solver(iparm, pM) {
    // get dimensions
    P    = Param_GetT(int, parm, "nbead", 1);
    PP   = P * P;
    PN   = P * N;
    PNPN = PN * PN;

    // pimd_flag, pimd_type

    std::string pimd_flag  = Param_GetT(std::string, parm, "pimd_flag", "#stag");
    std::string integrator = Param_GetT(std::string, parm, "integrator", "BAOAB");
    std::string traj_rst   = Param_GetT(std::string, parm, "traj_rst", "both");
    pimd_type              = PIMDTransformPolicy::_dict.at(pimd_flag);
    integ_type             = PIMDIntegratorPolicy::_dict.at(integrator);
    rst_type               = PIMDRestartPolicy::_dict.at(traj_rst);
    pThermo->init_alloc(PN);

    num_real beta = (pThermo->beta);
    bf2           = P / (beta * beta);
    betap         = beta / P;

    try {  // allocate arrays
        ALLOCATE_PTR_TO_VECTOR(nrs, PN);
        ALLOCATE_PTR_TO_VECTOR(nks, PN);
        ALLOCATE_PTR_TO_VECTOR(nps, PN);
        ALLOCATE_PTR_TO_VECTOR(nms, PN);
        ALLOCATE_PTR_TO_VECTOR(nfs, PN);
        ALLOCATE_PTR_TO_VECTOR(nfks, PN);
        ALLOCATE_PTR_TO_VECTOR(fk_spring, PN);

        ALLOCATE_PTR_TO_VECTOR(vpeses, P);
        ALLOCATE_PTR_TO_VECTOR(grads, PN);
        // ALLOCATE_PTR_TO_VECTOR(hesses, PN * N); // too large

        ALLOCATE_PTR_TO_VECTOR(masswgt, P);
        ALLOCATE_PTR_TO_VECTOR(bfwgt, P);
        ALLOCATE_PTR_TO_VECTOR(D2, PP);
        ALLOCATE_PTR_TO_VECTOR(Tran, PP);

        rc = new num_real[3], vc = new num_real[3], ac = new num_real[3];
        Lc = new num_real[3], Wc = new num_real[3], Mc = new num_real[3], Ac = new num_real[3];
        vectmp = new num_real[3], Mattmp = new num_real[9], Ic = new num_real[9], invIc = new num_real[9];
    } catch (std::bad_alloc& e) { LOG(FATAL) << e.what(); }

    save = name() + "_" + pForceField->tag;
}

PIMDTraj_Solver::~PIMDTraj_Solver() {
    delete[] rc, delete[] vc, delete[] ac;
    delete[] Lc, delete[] Wc, delete[] Mc, delete[] Ac;
    delete[] vectmp, delete[] Mattmp, delete[] Ic, delete[] invIc;
    delete pThermo;
};

int PIMDTraj_Solver::ff_calc1(const int& level) {  // multi-forcefield calculation at a fixed nr, np
    plFunction();
    all_K2X();  // From nks to nrs (\xi to r)

    // calc centroid r and saved in traditional md variables
    ARRAY_CLEAR(nr, N);
    for (int i = 0, idx = 0; i < P; ++i)
        for (int j = 0; j < N; ++j, ++idx) nr[j] += nrs[idx];
    for (int j = 0; j < N; ++j) nr[j] /= P;

    for (int i = 0, idx = 0; i < P; ++i, idx += N) {
        num_real* nri   = nrs + idx;
        num_real* npi   = nps + idx;
        num_real* vpesi = vpeses + i;
        num_real* gradi = grads + idx;
        num_real* nfi   = nfs + idx;
        num_real* hessi = hesses + 0;  ///< @todo don't use hessian

        pForceField->ForceField_npes(vpesi, gradi, hessi, nri, npi, level, N, itraj, istep);

        for (int j = 0; j < N; ++j) nfi[j] = gradi[j];  // a copy
    }
    all_FX2FK();  ///< transfrom f to fk
    return 0;
}

int PIMDTraj_Solver::spring_force() {
    plFunction();
    ARRAY_CLEAR(fk_spring, PN);
    switch (pimd_type) {
        case PIMDTransformPolicy::Primitive:
            for (int i = 0, idx = 0, idxp = N, idxm = PN - N; i < P; ++i) {
                if (i == 1) idxm = 0; // idxm = (idx - N) mod PN
                if (i == P - 1) idxp = 0; // idxp = (idx + N) mod PN
                for (int j = 0; j < N; ++j, ++idx, ++idxp, ++idxm) {
                    fk_spring[idx] = bf2 * bfwgt[i] * nms[idx] * (nks[idx] - nks[idxp] - nks[idxm]);
                }
            }
            break;
        case PIMDTransformPolicy::Staging:
            for (int i = 0, idx = 0; i < P; ++i) {
                for (int j = 0; j < N; ++j, ++idx) { fk_spring[idx] = bf2 * bfwgt[i] * nms[idx] * nks[idx]; }
            }
            break;
        case PIMDTransformPolicy::NormalMode: {
            for (int i = 0, idx = 0; i < P; ++i) {
                for (int j = 0; j < N; ++j, ++idx) { fk_spring[idx] = bf2 * bfwgt[i] * nms[idx] * nks[idx]; }
            }
            break;
        }
        case PIMDTransformPolicy::TEST: { // @TODO
            for (int i = 0, idx = 0; i < P; ++i) {
                for (int j = 0; j < N; ++j, ++idx) { fk_spring[idx] = bf2 * bfwgt[i] * nms[idx] * nks[idx]; }
            }
            break;
        }
        default:
            LOG(FATAL);
    }
    return 0;
}

int PIMDTraj_Solver::update_p(const num_real& dt_in) {  ///< overload with different size
    plFunction();
    for (int i = 0; i < PN; ++i) nps[i] -= nfks[i] * dt_in;
    return 0;
}

int PIMDTraj_Solver::update_p_harm(const num_real& dt_in) {  // for harmonic bead potential
    plFunction();
    spring_force();
    for (int i = 0, idx = 0; i < P; ++i) {
        for (int j = 0; j < N; ++j, ++idx) {
            /**
             * you can skip the first bead (P=0) in #stag & #norm, but adapted to #prim & #test flag.
             * for example, pimd for boson can only used in #prim formalism. so don't skip P=0 bead.
             * 
             * REPLY from WSH: fix the problem through a new function spring_force().
             */
            if (bfwgt[i] == 0) continue;
            nps[idx] -= fk_spring[idx] * dt_in;
        }
    }
    return 0;
}

int PIMDTraj_Solver::update_r(const num_real& dt_in) {  ///< overload with different size
    plFunction();
    for (int i = 0; i < PN; ++i) nks[i] += nps[i] / nms[i] * dt_in;
    return 0;
}

int PIMDTraj_Solver::caylay_update_half(const num_real& dt_in) {
    plFunction();
    for (int i = 0, idx = 0; i < P; ++i) {
        for (int j = 0; j < N; ++j, ++idx) {
            num_real nk_tmp_j =
                sqrt(1.0f / (4.0f + bf2 * bfwgt[i] * dt_in * dt_in)) * (2 * nks[idx] + nps[idx] / nms[idx] * dt_in);
            nps[idx] = sqrt(1.0f / (4.0f + bf2 * bfwgt[i] * dt_in * dt_in)) *
                       (-nms[idx] * bf2 * bfwgt[i] * nks[idx] * dt_in + 2 * nps[idx]);
            nks[idx] = nk_tmp_j;
        }
    }
    return 0;
}

int PIMDTraj_Solver::update_thermo(const num_real& dt_in) {  ///< overload with different size
    plFunction();
    num_real gamma_ref = sqrt(bf2);
    for (int i = 0, iN = 0; i < P; ++i, iN += N) {
        num_real *nki = nks + iN, *npi = nps + iN, *nmi = nms + iN;
        num_real gamma_opt = gamma_ref * (bfwgt[i] == 0 ? 1 : sqrt(bfwgt[i]));  // norm 0-f should use last-f
        pThermo->evolve(nki, npi, nmi, dt_in, N, iN, gamma_opt);
    }
    // pThermo->evolve(nks, nps, nms, dt_in, PN, 0);
    return 0;
}

int PIMDTraj_Solver::BAOAB(int& succ, int step) {
    plFunction();
    if (succ == 0) succ = update_p(halfdt);
    if (succ == 0) succ = update_p_harm(halfdt);
    if (succ == 0) succ = update_r(halfdt);
    if (succ == 0 && pThermo->dothermo(step)) succ = update_thermo(dt);
    if (succ == 0) update_r(halfdt);
    if (succ == 0) succ = ff_calc1(level);
    if (succ == 0) update_p(halfdt);
    if (succ == 0) succ = update_p_harm(halfdt);
    return succ;
}

int PIMDTraj_Solver::BCOCB(int& succ, int step) {
    plFunction();

    if (succ == 0) succ = update_p(halfdt);
    if (succ == 0) succ = caylay_update_half(dt);
    if (succ == 0 && pThermo->dothermo(step)) succ = update_thermo(dt);
    if (succ == 0) succ = caylay_update_half(dt);
    if (succ == 0) succ = ff_calc1(level);
    if (succ == 0) update_p(halfdt);
    return succ;
}

int PIMDTraj_Solver::traj(TCFnucl& tcfer, const int& PN) { return traj_Middle(tcfer, PN); };

int PIMDTraj_Solver::traj_Middle(TCFnucl& tcfer, const int& PN) {
    init(itraj);
    all_X2K();

    istep    = 0, ff_calc1(level);  // ff_calc1 need istep < 1 for some forcefield
    int succ = sampler(0, tcfer);
    for (int istep_dummy = 0; istep_dummy < nstep; ++istep_dummy) {
        if (succ == 0) istep = istep_dummy;  // otherwise froze istep
        if (succ == 0) {
            plFunction();
            switch (integ_type) {
                case PIMDIntegratorPolicy::BAOAB:
                    succ = BAOAB(succ, istep_dummy);
                    break;
                case PIMDIntegratorPolicy::BCOCB:
                    succ = BCOCB(succ, istep_dummy);
                    break;
                default:
                    LOG(FATAL);
            }
        }
        if (succ == 0) traj_property(dt);

        if ((istep_dummy + 1) % sstep == 0) {
            if (check_break(succ) != 0) break;  // check only when do sampling
            // sampling
            isamp = (istep_dummy + 1) / sstep;
            if (succ == 0) succ = sampler(isamp, tcfer);
        }
    }
    final(itraj);
    return 0;
}

int PIMDTraj_Solver::traj_property(const num_real& dt) {
    Khere = 0.0f, Vhere = 0.0f;
    for (int j = 0; j < PN; ++j) Khere += 0.5f * nps[j] * nps[j] / nms[j];
    for (int j = 0; j < P; ++j) Vhere += vpeses[j];
    for (int i = 0, idx = 0; i < P; ++i) {
        for (int j = 0; j < N; ++j, ++idx) {
            if (i == 0) continue;  // skip the first bead
            Vhere += 0.5 * bf2 * nms[idx] * nks[idx] * nks[idx];
        }
    }
    Hhere = Khere + Vhere;
    Lhere = Khere - Vhere;
    Shere += Lhere * dt;
    return 0;
}

int PIMDTraj_Solver::sampler(const int& isamp, TCFnucl& tcfer) {
    estimator(isamp, tcfer);
    ispec = pForceField->ForceField_spec(nr, np, nm, N);
    // tcfer.Count(isamp, ispec);
    return 0;
}

int PIMDTraj_Solver::estimator(const int& isamp, TCFnucl& tcfer) {
    num_real esti_V, esti_Kprim, esti_Kvir, esti_Kvir_c;

    // V
    esti_V = 0.0f;
    for (int i = 0; i < P; ++i) esti_V += vpeses[i];
    esti_V /= P;

    esti_Kprim  = 0.5f * N / betap;
    esti_Kvir   = 0.0f;
    esti_Kvir_c = 0.5f * N / pThermo->beta;
    for (int i = 0, idx = 0; i < P; ++i) {
        for (int j = 0; j < N; ++j, ++idx) {
            esti_Kprim -= 0.5f * bf2 * bfwgt[i] * nms[idx] * nks[idx] * nks[idx];  //@adapted to #prim
            esti_Kvir += 0.5f * nrs[idx] * nfs[idx] / P;
            esti_Kvir_c += 0.5f * (nrs[idx] - nr[j]) * nfs[idx] / P;  //@ bugs for TEST flag! why?
        }
    }

    if (isamp == 0) {
        ofs_ENER << FMT(8) << "flag" << FMT(8) << "time"  //
                 << FMT(8) << "N"                         //
                 << FMT(8) << "P"                         //
                 << FMT(8) << "beta"                      //
                 << FMT(8) << "esti_V"                    //
                 << FMT(8) << "esti_Kp"                   //
                 << FMT(8) << "esti_Kv"                   //
                 << FMT(8) << "esti_Kv_c"                 //
                 << std::endl;
    }

    num_real sampunit = dt * sstep * iou.time;
    ofs_ENER << FMT(8) << itraj << FMT(8) << isamp * sampunit  //
             << FMT(8) << N                                    //
             << FMT(8) << P                                    //
             << FMT(8) << pThermo->beta                        //
             << FMT(8) << esti_V                               //
             << FMT(8) << esti_Kprim                           //
             << FMT(8) << esti_Kvir                            //
             << FMT(8) << esti_Kvir_c                          //
             << std::endl;

    for (int idof = 0; idof < PN; ++idof) ofs_TRAJ << FMT(8) << nrs[idof];
    for (int idof = 0; idof < PN; ++idof) ofs_TRAJ << FMT(8) << nfs[idof];
    ofs_TRAJ << std::endl;

    return 0;
}

int PIMDTraj_Solver::run_impl() {
    num_real sampunit = dt * sstep * iou.time;
    init(-1);  // get occ0
    TCFnucl coll;
    {
        itraj = 0;  // global variables
        try {
            traj(coll, PN);
        } catch (std::runtime_error& e) {  // if some error, output currect results
            LOG(WARNING) << "runtime_error cause breakdown for traj=" << itraj;
        }
    }
    final(-1);
    return 0;
}

int PIMDTraj_Solver::run_parallel() {
    int nsave         = (ntraj > mpi_nprocs) ? ntraj / mpi_nprocs : 1;
    num_real sampunit = dt * sstep * iou.time;
    init(-1);  // get occ0
    TCFnucl coll;

    for (int isave = 0; isave < nsave; ++isave) {
        int eachstart = (isave * ntraj) / nsave, eachend = ((isave + 1) * ntraj) / nsave, istart, iend;
        mpi_range(eachstart, eachend, mpi_nprocs, mpi_rank, istart, iend);
        CHECK_EQ(ntraj % (nsave * mpi_nprocs), 0);
        LOG(INFO) << "During [" << eachstart << ", " << eachend << "), "
                  << "mpi-" << mpi_rank << " cycle in [" << istart << ", " << iend << ")";

        MPI_Barrier(MPI_COMM_WORLD);
        for (int icycle = istart; icycle < iend; ++icycle) {
            itraj = icycle;
            try {
                traj(coll, PN);
            } catch (std::runtime_error& e) {  // if some error, output currect results
                LOG(WARNING) << "runtime_error cause breakdown for traj=" << itraj;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    final(-1);
    return 0;
}

int PIMDTraj_Solver::init(const int& itraj) {  // initialize masswgt & Tran
    init_ofs(itraj);

    if (itraj < 0) {
        switch (pimd_type) {
            case PIMDTransformPolicy::Primitive:
                for (int i = 0; i < P; ++i) masswgt[i] = 1.0f;
                for (int i = 0; i < P; ++i) bfwgt[i] = 1.0f;
                break;
            case PIMDTransformPolicy::Staging:
                for (int i = 0; i < P; ++i) masswgt[i] = (i == 0) ? 1.0f : (i + 1.0f) / i;
                for (int i = 0; i < P; ++i) bfwgt[i] = (i == 0) ? 0.0f : 1.0f;  //< P=0 bead is free
                break;
            case PIMDTransformPolicy::NormalMode: {
                for (int i = 0; i < P; ++i) masswgt[i] = 1.0f;

                CHECK_EQ(P % 2, 0);
                int idxM = 0, idxT = 0;
                bfwgt[idxM++] = 0.0f;  //< P=0 bead is free
                for (int i = 0; i < P; ++i) Tran[idxT++] = sqrt(1.0f / P);
                for (int i = 1; i < P / 2; ++i) {
                    bfwgt[idxM++] = 2.0f * (1.0f - cos(phys::math::twopi * i / P));
                    for (int j = 0; j < P; ++j) Tran[idxT++] = cos(phys::math::twopi * i * j / P) * sqrt(2.0f / P);
                }
                bfwgt[idxM++] = 4.0f;
                for (int j = 0; j < P / 2; ++j) {
                    Tran[idxT++] = -sqrt(1.0f / P);
                    Tran[idxT++] = sqrt(1.0f / P);
                }
                for (int i = P / 2 + 1; i < P; ++i) {
                    bfwgt[idxM++] = 2.0f * (1.0f - cos(phys::math::twopi * i / P));
                    for (int j = 0; j < P; ++j) Tran[idxT++] = sin(phys::math::twopi * i * j / P) * sqrt(2.0f / P);
                }
                break;
            }
            case PIMDTransformPolicy::TEST: {
                CHECK_EQ(P % 2, 1);  // need test for odd bead (sum_ij D2ij !=0, why?)

                for (int i = 0; i < P; ++i) masswgt[i] = 1.0f;
                num_real pi2 = phys::math::pi * phys::math::pi;
                for (int i = 0, ik = 0; i < P; ++i) {
                    for (int k = 0; k < P; ++k, ++ik) {
                        if (i == k) {
                            D2[ik] = pi2 / 3.0f * (P * P - 1.0f) / (P * P);
                        } else {
                            int sign   = ((i - k) % 2 == 0) ? 1 : -1;
                            num_real s = std::sin((i - k) * phys::math::pi / P);
                            D2[ik]     = 2 * sign * pi2 * std::cos((i - k) * phys::math::pi / P) / (P * P * s * s);
                        }
                    }
                }
                EigenSolve(bfwgt, Tran, D2, P);
                bfwgt[0] = 0.0f;  // for safety
                break;
            }
            default:
                LOG(FATAL);
        }

        for (int i = 0, ibead_j = 0; i < P; ++i) {
            pForceField->ForceField_init(nr0, np0, nm, N, itraj);
            for (int j = 0; j < N; ++j, ++ibead_j) {
                nrs[ibead_j] = nr0[j], nps[ibead_j] = np0[j], nms[ibead_j] = masswgt[i] * nm[j];
            }
        }
    } else {
        for (int i = 0, ibead_j = 0; i < P; ++i) {
            pForceField->ForceField_init(nr0, np0, nm, N, itraj);
            for (int j = 0; j < N; ++j, ++ibead_j) {
                nrs[ibead_j] = nr0[j], nps[ibead_j] = np0[j], nms[ibead_j] = masswgt[i] * nm[j];
            }
        }
        if (FLAGS_r != "") rst_read(itraj);
    }
    return 0;
}

int PIMDTraj_Solver::rst_read(const int& traj_in) {
    switch  (rst_type) {
        case PIMDRestartPolicy::none:
        case PIMDRestartPolicy::out:
            hload(pContext_in, "position", nrs, N * P);
            hload(pContext_in, "momentum", nps, N * P);
            break;
        case PIMDRestartPolicy::in:
        case PIMDRestartPolicy::both:
            hload_extend(pContext_in, "position", nrs, N * P, ntraj, traj_in);
            hload_extend(pContext_in, "momentum", nps, N * P, ntraj, traj_in);
            break;
        default:
            LOG(FATAL);
    }
    return 0;
}

int PIMDTraj_Solver::rst_output(const int& traj_in) {
    switch  (rst_type) {
        case PIMDRestartPolicy::none:
        case PIMDRestartPolicy::in:
            hdump(pContext_out, "position", nrs, N * P);
            hdump(pContext_out, "momentum", nps, N * P);
            break;
        case PIMDRestartPolicy::out:
        case PIMDRestartPolicy::both:
            LOG(INFO) << "itraj:" << traj_in << " ntraj:" << ntraj << std::endl;
            hdump_extend(pContext_out, "position", nrs, N * P, ntraj, traj_in);
            hdump_extend(pContext_out, "momentum", nps, N * P, ntraj, traj_in);
            break;
        default:
            LOG(FATAL);
    }
    return 0;
}

int PIMDTraj_Solver::all_X2K() {
    switch (pimd_type) {
        case PIMDTransformPolicy::Primitive:
            for (int ibead_j = 0; ibead_j < PN; ++ibead_j) nks[ibead_j] = nrs[ibead_j];
            break;
        case PIMDTransformPolicy::Staging:
            // first bead
            for (int first_j = 0; first_j < N; ++first_j) nks[first_j] = nrs[first_j];
            if (P == 1) return 0;

            // last bead
            for (int first_j = 0, last_j = PN - N; first_j < N; ++first_j, ++last_j) {
                nks[last_j] = nrs[last_j] - nks[first_j];
            }

            // other beads (i-th bead), swap backward
            for (int i = P - 2, ibead_j = PN - N - 1, ipbead_j = PN - 1; i > 0; --i) {
                for (int first_j = N - 1; first_j > -1; --first_j, --ibead_j, --ipbead_j) {
                    nks[ibead_j] = nrs[ibead_j] - (nrs[ipbead_j] * i + nrs[first_j]) / (i + 1);
                }
            }
            break;
        case PIMDTransformPolicy::NormalMode:
        case PIMDTransformPolicy::TEST: {
            // num_real sqrt1overP = sqrt(1.0f / P);
            ARRAY_MATMUL(nks, Tran, nrs, P, P, N);
            // for (int i = 0; i < P; ++i) K_arr[i] *= sqrt1overP;
            break;
        }
        default:
            LOG(FATAL);
    }
    return 0;
}

int PIMDTraj_Solver::all_K2X() {  // @opted
    switch (pimd_type) {
        case PIMDTransformPolicy::Primitive:
            for (int ibead_j = 0; ibead_j < PN; ++ibead_j) nrs[ibead_j] = nks[ibead_j];
            break;
        case PIMDTransformPolicy::Staging:  // @fix bugs for P=2
            // first bead
            for (int first_j = 0; first_j < N; ++first_j) nrs[first_j] = nks[first_j];
            if (P == 1) return 0;

            // last bead
            for (int first_j = 0, last_j = PN - N; first_j < N; ++first_j, ++last_j) {
                nrs[last_j] = nks[last_j] + nks[first_j];
            }

            // other beads (i-th bead), swap backward
            for (int i = P - 2, ibead_j = PN - N - 1, ipbead_j = PN - 1; i > 0; --i) {
                for (int first_j = N - 1; first_j > -1; --first_j, --ibead_j, --ipbead_j) {
                    nrs[ibead_j] = nks[ibead_j] + (nrs[ipbead_j] * i + nks[first_j]) / (i + 1);
                }
            }
            break;
        case PIMDTransformPolicy::NormalMode:
        case PIMDTransformPolicy::TEST:
            ARRAY_MATMUL_TRANS1(nrs, Tran, nks, P, P, N);
            break;
        default:
            LOG(FATAL);
    }
    return 0;
}

int PIMDTraj_Solver::all_FX2FK() {
    switch (pimd_type) {
        case PIMDTransformPolicy::Primitive:
            for (int ibead_j = 0; ibead_j < PN; ++ibead_j) nfks[ibead_j] = nfs[ibead_j] / P;
            break;
        case PIMDTransformPolicy::Staging:
            // first bead (i = 0)
            for (int first_j = 0; first_j < N; ++first_j) {
                nfks[first_j] = 0;
                for (int i = 0, ibead_j = first_j; i < P; ++i, ibead_j += N) nfks[first_j] += nfs[ibead_j] / P;
            }
            // other beads (i = 1,...,P-1; im = i - 1)
            for (int i = 1, ibead_j = N, imbead_j = 0; i < P; ++i) {
                for (int j = 0; j < N; ++j, ++ibead_j, ++imbead_j) {
                    nfks[ibead_j] = nfs[ibead_j] / P + nfks[imbead_j] * (i - 1) / num_real(i);
                }
            }
            break;
        case PIMDTransformPolicy::NormalMode:
        case PIMDTransformPolicy::TEST:
            ARRAY_MATMUL(nfks, Tran, nfs, P, P, N);
            for (int ibead_j = 0; ibead_j < PN; ++ibead_j) nfks[ibead_j] /= P;
            break;
        default:
            LOG(FATAL);
    }
    return 0;
}

int PIMDTraj_Solver::rot_trans_corr(int Natom, num_real* m_in, num_real* x_in, num_real* p_in, num_real* F_in,
                                    bool cal_force) {
    Mt = 0.0f;
    ARRAY_CLEAR(rc, 3);
    // calculate total mass and rc
    for (int i = 0, idof = 0; i < Natom; ++i) {
        Mt += m_in[idof];
        for (int j = 0; j < 3; ++j, ++idof) rc[j] += x_in[idof] * m_in[idof];
    }
    for (int i = 0; i < 3; ++i) rc[i] /= Mt;
    // calculate relative position
    ARRAY_CLEAR(nr, N);
    for (int i = 0, idof = 0; i < Natom; ++i) {
        for (int j = 0; j < 3; ++j, ++idof) nr[idof] = x_in[idof] - rc[j];
    }
    // calculate inertia tensor
    ARRAY_CLEAR(Ic, 9);
    ARRAY_CLEAR(invIc, 9);
    for (int iatom = 0, idof = 0; iatom < Natom; ++iatom, idof += 3) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                Ic[i * 4] += m_in[idof] * nr[idof + j] * nr[idof + j];
                Ic[i * 3 + j] -= m_in[idof] * nr[idof + i] * nr[idof + j];
            }
        }
    }
    // cal pseudo inverse of I
    // int succ = pseudo_inv(3, Ic, invIc, vectmp, 1e-8);
    // if (succ != 0) LOG(FATAL) << succ;
    PseudoInverse(Ic, invIc, 3);
    // cal vc and angular momentum
    ARRAY_CLEAR(vc, 3);
    for (int iatom = 0, idof = 0; iatom < Natom; ++iatom) {
        for (int i = 0; i < 3; ++i, ++idof) vc[i] += p_in[idof];
    }
    for (int i = 0; i < 3; ++i) vc[i] /= Mt;

    ARRAY_CLEAR(np, N);
    for (int i = 0, idof = 0; i < Natom; ++i) {
        for (int j = 0; j < 3; ++j, ++idof) np[idof] = p_in[idof] / m_in[idof] - vc[j];
    }

    ARRAY_CLEAR(Lc, 3);
    for (int iatom = 0, idof = 0; iatom < Natom; ++iatom, idof += 3) {
        cross(nr + iatom * 3, np + iatom * 3, vectmp);
        for (int i = 0; i < 3; ++i) Lc[i] += vectmp[i] * m_in[idof];
    }
    // cal center angular velocity
    ARRAY_MATMUL(invIc, Lc, Wc, 3, 3, 1);
    // set tangintial velocity to zero
    for (int iatom = 0, idof = 0; iatom < Natom; ++iatom) {
        cross(Wc, nr + idof, vectmp);
        for (int i = 0; i < 3; ++i, ++idof) p_in[idof] = m_in[idof] * (np[idof] - vectmp[i]);
    }
    // calculate force correction
    if (cal_force == true) {
        ARRAY_CLEAR(ac, 3);
        for (int iatom = 0, idof = 0; iatom < Natom; ++iatom) {
            for (int i = 0; i < 3; ++i, ++idof) ac[i] += F_in[idof];
        }
        for (int i = 0; i < 3; ++i) ac[i] /= Mt;

        ARRAY_CLEAR(nf, N);
        for (int iatom = 0, idof = 0; iatom < Natom; ++iatom) {
            for (int i = 0; i < 3; ++i, ++idof) nf[idof] = F_in[idof] / m_in[idof] - ac[i];
        }
        // calculate centroid torque
        ARRAY_CLEAR(Mc, 3);
        for (int iatom = 0, idof = 0; iatom < Natom; ++iatom, idof += 3) {
            cross(nr + idof, nf + idof, vectmp);
            for (int i = 0; i < 3; ++i) Mc[i] += vectmp[i] * m_in[idof];
        }
        // calculate centroid angular acceleration
        ARRAY_MATMUL(invIc, Mc, Ac, 3, 3, 1);
        // set tangintial force to zero
        for (int iatom = 0, idof = 0; iatom < Natom; ++iatom) {
            cross(Ac, nr + idof, vectmp);
            for (int i = 0; i < 3; ++i, ++idof) F_in[idof] = m_in[idof] * (nf[idof] - vectmp[i]);
        }
    }
    return 0;
}

int PIMDTraj_Solver::pseudo_inv(int N, num_real* A, num_real* invA, num_real* vectmp, num_real eps) {
    int rank, info;
    ARRAY_EYE(invA, N);
    ARRAY_CLEAR(vectmp, N);
    info = LAPACKE_dgelss(LAPACK_ROW_MAJOR, N, N, N, A, N, invA, N, vectmp, eps, &rank);
    return info;
}

int PIMDTraj_Solver::cross(num_real* vec1, num_real* vec2, num_real* prod) {
    prod[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    prod[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    prod[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
    return 0;
}
