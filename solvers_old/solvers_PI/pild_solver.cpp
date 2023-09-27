#include "pild_solver.h"

using namespace ARRAY_EG;

PILD_Solver::PILD_Solver(Param iparm, Model* pM) : PIMDTraj_Solver(iparm, pM) {
    pimd_type = PIMDTransformPolicy::Staging;  // for all PILD methods

    pild_flag     = Param_GetT(std::string, iparm, "pild_flag", "#mPILD");
    pild_type     = PILDPolicy::_dict.at(pild_flag);
    constrain_rot = Param_GetT(bool, iparm, "constrain_rotation", false);
    integ_type    = PIMDIntegratorPolicy::BAOAB;
    gam_ad        = Param_GetT(num_real, iparm, "gamma_ad", 1);

    try {  // allocate arrays
        ALLOCATE_PTR_TO_VECTOR(M_therm, NN);
    } catch (std::bad_alloc& e) { LOG(FATAL) << e.what(); }

    std::string suffix = pild_flag;
    save               = save + "_" + PILD_Solver::name() + suffix + "_" + pForceField->tag;
}

PILD_Solver::~PILD_Solver(){};

int PILD_Solver::update_thermo(const num_real& dt_in) {  ///< overload with different size
    plFunction();
    num_real gamma_ref = sqrt(bf2);
    for (int i = 0, iN = 0; i < P; ++i, iN += N) {
        if (i == 0) continue;  // skip first bead in Thermostat
        num_real *nki = nks + iN, *npi = nps + iN, *nmi = nms + iN;
        num_real gamma_opt = gamma_ref * sqrt(bfwgt[i]);
        pThermo->evolve(nki, npi, nmi, dt_in, N, iN, gamma_opt);
    }
    return 0;
}

int PILD_Solver::update_p(const num_real& dt_in) {
    plFunction();
    switch (pild_type) {
        case PILDPolicy::PILD:
            // p_1 = p1 - M_therm*M^{-1}*Fk_1*dt/2
            for (int i = 0; i < N; ++i) np[i] = nfks[i] / nms[i] * dt_in;
            ARRAY_MATMUL(np, M_therm, np, N, N, 1);
            break;
        case PILDPolicy::mPILD:
            for (int i = 0; i < N; ++i) np[i] = nfks[i] * dt_in;
            break;
        default:
            LOG(FATAL);
            break;
    }

    for (int i = 0; i < N; ++i) nps[i] -= np[i];
    // p_i = p_i - Fk_i*dt/2    i = 2,...,P
    for (int i = N; i < PN; ++i) nps[i] -= nfks[i] * dt_in;
    // LOG(INFO) << "update p";
    return 0;
}

int PILD_Solver::estimator(const int& isamp, TCFnucl& tcfer) {
    num_real esti_V, esti_Kprim, esti_Kvir, esti_Eprim, esti_Evir, esti_Cprim, esti_Cvir;

    // V
    esti_V = 0.0f;
    for (int i = 0; i < P; ++i) esti_V += vpeses[i];
    esti_V /= P;

    esti_Kprim = 0.5f * N / betap;
    esti_Kvir  = 0.5f * N / pThermo->beta;
    for (int i = 0, idx = 0; i < P; ++i) {
        for (int j = 0; j < N; ++j, ++idx) {
            if (i != 0) esti_Kprim -= 0.5f * bf2 * bfwgt[i] * nms[idx] * nks[idx] * nks[idx];
            esti_Kvir += 0.5f * (nrs[idx] - nrs[j]) * nfs[idx] / P;
        }
    }

    esti_Eprim = esti_V + esti_Kprim;
    esti_Evir  = esti_V + esti_Kvir;

    esti_Cprim = 0.f;
    esti_Cvir  = 0.f;

    if (isamp == 0) {
        ofs_ENER << FMT(8) << "flag" << FMT(8) << "time"  //
                 << FMT(8) << "N"                         //
                 << FMT(8) << "P"                         //
                 << FMT(8) << "beta"                      //
                 << FMT(8) << "esti_V"                    //
                 << FMT(8) << "esti_Kp"                   //
                 << FMT(8) << "esti_Kv"                   //
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
             << std::endl;


    for (int idof = 0; idof < N; ++idof) ofs_TRAJ << FMT(8) << nks[idof];
    for (int idof = 0; idof < N; ++idof) ofs_TRAJ << FMT(8) << nps[idof];
    ofs_TRAJ << std::endl;

    return 0;
}

/* int PILD_Solver::run_parallel() {
    int nsave       = (ntraj > mpi_nprocs) ? ntraj / mpi_nprocs : 1;
    num_real sampunit = dt * sstep * iou.time;
    init(-1);  // get occ0
               //    ARRAY_SHOW(nm, 1, N);
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
            // coll.Clear();
            try {
                traj(coll, PN);
            } catch (std::runtime_error& e) {  // if some error, output currect results
                LOG(WARNING) << "runtime_error cause breakdown for traj=" << itraj;
            }
            // collsum.Amount(coll);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // collmpi.MPIAmount(collsum);
        // if (mpi_rank == 0) collmpi.report(save + "-cache" + std::to_string(isave), sampunit);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // collmpi.MPIAmount(collsum);
    // if (mpi_rank == 0) collmpi.report(save, sampunit);
    // collsum.report(save + "-mpi" + std::to_string(mpi_rank), sampunit);
    // solverofs_mpi.close();

    final(-1);
    return 0;
}
 */
int PILD_Solver::init(const int& itraj) {
    init_ofs(itraj);

    if (itraj < 0) {
        // M_eff = M_st * gam_ad
        // omega_ad = omega_P / sqrt(gam_ad)
        pThermo->set_gammal(1 / sqrt(gam_ad));

        for (int i = 0; i < P; ++i) masswgt[i] = (i == 0) ? 1.0f : (i + 1.0f) / i * gam_ad;
        for (int i = 0; i < P; ++i) bfwgt[i] = (i == 0) ? 0 : 1.0f / gam_ad;  // @please check it
        switch (pild_type) {
            case PILDPolicy::mPILD:
                ARRAY_EYE(M_therm, N);
                break;

            case PILDPolicy::PILD: {
                HighFive::File Mthermfile("Mtherm.h5", HighFive::File::ReadWrite);
                Mthermfile.getDataSet("Mtherm").read<num_real>(M_therm);
            } break;

            default:
                LOG(FATAL);
        }

        for (int i = 0, J = 0; i < P; ++i) {
            pForceField->ForceField_init(nr0, np0, nm, N, itraj);
            for (int j = 0; j < N; ++j, ++J) nrs[J] = nr0[j], nps[J] = np0[j], nms[J] = masswgt[i] * nm[j];
        }
    } else {
        for (int i = 0, J = 0; i < P; ++i) {
            pForceField->ForceField_init(nr0, np0, nm, N, itraj);
            for (int j = 0; j < N; ++j, ++J) nrs[J] = nr0[j], nps[J] = np0[j], nms[J] = masswgt[i] * nm[j];
        }
        if (FLAGS_r != "") { rst_read(itraj); }
    }
    return 0;
}



/**
 * @brief Algorithm to find least root mean square distance and corresponding rigid transition
 *  X1'=R*X1+r then X1'~X2
 *
 * @param Natom number of atoms
 * @param D number of dimension
 * @param X1 first position array with size N*D
 * @param X2 second position array with size N*D
 * @param m mass array, should be non negative and have at least one positive value (size = N)
 * @param R rotation matrix
 * @param r translation vector
 * @return least root mean square distance
 */
int Kabsch(int Natom, int D, num_real* X1, num_real* X2, num_real* m, num_real* R, num_real* r) {
    num_real rc1[D], rc2[D], C[D * D];
    Eigen::Map<EigX> MapX1(X1, Natom, D);
    Eigen::Map<EigX> MapX2(X2, Natom, D);
    Eigen::Map<EigX> MapR(R, Natom, D);
    Eigen::Map<Eigen::VectorXd> Mapm(m, Natom);
    Eigen::Map<EigX> MapC(C, D, D);
    num_real Mt = 0.0f;
    ARRAY_CLEAR(rc1, D);
    ARRAY_CLEAR(rc2, D);
    // calculate total mass and rc
    for (int i = 0, idof = 0; i < Natom; ++i) {
        Mt += m[idof];
        for (int j = 0; j < D; ++j, ++idof) {
            rc1[j] += X1[idof] * m[idof];
            rc2[j] += X2[idof] * m[idof];
        }
    }
    for (int i = 0; i < D; ++i) {
        rc1[i] /= Mt;
        rc2[i] /= Mt;
    }
    // move both X1 and X2 to center
    for (int i = 0, idof = 0; i < Natom; ++i) {
        for (int j = 0; j < D; ++j, ++idof) {
            X1[idof] -= rc1[j];
            X2[idof] -= rc2[j];
        }
    }
    // calculate C matrix
    MapC = MapX1.transpose() * Mapm.asDiagonal() / Mt * MapX2;
    // SVD decomposition
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(MapC, Eigen::ComputeFullU | Eigen::ComputeFullV);
    // calculate rotation matrix
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(D, D);
    if ((svd.matrixV() * svd.matrixU().transpose()).determinant() < 0) I(D, D) = -1;
    MapR = svd.matrixV() * I * svd.matrixU().transpose();
    // calculate translation vector
    ARRAY_CLEAR(r, D);
    for (int i = 0; i < D; ++i) {
        r[i] += rc2[i];
        for (int j = 0; j < D; ++j) r[i] -= MapR(i, j) * rc1[j];
    }
    return 0;
}
