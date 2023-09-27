#include "cmd_solver.h"

using namespace ARRAY_EG;

CentroidMD_Solver::CentroidMD_Solver(Param iparm, Model* pM) : PIMDTraj_Solver(iparm, pM) {
    pimd_type = PIMDTransformPolicy::NormalMode;  // for all centroid method

    CMD_flag               = Param_GetT(std::string, iparm, "cmd_flag", "#CMD");
    cmd_type               = CentroidMDPolicy::_dict.at(CMD_flag);
    std::string integrator = Param_GetT(std::string, parm, "integrator", "BCOCB");
    integ_type             = PIMDIntegratorPolicy::_dict.at(integrator);
    gam_ad                 = Param_GetT(num_real, iparm, "gamma_ad");  // Required, should be smaller than 1

    // delete switch scope, for optimized gammal can be obtained from bf2 & bfwgt, in CMD or TRPMD

    try {  // allocate arrays

    } catch (std::bad_alloc& e) { LOG(FATAL) << e.what(); }

    std::string suffix = CMD_flag;
    save               = save + "_" + name() + suffix + "_" + pForceField->tag;
}

CentroidMD_Solver::~CentroidMD_Solver(){};

int CentroidMD_Solver::update_thermo(const num_real& dt_in) {  ///< overload with different size
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

int CentroidMD_Solver::init(const int& itraj) {
    init_ofs(itraj);

    if (itraj < 0) {
        for (int i = 0; i < P; ++i) bfwgt[i] = 1.0f;
        for (int i = 0; i < P; ++i) masswgt[i] = 1.0f;

        CHECK_EQ(P % 2, 0);
        int idxM = 0, idxT = 0;
        bfwgt[idxM++] = 0.0f;
        for (int i = 0; i < P; ++i) Tran[idxT++] = sqrt(1.0f / P);
        for (int i = 1; i < P / 2; ++i) {
            bfwgt[idxM++] = 2.0f * (1.0f - cos(phys::math::twopi * i / P));  // @BUG
            for (int j = 0; j < P; ++j) Tran[idxT++] = cos(phys::math::twopi * i * j / P) * sqrt(2.0f / P);
        }

        bfwgt[idxM++] = 4.0f;
        for (int j = 0; j < P / 2; ++j) {
            Tran[idxT++] = -sqrt(1.0f / P);
            Tran[idxT++] = sqrt(1.0f / P);
        }

        for (int i = P / 2 + 1; i < P; ++i) {
            bfwgt[idxM++] = 2.0f * (1.0f - cos(phys::math::twopi * i / P));  // @BUG
            for (int j = 0; j < P; ++j) Tran[idxT++] = sin(phys::math::twopi * i * j / P) * sqrt(2.0f / P);
        }

        switch (cmd_type) {
            case CentroidMDPolicy::CMD:
                for (int i = 1; i < P; ++i) {
                    masswgt[i] = bfwgt[i] / (gam_ad * gam_ad);
                    bfwgt[i]   = gam_ad * gam_ad;
                }
                break;
            case CentroidMDPolicy::TRPMD:
                for (int i = 0; i < P; ++i) masswgt[i] = 1.0f;
                break;
            default:
                LOG(FATAL);
        }

        for (int i = 0, J = 0; i < P; ++i) {
            pForceField->ForceField_init(nr0, np0, nm, N, itraj);
            for (int j = 0; j < N; ++j, ++J) nrs[J] = nr0[j], nps[J] = np0[j], nms[J] = masswgt[i] * nm[j];
        }
    } else {
        // something wrong? why make
        for (int i = 0, J = 0; i < P; ++i) {
            pForceField->ForceField_init(nr0, np0, nm, N, itraj);
            for (int j = 0; j < N; ++j, ++J) nrs[J] = nr0[j], nps[J] = np0[j], nms[J] = masswgt[i] * nm[j];
        }
        if (FLAGS_r != "") rst_read(itraj);
    }
    return 0;
}
int CentroidMD_Solver::estimator(const int& isamp, TCFnucl& tcfer) {
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


    for (int idof = 0; idof < N; ++idof) ofs_TRAJ << FMT(8) << nks[idof] * Tran[0];
    switch (cmd_type)
    {
    case CentroidMDPolicy::CMD:
        for (int idof = 0; idof < N; ++idof) ofs_TRAJ << FMT(8) << nps[idof] * Tran[0];
        break;
    case CentroidMDPolicy::TRPMD:
        for (int idof = 0; idof < N; ++idof) ofs_TRAJ << FMT(8) << nps[idof];
        break;
    default:
        LOG(FATAL);
        break;
    }
    ofs_TRAJ << std::endl;


    return 0;
}

/* int CentroidMD_Solver::run_parallel() {
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
} */
