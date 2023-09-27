#include "ppimd_solver.h"

#include <sys/time.h>

#include <cstring>

#include "../../utils/mpi_utils.h"

using namespace ARRAY_EG;

PPIMDPARATraj_Solver::PPIMDPARATraj_Solver(Param iparm, Model* pM) : PIMDPARATraj_Solver(iparm, pM) {
    // n_term = Param_GetT(int, parm, "n_term", 3);
    n_term = 3;  // currently set to 3;
    CHECK_EQ(P % n_term, 0);
    pon = P / n_term;

    split_t1                   = Param_GetT(num_real, parm, "split_t", 1.0 / 3.0);
    split_v1                   = Param_GetT(num_real, parm, "split_v", 1.0 / 3.0);
    split_c1                   = Param_GetT(num_real, parm, "split_c", 1.0 / 3.0);
    std::string highorder_flag = Param_GetT(std::string, parm, "highorder_flag", "ppi");
    highorder_type             = PIMDPARAHighorderPolicy::_dict.at(highorder_flag);

    use_liquidne = Param_GetT(int, parm, "use_liquidne", 0);

    try {
        ALLOCATE_PTR_TO_VECTOR(split_t, P);
        ALLOCATE_PTR_TO_VECTOR(split_v, P);
        ALLOCATE_PTR_TO_VECTOR(split_c, P);
        ALLOCATE_PTR_TO_VECTOR(stag_a, P);
        ALLOCATE_PTR_TO_VECTOR(stag_b, P);
        ALLOCATE_PTR_TO_VECTOR(stag_c, P);
        ALLOCATE_PTR_TO_VECTOR(nhs, P);
        ALLOCATE_PTR_TO_VECTOR(nrcent, N);
    } catch (std::bad_alloc& e) { LOG(FATAL) << e.what(); }

    save = name() + "_" + pForceField->tag;
}

int PPIMDPARATraj_Solver::init(const int& itraj) {
    init_ofs(itraj);

    if (itraj < 0) {
        int pnbd = ceil((num_real) P / mpi_nprocs);
        for (int i = 0; i < mpi_nprocs; ++i) { mpi_range_array[i] = std::min((i + 1) * pnbd, P); }

        for (int i = 0; i < P; ++i) bfwgt[i] = 1.0;

        num_real t1 = split_t1, t2 = 1.0 - 2 * t1;
        num_real v1 = split_v1, v2 = (1.0 - v1) * 0.5;
        num_real c1 = split_c1, c2 = (1.0 - c1) * 0.5;

        kappa = ((2.0 - 3.0 * t2) * (1.0 - v1) * (1.0 - v1) + (3.0 * t2 * t2 - 1.0) * (1.0 - v1) + v1 * v1) / 24.0;

        num_real beta = pThermo->beta;
        // bf2       = pon / beta / beta;
        sum_V     = 0.0;
        sum_Kprim = 0.0;
        sum_Kvir  = 0.0;
        sum_F2    = 0.0;
        sum_F2v   = 0.0;
        sum_F2kp  = 0.0;
        sum_F2kv  = 0.0;

        num_real t_series[3] = {t1, t1, t2};
        num_real v_series[3] = {v1, v2, v2};
        num_real c_series[3] = {c1, c2, c2};

        for (int idx = 0; idx < P; idx += n_term) {
            num_real* t = split_t + idx;
            num_real* v = split_v + idx;
            num_real* c = split_c + idx;
            memcpy(t, t_series, sizeof(t_series));
            memcpy(v, v_series, sizeof(v_series));
            memcpy(c, c_series, sizeof(c_series));
        }

        switch (pimd_type) {
            case PIMDPARATransformPolicy::Primitive:  // to be implimented
                LOG(FATAL);
                break;
            case PIMDPARATransformPolicy::Staging: {
                num_real stag_q[P];

                stag_a[0] = 0.0, stag_b[0] = 1.0, stag_c[0] = 1.0 / split_t1, stag_q[0] = 0.0;
                for (int i = 1; i < P; ++i) {
                    stag_c[i]  = 1.0 / split_t[(i + 1) % P];
                    stag_q[i]  = stag_c[i] + stag_c[i - 1] - stag_c[i - 1] * stag_a[i - 1];
                    masswgt[i] = stag_q[i] / n_term;
                    stag_a[i]  = stag_c[i] / stag_q[i];
                    stag_b[i]  = stag_c[i - 1] * stag_b[i - 1] / stag_q[i];
                }
                masswgt[0] = 1.0;

                break;
            }
            case PIMDPARATransformPolicy::NormalMode:  // to be implemented
                LOG(FATAL);
                break;
            default:
                LOG(FATAL);
        }

        switch (highorder_type) {
            case PIMDPARAHighorderPolicy::PPI:
                fppi = beta * beta * beta * kappa / pon / pon;
                break;
            case PIMDPARAHighorderPolicy::FD4T3V:
                fppi = beta * beta * kappa / pon / pon;
                break;
            default:
                LOG(FATAL);
        }
        if (!use_liquidne) {
            for (int i = 0, ibead_j = 0; i < P; ++i) {
                pForceField->ForceField_init(nr0, np0, nm, N, itraj);
                for (int j = 0; j < N; ++j, ++ibead_j) {
                    nrs[ibead_j] = nr0[j], nps[ibead_j] = np0[j], nms[ibead_j] = masswgt[i] * nm[j];
                    // nrs[ibead_j] = nr0[j], nps[ibead_j] = 0.0, nms[ibead_j] = masswgt[i] * nm[j];
                }
            }
        } else {
            pForceField->ForceField_init(nrs, nps, nm, PN, itraj);
            for (int i = 0, ibead_j = 0; i < P; ++i)
                for (int j = 0; j < N; ++j, ++ibead_j) {
                    // nrs[ibead_j] /= phys::au_2_ang * 1e-10L;
                    nms[ibead_j] = masswgt[i] * nm[j];
                }
        }
    } else {
        if (!use_liquidne) {
            for (int i = 0, ibead_j = 0; i < P; ++i) {
                pForceField->ForceField_init(nr0, np0, nm, N, itraj);
                for (int j = 0; j < N; ++j, ++ibead_j) {
                    nrs[ibead_j] = nr0[j], nps[ibead_j] = np0[j], nms[ibead_j] = masswgt[i] * nm[j];
                    // nrs[ibead_j] = nr0[j], nps[ibead_j] = 0.0, nms[ibead_j] = masswgt[i] * nm[j];
                }
            }
        } else {
            pForceField->ForceField_init(nrs, nps, nm, PN, itraj);
            for (int i = 0, ibead_j = 0; i < P; ++i)
                for (int j = 0; j < N; ++j, ++ibead_j) {
                    // nrs[ibead_j] /= phys::au_2_ang * 1e-10L;
                    nms[ibead_j] = masswgt[i] * nm[j];
                }
        }

        if (FLAGS_r != "") rst_read(itraj);
    }

    return 0;
}

int PPIMDPARATraj_Solver::all_X2K() {
    switch (pimd_type) {
        case PIMDPARATransformPolicy::Primitive:
            LOG(FATAL);
            break;
        case PIMDPARATransformPolicy::Staging:
            // first bead
            for (int first_j = 0; first_j < N; ++first_j) nks[first_j] = nrs[first_j];

            // other beads
            for (int i = P - 2, ibead_j = PN - N - 1, ipbead_j = PN - 1; i > 0; --i) {
                for (int first_j = N - 1; first_j > -1; --first_j, --ibead_j, --ipbead_j) {
                    nks[ibead_j] = nrs[ibead_j] - stag_a[i] * nrs[ipbead_j] - stag_b[i] * nrs[first_j];
                }
            }
            // last bead
            for (int first_j = 0, last_j = PN - N; first_j < N; ++first_j, ++last_j) {
                nks[last_j] = nrs[last_j] - (stag_a[P - 1] + stag_b[P - 1]) * nrs[first_j];
            }

            break;
        case PIMDPARATransformPolicy::NormalMode:
            LOG(FATAL);
            break;
        default:
            LOG(FATAL);
    }
    return 0;
}

int PPIMDPARATraj_Solver::all_K2X() {
    switch (pimd_type) {
        case PIMDPARATransformPolicy::Primitive:
            LOG(FATAL);
            break;
        case PIMDPARATransformPolicy::Staging:
            // first bead
            for (int first_j = 0; first_j < N; ++first_j) nrs[first_j] = nks[first_j];

            // last bead
            for (int first_j = 0, last_j = PN - N; first_j < N; ++first_j, ++last_j) {
                nrs[last_j] = nks[last_j] + nks[first_j];
            }

            // other beads
            for (int i = P - 2, ibead_j = PN - N - 1, ipbead_j = PN - 1; i > 0; --i) {
                for (int first_j = N - 1; first_j > -1; --first_j, --ibead_j, --ipbead_j) {
                    nrs[ibead_j] = stag_a[i] * nrs[ipbead_j] + stag_b[i] * nks[first_j] + nks[ibead_j];
                }
            }
            break;
        case PIMDPARATransformPolicy::NormalMode:
            LOG(FATAL);
            break;
        default:
            LOG(FATAL);
    }
    return 0;
}

int PPIMDPARATraj_Solver::all_FX2FK() {
    switch (pimd_type) {
        case PIMDPARATransformPolicy::Primitive:
            LOG(FATAL);
            break;
        case PIMDPARATransformPolicy::Staging:
            switch (highorder_type) {
                case PIMDPARAHighorderPolicy::PPI:
                    if (mpi_isroot) {
                        for (int first_j = 0; first_j < N; ++first_j) {
                            nfks[first_j] = split_v[0] * nfs[first_j] / pon;
                        }

                        for (int i = 1, ibead_j = N, imbead_j = 0; i < P; ++i) {
                            for (int first_j = 0; first_j < N; ++first_j, ++ibead_j, ++imbead_j) {
                                nfks[ibead_j] = split_v[i] * nfs[ibead_j] / pon + stag_a[i - 1] * nfks[imbead_j];
                                nfks[first_j] += split_v[i] * nfs[ibead_j] / pon;
                            }
                        }
                    }
                    break;
                case PIMDPARAHighorderPolicy::FD4T3V:
                    cal_FDhessian();
                    if (mpi_isroot) {
                        for (int first_j = 0; first_j < N; ++first_j) {
                            nfks[first_j] =
                                split_v[0] * nfs[first_j] / pon + 2.0 * fppi * split_c[0] * nhs[first_j] / pon;
                        }

                        for (int i = 1, ibead_j = N, imbead_j = 0; i < P; ++i) {
                            for (int first_j = 0; first_j < N; ++first_j, ++ibead_j, ++imbead_j) {
                                num_real ftmp =
                                    split_v[i] * nfs[ibead_j] / pon + 2.0 * fppi * split_c[i] * nhs[ibead_j] / pon;
                                nfks[ibead_j] = ftmp + stag_a[i - 1] * nfks[imbead_j];
                                nfks[first_j] += ftmp;
                            }
                        }
                    }
                    break;
                default:
                    LOG(FATAL);
            }
            break;
        case PIMDPARATransformPolicy::NormalMode:
            LOG(FATAL);
            break;
        default:
            LOG(FATAL);
    }
    return 0;
}

int PPIMDPARATraj_Solver::ppiestimator(const int& isamp, TCFnucl& tcfer) {
    // num_real esti_V, esti_Kprim, esti_Kvir;//, esti_Eprim, esti_Evir;
    // num_real esti_F2, esti_Bf2;//, esti_F2v, esti_F2kp, esti_F2kv;
    num_real esti_Vppi, esti_Kpppi, esti_Kvppi;
    num_real head_F2kp, head_F2kv;

    num_real beta = pThermo->beta;

    esti_V     = 0.0;
    esti_Kprim = 0.0;
    esti_Kvir  = 0.0;
    esti_F2    = 0.0;
    esti_Bf2   = 0.0;

    cal_FDhessian();

    if (mpi_isroot) {
        for (int i = 0; i < P; ++i) esti_V += split_v[i] * vpeses[i] / pon;
        // for (int i = 0; i < P; ++i) esti_V += vpeses[i] / P;

        for (int i = 0, ibead_j = 0; i < P; ++i) {
            for (int j = 0; j < N; ++j, ++ibead_j) {
                // if (i) esti_Kprim += 0.5 * bf2 * bfwgt[i] * stag_c[i] * nms[ibead_j] * nks[ibead_j] * nks[ibead_j] /
                // n_term;
                if (i) esti_Kprim += 0.5 * bf2 * bfwgt[i] * nms[ibead_j] * nks[ibead_j] * nks[ibead_j];
                esti_Kvir +=
                    0.5 * split_v[i] * (nrs[ibead_j] - nr[j]) * nfs[ibead_j] / pon;  // caution of the sign of nfs
                esti_Bf2 += split_c[i] * nhs[ibead_j] * (nrs[ibead_j] - nr[j]) / beta / pon;
                esti_F2 += split_c[i] * nfs[ibead_j] * nfs[ibead_j] / nm[j] / pon;
            }
        }

        sum_V += esti_V;
        sum_Kprim += esti_Kprim;
        sum_Kvir += esti_Kvir;

        sum_F2 += esti_F2;
        sum_F2v += esti_V * esti_F2;
        sum_F2kp += esti_Kprim * esti_F2;
        sum_F2kv += esti_Kvir * esti_F2 - esti_Bf2;
        head_F2kp  = 0.5 * N / betap * sum_F2 - sum_F2kp;
        head_F2kv  = 0.5 * N / beta * sum_F2 + sum_F2kv;
        esti_V     = sum_V / (isamp + 1);
        esti_Kprim = 0.5 * N / betap - sum_Kprim / (isamp + 1);
        esti_Kvir  = 0.5 * N / beta + sum_Kvir / (isamp + 1);
        esti_Vppi  = esti_V + fppi * (esti_V * sum_F2 - sum_F2v + 2.0 * sum_F2 / beta) / (isamp + 1);
        esti_Kpppi = esti_Kprim + fppi * (esti_Kprim * sum_F2 - head_F2kp + sum_F2 / beta) / (isamp + 1);
        esti_Kvppi = esti_Kvir + fppi * (esti_Kvir * sum_F2 - head_F2kv + sum_F2 / beta) / (isamp + 1);

        num_real sampunit = dt * sstep * iou.time;
        if (!isamp) {
            ofs_ENER << FMT(8) << "itraj" << FMT(8) << "time"  //
                     << FMT(8) << "N"                          //
                     << FMT(8) << "P"                          //
                     << FMT(8) << "beta"                       //
                     << FMT(8) << "esti_V"                     //
                     << FMT(8) << "esti_Kprim"                 //
                     << FMT(8) << "esti_Kvir"                  //
                     << FMT(8) << "esti_Vppi"                  //
                     << FMT(8) << "esti_Kpppi"                 //
                     << FMT(8) << "esti_Kvppi"                 //
                     << FMT(8) << "esti_F2"                    //
                     << FMT(8) << "esti_Bf2"                   //
                     << FMT(8) << "sum_F2"                     //
                     << FMT(8) << "sum_F2v"                    //
                     << FMT(8) << "sum_F2kp"                   //
                     << FMT(8) << "sum_F2kv"                   //
                     << FMT(8) << "head_F2kp"                  //
                     << FMT(8) << "head_F2kv"                  //
                     << std::endl;
        }
        ofs_ENER << FMT(8) << itraj << FMT(8) << isamp * sampunit  //
                 << FMT(8) << N                                    //
                 << FMT(8) << P                                    //
                 << FMT(8) << pThermo->beta                        //
                 << FMT(8) << esti_V                               //
                 << FMT(8) << esti_Kprim                           //
                 << FMT(8) << esti_Kvir                            //
                 << FMT(8) << esti_Vppi                            //
                 << FMT(8) << esti_Kpppi                           //
                 << FMT(8) << esti_Kvppi                           //
                 << FMT(8) << esti_F2                              //
                 << FMT(8) << esti_Bf2                             //
                 << FMT(8) << sum_F2                               //
                 << FMT(8) << sum_F2v                              //
                 << FMT(8) << sum_F2kp                             //
                 << FMT(8) << sum_F2kv                             //
                 << FMT(8) << head_F2kp                            //
                 << FMT(8) << head_F2kv                            //
                 << std::endl;

        for (int idof = 0; idof < PN; ++idof) ofs_TRAJ << FMT(8) << nrs[idof];
        for (int idof = 0; idof < PN; ++idof) ofs_TRAJ << FMT(8) << nfs[idof];
        ofs_TRAJ << std::endl;
    }

    return 0;
}

int PPIMDPARATraj_Solver::fdestimator(const int& isamp, TCFnucl& tcfer) {
    // num_real esti_V, esti_Kprim, esti_Kvir;//, esti_Eprim, esti_Evir;
    // num_real esti_F2, esti_Bf2;//, esti_F2v, esti_F2kp, esti_F2kv;
    // num_real esti_Vppi, esti_Kpppi, esti_Kvppi;
    // num_real head_F2kp, head_F2kv;

    num_real beta = pThermo->beta;

    esti_V     = 0.0;
    esti_Kprim = 0.0;
    esti_Kvir  = 0.0;
    esti_F2    = 0.0;
    esti_Bf2   = 0.0;

    if (mpi_isroot) {
        for (int i = 0; i < P; ++i) esti_V += split_v[i] * vpeses[i] / pon;
        // for (int i = 0; i < P; ++i) esti_V += vpeses[i] / P;

        for (int i = 0, ibead_j = 0; i < P; ++i) {
            for (int j = 0; j < N; ++j, ++ibead_j) {
                if (i) esti_Kprim += 0.5 * bf2 * bfwgt[i] * nms[ibead_j] * nks[ibead_j] * nks[ibead_j];
                esti_Kvir += (nrs[ibead_j] - nr[j]) *
                             (0.5 * split_v[i] * nfs[ibead_j] + fppi * split_c[i] * nhs[ibead_j]) /
                             pon;  // caution of the sign of nfs
                // esti_Bf2  += split_c[i] * nhs[ibead_j] * (nrs[ibead_j] - nr[j]) / beta / pon;
                esti_F2 += fppi * split_c[i] * nfs[ibead_j] * nfs[ibead_j] / nm[j] / pon;
            }
        }

        esti_V += 2.0 * esti_F2;
        esti_Kprim = 0.5 * N / betap - esti_Kprim + esti_F2;
        esti_Kvir += 0.5 * N / beta + esti_F2;

        sum_V += esti_V + 2.0 * esti_F2;

        // sum_V      += esti_V;
        // sum_Kprim  += esti_Kprim;
        // sum_Kvir   += esti_Kvir;

        // sum_F2     += esti_F2;
        // sum_F2v    += esti_V * esti_F2;
        // sum_F2kp   += esti_Kprim * esti_F2;
        // sum_F2kv   += esti_Kvir * esti_F2 - esti_Bf2;
        // head_F2kp  = 0.5 * N / betap * sum_F2 - sum_F2kp;
        // head_F2kv  = 0.5 * N / beta * sum_F2 + sum_F2kv;
        // esti_V     = sum_V / (isamp+1);
        // esti_Kprim = 0.5 * N / betap - sum_Kprim / (isamp+1);
        // esti_Kvir  = 0.5 * N / beta + sum_Kvir / (isamp+1);
        // esti_Vppi  = esti_V + fppi * (esti_V * sum_F2 - sum_F2v + 2.0 * sum_F2 / beta) / (isamp+1);
        // esti_Kpppi = esti_Kprim + fppi * (esti_Kprim * sum_F2 - head_F2kp + sum_F2 / beta) / (isamp+1);
        // esti_Kvppi = esti_Kvir + fppi * (esti_Kvir * sum_F2 - head_F2kv + sum_F2 / beta) / (isamp+1);

        // sum_F2kv   += sum_V + sum_Kprim + sum_Kvir + esti_V + esti_Kprim + esti_Kvir + esti_Vppi + esti_Kpppi +
        // esti_Kvppi; sum_F2kv   += sum_F2 + sum_F2v + sum_F2kp + esti_Bf2 + esti_F2 + beta;

        num_real sampunit = dt * sstep * iou.time;
        if (!isamp) {
            ofs_ENER << FMT(8) << "itraj" << FMT(8) << "time"  //
                     << FMT(8) << "N"                          //
                     << FMT(8) << "P"                          //
                     << FMT(8) << "beta"                       //
                     << FMT(8) << "esti_V"                     //
                     << FMT(8) << "esti_Kprim"                 //
                     << FMT(8)
                     << "esti_Kvir"  //
                     //<< FMT(8) << "esti_Vppi"                            //
                     //<< FMT(8) << "esti_Kpppi"                           //
                     //<< FMT(8) << "esti_Kvppi"                           //
                     << FMT(8)
                     << "esti_F2"  //
                     //<< FMT(8) << "esti_Bf2"                           //
                     //<< FMT(8) << "sum_F2"                           //
                     //<< FMT(8) << "sum_F2v"                           //
                     //<< FMT(8) << "sum_F2kp"                           //
                     //<< FMT(8) << "sum_F2kv"                           //
                     << std::endl;
        }
        ofs_ENER << FMT(8) << itraj << FMT(8) << isamp * sampunit  //
                 << FMT(8) << N                                    //
                 << FMT(8) << P                                    //
                 << FMT(8) << pThermo->beta                        //
                 << FMT(8) << esti_V                               //
                 << FMT(8) << esti_Kprim                           //
                 << FMT(8)
                 << esti_Kvir  //
                 //<< FMT(8) << esti_Vppi                            //
                 //<< FMT(8) << esti_Kpppi                           //
                 //<< FMT(8) << esti_Kvppi                           //
                 << FMT(8)
                 << esti_F2  //
                 //<< FMT(8) << esti_Bf2                           //
                 //<< FMT(8) << sum_F2                           //
                 //<< FMT(8) << sum_F2v                           //
                 //<< FMT(8) << sum_F2kp                           //
                 //<< FMT(8) << sum_F2kv                           //
                 << std::endl;

        for (int idof = 0; idof < PN; ++idof) ofs_TRAJ << FMT(8) << nrs[idof];
        for (int idof = 0; idof < PN; ++idof) ofs_TRAJ << FMT(8) << nfs[idof];
        ofs_TRAJ << std::endl;
    }

    return 0;
}

int PPIMDPARATraj_Solver::estimator(const int& isamp, TCFnucl& tcfer) {
    LOG(INFO) << vpeses[0] << ' ' << vpeses[1] << ' ' << vpeses[2];
    switch (highorder_type) {
        case PIMDPARAHighorderPolicy::PPI:
            ppiestimator(isamp, tcfer);
            break;
        case PIMDPARAHighorderPolicy::FD4T3V:
            fdestimator(isamp, tcfer);
            break;
        default:
            LOG(FATAL);
    }
    return 0;
}

int PPIMDPARATraj_Solver::cal_cent() {
    ARRAY_CLEAR(nr, N);

    for (int i = 0, ibead_j = 0; i < P; ++i) {
        for (int j = 0; j < N; ++j, ++ibead_j) { nr[j] += nrs[ibead_j] / P; }
    }
    return 0;
}

int PPIMDPARATraj_Solver::ff_calc1(const int& level) {  // multi-forcefield calculation at a fixed nr, np
    plFunction();
    if (mpi_isroot) {
        all_K2X();  // From nks to nrs (\xi to r)

        // calc centroid r and saved in traditional md variables
        ARRAY_CLEAR(nr, N);
        cal_cent();
    }
    mpiSendx(nrs);

    ff_calc2(nrs, nps, vpeses, nfs);
    mpiRecvf(nfs, N);
    mpiRecvf(vpeses, 1);
    // if (mpi_isroot) {
    all_FX2FK();  ///< transfrom f to fk
    //}
    return 0;
}

int PPIMDPARATraj_Solver::ff_calc2(num_real* trs, num_real* tps, num_real* tvps, num_real* tfs, const int& level) {
    if (!use_liquidne) {
        int pstart = 0, pend = P;
        if (para_type == ParaPolicy::calc_para) {
            if (mpi_isroot) {
                pstart = 0;
                pend   = mpi_range_array[0];
            } else {
                pstart = mpi_range_array[mpi_rank - 1];
                pend   = mpi_range_array[mpi_rank];
            }
        }

        for (int i = pstart, idx = pstart * N; i < pend; ++i, idx += N) {
            num_real* nri   = trs + idx;
            num_real* npi   = tps + idx;
            num_real* vpesi = tvps + i;
            num_real* gradi = grads + idx;
            num_real* nfi   = tfs + idx;
            num_real* hessi = hesses + 0;  ///< @todo don't use hessian

            pForceField->ForceField_npes(vpesi, gradi, hessi, nri, npi, level, N, itraj, istep);

            for (int j = 0; j < N; ++j) nfi[j] = gradi[j];  // a copy
        }
    } else {
        pForceField->ForceField_npes(tvps, tfs, hesses, trs, tps, level, PN, itraj, istep);

        // for (int i = 0; i < PN; ++i) {
        //    trs[i] /= phys::au_2_ang * 1e-10L;
        //    tfs[i] = -tfs[i] * phys::au_2_ang * 1e-10L;
        //}
    }
    return 0;
}

int PPIMDPARATraj_Solver::cal_FDhessian() {
    num_real u_temp[PN], nrs_temp1[PN], nrs_temp2[PN], nfs_temp1[PN], nfs_temp2[PN];
    num_real vp_temp[P];
    num_real delta, sum_f2m2;
    num_real eps = 1.e-15;

    if (mpi_isroot) {
        ARRAY_CLEAR(u_temp, PN);
        sum_f2m2 = 0.0;
        for (int i = 0, ibead_j = 0; i < P; ++i) {
            for (int j = 0; j < N; ++j, ++ibead_j) {
                u_temp[ibead_j] = nfs[ibead_j] / nm[j];
                sum_f2m2 += nfs[ibead_j] * nfs[ibead_j] / nm[j] / nm[j];
            }
        }

        delta = eps / sqrt(sum_f2m2 / PN);

        for (int i = 0, ibead_j = 0; i < P; ++i) {
            for (int j = 0; j < N; ++j, ++ibead_j) {
                nrs_temp1[ibead_j] = nrs[ibead_j] + delta * u_temp[ibead_j];
                nrs_temp2[ibead_j] = nrs[ibead_j] - delta * u_temp[ibead_j];
            }
        }
    }

    mpiSendx(nrs_temp1);

    ff_calc2(nrs_temp1, nps, vp_temp, nfs_temp1);

    mpiRecvf(nfs_temp1, N);

    // LOG(INFO) << "vp_temp1: " << vp_temp[0] << " " << vp_temp[2] << " " << vp_temp[4] << " " << vp_temp[6];

    mpiSendx(nrs_temp2);

    ff_calc2(nrs_temp2, nps, vp_temp, nfs_temp2);

    mpiRecvf(nfs_temp2, N);

    // LOG(INFO) << "vp_temp2: " << vp_temp[0] << " " << vp_temp[2] << " " << vp_temp[4] << " " << vp_temp[6];

    if (mpi_isroot) {
        for (int i = 0, ibead_j = 0; i < P; ++i) {
            for (int j = 0; j < N; ++j, ++ibead_j) {
                nhs[ibead_j] = 0.5 * (nfs_temp1[ibead_j] - nfs_temp2[ibead_j]) / delta;
            }
        }
    }

    return 0;
}
