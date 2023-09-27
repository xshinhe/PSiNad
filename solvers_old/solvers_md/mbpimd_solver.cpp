#include "mbpimd_solver.h"

MBPIMDTraj_Solver::MBPIMDTraj_Solver(Param iparm, Model* pM) : PIMDTraj_Solver(iparm, pM) {
    // reset flags
    pimd_type  = PIMDTransformPolicy::Primitive;
    integ_type = PIMDIntegratorPolicy::BAOAB;

    try {  // allocate arrays
        ALLOCATE_PTR_TO_VECTOR(VHO, Natom);
        ALLOCATE_PTR_TO_VECTOR(fV, Natom);
        ALLOCATE_PTR_TO_VECTOR(fE, Natom * Natom);
        ALLOCATE_PTR_TO_VECTOR(dV_spring, Natom * PN);
        ALLOCATE_PTR_TO_VECTOR(dE_spring, Natom * Natom * PN);
    } catch (std::bad_alloc& e) { LOG(FATAL) << e.what(); }

    save = name() + "_" + pForceField->tag;
}

MBPIMDTraj_Solver::~MBPIMDTraj_Solver(){};

int MBPIMDTraj_Solver::spring_force() {
    int PendN     = (P - 1) * N;
    num_real* bd1 = nrs;          // the first bead
    num_real* bdP = nrs + PendN;  // the last bead
    /**
     * @note calculate fE factor and dE_spring
     *
     *     fE(a0, a1) = spring {a0_1~a0_2~~a0_P:ax_1~~:a1_1~a1_2~~a1_P}
     *                          |__________________________________|
     */
    for (int ax = 0; ax < Natom; ++ax) {
        for (int a1 = ax; a1 < Natom; ++a1) {
            int a0   = a1 - ax;
            int a0a1 = a0 * Natom + a1;
            if (ax == 0) {  // directly
                fE[a0a1] = 0.0f;
                for (int d = 0, a1j = a1 * Ndim; d < Ndim; ++d, ++a1j) {
                    for (int i = 0, ip = i + 1; i < P; ++i, ++ip %= P) {
                        num_real Diff = nrs[i * N + a1j] - nrs[ip * N + a1j];
                        fE[a0a1] += 0.5f * nm[a1j] * bf2 * Diff * Diff;
                    }
                }
                fE[a0a1] = std::exp(-pThermo->beta * fE[a0a1]);
            } else {  // or by iterative scheme
                int ah   = a0 + 1;
                int a0a0 = a0 * Natom + a0;
                int aha1 = ah * Natom + a1;

                fE[a0a1] = fE[a0a0] * fE[aha1];
                for (int d = 0, a0j = a0 * Ndim, ahj = ah * Ndim, a1j = a1 * Ndim;  //
                     d < Ndim; ++d, ++a0j, ++ahj, ++a1j) {
                    num_real corr = (bd1[a0j] - bd1[ahj]) * (bdP[a0j] - bdP[a1j]);  // cancel factor 2
                    fE[a0a1] *= std::exp(-pThermo->beta * nm[a1j] * bf2 * corr);    // cancel factor 0.5
                }
            }
            num_real* dE01 = dE_spring + a0a1 * PN;
            for (int i = 0; i < P; ++i) {
                for (int b = a0; b <= a1; ++b) {              // loop with atom b
                    int bn  = (b + 1 - a0) % (ax + 1) + a0;   // next atom index
                    int bp  = (b + ax - a0) % (ax + 1) + a0;  // previous atom index
                    int ib  = i * N + b * Ndim;               // i-th bead b-th atom
                    int inb = (i < P - 1) ? ib + N : bn * Ndim;
                    int ipb = (i > 0) ? ib - N : PendN + bp * Ndim;
                    for (int d = 0, ibj = ib, inbj = inb, ipbj = ipb;  //
                         d < Ndim; ++d, ++ibj, ++inbj, ++ipbj) {
                        dE01[ibj] = nm[ibj] * bf2 * (2 * nrs[ibj] - nrs[ipbj] - nrs[inbj]);
                    }
                }
            }
        }
    }
    /**
     * @note calculate fV factor & dV_spring
     */
    for (int a = 0; a < N; ++a) {
        fV[a] = fE[a];
        for (int k = 0, kpa = Natom + a; k < a; ++k, kpa += Natom) fV[a] += fV[k] * fE[kpa];
        fV[a] /= (a + 1);

        num_real* dVa  = dV_spring + a * PN;
        num_real* dE0a = dE_spring + a * PN;
        for (int J = 0; J < PN; ++J) dVa[J] = dE0a[J] * fE[a];
        for (int k = 0, kpa = Natom + a; k < a; ++k, kpa += Natom) {
            num_real* dVk   = dV_spring + k * PN;
            num_real* dEkpa = dE_spring + kpa * PN;
            num_real fVfE   = fV[k] * fE[kpa];
            for (int J = 0; J < PN; ++J) dVa[J] += (dVk[J] + dEkpa[J]) * fVfE;
        }
        num_real wa = (a + 1) * fE[a];
        for (int J = 0; J < PN; ++J) dVa[J] /= wa;
    }
    return 0;
}

int MBPIMDTraj_Solver::update_p_harm(const num_real& dt_in) {  // used in BAOAB
    plFunction();
    spring_force();
    num_real* force_spring_mb = dV_spring + (Natom - 1) * PN;
    for (int i = 0, idx = 0; i < P; ++i) {
        for (int j = 0; j < N; ++j, ++idx) nps[idx] -= force_spring_mb[idx] * dt_in;
    }
    return 0;
}

int MBPIMDTraj_Solver::estimator(const int& isamp, TCFnucl& tcfer) {  // used in BAOAB
    num_real esti_V, esti_Kprim;
    // V
    esti_V = 0.0f;
    for (int i = 0; i < P; ++i) esti_V += vpeses[i];
    esti_V /= P;

    // Kprim
    esti_Kprim = 0.5f * N / betap;
    for (int a = 0; a < N; ++a) {
        VHO[a] = std::log(fE[a]) / pThermo->beta * fE[a];
        for (int k = 0, kpa = Natom + a; k < a; ++k, kpa += Natom) {
            VHO[a] += (VHO[k] + std::log(fE[kpa]) / pThermo->beta) * fV[k] * fE[kpa];
        }
        VHO[a] /= (a + 1) * fV[a];
    }
    esti_Kprim += VHO[Natom - 1];

    num_real sampunit = dt * sstep * iou.time;
    ofs_ENER << FMT(8) << itraj << FMT(8) << isamp * sampunit  //
             << FMT(8) << Ndim                                 //
             << FMT(8) << N                                    //
             << FMT(8) << P                                    //
             << FMT(8) << pThermo->beta                        //
             << FMT(8) << esti_V                               //
             << FMT(8) << esti_Kprim                           //
             << std::endl;

    for (int idof = 0; idof < PN; ++idof) ofs_TRAJ << FMT(8) << nrs[idof];
    for (int idof = 0; idof < PN; ++idof) ofs_TRAJ << FMT(8) << nfs[idof];
    ofs_TRAJ << std::endl;

    return 0;
}
