#include "solver_scf.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../thirdpart/Eigen/Dense"
#include "../utils/elements.h"

#ifdef USE_LIBCINT
#include "libcint/cint.h"
// #include "libcint/cint_funcs.h"
#endif  // USE_LIBCINT

#include "simple_integral_gto.h"

using namespace ARRAY_EG;

SCF_Solver::SCF_Solver(const Param& iparm, Model* pM) : Solver(iparm, pM) {
    plFunction();

    //  // add it in future, and in
    /**
     * @todo Add QMolecule_Model class in future
     *
     *  ```
     *      pMolecule = dynamic_cast<QMolecule_Model*>(pM);
     *  ```
     *  now read atoms directly
     *
     */

    Natom   = Param_GetT(int, parm, "Natom", 0);
    Q       = Param_GetT(int, parm, "Q", 0);
    S       = Param_GetT(int, parm, "S", 1);
    maxiter = Param_GetT(int, parm, "maxiter", 50);
    prec    = Param_GetT(double, parm, "prec", 1e-8);

    std::string xyzfile = Param_GetT(std::string, parm, "xyzfile", "example.xyz");
    std::string bsname  = Param_GetT(std::string, parm, "bsname", "6-311g");

    atoms = new Atom[Natom];

    int Natom_comp;
    std::string stmp;
    double dtmp;

    std::ifstream ifs(xyzfile);

    ifs >> Natom_comp;  ///< read Natom_comp in first line in xyzfile
    CHECK_EQ(Natom, Natom_comp);
    getline(ifs, stmp);  ///< read break of first line
    getline(ifs, stmp);  ///< read second line as comment in xyzfile

    Nelec = 0;
    for (int i = 0, idx = 0; i < Natom; ++i) {          ///< loop for reading atoms
        ifs >> stmp;                                    ///< reading labels
        atoms[i].znum = Elements::ELEMENTS_ZNUM(stmp);  // parse znum
        Nelec += atoms[i].znum;

        /** build atom_count information */
        if (atom_count.count(atoms[i].znum) == 0) {
            atom_count[atoms[i].znum] = 1;  // add this kind of atom
        } else {
            atom_count[atoms[i].znum] += 1;  // counting number++
        }

        ifs >> atoms[i].x >> atoms[i].y >> atoms[i].z;  ///< read cartesian coordinates

        atoms[i].x /= phys::au_2_ang;
        atoms[i].y /= phys::au_2_ang;
        atoms[i].z /= phys::au_2_ang;
    }
    Nelec -= Q;
    Nocc = (Nelec + 1) / 2;

    /* build basis_list */
    Nb = 0;
    for (auto i : atom_count) {  ///< i.first is znum, i.second is No. of atoms of znum
        int current_idx      = basis_list.size();
        atom_bslink[i.first] = current_idx;  ///< insert (key,val) bindings
        basis_list.push_back(Atomic_BasisSet(i.first, bsname));
        Nb += i.second * basis_list[current_idx].nb;
    }
    NNb  = Nb * Nb;
    NNNb = Nb * Nb * Nb;

    ppbasis_helper = new Atomic_BasisSet*[Natom];  ///< only allocate 2nd pointer, don't allocate 1st pointers
    for (int i = 0; i < Natom; ++i) {
        ppbasis_helper[i] = &basis_list[atom_bslink[atoms[i].znum]];  ///< point pointers to (global) basis_list
    }

    S_AO   = new double[NNb];        ///< overlap integral in AO basis
    H_AO   = new double[NNb];        ///< single electron integral in AO basis
    ERI_AO = new double[NNb * NNb];  ///< double electrons integral in AO basis

    matX = new double[NNb];
    matY = new double[NNb];

    G_AO = new double[NNb];
    ARRAY_CLEAR(G_AO, NNb);
    J_AO = new double[NNb];
    ARRAY_CLEAR(J_AO, NNb);
    K_AO = new double[NNb];
    ARRAY_CLEAR(K_AO, NNb);

    F_AO = new double[NNb];
    ARRAY_CLEAR(F_AO, NNb);
    C_AO = new double[NNb];
    ARRAY_CLEAR(C_AO, NNb);
    P_AO = new double[NNb];
    ARRAY_CLEAR(P_AO, NNb);
    P_AO_last = new double[NNb];
    ARRAY_CLEAR(P_AO_last, NNb);

    F_MO = new double[NNb];
    ARRAY_CLEAR(F_MO, NNb);
    C_MO = new double[NNb];
    ARRAY_CLEAR(C_MO, NNb);
    P_MO = new double[NNb];
    ARRAY_CLEAR(P_MO, NNb);
    E_MO = new double[Nb];
    ARRAY_CLEAR(E_MO, Nb);

    P_NO = new double[NNb];
    ARRAY_CLEAR(P_NO, NNb);
}

SCF_Solver::~SCF_Solver() {
    delete[] atoms;
    delete[] S_AO, delete[] H_AO, delete[] ERI_AO;
    delete[] matX, delete[] matY;
    delete[] G_AO, delete[] J_AO, delete[] K_AO;
    delete[] F_AO, delete[] C_AO, delete[] P_AO, delete[] P_AO_last;
    delete[] F_MO, delete[] C_MO, delete[] P_MO, delete[] E_MO;
    delete[] P_NO;
    delete[] ppbasis_helper;
}

int SCF_Solver::print_info() {
    plFunction();
    for (int a1 = 0, B1 = 0; a1 < Natom; ++a1) {
        Atom& A1              = atoms[a1];            ///< atom1
        Atomic_BasisSet& ABS1 = *ppbasis_helper[a1];  ///< Atomic_BasisSet of atom1
        /* Get information in Atomic_BasisSet */
        int& nb1        = ABS1.nb;
        double* alphas1 = ABS1.alphas;
        double* coeffs1 = ABS1.coeffs;
        Qnum* qnums1    = ABS1.qnums;
        int* conts1     = ABS1.nconts;

        LOG(INFO) << "Print for " << a1 << "-th Atom: " << Elements::ELEMENTS_NAME(A1.znum);
        LOG(INFO) << "\tCoordinate:" << FMT(8) << A1.x << ", " << FMT(8) << A1.y << ", " << FMT(8) << A1.z;
        LOG(INFO) << "\tNo. of basis = " << nb1;

        for (int b1 = 0; b1 < nb1; ++b1, ++B1) {
            LOG(INFO) << "\tPrint for " << b1 << "-th Basis for " << a1 << "-th Atom: (total " << B1 << "/" << Nb
                      << ")";
            LOG(INFO) << "\t\tQnum:" << FMT(1) << qnums1[b1].L << ", " << FMT(1) << qnums1[b1].M << ", " << FMT(1)
                      << qnums1[b1].N;
            for (int i = 0; i < conts1[b1]; ++i) {
                LOG(INFO) << "\t\t\talpha = " << FMT(8) << alphas1[i] << ", coeffs = " << FMT(8) << coeffs1[i];
            }
            alphas1 += conts1[b1];
            coeffs1 += conts1[b1];
        }
    }
    return 0;
}

int SCF_Solver::call_ao_integral() {
    Eigen::Map<EigMXr> MapS_AO(S_AO, Nb, Nb);
    Eigen::Map<EigMXr> MapH_AO(H_AO, Nb, Nb);
    Eigen::Map<EigMXr> MapERI_AO(ERI_AO, NNb, NNb);

#ifdef USE_LIBCINT
    call_ao_integral_libcint();
#else
    call_ao_integral_simple();
#endif  // USE_LIBCINT

    LOG(WARNING) << "S\n" << MapS_AO;
    LOG(WARNING) << "H\n" << MapH_AO;
    LOG(WARNING) << "G\n" << MapERI_AO;
    return 0;
}

#ifdef USE_LIBCINT
int SCF_Solver::call_ao_integral_libcint() {
    int natm = Natom;
    int nbas = 0;
    for (int a1 = 0; a1 < Natom; ++a1) {
        Atomic_BasisSet& ABS1 = *ppbasis_helper[a1];  ///< Atomic_BasisSet of atom1
        for (int ib = 0; ib < ABS1.rdata_angl.size(); ++ib) nbas += ABS1.rdata_angl[ib].size();
    }

    // ATM_SLOTS = 6; BAS_SLOTS = 8;
    int* atm    = (int*) malloc(sizeof(int) * natm * ATM_SLOTS);
    int* bas    = (int*) malloc(sizeof(int) * nbas * BAS_SLOTS);
    double* env = (double*) malloc(sizeof(double) * 10000);

    int* bas_indx = (int*) malloc(sizeof(int) * nbas);

    int off = PTR_ENV_START;  // = 20
    for (int iatom = 0, idxR = 0; iatom < Natom; ++iatom) {
        atm[CHARGE_OF + ATM_SLOTS * iatom] = atoms[iatom].znum;
        atm[PTR_COORD + ATM_SLOTS * iatom] = off;
        env[off++]                         = atoms[iatom].x;
        env[off++]                         = atoms[iatom].y;
        env[off++]                         = atoms[iatom].z;
    }

    for (int a1 = 0, ibas = 0; a1 < Natom; ++a1) {
        Atomic_BasisSet& ABS1 = *ppbasis_helper[a1];  ///< Atomic_BasisSet of atom1
        for (int ib = 0; ib < ABS1.rdata_angl.size(); ++ib) {
            int nrow = ABS1.rdata_nrow[ib], ncol = ABS1.rdata_ncol[ib];

            LOG(WARNING) << nrow << "; " << ncol;

            double* data = ABS1.rdata.data() + ABS1.rdata_posi[ib];
            if (ABS1.rdata_angl[ib].size() == 1) {
                bas[ATOM_OF + BAS_SLOTS * ibas]  = a1;
                bas[ANG_OF + BAS_SLOTS * ibas]   = angular_map.at(ABS1.rdata_angl[ib][0]);
                bas[NPRIM_OF + BAS_SLOTS * ibas] = nrow;
                bas[NCTR_OF + BAS_SLOTS * ibas]  = ncol - 1;

                bas[PTR_EXP + BAS_SLOTS * ibas] = off;
                for (int i = 0; i < nrow; ++i, ++off) env[off] = data[i * ncol + 0];

                bas[PTR_COEFF + BAS_SLOTS * ibas] = off;
                for (int j = 1; j < ncol; ++j)
                    for (int i = 0; i < nrow; ++i, ++off) {
                        env[off] = data[i * ncol + j] * CINTgto_norm(bas[ANG_OF + BAS_SLOTS * ibas],
                                                                     env[bas[PTR_EXP + BAS_SLOTS * ibas] + i]);
                    }
                ibas++;
            } else if (ABS1.rdata_angl[ib].size() == ncol - 1) {
                for (int ich = 0; ich < ABS1.rdata_angl[ib].size(); ++ich) {
                    bas[ATOM_OF + BAS_SLOTS * ibas]  = a1;
                    bas[ANG_OF + BAS_SLOTS * ibas]   = angular_map.at(ABS1.rdata_angl[ib][ich]);
                    bas[NPRIM_OF + BAS_SLOTS * ibas] = nrow;
                    bas[NCTR_OF + BAS_SLOTS * ibas]  = 1;

                    bas[PTR_EXP + BAS_SLOTS * ibas] = off;
                    for (int i = 0; i < nrow; ++i, ++off) env[off] = data[i * ncol + 0];

                    bas[PTR_COEFF + BAS_SLOTS * ibas] = off;
                    for (int i = 0; i < nrow; ++i, ++off)
                        env[off] = data[i * ncol + ich] *
                                   CINTgto_norm(bas[ANG_OF + BAS_SLOTS * ibas], env[bas[PTR_EXP + BAS_SLOTS * ibas]]);
                    ibas++;
                }
            } else {
                LOG(FATAL);
            }
        }
    }

    int shls[4];
    double* buf;
    // indexing and normalization
    for (int i = 0, cnt = 0; i < nbas; ++i) {
        bas_indx[i] = cnt;
        int di      = CINTcgto_cart(i, bas);
        shls[0] = i, shls[1] = i;

        buf = (double*) malloc(sizeof(double) * di * di);
        cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env);

        LOG(WARNING) << di << " " << di << ": " << buf[0];

        // not suitable for DZ-basis
        for (int ip = 0; ip < bas[NPRIM_OF + BAS_SLOTS * i]; ++ip) {
            env[bas[PTR_COEFF + BAS_SLOTS * i] + ip] /= sqrt(buf[0]);
        }

        free(buf);
        cnt += di;
    }

    CINTOpt* opt = NULL;
    cint2e_cart_optimizer(&opt, atm, natm, bas, nbas, env);
    for (int i = 0; i < nbas; ++i) {
        int di = CINTcgto_cart(i, bas);
        int i0 = bas_indx[i];
        for (int j = 0; j < nbas; ++j) {
            int dj = CINTcgto_cart(j, bas);
            int j0 = bas_indx[j];

            shls[0] = i, shls[1] = j;

            buf = (double*) malloc(sizeof(double) * di * dj);

            // overlap of i-bas & j-bas
            cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env);
            for (int ii = 0, ix = bas_indx[i], ibuf = 0; ii < di; ++ii, ++ix) {
                for (int jj = 0, jx = bas_indx[j]; jj < dj; ++jj, ++jx, ++ibuf) S_AO[ix + jx * Nb] = buf[ibuf];
            }

            // kinetic part
            cint1e_kin_cart(buf, shls, atm, natm, bas, nbas, env);
            for (int ii = 0, ix = bas_indx[i], ibuf = 0; ii < di; ++ii, ++ix) {
                for (int jj = 0, jx = bas_indx[j]; jj < dj; ++jj, ++jx, ++ibuf) H_AO[ix + jx * Nb] = buf[ibuf];
            }

            // Vnuclear part
            cint1e_nuc_cart(buf, shls, atm, natm, bas, nbas, env);
            for (int ii = 0, ix = bas_indx[i], ibuf = 0; ii < di; ++ii, ++ix) {
                for (int jj = 0, jx = bas_indx[j]; jj < dj; ++jj, ++jx, ++ibuf) H_AO[ix + jx * Nb] += buf[ibuf];
            }
            free(buf);

            for (int k = 0; k < nbas; ++k) {
                int dk = CINTcgto_cart(k, bas);
                int k0 = bas_indx[k];
                for (int l = 0; l < nbas; ++l) {
                    int dl = CINTcgto_cart(l, bas);
                    int l0 = bas_indx[l];

                    shls[2] = k, shls[3] = l;
                    buf = (double*) malloc(sizeof(double) * di * dj * dk * dl);
                    cint2e_cart(buf, shls, atm, natm, bas, nbas, env, opt);
                    for (int ii = 0, ix = bas_indx[i], ibuf = 0; ii < di; ++ii, ++ix) {
                        for (int jj = 0, jx = bas_indx[j]; jj < dj; ++jj, ++jx) {
                            for (int kk = 0, kx = bas_indx[k]; kk < dk; ++kk, ++kx) {
                                for (int ll = 0, lx = bas_indx[l]; ll < dl; ++ll, ++lx, ++ibuf) {
                                    ERI_AO[ix + jx * Nb + kx * NNb + lx * NNNb] = buf[ibuf];
                                }
                            }
                        }
                    }
                    free(buf);
                }
            }
        }
    }
    CINTdel_optimizer(&opt);
    free(atm);
    free(bas);
    free(env);
    free(bas_indx);
    return 0;
}
#endif  // USE_LIBCINT

int SCF_Solver::call_ao_integral_simple() {  ///< a simple realization of call_ao_integral() based on
                                             ///< simple_integral_gto.h
    plFunction();

    for (int a1 = 0, B1 = 0; a1 < Natom; ++a1) {
        Atom& A1              = atoms[a1];            ///< atom1
        Atomic_BasisSet& ABS1 = *ppbasis_helper[a1];  ///< Atomic_BasisSet of atom1
        /* Get information in Atomic_BasisSet */
        int& nb1        = ABS1.nb;
        double* alphas1 = ABS1.alphas;
        double* coeffs1 = ABS1.coeffs;
        Qnum* qnums1    = ABS1.qnums;
        int* conts1     = ABS1.nconts;

        for (int b1 = 0; b1 < nb1; ++b1, ++B1) {  ///< loop in Atomic_BasisSet of atom1
            int& nc1 = conts1[b1];
            Qnum& q1 = qnums1[b1];

            for (int a2 = 0, B2 = 0; a2 < Natom; ++a2) {  ///< Atomic_BasisSet of atom2
                Atom& A2              = atoms[a2];
                Atomic_BasisSet& ABS2 = *ppbasis_helper[a2];  ///< Atomic_BasisSet of atom1

                int& nb2        = ABS2.nb;
                double* alphas2 = ABS2.alphas;
                double* coeffs2 = ABS2.coeffs;
                Qnum* qnums2    = ABS2.qnums;
                int* conts2     = ABS2.nconts;
                for (int b2 = 0; b2 < nb2; ++b2, ++B2) {  ///< loop in Atomic_BasisSet of atom2
                    int& nc2 = conts2[b2];
                    Qnum& q2 = qnums2[b2];

                    if (B1 <= B2) {  // only do integral for B1 <= B2, otherwise it copy values
                        // LOG(WARNING) << "S-int: " << B1 << ", " << B2;

                        S_AO[B1 * Nb + B2] = IntecGTO_S(  ///< S-integral
                            coeffs1, alphas1, A1.x, A1.y, A1.z, q1.L, q1.M, q1.N, nc1, coeffs2, alphas2, A2.x, A2.y,
                            A2.z, q2.L, q2.M, q2.N, nc2);

                        H_AO[B1 * Nb + B2] = IntecGTO_T(  ///< T-integral in single-part integral
                            coeffs1, alphas1, A1.x, A1.y, A1.z, q1.L, q1.M, q1.N, nc1, coeffs2, alphas2, A2.x, A2.y,
                            A2.z, q2.L, q2.M, q2.N, nc2);
                    } else { /* copy values when B1 > B2 */
                        S_AO[B1 * Nb + B2] = S_AO[B2 * Nb + B1];
                        H_AO[B1 * Nb + B2] = H_AO[B2 * Nb + B1];
                    }

                    for (int a3 = 0, B3 = 0; a3 < Natom; ++a3) {  ///< loop in Atomic_BasisSet of atom3
                        Atom& A3              = atoms[a3];
                        Atomic_BasisSet& ABS3 = *ppbasis_helper[a3];  ///< Atomic_BasisSet of atom1

                        int& nb3        = ABS3.nb;
                        double* alphas3 = ABS3.alphas;
                        double* coeffs3 = ABS3.coeffs;
                        Qnum* qnums3    = ABS3.qnums;
                        int* conts3     = ABS3.nconts;

                        if (B1 <= B2) {
                            H_AO[B1 * Nb + B2] += IntecGTO_V(  ///< V-integral in single-part integral
                                coeffs1, alphas1, A1.x, A1.y, A1.z, q1.L, q1.M, q1.N, nc1, coeffs2, alphas2, A2.x, A2.y,
                                A2.z, q2.L, q2.M, q2.N, nc2, A3.x, A3.y, A3.z, A3.znum);
                        }

                        for (int b3 = 0; b3 < nb3; ++b3, ++B3) {  ///< loop in Atomic_BasisSet of atom3
                            int& nc3 = conts3[b3];
                            Qnum& q3 = qnums3[b3];

                            for (int a4 = 0, B4 = 0; a4 < Natom; ++a4) {  ///< loop in Atomic_BasisSet of atom4
                                Atom& A4              = atoms[a4];
                                Atomic_BasisSet& ABS4 = *ppbasis_helper[a4];  ///< Atomic_BasisSet of atom1

                                int nb4         = ABS4.nb;
                                double* alphas4 = ABS4.alphas;
                                double* coeffs4 = ABS4.coeffs;
                                Qnum* qnums4    = ABS4.qnums;
                                int* conts4     = ABS4.nconts;
                                for (int b4 = 0; b4 < nb4; ++b4, ++B4) {  ///< loop in Atomic_BasisSet of atom4
                                    int& nc4 = conts4[b4];
                                    Qnum& q4 = qnums4[b4];

                                    if (B1 <= B2 && B3 <= B4 && B2 * (B2 + 1) + 2 * B1 <= B4 * (B4 + 1) + 2 * B3) {
                                        ERI_AO[B1 * NNNb + B2 * NNb + B3 * Nb + B4] =
                                            IntecGTO_ERI(coeffs1, alphas1, A1.x, A1.y, A1.z, q1.L, q1.M, q1.N, nc1,
                                                         coeffs2, alphas2, A2.x, A2.y, A2.z, q2.L, q2.M, q2.N, nc2,
                                                         coeffs3, alphas3, A3.x, A3.y, A3.z, q3.L, q3.M, q3.N, nc3,
                                                         coeffs4, alphas4, A4.x, A4.y, A4.z, q4.L, q4.M, q4.N, nc4);

                                        // LOG(WARNING) << "(" << B1 << "," << B2 << "," << B3 << "," << B4 << "):   "
                                        // <<  ERI_AO[B1*NNNb + B2*NNb + B3*Nb + B4] << std::endl;
                                        // ERI_AO copies
                                        ERI_AO[B2 * NNNb + B1 * NNb + B3 * Nb + B4] =
                                            ERI_AO[B1 * NNNb + B2 * NNb + B3 * Nb + B4];
                                        ERI_AO[B1 * NNNb + B2 * NNb + B4 * Nb + B3] =
                                            ERI_AO[B1 * NNNb + B2 * NNb + B3 * Nb + B4];
                                        ERI_AO[B2 * NNNb + B1 * NNb + B4 * Nb + B3] =
                                            ERI_AO[B1 * NNNb + B2 * NNb + B3 * Nb + B4];
                                        ERI_AO[B3 * NNNb + B4 * NNb + B1 * Nb + B2] =
                                            ERI_AO[B1 * NNNb + B2 * NNb + B3 * Nb + B4];
                                        ERI_AO[B3 * NNNb + B4 * NNb + B2 * Nb + B1] =
                                            ERI_AO[B1 * NNNb + B2 * NNb + B3 * Nb + B4];
                                        ERI_AO[B4 * NNNb + B3 * NNb + B1 * Nb + B2] =
                                            ERI_AO[B1 * NNNb + B2 * NNb + B3 * Nb + B4];
                                        ERI_AO[B4 * NNNb + B3 * NNb + B2 * Nb + B1] =
                                            ERI_AO[B1 * NNNb + B2 * NNb + B3 * Nb + B4];
                                    }
                                    alphas4 += nc4;  ///< shift of pointer
                                    coeffs4 += nc4;  ///< shift of pointer
                                }
                            }
                            alphas3 += nc3;  ///< shift of pointer
                            coeffs3 += nc3;  ///< shift of pointer
                        }
                    }
                    alphas2 += nc2;  ///< shift of pointer
                    coeffs2 += nc2;  ///< shift of pointer
                }
            }
            alphas1 += nc1;  ///< shift of pointer
            coeffs1 += nc1;  ///< shift of pointer
        }
    }

    return 0;
}

int SCF_Solver::run_impl() {
    plFunction();

    Eigen::Map<EigMXr> MapP_AO(P_AO, Nb, Nb);
    Eigen::Map<EigMXr> MapP_AO_last(P_AO_last, Nb, Nb);

    print_info();
    call_ao_integral();
    solve_orthonormalize();
    initial_guess();

    bool converged = false;
    int cnt        = 1;
    while (!converged && cnt++ <= maxiter) {
        iter(cnt);
        analysis(cnt);

        /* check convergence */
        converged    = (std::abs(E - E_last) < prec) && ((MapP_AO - MapP_AO_last).norm() < prec);
        E_last       = E;
        MapP_AO_last = MapP_AO;
    }
    return 0;
}

int SCF_Solver::iter(const int& cnt) { return iter_simple(cnt); }

int SCF_Solver::iter_simple(const int& cnt) {  ///< a simple realization of iter()
    build_Fock();                              // from P build F
    solve_Fock();                              // solve eigen problem for F
    sync_mo2ao();                              // from solution rebuild P
    return 0;
}

int SCF_Solver::iter_diis(const int& cnt) {  ///< @todo
    LOG(FATAL);
    return 0;
}

int SCF_Solver::solve_orthonormalize() {
    plFunction();

    Eigen::Map<EigMXr> MapS_AO(S_AO, Nb, Nb);
    Eigen::Map<EigMXr> MapX(matX, Nb, Nb);
    Eigen::Map<EigMXr> MapY(matY, Nb, Nb);
    Eigen::SelfAdjointEigenSolver<EigMXr> eig(MapS_AO);

    auto e = eig.eigenvalues();
    auto T = eig.eigenvectors();
    auto L = (e.array()).sqrt();

    MapX = T * (1 / L).matrix().asDiagonal() * T.transpose();
    MapY = T * L.matrix().asDiagonal() * T.transpose();

    return 0;
}

int SCF_Solver::initial_guess() {
    plFunction();

    // from H we guess a C
    Eigen::Map<EigMXr> MapP_NO(P_NO, Nb, Nb);
    Eigen::Map<EigMXr> MapF_MO(F_MO, Nb, Nb);
    Eigen::Map<EigMXr> MapH_AO(H_AO, Nb, Nb);
    Eigen::Map<EigMXr> MapX(matX, Nb, Nb);

    /* initalize P_NO: only once!!! */
    MapP_NO = EigMXr::Zero(Nb, Nb);
    for (int i = 0; i < Nocc; ++i) MapP_NO(i, i) = 1;

    /* guess initial Fock problem */
    MapF_MO = MapX.transpose() * MapH_AO * MapX;  ///< instead of build_Fock
    solve_Fock();

    /* further guess initial P_MO & P_AO */
    sync_mo2ao();
    analysis(0);
    E_last = E;
    return 0;
}

int SCF_Solver::build_Fock() {  ///< @todo to be optimized
    plFunction();

    Eigen::Map<EigMXr> MapG_AO(G_AO, Nb, Nb);
    Eigen::Map<EigMXr> MapJ_AO(J_AO, Nb, Nb);
    Eigen::Map<EigMXr> MapK_AO(K_AO, Nb, Nb);
    Eigen::Map<EigMXr> MapP_AO(P_AO, Nb, Nb);
    Eigen::Map<EigMXr> MapX(matX, Nb, Nb);
    Eigen::Map<EigMXr> MapH_AO(H_AO, Nb, Nb);
    Eigen::Map<EigMXr> MapF_AO(F_AO, Nb, Nb);
    Eigen::Map<EigMXr> MapF_MO(F_MO, Nb, Nb);

    for (int B1 = 0; B1 < Nb; ++B1) {
        for (int B2 = 0; B2 < Nb; ++B2) {
            MapF_AO(B1, B2) = MapH_AO(B1, B2);
            MapG_AO(B1, B2) = 0.0f;
            MapJ_AO(B1, B2) = 0.0f;
            MapK_AO(B1, B2) = 0.0f;
            for (int B3 = 0; B3 < Nb; ++B3) {
                for (int B4 = 0; B4 < Nb; ++B4) {
                    MapJ_AO(B1, B2) += ERI_AO[B1 * NNNb + B2 * NNb + B3 * Nb + B4] * MapP_AO(B4, B3);
                    MapK_AO(B1, B2) += ERI_AO[B1 * NNNb + B4 * NNb + B3 * Nb + B2] * MapP_AO(B4, B3);
                    MapG_AO(B1, B2) += (2 * ERI_AO[B1 * NNNb + B2 * NNb + B3 * Nb + B4] -
                                        ERI_AO[B1 * NNNb + B4 * NNb + B3 * Nb + B2]) *
                                       MapP_AO(B4, B3);
                }
            }
            MapF_AO(B1, B2) += MapG_AO(B1, B2);
        }
    }
    MapF_MO = MapX.transpose() * MapF_AO * MapX;
    return 0;
}

int SCF_Solver::solve_Fock() {
    plFunction();
    Eigen::Map<EigMXr> MapE_MO(E_MO, Nb, 1);
    Eigen::Map<EigMXr> MapC_MO(C_MO, Nb, Nb);
    Eigen::Map<EigMXr> MapF_MO(F_MO, Nb, Nb);
    Eigen::SelfAdjointEigenSolver<EigMXr> eig(MapF_MO);

    MapE_MO = eig.eigenvalues();
    MapC_MO = eig.eigenvectors();

    return 0;
}

int SCF_Solver::sync_mo2ao() {
    plFunction();

    Eigen::Map<EigMXr> MapX(matX, Nb, Nb);
    Eigen::Map<EigMXr> MapC_MO(C_MO, Nb, Nb);
    Eigen::Map<EigMXr> MapC_AO(C_AO, Nb, Nb);
    Eigen::Map<EigMXr> MapP_NO(P_NO, Nb, Nb);
    Eigen::Map<EigMXr> MapP_MO(P_MO, Nb, Nb);
    Eigen::Map<EigMXr> MapP_AO(P_AO, Nb, Nb);

    MapC_AO = MapX * MapC_MO;
    MapP_MO = MapC_MO * MapP_NO * MapC_MO.transpose();
    MapP_AO = MapC_AO * MapP_NO * MapC_AO.transpose();
    return 0;
}

int SCF_Solver::analysis(const int& cnt) {
    Eigen::Map<EigMXr> MapP_AO(P_AO, Nb, Nb);
    Eigen::Map<EigMXr> MapH_AO(H_AO, Nb, Nb);

    E = E_nucl_repul();  ///< get nuclear repulsion energy

    /* count once single body interaction & once double body interaction */
    for (int k = 0; k < Nocc; ++k) E += E_MO[k];
    /* count once single body interaction again (for each MO has 2-electrons) */
    E += (MapP_AO * MapH_AO).trace();

    scf_report(cnt);
    return 0;
}

int SCF_Solver::scf_report(const int& cnt) {
    Eigen::Map<EigMXr> MapC_MO(C_MO, Nb, Nb);
    Eigen::Map<EigMXr> MapC_AO(C_AO, Nb, Nb);
    Eigen::Map<EigMXr> MapE_MO(E_MO, 1, Nb);

    LOG(INFO) << "####################" << FMT(1) << cnt << "  ###################";
    LOG(INFO) << "E_MO" << std::endl << MapE_MO;
    LOG(INFO) << "C_AO" << std::endl << MapC_AO;
    LOG(INFO) << "Energy " << FMT(8) << E << "    "
              << "Convergence: " << FMT(8) << E_last - E;
    return 0;
}

double SCF_Solver::E_nucl_repul() {
    plFunction();

    double sum = 0;
    for (int a1 = 0; a1 < Natom; ++a1) {
        Atom& A1 = atoms[a1];
        for (int a2 = a1 + 1; a2 < Natom; ++a2) {
            Atom& A2 = atoms[a2];
            sum += A1.znum * A2.znum /
                   sqrt((A1.x - A2.x) * (A1.x - A2.x) + (A1.y - A2.y) * (A1.y - A2.y) + (A1.z - A2.z) * (A1.z - A2.z));
        }
    }
    return sum;
}
