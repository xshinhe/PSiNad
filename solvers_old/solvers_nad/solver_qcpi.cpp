#include "solver_qcpi.h"

#include "../utils/definitions.h"
#include "../utils/nad_utils.h"

QCPI_Solver::QCPI_Solver(Param iparm, Model* pM) : NadTraj_Solver(iparm, pM) {
    mem   = Param_GetT(int, parm, "mem", 0);
    niter = Param_GetT(int, parm, "niter", 10);

    std::string suffix = std::to_string(mem);
    save               = save + "_" + QCPI_Solver::name() + suffix + "_" + pForceField->tag;

    Lx = mem + 1;
    Ly = Lx + 1;
    Nx = (int) pow(FF, Lx);
    Ny = (int) pow(FF, Ly);

    LOG(WARNING) << "Using Ntraj = " << Ny;

    try {
        ALLOCATE_PTR_TO_VECTOR(idx_arr, Ly);
        ALLOCATE_PTR_TO_VECTOR(nrs, Ny * N);
        ALLOCATE_PTR_TO_VECTOR(nps, Ny * N);
        ALLOCATE_PTR_TO_VECTOR(nrs_copy, Ny * N);
        ALLOCATE_PTR_TO_VECTOR(nps_copy, Ny * N);

        ALLOCATE_PTR_TO_VECTOR(nx, Nx);
        ALLOCATE_PTR_TO_VECTOR(nx_copy, Nx);
        ALLOCATE_PTR_TO_VECTOR(ny, Ny);
        ALLOCATE_PTR_TO_VECTOR(ny_copy, Ny);
    } catch (std::bad_alloc& e) { LOG(FATAL) << e.what(); }

    save = name() + "_" + pForceField->tag;
}

QCPI_Solver::~QCPI_Solver(){};

int QCPI_Solver::propagate_tensor() {
    plFunction();

    // first step: calculate Y tensor
    memset(ny_copy, 0, Ny * sizeof(num_complex));
    if (Ly_now < Ly) {  // exactly calculate Y tensor
        plScope("expay");

        // declaration of sizes
        int L_prev = Ly_now;
        int L_next = L_prev + 1;
        int N_prev = (int) pow(FF, L_prev);
        int N_next = (int) pow(FF, L_next);

        for (int i_prev = 0; i_prev < N_prev; ++i_prev) {
            num_real* nr_prev = nrs + i_prev * N;
            num_real* np_prev = nps + i_prev * N;
            get_index(idx_arr, i_prev, FF, L_prev);

            int prev = idx_arr[L_prev - 1];
            for (int next = 0; next < FF; ++next) {
                idx_arr[L_next - 1] = next;
                int i_next          = get_Nnum(idx_arr, FF, L_next);
                num_real* nr_next   = nrs_copy + i_next * N;
                num_real* np_next   = nps_copy + i_next * N;
                ny_copy[i_next]     = propagate_nuc(nr_next, np_next, nr_prev, np_prev, prev, next, dt);
            }
        }
        Ly_now++;
    } else {  // propagate Y tensor by random choice of memory tails
        plScope("propy");

        for (int i_next = 0; i_next < Ny; ++i_next) {
            get_index(idx_arr, i_next, FF, Ly);
            int next = idx_arr[Ly - 1];

            for (int i = Ly - 1; i > 0; --i) idx_arr[i] = idx_arr[i - 1];
            rand_catalog(idx_arr);
            int i_prev = get_Nnum(idx_arr, FF, Ly);
            int prev   = idx_arr[Ly - 1];

            num_real* nr_next = nrs_copy + i_next * N;
            num_real* np_next = nps_copy + i_next * N;
            num_real* nr_prev = nrs + i_prev * N;
            num_real* np_prev = nps + i_prev * N;
            ny_copy[i_next]   = propagate_nuc(nr_next, np_next, nr_prev, np_prev, prev, next, dt);
        }
    }
    for (int iy = 0, J = 0; iy < Ny; ++iy) ny[iy] = ny_copy[iy];
    for (int iy = 0, J = 0; iy < Ny; ++iy) {
        for (int j = 0; j < N; ++j, ++J) { nrs[J] = nrs_copy[J], nps[J] = nps_copy[J]; }
    }

    // second step: calculate X tensor
    memset(nx_copy, 0, Nx * sizeof(num_complex));
    if (Lx_now < Lx) {  // copy expand X
        plScope("expax");

        // declaration of sizes
        int L_prev = Lx_now;
        int L_next = L_prev + 1;
        int N_prev = (int) pow(FF, L_prev);
        int N_next = (int) pow(FF, L_next);

        for (int i_prev = 0; i_prev < N_prev; ++i_prev) {
            get_index(idx_arr, i_prev, FF, L_prev);
            int prev = idx_arr[L_prev - 1];
            for (int next = 0; next < FF; ++next) {
                idx_arr[L_next - 1] = next;
                int i_next          = get_Nnum(idx_arr, FF, L_next);
                nx_copy[i_next]     = ny[i_next] * nx[i_prev];
            }
        }
        Lx_now++;
    } else {  // propagate X
        plScope("propx");

        int* xidx_arr = idx_arr + 1;
        for (int i_next = 0; i_next < Nx; ++i_next) {
            get_index(xidx_arr, i_next, FF, Lx);
            for (int iend = 0; iend < FF; ++iend) {
                idx_arr[0] = iend;
                int i_prev = get_Nnum(idx_arr, FF, Lx);
                int iy     = get_Nnum(idx_arr, FF, Ly);
                nx_copy[i_next] += nx[i_prev] * ny[iy];
            }
        }
    }
    for (int ix = 0; ix < Nx; ++ix) nx[ix] = nx_copy[ix];

    return 0;
}

num_complex QCPI_Solver::propagate_nuc(num_real* nr2_trj, num_real* np2_trj,  // updated traj
                                       num_real* nr1_trj, num_real* np1_trj,  // old traj
                                       const int& prev, const int& next, const num_real& dt) {
    plFunction();

    num_real ddt = dt / (2 * niter);
    int fp1 = prev % F, fp2 = prev / F;
    int fn1 = next % F, fn2 = next / F;
    int Fadd1 = F + 1;
    num_real* ForceMat;
    switch (rep_type) {
        case representation::diabatic:
            ForceMat = dV;  // using diabatic force matrix
            break;
        case representation::adiabatic:
        case representation::onthefly:
            ForceMat = dE;  // using adiabatic force matrix
            break;
    }

    // copy to rho, nr, np & call Parent NadTraj_Solver's ff_calc1
    for (int i = 0; i < FF; ++i) { rho[i] = (i == prev) ? phys::math::iu : phys::math::iz; }
    for (int j = 0; j < N; ++j) nr[j] = nr1_trj[j], np[j] = np1_trj[j];

    // clang-format off
    ff_calc1(level);
    for (int j = 0, idxdV1 = fp1 * Fadd1, idxdV2 = fp2 * Fadd1;
        j < N; ++j, idxdV1 += FF, idxdV2 += FF) {
        nf[j] = grad[j] + 0.5f * (ForceMat[idxdV1] + ForceMat[idxdV2]);
    }
    // clang-format on
    for (int i = 0; i < niter; ++i) {
        plScope("1");
        for (int j = 0; j < N; ++j) np[j] -= nf[j] * 0.5f * ddt;
        for (int j = 0; j < N; ++j) nr[j] += np[j] / nm[j] * ddt;
        ff_calc1(level);

        solve_Ut(U, S, L, T, E, ddt, rep_type, N, F, workr, workc);
        update_rho(rho, U, N, F, workr, workc);

        // clang-format off
        for (int j = 0, idxdV1 = fp1 * Fadd1, idxdV2 = fp2 * Fadd1;
            j < N; ++j, idxdV1 += FF, idxdV2 += FF) {
            nf[j] = grad[j] + 0.5f * (ForceMat[idxdV1] + ForceMat[idxdV2]);
        }
        // clang-format on
        for (int j = 0; j < N; ++j) np[j] -= nf[j] * 0.5f * ddt;
    }
    // clang-format off
    for (int j = 0, idxdV1 = fn1 * Fadd1, idxdV2 = fn2 * Fadd1;
        j < N; ++j, idxdV1 += FF, idxdV2 += FF) {
        nf[j] = grad[j] + 0.5f * (ForceMat[idxdV1] + ForceMat[idxdV2]);
    }
    // clang-format on
    for (int i = 0; i < niter; ++i) {
        plScope("2");

        for (int j = 0; j < N; ++j) np[j] -= nf[j] * 0.5f * ddt;
        for (int j = 0; j < N; ++j) nr[j] += np[j] / nm[j] * ddt;
        ff_calc1(level);

        solve_Ut(U, S, L, T, E, ddt, rep_type, N, F, workr, workc);
        update_rho(rho, U, N, F, workr, workc);

        // clang-format off
        for (int j = 0, idxdV1 = fn1 * Fadd1, idxdV2 = fn2 * Fadd1;
            j < N; ++j, idxdV1 += FF, idxdV2 += FF) {
            nf[j] = grad[j] + 0.5f * (ForceMat[idxdV1] + ForceMat[idxdV2]);
        }
        // clang-format on
        for (int j = 0; j < N; ++j) np[j] -= nf[j] * 0.5f * ddt;
    }

    // copy back
    for (int j = 0; j < N; ++j) nr2_trj[j] = nr[j], np2_trj[j] = np[j];
    return rho[next];
};

int QCPI_Solver::init(const int& itraj) {
    if (itraj < 0) {
        int tmpi;
        rand_catalog(&tmpi, 1, true, 0, FF - 1);
    }

    NadTraj_Solver::init(itraj);
    memset(nx, 0, Nx * sizeof(num_complex));
    memset(ny, 0, Ny * sizeof(num_complex));

    Lx_now = 1, Ly_now = 1;
    for (int i = 0, J = 0; i < FF; ++i) {
        nx[i] = (i % F == occ0 && i / F == occ0) ? phys::math::iu : phys::math::iz;
        for (int j = 0; j < N; ++j, ++J) nrs[J] = nr[j], nps[J] = np[j];
    }
    return 0;
};

int QCPI_Solver::kernel0(num_complex* rhox, const int& F) {
    memset(rhox, 0, FF * sizeof(num_complex));
    int Nx_now = (int) pow(FF, Lx_now);
    for (int ix = 0; ix < Nx_now; ++ix) {
        get_index(idx_arr, ix, FF, Lx);
        rhox[idx_arr[Lx_now - 1]] += nx[ix];
    }
    return 0;
};
int QCPI_Solver::kernelt(num_complex* rhox, const int& F) { return kernel0(rhox, F); }

int QCPI_Solver::traj(NAD_TCFer& tcfer, const int& N, const int& F) {
    plFunction();
    init(itraj);
    int succ = sampler(0, tcfer);
    for (int istep_dummy = 0; istep_dummy < nstep; ++istep_dummy) {
        propagate_tensor();
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
