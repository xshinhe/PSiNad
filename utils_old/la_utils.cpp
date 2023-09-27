#include "la_utils.h"

#ifdef ARRAY_USE_MKL

void EigenSolve(num_real* E, num_real* T, num_real* A, const int& N) {
    plFunction();
    plScope("MKL");
    int NN = N * N;
    for (int i = 0; i < NN; ++i) T[i] = A[i];
    int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', N, T, N, E);
    if (info != 0) LOG(FATAL) << "LAPACK FAILED";
}

void EigenSolve(num_real* E, num_complex* T, num_complex* A, const int& N) {
    int NN = N * N;
    for (int i = 0; i < NN; ++i) T[i] = A[i];
    int info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U', N, T, N, E);
    if (info != 0) LOG(FATAL) << "LAPACK FAILED";
}

void EigenSolve_zgeev(num_complex* E, num_complex* Tl, num_complex* Tr, num_complex* A, const int& N) {
    int info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'V', 'V', N, A, N, E, Tl, N, Tr, N);
    if (info != 0) LOG(FATAL) << "LAPACK FAILED";
}

#elif defined(ARRAY_USE_EIGEN)
#include "../thirdpart/Eigen/Dense"

void LinearSolve(num_real* x, num_real* A, num_real* b, const int& N) {
    MapMXr Mapx(x, N, 1);
    MapMXr MapA(A, N, N);
    MapMXr Mapb(b, N, 1);
    Mapx = MapA.householderQr().solve(Mapb);
}

void EigenSolve(num_real* E, num_real* T, num_real* A, const int& N) {
    // int NN=N*N; ARRAY_1_COPY(T, A, NN);
    plFunction();

    MapMXr MapE(E, N, 1);
    MapMXr MapT(T, N, N);
    MapMXr MapA(A, N, N);
    {
        plScope("pureES");
        Eigen::SelfAdjointEigenSolver<EigMXr> eig(MapA);
        MapE = eig.eigenvalues().real();
        MapT = eig.eigenvectors().real();
    }
}

void EigenSolve(num_real* E, num_complex* T, num_complex* A, const int& N) {
    // int NN=N*N; ARRAY_1_COPY(T, A, NN);
    MapMXr MapE(E, N, 1);
    MapMXc MapT(T, N, N);
    MapMXc MapA(A, N, N);
    Eigen::SelfAdjointEigenSolver<EigMXc> eig(MapA);
    MapE = eig.eigenvalues().real();
    MapT = eig.eigenvectors();
}

void EigenSolve(num_complex* E, num_complex* T, num_complex* A, const int& N) {
    // int NN=N*N; ARRAY_1_COPY(T, A, NN);
    MapMXc MapE(E, N, 1);
    MapMXc MapT(T, N, N);
    MapMXc MapA(A, N, N);
    Eigen::ComplexEigenSolver<EigMXc> eig(MapA);
    MapE = eig.eigenvalues();
    MapT = eig.eigenvectors();
}

void PseudoInverse(num_real* A, num_real* invA, const int& N, num_real e) {
    MapMXr MapA(A, N, N);
    MapMXr MapInvA(invA, N, N);
    MapInvA = MapA.completeOrthogonalDecomposition().pseudoInverse();
    /* Eigen::JacobiSVD<EigMXr> svd = MapA.jacobiSvd(Eigen::ComputeFullU|Eigen::ComputeFullV);
    EigMXr invS = EigMXr::Zero(N, N);
    for (int i = 0; i < N; ++i) {
        if (svd.singularValues()(i) > e || svd.singularValues()(i) < -e) {
            invS(i, i) = 1 / svd.singularValues()(i);
        }
    }
    MapInvA = svd.matrixV()*invS*svd.matrixU().transpose(); */
}
#endif
