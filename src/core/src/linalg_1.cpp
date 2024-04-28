#include <cmath>
#include <complex>

#include "Eigen/Dense"
#include "Eigen/QR"
#include "kids/Types.h"
#include "kids/linalg.h"

#define EigMajor Eigen::RowMajor

namespace PROJECT_NS {

template <class T>
using EigVX = Eigen::Matrix<T, Eigen::Dynamic, 1, EigMajor>;

template <class T>
using EigMX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, EigMajor>;

template <class T>
using EigAX = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, EigMajor>;

using EigVXr = EigVX<kids_real>;
using EigVXc = EigVX<kids_complex>;
using EigMXr = EigMX<kids_real>;
using EigMXc = EigMX<kids_complex>;
using EigAXr = EigMX<kids_real>;
using EigAXc = EigMX<kids_complex>;
using MapVXr = Eigen::Map<EigVXr>;
using MapVXc = Eigen::Map<EigVXc>;
using MapMXr = Eigen::Map<EigMXr>;
using MapMXc = Eigen::Map<EigMXc>;
using MapAXr = Eigen::Map<EigAXr>;
using MapAXc = Eigen::Map<EigAXc>;

void LinearSolve(kids_real* x, kids_real* A, kids_real* b, size_t N) {
    MapMXr Mapx(x, N, 1);
    MapMXr MapA(A, N, N);
    MapMXr Mapb(b, N, 1);
    Mapx = MapA.householderQr().solve(Mapb);
}

void EigenSolve(kids_real* E, kids_real* T, kids_real* A, size_t N) {
    MapMXr                                MapE(E, N, 1);
    MapMXr                                MapT(T, N, N);
    MapMXr                                MapA(A, N, N);
    Eigen::SelfAdjointEigenSolver<EigMXr> eig(MapA);
    MapE = eig.eigenvalues().real();
    MapT = eig.eigenvectors().real();
}

void EigenSolve(kids_real* E, kids_complex* T, kids_complex* A, size_t N) {
    MapMXr                                MapE(E, N, 1);
    MapMXc                                MapT(T, N, N);
    MapMXc                                MapA(A, N, N);
    Eigen::SelfAdjointEigenSolver<EigMXc> eig(MapA);
    MapE = eig.eigenvalues().real();
    MapT = eig.eigenvectors();
}

void EigenSolve(kids_complex* E, kids_complex* T, kids_complex* A, size_t N) {
    MapMXc                            MapE(E, N, 1);
    MapMXc                            MapT(T, N, N);
    MapMXc                            MapA(A, N, N);
    Eigen::ComplexEigenSolver<EigMXc> eig(MapA);
    MapE = eig.eigenvalues();
    MapT = eig.eigenvectors();
}

void PseudoInverse(kids_real* A, kids_real* invA, size_t N, kids_real e) {
    MapMXr MapA(A, N, N);
    MapMXr MapInvA(invA, N, N);
    MapInvA = MapA.completeOrthogonalDecomposition().pseudoInverse();
    // Eigen::JacobiSVD<EigMXr> svd = MapA.jacobiSvd(Eigen::ComputeFullU|Eigen::ComputeFullV);
    // EigMXr invS = EigMXr::Zero(N, N);
    // for (int i = 0; i < N; ++i) {
    //     if (svd.singularValues()(i) > e || svd.singularValues()(i) < -e) {
    //         invS(i, i) = 1 / svd.singularValues()(i);
    //     }
    // }
    // MapInvA = svd.matrixV()*invS*svd.matrixU().transpose();
}

void ARRAY_INV_MAT(kids_real* invA, kids_real* A, size_t N) {
    MapMXr Map_invA(invA, N, N);
    MapMXr Map_A(A, N, N);
    Map_invA = Map_A.lu().inverse();
}

void ARRAY_INV_MAT(kids_complex* invA, kids_complex* A, size_t N) {
    MapMXc Map_invA(invA, N, N);
    MapMXc Map_A(A, N, N);
    Map_invA = Map_A.lu().inverse();
}

void ARRAY_EXP_MAT_GENERAL(kids_complex* expkA, kids_complex* A, kids_complex k, size_t N) {
    // MapMXc Map_A(A, N, N);
    // MapMXc Map_expkA(expkA, N, N);
    // auto eigr = Eigen::ComplexEigenSolver<EigMXc>(Map_A);
    // auto Er   = eigr.eigenvalues();
    // auto Vr   = eigr.eigenvectors();
    // auto eigl = Eigen::ComplexEigenSolver<EigMXc>(Map_A.adjoint());
    // auto El   = eigl.eigenvalues();
    // auto Vl   = eigl.eigenvectors();
    // auto Slr  = (Vl.adjoint() * Vr).diagonal();
    // Map_expkA = Vr * ((k * Er.array()).exp() / Slr.array()).matrix().asDiagonal() * Vl.adjoint();

    kids_complex* Vr_ptr   = new kids_complex[N * N];
    kids_complex* Vl_ptr   = new kids_complex[N * N];
    kids_complex* At_ptr   = new kids_complex[N * N];
    kids_complex* Slr_ptr  = new kids_complex[N * N];
    kids_complex* lamb_ptr = new kids_complex[N];

    for (int i = 0; i < N * N; ++i) At_ptr[i] = A[i];
    ARRAY_TRANSPOSE(At_ptr, N, N);
    EigenSolve(lamb_ptr, Vl_ptr, At_ptr, N);
    EigenSolve(lamb_ptr, Vr_ptr, A, N);
    ARRAY_MATMUL_TRANS1(Slr_ptr, Vl_ptr, Vr_ptr, N, N, N);
    for (int i = 0, ii = 0; i < N; ++i, ii += (N + 1)) {  //
        lamb_ptr[i] = std::exp(k * lamb_ptr[i]) / Slr_ptr[ii];
    }
    ARRAY_MATMUL3_TRANS2(expkA, Vr_ptr, lamb_ptr, Vl_ptr, N, N, 0, N);

    delete[] Vr_ptr;
    delete[] Vl_ptr;
    delete[] At_ptr;
    delete[] Slr_ptr;
    delete[] lamb_ptr;
}

void ARRAY_CORRECT_U(kids_complex* U, size_t N) {
    // MapMXc Map_U(U, N, N);
    // auto   eigr = Eigen::SelfAdjointEigenSolver<EigMXc>(-0.5e0 * kids_complex(0, 1) * (Map_U - Map_U.adjoint()));
    // auto   Er   = eigr.eigenvalues().real();
    // auto   Vr   = eigr.eigenvectors();
    // Map_U       = Vr * ((kids_complex(0, 1) * Er.array()).exp()).matrix().asDiagonal() * Vr.adjoint();

    kids_complex* V_ptr     = new kids_complex[N * N];
    kids_complex* A_ptr     = new kids_complex[N * N];
    kids_real*    lamb_ptr  = new kids_real[N];
    kids_complex* lamb2_ptr = new kids_complex[N];

    for (int i = 0, ik = 0; i < N; ++i) {
        for (int k = 0, ki = i; k < N; ++k, ++ik, ki += N) {
            A_ptr[ik] = -0.5e0 * kids_complex(0, 1) * (U[ik] - std::conj(U[ki]));
        }
    }
    EigenSolve(lamb_ptr, V_ptr, A_ptr, N);
    for (int i = 0; i < N; ++i) lamb2_ptr[i] = std::exp(kids_complex(0, 1) * lamb_ptr[i]);
    ARRAY_MATMUL3_TRANS2(U, V_ptr, lamb2_ptr, V_ptr, N, N, 0, N);
    delete[] V_ptr;
    delete[] A_ptr;
    delete[] lamb_ptr;
}

};  // namespace PROJECT_NS
