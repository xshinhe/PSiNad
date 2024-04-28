#include "kids/linalg.h"

#include <cmath>
#include <complex>

#include "Eigen/Dense"
#include "Eigen/QR"
#include "kids/Types.h"

#define EigMajor Eigen::RowMajor

namespace PROJECT_NS {

template <class T>
using EigVX = Eigen::Matrix<T, Eigen::Dynamic, 1, EigMajor>;

template <class T>
using EigMX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, EigMajor>;

template <class T>
using EigAX = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, EigMajor>;

typedef EigVX<kids_real> EigVXr;
typedef EigVX<kids_complex> EigVXc;
typedef EigMX<kids_real> EigMXr;
typedef EigMX<kids_complex> EigMXc;
typedef EigMX<kids_real> EigAXr;
typedef EigMX<kids_complex> EigAXc;
typedef Eigen::Map<EigVXr> MapVXr;
typedef Eigen::Map<EigVXc> MapVXc;
typedef Eigen::Map<EigMXr> MapMXr;
typedef Eigen::Map<EigMXc> MapMXc;
typedef Eigen::Map<EigAXr> MapAXr;
typedef Eigen::Map<EigAXc> MapAXc;

bool ARRAY_ISFINITE(kids_real* A, size_t n) {
    for (int i = 0; i < n; ++i)
        if (!std::isfinite((A[i]))) return false;
    return true;
}

bool ARRAY_ISFINITE(kids_complex* A, size_t n) {
    for (int i = 0; i < n; ++i)
        if (!std::isfinite(std::abs(A[i]))) return false;
    return true;
}

void ARRAY_CLEAR(kids_int* A, size_t N) { memset(A, 0, N * sizeof(kids_int)); }

void ARRAY_CLEAR(kids_real* A, size_t N) { memset(A, 0, N * sizeof(kids_real)); }

void ARRAY_CLEAR(kids_complex* A, size_t N) { memset(A, 0, N * sizeof(kids_complex)); }

void ARRAY_MATMUL(kids_real* A, kids_real* B, kids_real* C, size_t N1, size_t N2, size_t N3) {
    MapMXr MapA(A, N1, N3);
    MapMXr MapB(B, N1, N2);
    MapMXr MapC(C, N2, N3);
    MapA = (MapB * MapC).eval();
}

void ARRAY_MATMUL(kids_complex* A, kids_complex* B, kids_complex* C, size_t N1, size_t N2, size_t N3) {
    MapMXc MapA(A, N1, N3);
    MapMXc MapB(B, N1, N2);
    MapMXc MapC(C, N2, N3);
    MapA = (MapB * MapC).eval();
}

void ARRAY_MATMUL(kids_complex* A, kids_real* B, kids_complex* C, size_t N1, size_t N2, size_t N3) {
    MapMXc MapA(A, N1, N3);
    MapMXr MapB(B, N1, N2);
    MapMXc MapC(C, N2, N3);
    MapA = (MapB * MapC).eval();
}

void ARRAY_MATMUL(kids_complex* A, kids_complex* B, kids_real* C, size_t N1, size_t N2, size_t N3) {
    MapMXc MapA(A, N1, N3);
    MapMXc MapB(B, N1, N2);
    MapMXr MapC(C, N2, N3);
    MapA = (MapB * MapC).eval();
}

void ARRAY_MATMUL_TRANS1(kids_real* A, kids_real* B, kids_real* C, size_t N1, size_t N2, size_t N3) {
    MapMXr MapA(A, N1, N3);
    MapMXr MapB(B, N2, N1);
    MapMXr MapC(C, N2, N3);
    MapA = (MapB.adjoint() * MapC).eval();
}

void ARRAY_MATMUL_TRANS1(kids_complex* A, kids_complex* B, kids_complex* C, size_t N1, size_t N2, size_t N3) {
    MapMXc MapA(A, N1, N3);
    MapMXc MapB(B, N2, N1);
    MapMXc MapC(C, N2, N3);
    MapA = (MapB.adjoint() * MapC).eval();
}

void ARRAY_MATMUL_TRANS1(kids_complex* A, kids_real* B, kids_complex* C, size_t N1, size_t N2, size_t N3) {
    MapMXc MapA(A, N1, N3);
    MapMXr MapB(B, N2, N1);
    MapMXc MapC(C, N2, N3);
    MapA = (MapB.adjoint() * MapC).eval();
}

void ARRAY_MATMUL_TRANS1(kids_complex* A, kids_complex* B, kids_real* C, size_t N1, size_t N2, size_t N3) {
    MapMXc MapA(A, N1, N3);
    MapMXc MapB(B, N2, N1);
    MapMXr MapC(C, N2, N3);
    MapA = (MapB.adjoint() * MapC).eval();
}

void ARRAY_MATMUL_TRANS2(kids_real* A, kids_real* B, kids_real* C, size_t N1, size_t N2, size_t N3) {
    MapMXr MapA(A, N1, N3);
    MapMXr MapB(B, N1, N2);
    MapMXr MapC(C, N3, N2);
    MapA = (MapB * MapC.adjoint()).eval();
}

void ARRAY_MATMUL_TRANS2(kids_complex* A, kids_complex* B, kids_complex* C, size_t N1, size_t N2, size_t N3) {
    MapMXc MapA(A, N1, N3);
    MapMXc MapB(B, N1, N2);
    MapMXc MapC(C, N3, N2);
    MapA = (MapB * MapC.adjoint()).eval();
}

void ARRAY_MATMUL_TRANS2(kids_complex* A, kids_real* B, kids_complex* C, size_t N1, size_t N2, size_t N3) {
    MapMXc MapA(A, N1, N3);
    MapMXr MapB(B, N1, N2);
    MapMXc MapC(C, N3, N2);
    MapA = (MapB * MapC.adjoint()).eval();
}

void ARRAY_MATMUL_TRANS2(kids_complex* A, kids_complex* B, kids_real* C, size_t N1, size_t N2, size_t N3) {
    MapMXc MapA(A, N1, N3);
    MapMXc MapB(B, N1, N2);
    MapMXr MapC(C, N3, N2);
    MapA = (MapB * MapC.adjoint()).eval();
}

void ARRAY_OUTER_TRANS2(kids_real* A, kids_real* B, kids_real* C, size_t N1, size_t N2) {
    ARRAY_MATMUL_TRANS2(A, B, C, N1, 1, N2);
}

void ARRAY_OUTER_TRANS2(kids_complex* A, kids_complex* B, kids_complex* C, size_t N1, size_t N2) {
    ARRAY_MATMUL_TRANS2(A, B, C, N1, 1, N2);
}


void ARRAY_MATMUL3_TRANS1(kids_real* A, kids_real* B, kids_real* C, kids_real* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3) {
    if (N0 == 0) {
        MapMXr MapA(A, N1, N3);
        MapMXr MapB(B, N2, N1);
        MapMXr MapC(C, N2, 1);
        MapMXr MapD(D, N2, N3);
        MapA = (MapB.adjoint() * (MapC.asDiagonal() * MapD)).eval();
    } else {  // N0 == N2
        MapMXr MapA(A, N1, N3);
        MapMXr MapB(B, N2, N1);
        MapMXr MapC(C, N2, N2);
        MapMXr MapD(D, N2, N3);
        MapA = (MapB.adjoint() * (MapC * MapD)).eval();
    }
}

void ARRAY_MATMUL3_TRANS1(kids_complex* A, kids_complex* B, kids_complex* C, kids_complex* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3) {
    if (N0 == 0) {
        MapMXc MapA(A, N1, N3);
        MapMXc MapB(B, N2, N1);
        MapMXc MapC(C, N2, 1);
        MapMXc MapD(D, N2, N3);
        MapA = (MapB.adjoint() * (MapC.asDiagonal() * MapD)).eval();
    } else {  // N0 == N2
        MapMXc MapA(A, N1, N3);
        MapMXc MapB(B, N2, N1);
        MapMXc MapC(C, N2, N2);
        MapMXc MapD(D, N2, N3);
        MapA = (MapB.adjoint() * (MapC * MapD)).eval();
    }
}

void ARRAY_MATMUL3_TRANS1(kids_complex* A, kids_real* B, kids_complex* C, kids_real* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3) {
    if (N0 == 0) {
        MapMXc MapA(A, N1, N3);
        MapMXr MapB(B, N2, N1);
        MapMXc MapC(C, N2, 1);
        MapMXr MapD(D, N2, N3);
        MapA = (MapB.adjoint() * (MapC.asDiagonal() * MapD)).eval();
    } else {  // N0 == N2
        MapMXc MapA(A, N1, N3);
        MapMXr MapB(B, N2, N1);
        MapMXc MapC(C, N2, N2);
        MapMXr MapD(D, N2, N3);
        MapA = (MapB.adjoint() * (MapC * MapD)).eval();
    }
}

void ARRAY_MATMUL3_TRANS1(kids_complex* A, kids_complex* B, kids_real* C, kids_complex* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3) {
    if (N0 == 0) {
        MapMXc MapA(A, N1, N3);
        MapMXc MapB(B, N2, N1);
        MapMXr MapC(C, N2, 1);
        MapMXc MapD(D, N2, N3);
        MapA = (MapB.adjoint() * (MapC.asDiagonal() * MapD)).eval();
    } else {  // N0 == N2
        MapMXc MapA(A, N1, N3);
        MapMXc MapB(B, N2, N1);
        MapMXr MapC(C, N2, N2);
        MapMXc MapD(D, N2, N3);
        MapA = (MapB.adjoint() * (MapC * MapD)).eval();
    }
}

void ARRAY_MATMUL3_TRANS2(kids_real* A, kids_real* B, kids_real* C, kids_real* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3) {
    if (N0 == 0) {
        MapMXr MapA(A, N1, N3);
        MapMXr MapB(B, N1, N2);
        MapMXr MapC(C, N2, 1);
        MapMXr MapD(D, N3, N2);
        MapA = (MapB * (MapC.asDiagonal() * MapD.adjoint())).eval();
    } else {  // N0 == N2
        MapMXr MapA(A, N1, N3);
        MapMXr MapB(B, N1, N2);
        MapMXr MapC(C, N2, N2);
        MapMXr MapD(D, N3, N2);
        MapA = (MapB * (MapC * MapD.adjoint())).eval();
    }
}

void ARRAY_MATMUL3_TRANS2(kids_complex* A, kids_complex* B, kids_complex* C, kids_complex* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3) {
    if (N0 == 0) {
        MapMXc MapA(A, N1, N3);
        MapMXc MapB(B, N1, N2);
        MapMXc MapC(C, N2, 1);
        MapMXc MapD(D, N3, N2);
        MapA = (MapB * (MapC.asDiagonal() * MapD.adjoint())).eval();
    } else {  // N0 == N2
        MapMXc MapA(A, N1, N3);
        MapMXc MapB(B, N1, N2);
        MapMXc MapC(C, N2, N2);
        MapMXc MapD(D, N3, N2);
        MapA = (MapB * (MapC * MapD.adjoint())).eval();
    }
}

void ARRAY_MATMUL3_TRANS2(kids_complex* A, kids_real* B, kids_complex* C, kids_real* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3) {
    if (N0 == 0) {
        MapMXc MapA(A, N1, N3);
        MapMXr MapB(B, N1, N2);
        MapMXc MapC(C, N2, 1);
        MapMXr MapD(D, N3, N2);
        MapA = (MapB * (MapC.asDiagonal() * MapD.adjoint())).eval();
    } else {  // N0 == N2
        MapMXc MapA(A, N1, N3);
        MapMXr MapB(B, N1, N2);
        MapMXc MapC(C, N2, N2);
        MapMXr MapD(D, N3, N2);
        MapA = (MapB * (MapC * MapD.adjoint())).eval();
    }
}

void ARRAY_MATMUL3_TRANS2(kids_complex* A, kids_complex* B, kids_real* C, kids_complex* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3) {
    if (N0 == 0) {
        MapMXc MapA(A, N1, N3);
        MapMXc MapB(B, N1, N2);
        MapMXr MapC(C, N2, 1);
        MapMXc MapD(D, N3, N2);
        MapA = (MapB * (MapC.asDiagonal() * MapD.adjoint())).eval();
    } else {  // N0 == N2
        MapMXc MapA(A, N1, N3);
        MapMXc MapB(B, N1, N2);
        MapMXr MapC(C, N2, N2);
        MapMXc MapD(D, N3, N2);
        MapA = (MapB * (MapC * MapD.adjoint())).eval();
    }
}

kids_real ARRAY_TRACE2(kids_real* B, kids_real* C, size_t N1, size_t N2) {
    MapMXr MapB(B, N1, N2);
    MapMXr MapC(C, N2, N1);
    auto res = (MapB.array() * (MapC.transpose()).array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2(kids_complex* B, kids_complex* C, size_t N1, size_t N2) {
    MapMXc MapB(B, N1, N2);
    MapMXc MapC(C, N2, N1);
    auto res = (MapB.array() * (MapC.transpose()).array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2(kids_complex* B, kids_real* C, size_t N1, size_t N2) {
    MapMXc MapB(B, N1, N2);
    MapMXr MapC(C, N2, N1);
    auto res = (MapB.array() * (MapC.transpose()).array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2(kids_real* B, kids_complex* C, size_t N1, size_t N2) {
    MapMXr MapB(B, N1, N2);
    MapMXc MapC(C, N2, N1);
    auto res = (MapB.array() * (MapC.transpose()).array()).sum();
    return res;
}

kids_real ARRAY_TRACE2_DIAG(kids_real* B, kids_real* C, size_t N1, size_t N2) {
    MapMXr MapB(B, N1, N2);
    MapMXr MapC(C, N2, N1);
    auto res = (MapB.diagonal().array() * MapC.diagonal().array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2_DIAG(kids_complex* B, kids_complex* C, size_t N1, size_t N2) {
    MapMXc MapB(B, N1, N2);
    MapMXc MapC(C, N2, N1);
    auto res = (MapB.diagonal().array() * MapC.diagonal().array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2_DIAG(kids_complex* B, kids_real* C, size_t N1, size_t N2) {
    MapMXc MapB(B, N1, N2);
    MapMXr MapC(C, N2, N1);
    auto res = (MapB.diagonal().array() * MapC.diagonal().array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2_DIAG(kids_real* B, kids_complex* C, size_t N1, size_t N2) {
    MapMXr MapB(B, N1, N2);
    MapMXc MapC(C, N2, N1);
    auto res = (MapB.diagonal().array() * MapC.diagonal().array()).sum();
    return res;
}

kids_real ARRAY_TRACE2_OFFD(kids_real* B, kids_real* C, size_t N1, size_t N2) {
    MapMXr MapB(B, N1, N2);
    MapMXr MapC(C, N2, N1);
    auto res = (MapB.array() * (MapC.transpose()).array()).sum()  //
               - (MapB.diagonal().array() * MapC.diagonal().array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2_OFFD(kids_complex* B, kids_complex* C, size_t N1, size_t N2) {
    MapMXc MapB(B, N1, N2);
    MapMXc MapC(C, N2, N1);
    auto res = (MapB.array() * (MapC.transpose()).array()).sum()  //
               - (MapB.diagonal().array() * MapC.diagonal().array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2_OFFD(kids_complex* B, kids_real* C, size_t N1, size_t N2) {
    MapMXc MapB(B, N1, N2);
    MapMXr MapC(C, N2, N1);
    auto res = (MapB.array() * (MapC.transpose()).array()).sum()  //
               - (MapB.diagonal().array() * MapC.diagonal().array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2_OFFD(kids_real* B, kids_complex* C, size_t N1, size_t N2) {
    MapMXr MapB(B, N1, N2);
    MapMXc MapC(C, N2, N1);
    auto res = (MapB.array() * (MapC.transpose()).array()).sum()  //
               - (MapB.diagonal().array() * MapC.diagonal().array()).sum();
    return res;
}

kids_real ARRAY_INNER_TRANS1(kids_real* B, kids_real* C, size_t N1) {
    MapMXr MapB(B, N1, 1);
    MapMXr MapC(C, N1, 1);
    return (MapB.adjoint() * MapC).sum();
}

kids_complex ARRAY_INNER_TRANS1(kids_complex* B, kids_complex* C, size_t N1) {
    MapMXc MapB(B, N1, 1);
    MapMXc MapC(C, N1, 1);
    return (MapB.adjoint() * MapC).sum();
}

kids_complex ARRAY_INNER_TRANS1(kids_complex* B, kids_real* C, size_t N1) {
    MapMXc MapB(B, N1, 1);
    MapMXr MapC(C, N1, 1);
    return (MapB.adjoint() * MapC).sum();
}

kids_complex ARRAY_INNER_TRANS1(kids_real* B, kids_complex* C, size_t N1) {
    MapMXr MapB(B, N1, 1);
    MapMXc MapC(C, N1, 1);
    return (MapB.adjoint() * MapC).sum();
}

kids_real ARRAY_INNER_VMV_TRANS1(kids_real* B, kids_real* C, kids_real* D, size_t N1, size_t N2) {
    MapMXr MapB(B, N1, 1);
    MapMXr MapC(C, N1, N2);
    MapMXr MapD(D, N2, 1);
    return (MapB.adjoint() * MapC * MapD).sum();
}

kids_complex ARRAY_INNER_VMV_TRANS1(kids_complex* B, kids_complex* C, kids_complex* D, size_t N1, size_t N2) {
    MapMXc MapB(B, N1, 1);
    MapMXc MapC(C, N1, N2);
    MapMXc MapD(D, N2, 1);
    return (MapB.adjoint() * MapC * MapD).sum();
}

kids_complex ARRAY_INNER_VMV_TRANS1(kids_complex* B, kids_real* C, kids_complex* D, size_t N1, size_t N2) {
    MapMXc MapB(B, N1, 1);
    MapMXr MapC(C, N1, N2);
    MapMXc MapD(D, N2, 1);
    return (MapB.adjoint() * MapC * MapD).sum();
}

kids_complex ARRAY_INNER_VMV_TRANS1(kids_real* B, kids_complex* C, kids_real* D, size_t N1, size_t N2) {
    MapMXr MapB(B, N1, 1);
    MapMXc MapC(C, N1, N2);
    MapMXr MapD(D, N2, 1);
    return (MapB.adjoint() * MapC * MapD).sum();
}

void ARRAY_EYE(kids_real* A, size_t n) {
    MapMXr MapA(A, n, n);
    MapA = EigMX<kids_real>::Identity(n, n);
}

void ARRAY_EYE(kids_complex* A, size_t n) {
    MapMXc MapA(A, n, n);
    MapA = EigMX<kids_complex>::Identity(n, n);
}

void ARRAY_MAT_DIAG(kids_real* A, kids_real* B, size_t N1) {
    MapMXr MapA(A, N1, N1);
    MapMXr MapB(B, N1, N1);
    ARRAY_CLEAR(A, N1 * N1);
    MapA.diagonal() = MapB.diagonal();
}

void ARRAY_MAT_DIAG(kids_complex* A, kids_complex* B, size_t N1) {
    MapMXc MapA(A, N1, N1);
    MapMXc MapB(B, N1, N1);
    ARRAY_CLEAR(A, N1 * N1);
    MapA.diagonal() = MapB.diagonal();
}

void ARRAY_MAT_DIAG(kids_complex* A, kids_real* B, size_t N1) {
    MapMXc MapA(A, N1, N1);
    MapMXr MapB(B, N1, N1);
    ARRAY_CLEAR(A, N1 * N1);
    MapA.diagonal() = MapB.diagonal();
}

void ARRAY_MAT_OFFD(kids_real* A, kids_real* B, size_t N1) {
    MapMXr MapA(A, N1, N1);
    MapMXr MapB(B, N1, N1);
    MapA                    = MapB;
    MapA.diagonal().array() = kids_real(0);
}

void ARRAY_MAT_OFFD(kids_complex* A, kids_complex* B, size_t N1) {
    MapMXc MapA(A, N1, N1);
    MapMXc MapB(B, N1, N1);
    MapA                    = MapB;
    MapA.diagonal().array() = kids_complex(0, 0);
}

void ARRAY_MAT_OFFD(kids_complex* A, kids_real* B, size_t N1) {
    MapMXc MapA(A, N1, N1);
    MapMXr MapB(B, N1, N1);
    MapA                    = MapB;
    MapA.diagonal().array() = kids_complex(0, 0);
}

/*===========================================================
=            realize interface of linear algebra            =
===========================================================*/

void LinearSolve(kids_real* x, kids_real* A, kids_real* b, size_t N) {
    MapMXr Mapx(x, N, 1);
    MapMXr MapA(A, N, N);
    MapMXr Mapb(b, N, 1);
    Mapx = MapA.householderQr().solve(Mapb);
}

void EigenSolve(kids_real* E, kids_real* T, kids_real* A, size_t N) {
    MapMXr MapE(E, N, 1);
    MapMXr MapT(T, N, N);
    MapMXr MapA(A, N, N);
    Eigen::SelfAdjointEigenSolver<EigMXr> eig(MapA);
    MapE = eig.eigenvalues().real();
    MapT = eig.eigenvectors().real();
}

void EigenSolve(kids_real* E, kids_complex* T, kids_complex* A, size_t N) {
    MapMXr MapE(E, N, 1);
    MapMXc MapT(T, N, N);
    MapMXc MapA(A, N, N);
    Eigen::SelfAdjointEigenSolver<EigMXc> eig(MapA);
    MapE = eig.eigenvalues().real();
    MapT = eig.eigenvectors();
}

void EigenSolve(kids_complex* E, kids_complex* T, kids_complex* A, size_t N) {
    MapMXc MapE(E, N, 1);
    MapMXc MapT(T, N, N);
    MapMXc MapA(A, N, N);
    Eigen::ComplexEigenSolver<EigMXc> eig(MapA);
    MapE = eig.eigenvalues();
    MapT = eig.eigenvectors();
}

void PseudoInverse(kids_real* A, kids_real* invA, size_t N, kids_real e) {
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
    MapMXc Map_A(A, N, N);
    MapMXc Map_expkA(expkA, N, N);
    auto eigr = Eigen::ComplexEigenSolver<EigMXc>(Map_A);
    auto Er   = eigr.eigenvalues();
    auto Vr   = eigr.eigenvectors();
    auto eigl = Eigen::ComplexEigenSolver<EigMXc>(Map_A.adjoint());
    auto El   = eigl.eigenvalues();
    auto Vl   = eigl.eigenvectors();
    auto Slr  = (Vl.adjoint() * Vr).diagonal();
    Map_expkA = Vr * ((k * Er.array()).exp() / Slr.array()).matrix().asDiagonal() * Vl.adjoint();
}

void ARRAY_CORRECT_U(kids_complex* U, size_t N) {
    MapMXc Map_U(U, N, N);
    auto eigr = Eigen::SelfAdjointEigenSolver<EigMXc>(-0.5e0 * kids_complex(0, 1) * (Map_U - Map_U.adjoint()));
    auto Er   = eigr.eigenvalues().real();
    auto Vr   = eigr.eigenvectors();
    Map_U     = Vr * ((kids_complex(0, 1) * Er.array()).exp()).matrix().asDiagonal() * Vr.adjoint();
}

void ARRAY_TRANSPOSE(kids_real* A, size_t N1, size_t N2) {
    MapMXr Map_A(A, N1, N2);
    MapMXr Map_Anew(A, N2, N1);
    Map_Anew = Map_A.adjoint().eval();
}

void ARRAY_TRANSPOSE(kids_complex* A, size_t N1, size_t N2) {
    MapMXc Map_A(A, N1, N2);
    MapMXc Map_Anew(A, N2, N1);
    Map_Anew = Map_A.adjoint().eval();
}

};  // namespace PROJECT_NS
