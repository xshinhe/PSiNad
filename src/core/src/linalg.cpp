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

static void ARRAY_MATMUL_UNIVERSAL(kids_real* A, kids_real* B, kids_real* C, size_t N1, size_t N2, size_t N3,
                                   bool ad1 = false, bool ad2 = false) {
    size_t NB1 = (N2 == 0) ? ((!ad1) ? N1 : N3) : ((!ad1) ? N1 : N2);
    size_t NB2 = (N2 == 0) ? ((!ad1) ? N3 : N1) : ((!ad1) ? N2 : N1);
    size_t NC1 = (N2 == 0) ? N3 : ((!ad2) ? N2 : N3);
    size_t NC2 = (N2 == 0) ? 1 : ((!ad2) ? N3 : N2);
    MapMXr MapA(A, N1, N3);
    MapMXr MapB(B, NB1, NB2);
    MapMXr MapC(C, NC1, NC2);
    auto   M1 = (!ad1) ? MapB.eval() : MapB.adjoint();
    auto   M2 = (N2 == 0) ? MapC.asDiagonal() : ((!ad2) ? MapC.eval() : MapC.adjoint());
    MapA      = M1 * M2;
}

static void ARRAY_MATMUL_UNIVERSAL(kids_complex* A, kids_complex* B, kids_complex* C,  //
                                   size_t N1, size_t N2, size_t N3, bool ad1 = false, bool ad2 = false) {
    size_t NB1 = (N2 == 0) ? ((!ad1) ? N1 : N3) : ((!ad1) ? N1 : N2);
    size_t NB2 = (N2 == 0) ? ((!ad1) ? N3 : N1) : ((!ad1) ? N2 : N1);
    size_t NC1 = (N2 == 0) ? N3 : ((!ad2) ? N2 : N3);
    size_t NC2 = (N2 == 0) ? 1 : ((!ad2) ? N3 : N2);
    MapMXc MapA(A, N1, N3);
    MapMXc MapB(B, NB1, NB2);
    MapMXc MapC(C, NC1, NC2);
    auto   M1 = (!ad1) ? MapB.eval() : MapB.adjoint();
    auto   M2 = (N2 == 0) ? MapC.asDiagonal() : ((!ad2) ? MapC.eval() : MapC.adjoint());
    MapA      = M1 * M2;
}

static void ARRAY_MATMUL_UNIVERSAL(kids_complex* A, kids_real* B, kids_complex* C,  //
                                   size_t N1, size_t N2, size_t N3, bool ad1 = false, bool ad2 = false) {
    size_t NB1 = (N2 == 0) ? ((!ad1) ? N1 : N3) : ((!ad1) ? N1 : N2);
    size_t NB2 = (N2 == 0) ? ((!ad1) ? N3 : N1) : ((!ad1) ? N2 : N1);
    size_t NC1 = (N2 == 0) ? N3 : ((!ad2) ? N2 : N3);
    size_t NC2 = (N2 == 0) ? 1 : ((!ad2) ? N3 : N2);
    MapMXc MapA(A, N1, N3);
    MapMXr MapB(B, NB1, NB2);
    MapMXc MapC(C, NC1, NC2);
    auto   M1 = (!ad1) ? MapB.eval() : MapB.adjoint();
    auto   M2 = (N2 == 0) ? MapC.asDiagonal() : ((!ad2) ? MapC.eval() : MapC.adjoint());
    MapA      = M1 * M2;
}

static void ARRAY_MATMUL_UNIVERSAL(kids_complex* A, kids_complex* B, kids_real* C,  //
                                   size_t N1, size_t N2, size_t N3, bool ad1 = false, bool ad2 = false) {
    size_t NB1 = (N2 == 0) ? ((!ad1) ? N1 : N3) : ((!ad1) ? N1 : N2);
    size_t NB2 = (N2 == 0) ? ((!ad1) ? N3 : N1) : ((!ad1) ? N2 : N1);
    size_t NC1 = (N2 == 0) ? N3 : ((!ad2) ? N2 : N3);
    size_t NC2 = (N2 == 0) ? 1 : ((!ad2) ? N3 : N2);
    MapMXc MapA(A, N1, N3);
    MapMXc MapB(B, NB1, NB2);
    MapMXr MapC(C, NC1, NC2);
    auto   M1 = (!ad1) ? MapB.eval() : MapB.adjoint();
    auto   M2 = (N2 == 0) ? MapC.asDiagonal() : ((!ad2) ? MapC.eval() : MapC.adjoint());
    MapA      = M1 * M2;
}

void ARRAY_MATMUL(kids_real* A, kids_real* B, kids_real* C, size_t N1, size_t N2, size_t N3) {
    ARRAY_MATMUL_UNIVERSAL(A, B, C, N1, N2, N3, false, false);
}

void ARRAY_MATMUL(kids_complex* A, kids_complex* B, kids_complex* C, size_t N1, size_t N2, size_t N3) {
    ARRAY_MATMUL_UNIVERSAL(A, B, C, N1, N2, N3, false, false);
}

void ARRAY_MATMUL(kids_complex* A, kids_real* B, kids_complex* C, size_t N1, size_t N2, size_t N3) {
    ARRAY_MATMUL_UNIVERSAL(A, B, C, N1, N2, N3, false, false);
}

void ARRAY_MATMUL(kids_complex* A, kids_complex* B, kids_real* C, size_t N1, size_t N2, size_t N3) {
    ARRAY_MATMUL_UNIVERSAL(A, B, C, N1, N2, N3, false, false);
}

void ARRAY_MATMUL_TRANS1(kids_real* A, kids_real* B, kids_real* C, size_t N1, size_t N2, size_t N3) {
    ARRAY_MATMUL_UNIVERSAL(A, B, C, N1, N2, N3, true, false);
}

void ARRAY_MATMUL_TRANS1(kids_complex* A, kids_complex* B, kids_complex* C, size_t N1, size_t N2, size_t N3) {
    ARRAY_MATMUL_UNIVERSAL(A, B, C, N1, N2, N3, true, false);
}

void ARRAY_MATMUL_TRANS1(kids_complex* A, kids_real* B, kids_complex* C, size_t N1, size_t N2, size_t N3) {
    ARRAY_MATMUL_UNIVERSAL(A, B, C, N1, N2, N3, true, false);
}

void ARRAY_MATMUL_TRANS1(kids_complex* A, kids_complex* B, kids_real* C, size_t N1, size_t N2, size_t N3) {
    ARRAY_MATMUL_UNIVERSAL(A, B, C, N1, N2, N3, true, false);
}

void ARRAY_MATMUL_TRANS2(kids_real* A, kids_real* B, kids_real* C, size_t N1, size_t N2, size_t N3) {
    ARRAY_MATMUL_UNIVERSAL(A, B, C, N1, N2, N3, false, true);
}

void ARRAY_MATMUL_TRANS2(kids_complex* A, kids_complex* B, kids_complex* C, size_t N1, size_t N2, size_t N3) {
    ARRAY_MATMUL_UNIVERSAL(A, B, C, N1, N2, N3, false, true);
}

void ARRAY_MATMUL_TRANS2(kids_complex* A, kids_real* B, kids_complex* C, size_t N1, size_t N2, size_t N3) {
    ARRAY_MATMUL_UNIVERSAL(A, B, C, N1, N2, N3, false, true);
}

void ARRAY_MATMUL_TRANS2(kids_complex* A, kids_complex* B, kids_real* C, size_t N1, size_t N2, size_t N3) {
    ARRAY_MATMUL_UNIVERSAL(A, B, C, N1, N2, N3, false, true);
}

void ARRAY_OUTER_TRANS2(kids_real* A, kids_real* B, kids_real* C, size_t N1, size_t N2) {
    ARRAY_MATMUL_TRANS2(A, B, C, N1, 1, N2);
}

void ARRAY_OUTER_TRANS2(kids_complex* A, kids_complex* B, kids_complex* C, size_t N1, size_t N2) {
    ARRAY_MATMUL_TRANS2(A, B, C, N1, 1, N2);
}

void ARRAY_MATMUL3_TRANS1(kids_real* A, kids_real* B, kids_real* C, kids_real* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3) {
    // MapMXr MapA(A, N1, N3);
    // MapMXr MapB(B, N2, N1);
    // MapMXr MapD(D, N2, N3);
    if (N0 == 0) {  // assuming N1 == N2 == N3
        ARRAY_MATMUL_TRANS1(A, B, C, N1, N0, N2);
        ARRAY_MATMUL(A, A, D, N1, N2, N3);
        // MapMXr MapC(C, N2, 1);
        // MapA = (MapB.adjoint() * (MapC.asDiagonal() * MapD)).eval();
    } else {  // assuming N1 == N2 == N0 == N3
        ARRAY_MATMUL_TRANS1(A, B, C, N1, N2, N2);
        ARRAY_MATMUL(A, A, D, N1, N2, N3);
        // MapMXr MapC(C, N2, N2);
        // MapA = (MapB.adjoint() * (MapC * MapD)).eval();
    }
}

void ARRAY_MATMUL3_TRANS1(kids_complex* A, kids_complex* B, kids_complex* C, kids_complex* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3) {
    // MapMXc MapA(A, N1, N3);
    // MapMXc MapB(B, N2, N1);
    // MapMXc MapD(D, N2, N3);
    if (N0 == 0) {  // assuming N1 == N2 == N3
        ARRAY_MATMUL_TRANS1(A, B, C, N1, N0, N2);
        ARRAY_MATMUL(A, A, D, N1, N2, N3);
        // MapMXc MapC(C, N2, 1);
        // MapA = (MapB.adjoint() * (MapC.asDiagonal() * MapD)).eval();
    } else {  // assuming N1 == N2 == N0 == N3
        ARRAY_MATMUL_TRANS1(A, B, C, N1, N2, N2);
        ARRAY_MATMUL(A, A, D, N1, N2, N3);
        // MapMXc MapC(C, N2, N2);
        // MapA = (MapB.adjoint() * (MapC * MapD)).eval();
    }
}

void ARRAY_MATMUL3_TRANS1(kids_complex* A, kids_real* B, kids_complex* C, kids_real* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3) {
    // MapMXc MapA(A, N1, N3);
    // MapMXr MapB(B, N2, N1);
    // MapMXr MapD(D, N2, N3);
    if (N0 == 0) {  // assuming N1 == N2 == N3
        ARRAY_MATMUL_TRANS1(A, B, C, N1, N0, N2);
        ARRAY_MATMUL(A, A, D, N1, N2, N3);
        // MapMXc MapC(C, N2, 1);
        // MapA = (MapB.adjoint() * (MapC.asDiagonal() * MapD)).eval();
    } else {  // assuming N1 == N2 == N0 == N3
        ARRAY_MATMUL_TRANS1(A, B, C, N1, N2, N2);
        ARRAY_MATMUL(A, A, D, N1, N2, N3);
        // MapMXc MapC(C, N2, N2);
        // MapA = (MapB.adjoint() * (MapC * MapD)).eval();
    }
}

void ARRAY_MATMUL3_TRANS1(kids_complex* A, kids_complex* B, kids_real* C, kids_complex* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3) {
    // MapMXc MapA(A, N1, N3);
    // MapMXc MapB(B, N2, N1);
    // MapMXc MapD(D, N2, N3);
    if (N0 == 0) {  // assuming N1 == N2 == N3
        ARRAY_MATMUL_TRANS1(A, B, C, N1, N0, N2);
        ARRAY_MATMUL(A, A, D, N1, N2, N3);
        // MapMXr MapC(C, N2, 1);
        // MapA = (MapB.adjoint() * (MapC.asDiagonal() * MapD)).eval();
    } else {  // assuming N1 == N2 == N0 == N3
        ARRAY_MATMUL_TRANS1(A, B, C, N1, N2, N2);
        ARRAY_MATMUL(A, A, D, N1, N2, N3);
        // MapMXr MapC(C, N2, N2);
        // MapA = (MapB.adjoint() * (MapC * MapD)).eval();
    }
}

void ARRAY_MATMUL3_TRANS2(kids_real* A, kids_real* B, kids_real* C, kids_real* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3) {
    // MapMXr MapA(A, N1, N3);
    // MapMXr MapB(B, N1, N2);
    // MapMXr MapD(D, N3, N2);
    if (N0 == 0) {  // assuming N1 == N2 == N3
        ARRAY_MATMUL(A, B, C, N1, N0, N2);
        ARRAY_MATMUL_TRANS2(A, A, D, N1, N2, N3);
        // MapMXr MapC(C, N2, 1);
        // MapA = (MapB * (MapC.asDiagonal() * MapD.adjoint())).eval();
    } else {  // assuming N1 == N2 == N0 == N3
        ARRAY_MATMUL(A, B, C, N1, N2, N2);
        ARRAY_MATMUL_TRANS2(A, A, D, N1, N2, N3);
        // MapMXr MapC(C, N2, N2);
        // MapA = (MapB * (MapC * MapD.adjoint())).eval();
    }
}

void ARRAY_MATMUL3_TRANS2(kids_complex* A, kids_complex* B, kids_complex* C, kids_complex* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3) {
    // MapMXc MapA(A, N1, N3);
    // MapMXc MapB(B, N1, N2);
    // MapMXc MapD(D, N3, N2);
    if (N0 == 0) {  // assuming N1 == N2 == N3
        ARRAY_MATMUL(A, B, C, N1, N0, N2);
        ARRAY_MATMUL_TRANS2(A, A, D, N1, N2, N3);
        // MapMXc MapC(C, N2, 1);
        // MapA = (MapB * (MapC.asDiagonal() * MapD.adjoint())).eval();
    } else {  // assuming N1 == N2 == N0 == N3
        ARRAY_MATMUL(A, B, C, N1, N2, N2);
        ARRAY_MATMUL_TRANS2(A, A, D, N1, N2, N3);
        // MapMXc MapC(C, N2, N2);
        // MapA = (MapB * (MapC * MapD.adjoint())).eval();
    }
}

void ARRAY_MATMUL3_TRANS2(kids_complex* A, kids_real* B, kids_complex* C, kids_real* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3) {
    // MapMXc MapA(A, N1, N3);
    // MapMXr MapB(B, N1, N2);
    // MapMXr MapD(D, N3, N2);
    if (N0 == 0) {  // assuming N1 == N2 == N3
        ARRAY_MATMUL(A, B, C, N1, N0, N2);
        ARRAY_MATMUL_TRANS2(A, A, D, N1, N2, N3);
        // MapMXc MapC(C, N2, 1);
        // MapA = (MapB * (MapC.asDiagonal() * MapD.adjoint())).eval();
    } else {  // assuming N1 == N2 == N0 == N3
        ARRAY_MATMUL(A, B, C, N1, N2, N2);
        ARRAY_MATMUL_TRANS2(A, A, D, N1, N2, N3);
        // MapMXc MapC(C, N2, N2);
        // MapA = (MapB * (MapC * MapD.adjoint())).eval();
    }
}

void ARRAY_MATMUL3_TRANS2(kids_complex* A, kids_complex* B, kids_real* C, kids_complex* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3) {
    // MapMXc MapA(A, N1, N3);
    // MapMXc MapB(B, N1, N2);
    // MapMXc MapD(D, N3, N2);
    if (N0 == 0) {  // assuming N1 == N2 == N3
        ARRAY_MATMUL(A, B, C, N1, N0, N2);
        ARRAY_MATMUL_TRANS2(A, A, D, N1, N2, N3);
        // MapMXr MapC(C, N2, 1);
        // MapA = (MapB * (MapC.asDiagonal() * MapD.adjoint())).eval();
    } else {  // assuming N1 == N2 == N0 == N3
        ARRAY_MATMUL(A, B, C, N1, N2, N2);
        ARRAY_MATMUL_TRANS2(A, A, D, N1, N2, N3);
        // MapMXr MapC(C, N2, N2);
        // MapA = (MapB * (MapC * MapD.adjoint())).eval();
    }
}

kids_real ARRAY_TRACE2(kids_real* B, kids_real* C, size_t N1, size_t N2) {
    MapMXr MapB(B, N1, N2);
    MapMXr MapC(C, N2, N1);
    auto   res = (MapB.array() * (MapC.transpose()).array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2(kids_complex* B, kids_complex* C, size_t N1, size_t N2) {
    MapMXc MapB(B, N1, N2);
    MapMXc MapC(C, N2, N1);
    auto   res = (MapB.array() * (MapC.transpose()).array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2(kids_complex* B, kids_real* C, size_t N1, size_t N2) {
    MapMXc MapB(B, N1, N2);
    MapMXr MapC(C, N2, N1);
    auto   res = (MapB.array() * (MapC.transpose()).array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2(kids_real* B, kids_complex* C, size_t N1, size_t N2) {
    MapMXr MapB(B, N1, N2);
    MapMXc MapC(C, N2, N1);
    auto   res = (MapB.array() * (MapC.transpose()).array()).sum();
    return res;
}

kids_real ARRAY_TRACE2_DIAG(kids_real* B, kids_real* C, size_t N1, size_t N2) {
    MapMXr MapB(B, N1, N2);
    MapMXr MapC(C, N2, N1);
    auto   res = (MapB.diagonal().array() * MapC.diagonal().array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2_DIAG(kids_complex* B, kids_complex* C, size_t N1, size_t N2) {
    MapMXc MapB(B, N1, N2);
    MapMXc MapC(C, N2, N1);
    auto   res = (MapB.diagonal().array() * MapC.diagonal().array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2_DIAG(kids_complex* B, kids_real* C, size_t N1, size_t N2) {
    MapMXc MapB(B, N1, N2);
    MapMXr MapC(C, N2, N1);
    auto   res = (MapB.diagonal().array() * MapC.diagonal().array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2_DIAG(kids_real* B, kids_complex* C, size_t N1, size_t N2) {
    MapMXr MapB(B, N1, N2);
    MapMXc MapC(C, N2, N1);
    auto   res = (MapB.diagonal().array() * MapC.diagonal().array()).sum();
    return res;
}

kids_real ARRAY_TRACE2_OFFD(kids_real* B, kids_real* C, size_t N1, size_t N2) {
    MapMXr MapB(B, N1, N2);
    MapMXr MapC(C, N2, N1);
    auto   res = (MapB.array() * (MapC.transpose()).array()).sum()  //
               - (MapB.diagonal().array() * MapC.diagonal().array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2_OFFD(kids_complex* B, kids_complex* C, size_t N1, size_t N2) {
    MapMXc MapB(B, N1, N2);
    MapMXc MapC(C, N2, N1);
    auto   res = (MapB.array() * (MapC.transpose()).array()).sum()  //
               - (MapB.diagonal().array() * MapC.diagonal().array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2_OFFD(kids_complex* B, kids_real* C, size_t N1, size_t N2) {
    MapMXc MapB(B, N1, N2);
    MapMXr MapC(C, N2, N1);
    auto   res = (MapB.array() * (MapC.transpose()).array()).sum()  //
               - (MapB.diagonal().array() * MapC.diagonal().array()).sum();
    return res;
}

kids_complex ARRAY_TRACE2_OFFD(kids_real* B, kids_complex* C, size_t N1, size_t N2) {
    MapMXr MapB(B, N1, N2);
    MapMXc MapC(C, N2, N1);
    auto   res = (MapB.array() * (MapC.transpose()).array()).sum()  //
               - (MapB.diagonal().array() * MapC.diagonal().array()).sum();
    return res;
}

kids_real ARRAY_INNER_TRANS1(kids_real* B, kids_real* C, size_t N1) {
    kids_real res;
    ARRAY_MATMUL_TRANS1(&res, B, C, 1, N1, 1);
    return res;
}

kids_complex ARRAY_INNER_TRANS1(kids_complex* B, kids_complex* C, size_t N1) {
    kids_complex res;
    ARRAY_MATMUL_TRANS1(&res, B, C, 1, N1, 1);
    return res;
}

kids_complex ARRAY_INNER_TRANS1(kids_complex* B, kids_real* C, size_t N1) {
    kids_complex res;
    ARRAY_MATMUL_TRANS1(&res, B, C, 1, N1, 1);
    return res;
}

kids_complex ARRAY_INNER_TRANS1(kids_real* B, kids_complex* C, size_t N1) {
    kids_complex res;
    ARRAY_MATMUL_TRANS1(&res, B, C, 1, N1, 1);
    return res;
}

kids_real ARRAY_INNER_VMV_TRANS1(kids_real* B, kids_real* C, kids_real* D, size_t N1, size_t N2) {
    kids_real res;
    ARRAY_MATMUL3_TRANS1(&res, B, C, D, 1, N1, N2, 1);
    return res;
}

kids_complex ARRAY_INNER_VMV_TRANS1(kids_complex* B, kids_complex* C, kids_complex* D, size_t N1, size_t N2) {
    kids_complex res;
    ARRAY_MATMUL3_TRANS1(&res, B, C, D, 1, N1, N2, 1);
    return res;
}

kids_complex ARRAY_INNER_VMV_TRANS1(kids_complex* B, kids_real* C, kids_complex* D, size_t N1, size_t N2) {
    kids_complex res;
    ARRAY_MATMUL3_TRANS1(&res, B, C, D, 1, N1, N2, 1);
    return res;
}

kids_complex ARRAY_INNER_VMV_TRANS1(kids_real* B, kids_complex* C, kids_real* D, size_t N1, size_t N2) {
    kids_complex res;
    ARRAY_MATMUL3_TRANS1(&res, B, C, D, 1, N1, N2, 1);
    return res;
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
