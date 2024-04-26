#ifndef Array_Eigen_H
#define Array_Eigen_H

#include <Eigen/Dense>

namespace ARRAY_EG {

template <class T>
using EigMX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <class T>
void ARRAY_CLEAR(T* A, size_t N) {
    memset(A, 0, N * sizeof(T));
}

template <class TA, class TB, class TC>
void ARRAY_MATMUL(TA* A, TB* B, TC* C, size_t N1, size_t N2, size_t N3) {
    Eigen::Map<EigMX<TA>> MapA(A, N1, N3);
    Eigen::Map<EigMX<TB>> MapB(B, N1, N2);
    Eigen::Map<EigMX<TC>> MapC(C, N2, N3);
    MapA = MapB * MapC;
}

template <class TA, class TB, class TC>
void ARRAY_MATMUL_TRANS1(TA* A, TB* B, TC* C, size_t N1, size_t N2, size_t N3) {
    Eigen::Map<EigMX<TA>> MapA(A, N1, N3);
    Eigen::Map<EigMX<TB>> MapB(B, N2, N1);
    Eigen::Map<EigMX<TC>> MapC(C, N2, N3);
    MapA = MapB.adjoint() * MapC;
}

template <class TA, class TB, class TC>
void ARRAY_MATMUL_TRANS2(TA* A, TB* B, TC* C, size_t N1, size_t N2, size_t N3) {
    Eigen::Map<EigMX<TA>> MapA(A, N1, N3);
    Eigen::Map<EigMX<TB>> MapB(B, N1, N2);
    Eigen::Map<EigMX<TC>> MapC(C, N3, N2);
    MapA = MapB * MapC.adjoint();
}

template <class T>
void ARRAY_OUTER_CONJ2(T* A, T* B, T* C, size_t N1, size_t N2) {
    ARRAY_MATMUL_TRANS2(A, B, C, N1, 1, N2);
}


template <class TA, class T, class TC>
void ARRAY_MATMUL3_TRANS1(TA* A, T* B, TC* C, T* D, size_t N1, size_t N2, size_t N0, size_t N3) {
    if (N0 == 0) {
        Eigen::Map<EigMX<TA>> MapA(A, N1, N3);
        Eigen::Map<EigMX<T>>  MapB(B, N2, N1);
        Eigen::Map<EigMX<TC>> MapC(C, N2, 1);
        Eigen::Map<EigMX<T>>  MapD(D, N2, N3);
        MapA = MapB.adjoint() * (MapC.asDiagonal() * MapD);
    } else {  // N0 == N2
        Eigen::Map<EigMX<TA>> MapA(A, N1, N3);
        Eigen::Map<EigMX<T>>  MapB(B, N2, N1);
        Eigen::Map<EigMX<TC>> MapC(C, N2, N2);
        Eigen::Map<EigMX<T>>  MapD(D, N2, N3);
        MapA = MapB.adjoint() * (MapC * MapD);
    }
}

template <class TA, class T, class TC>
void ARRAY_MATMUL3_TRANS2(TA* A, T* B, TC* C, T* D, size_t N1, size_t N2, size_t N0, size_t N3) {
    if (N0 == 0) {
        Eigen::Map<EigMX<TA>> MapA(A, N1, N3);
        Eigen::Map<EigMX<T>>  MapB(B, N1, N2);
        Eigen::Map<EigMX<TC>> MapC(C, N2, 1);
        Eigen::Map<EigMX<T>>  MapD(D, N3, N2);
        MapA = MapB * (MapC.asDiagonal() * MapD.adjoint());
    } else {  // N0 == N2
        Eigen::Map<EigMX<TA>> MapA(A, N1, N3);
        Eigen::Map<EigMX<T>>  MapB(B, N1, N2);
        Eigen::Map<EigMX<TC>> MapC(C, N2, N2);
        Eigen::Map<EigMX<T>>  MapD(D, N3, N2);
        MapA = MapB * (MapC * MapD.adjoint());
    }
}

template <class TB, class TC>
TB ARRAY_TRACE2(TB* B, TC* C, size_t N1, size_t N2) {
    Eigen::Map<EigMX<TB>> MapB(B, N1, N2);
    Eigen::Map<EigMX<TC>> MapC(C, N2, N1);
    TB                    res = (MapB.array() * (MapC.transpose()).array()).sum();
    return res;
}

template <class T>
void ARRAY_EYE(T* A, size_t n) {
    Eigen::Map<EigMX<T>> MapA(A, n, n);
    MapA = EigMX<T>::Identity(n, n);
}

};  // namespace ARRAY_EG

#endif  // Array_Eigen_H
