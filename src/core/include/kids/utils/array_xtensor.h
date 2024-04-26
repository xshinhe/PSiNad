#ifndef Array_XTENSOR_H
#define Array_XTENSOR_H

#include <array>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xarray.hpp>

namespace ARRAY_XT {

template <class T>
void ARRAY_CLEAR(T* A, size_t N) {
    memset(A, 0, N * sizeof(T));
}

template <class TA, class TB, class TC>
void ARRAY_MATMUL(TA* A, TB* B, TC* C, size_t N1, size_t N2, size_t N3) {
    std::array<size_t, 2> shapeA = {N1, N3};
    std::array<size_t, 2> shapeB = {N1, N2};
    std::array<size_t, 2> shapeC = {N2, N3};
    auto                  a      = xt::adapt(A, N1 * N3, xt::no_ownership(), shapeA);
    auto                  b      = xt::adapt(B, N1 * N2, xt::no_ownership(), shapeB);
    auto                  c      = xt::adapt(C, N2 * N3, xt::no_ownership(), shapeC);
    a                            = xt::linalg::dot(b, c);
}

template <class TA, class TC>
void ARRAY_MATMUL_TRANS1(TA* A, num_real* B, TC* C, size_t N1, size_t N2, size_t N3) {
    std::array<size_t, 2> shapeA = {N1, N3};
    std::array<size_t, 2> shapeB = {N2, N1};
    std::array<size_t, 2> shapeC = {N2, N3};
    auto                  a      = xt::adapt(A, N1 * N3, xt::no_ownership(), shapeA);
    auto                  b      = xt::adapt(B, N1 * N2, xt::no_ownership(), shapeB);
    auto                  c      = xt::adapt(C, N2 * N3, xt::no_ownership(), shapeC);
    a                            = xt::linalg::dot(xt::transpose(b), c);
}

template <class TA, class TC>
void ARRAY_MATMUL_TRANS1(TA* A, num_complex* B, TC* C, size_t N1, size_t N2, size_t N3) {
    std::array<size_t, 2> shapeA = {N1, N3};
    std::array<size_t, 2> shapeB = {N2, N1};
    std::array<size_t, 2> shapeC = {N2, N3};
    auto                  a      = xt::adapt(A, N1 * N3, xt::no_ownership(), shapeA);
    auto                  b      = xt::adapt(B, N1 * N2, xt::no_ownership(), shapeB);
    auto                  c      = xt::adapt(C, N2 * N3, xt::no_ownership(), shapeC);
    a                            = xt::linalg::dot(xt::transpose(xt::conj(b)), c);
}

template <class TA, class TB>
void ARRAY_MATMUL_TRANS2(TA* A, TB* B, num_real* C, size_t N1, size_t N2, size_t N3) {
    std::array<size_t, 2> shapeA = {N1, N3};
    std::array<size_t, 2> shapeB = {N1, N2};
    std::array<size_t, 2> shapeC = {N3, N2};
    auto                  a      = xt::adapt(A, N1 * N3, xt::no_ownership(), shapeA);
    auto                  b      = xt::adapt(B, N1 * N2, xt::no_ownership(), shapeB);
    auto                  c      = xt::adapt(C, N2 * N3, xt::no_ownership(), shapeC);
    a                            = xt::linalg::dot(b, xt::transpose(c));
}

template <class TA, class TB>
void ARRAY_MATMUL_TRANS2(TA* A, TB* B, num_complex* C, size_t N1, size_t N2, size_t N3) {
    std::array<size_t, 2> shapeA = {N1, N3};
    std::array<size_t, 2> shapeB = {N1, N2};
    std::array<size_t, 2> shapeC = {N3, N2};
    auto                  a      = xt::adapt(A, N1 * N3, xt::no_ownership(), shapeA);
    auto                  b      = xt::adapt(B, N1 * N2, xt::no_ownership(), shapeB);
    auto                  c      = xt::adapt(C, N2 * N3, xt::no_ownership(), shapeC);
    a                            = xt::linalg::dot(b, xt::transpose(xt::conj(c)));
}

template <class T>
void ARRAY_OUTER_CONJ2(T* A, T* B, T* C, size_t N1, size_t N2) {
    ARRAY_MATMUL_TRANS2(A, B, C, N1, 1, N2);
}

template <class TA, class TC>
void ARRAY_MATMUL3_TRANS1(TA* A, num_real* B, TC* C, num_real* D, size_t N1, size_t N2, size_t N0, size_t N3) {
    if (N0 == 0) {
        std::array<size_t, 2> shapeA = {N1, N3};
        std::array<size_t, 2> shapeB = {N2, N1};
        std::array<size_t, 2> shapeC = {N2, 1};
        std::array<size_t, 2> shapeD = {N2, N3};
        auto                  a      = xt::adapt(A, N1 * N3, xt::no_ownership(), shapeA);
        auto                  b      = xt::adapt(B, N1 * N2, xt::no_ownership(), shapeB);
        auto                  c      = xt::adapt(C, N2 * 1, xt::no_ownership(), shapeC);
        auto                  d      = xt::adapt(D, N2 * N3, xt::no_ownership(), shapeD);

        a = xt::linalg::dot(xt::transpose(b), xt::linalg::dot(xt::diag(c), d));

    } else {  // N0 == N2
        std::array<size_t, 2> shapeA = {N1, N3};
        std::array<size_t, 2> shapeB = {N2, N1};
        std::array<size_t, 2> shapeC = {N2, N2};
        std::array<size_t, 2> shapeD = {N2, N3};
        auto                  a      = xt::adapt(A, N1 * N3, xt::no_ownership(), shapeA);
        auto                  b      = xt::adapt(B, N1 * N2, xt::no_ownership(), shapeB);
        auto                  c      = xt::adapt(C, N2 * N2, xt::no_ownership(), shapeC);
        auto                  d      = xt::adapt(D, N2 * N3, xt::no_ownership(), shapeD);

        a = xt::linalg::dot(xt::transpose(b), xt::linalg::dot(c, d));
    }
}

template <class TA, class TC>
void ARRAY_MATMUL3_TRANS1(TA* A, num_complex* B, TC* C, num_complex* D, size_t N1, size_t N2, size_t N0, size_t N3) {
    if (N0 == 0) {
        std::array<size_t, 2> shapeA = {N1, N3};
        std::array<size_t, 2> shapeB = {N2, N1};
        std::array<size_t, 2> shapeC = {N2, 1};
        std::array<size_t, 2> shapeD = {N2, N3};
        auto                  a      = xt::adapt(A, N1 * N3, xt::no_ownership(), shapeA);
        auto                  b      = xt::adapt(B, N1 * N2, xt::no_ownership(), shapeB);
        auto                  c      = xt::adapt(C, N2 * 1, xt::no_ownership(), shapeC);
        auto                  d      = xt::adapt(D, N2 * N3, xt::no_ownership(), shapeD);

        a = xt::linalg::dot(xt::transpose(xt::conj(b)), xt::linalg::dot(xt::diag(c), d));

    } else {  // N0 == N2
        std::array<size_t, 2> shapeA = {N1, N3};
        std::array<size_t, 2> shapeB = {N2, N1};
        std::array<size_t, 2> shapeC = {N2, N2};
        std::array<size_t, 2> shapeD = {N2, N3};
        auto                  a      = xt::adapt(A, N1 * N3, xt::no_ownership(), shapeA);
        auto                  b      = xt::adapt(B, N1 * N2, xt::no_ownership(), shapeB);
        auto                  c      = xt::adapt(C, N2 * N2, xt::no_ownership(), shapeC);
        auto                  d      = xt::adapt(D, N2 * N3, xt::no_ownership(), shapeD);

        a = xt::linalg::dot(xt::transpose(xt::conj(b)), xt::linalg::dot(c, d));
    }
}

template <class TA, class TC>
void ARRAY_MATMUL3_TRANS2(TA* A, num_real* B, TC* C, num_real* D, size_t N1, size_t N2, size_t N0, size_t N3) {
    if (N0 == 0) {
        std::array<size_t, 2> shapeA = {N1, N3};
        std::array<size_t, 2> shapeB = {N1, N2};
        std::array<size_t, 2> shapeC = {N2, 1};
        std::array<size_t, 2> shapeD = {N3, N2};
        auto                  a      = xt::adapt(A, N1 * N3, xt::no_ownership(), shapeA);
        auto                  b      = xt::adapt(B, N1 * N2, xt::no_ownership(), shapeB);
        auto                  c      = xt::adapt(C, N2 * 1, xt::no_ownership(), shapeC);
        auto                  d      = xt::adapt(D, N2 * N3, xt::no_ownership(), shapeD);

        a = xt::linalg::dot(b, xt::linalg::dot(xt::diag(c), xt::transpose(d)));
    } else {  // N0 == N2
        std::array<size_t, 2> shapeA = {N1, N3};
        std::array<size_t, 2> shapeB = {N1, N2};
        std::array<size_t, 2> shapeC = {N2, N2};
        std::array<size_t, 2> shapeD = {N3, N2};
        auto                  a      = xt::adapt(A, N1 * N3, xt::no_ownership(), shapeA);
        auto                  b      = xt::adapt(B, N1 * N2, xt::no_ownership(), shapeB);
        auto                  c      = xt::adapt(C, N2 * N2, xt::no_ownership(), shapeC);
        auto                  d      = xt::adapt(D, N2 * N3, xt::no_ownership(), shapeD);

        a = xt::linalg::dot(b, xt::linalg::dot(c, xt::transpose(d)));
    }
}

template <class TA, class TC>
void ARRAY_MATMUL3_TRANS2(TA* A, num_complex* B, TC* C, num_complex* D, size_t N1, size_t N2, size_t N0, size_t N3) {
    if (N0 == 0) {
        std::array<size_t, 2> shapeA = {N1, N3};
        std::array<size_t, 2> shapeB = {N1, N2};
        std::array<size_t, 2> shapeC = {N2, 1};
        std::array<size_t, 2> shapeD = {N3, N2};
        auto                  a      = xt::adapt(A, N1 * N3, xt::no_ownership(), shapeA);
        auto                  b      = xt::adapt(B, N1 * N2, xt::no_ownership(), shapeB);
        auto                  c      = xt::adapt(C, N2 * 1, xt::no_ownership(), shapeC);
        auto                  d      = xt::adapt(D, N2 * N3, xt::no_ownership(), shapeD);

        a = xt::linalg::dot(b, xt::linalg::dot(xt::diag(c), xt::transpose(xt::conj(d))));
    } else {  // N0 == N2
        std::array<size_t, 2> shapeA = {N1, N3};
        std::array<size_t, 2> shapeB = {N1, N2};
        std::array<size_t, 2> shapeC = {N2, N2};
        std::array<size_t, 2> shapeD = {N3, N2};
        auto                  a      = xt::adapt(A, N1 * N3, xt::no_ownership(), shapeA);
        auto                  b      = xt::adapt(B, N1 * N2, xt::no_ownership(), shapeB);
        auto                  c      = xt::adapt(C, N2 * N2, xt::no_ownership(), shapeC);
        auto                  d      = xt::adapt(D, N2 * N3, xt::no_ownership(), shapeD);

        a = xt::linalg::dot(b, xt::linalg::dot(c, xt::transpose(xt::conj(d))));
    }
}

template <class TB, class TC>
TB ARRAY_TRACE2(TB* B, TC* C, size_t N1, size_t N2) {
    std::array<size_t, 2> shapeB = {N1, N2};
    std::array<size_t, 2> shapeC = {N2, N1};
    auto                  b      = xt::adapt(B, N1 * N2, xt::no_ownership(), shapeB);
    auto                  c      = xt::adapt(C, N2 * N1, xt::no_ownership(), shapeC);

    TB res = xt::sum(b * xt::transpose(c))();
    return res;
}

template <class T>
void ARRAY_EYE(T* A, size_t n) {
    std::array<size_t, 2> shapeA = {n, n};
    auto                  a      = xt::adapt(A, n * n, xt::no_ownership(), shapeA);
    a                            = xt::eye<T>(n);
}


};  // namespace ARRAY_XT

#endif  // Array_XTENSOR_H
