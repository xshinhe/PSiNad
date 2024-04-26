/**
 * @file array_macro.h
 * @author xshinhe
 * @version 1.0
 * @date 2019-01-01
 * @brief utils for manipulation of array
 * @details
 *
 *      https://en.cppreference.com/w/cpp/header/type_traits
 */

#ifndef Array_Macro_H
#define Array_Macro_H
#include <type_traits>

#define ARRAY_0_OP1(_A, _op1, _C, _n)                        \
    ({                                                       \
        for (int _i = 0; _i < (_n); ++_i) (_A)[_i] _op1(_C); \
    })

#define ARRAY_1_OP1(_A, _op1, _C, _n)                            \
    ({                                                           \
        for (int _i = 0; _i < (_n); ++_i) (_A)[_i] _op1(_C)[_i]; \
    })

#define ARRAY_0_OP2(_A, _B, _op2, _C, _n)                               \
    ({                                                                  \
        for (int _i = 0; _i < (_n); ++_i) (_A)[_i] = (_B)[_i] _op2(_C); \
    })

#define ARRAY_1_OP2(_A, _B, _op2, _C, _n)                                   \
    ({                                                                      \
        for (int _i = 0; _i < (_n); ++_i) (_A)[_i] = (_B)[_i] _op2(_C)[_i]; \
    })

#define ARRAY_CLEAR(_A, _n)         (memset((_A), 0, (_n) * sizeof(*(_A))))
#define ARRAY_0_COPY(_A, _C, _n)    (ARRAY_0_OP1(_A, =, _C, _n))
#define ARRAY_1_COPY(_A, _C, _n)    (ARRAY_1_OP1(_A, =, _C, _n))
#define ARRAY_1_ADD(_A, _B, _C, _n) (ARRAY_1_OP2(_A, _B, +, _C, _n))
#define ARRAY_1_MIN(_A, _B, _C, _n) (ARRAY_1_OP2(_A, _B, -, _C, _n))
#define ARRAY_1_MUL(_A, _B, _C, _n) (ARRAY_1_OP2(_A, _B, *, _C, _n))
#define ARRAY_1_DIV(_A, _B, _C, _n) (ARRAY_1_OP2(_A, _B, /, _C, _n))

// @noused
#define ARRAY_SUM(_A, _n)                                            \
    ({                                                               \
        auto _res = (std::remove_reference<decltype(*_A)>::type)(0); \
        for (int _i = 0; _i < (_n); ++_i) _res += (_A)[_i];          \
        _res;                                                        \
    })

#define ARRAY_0_FUN(_A, _fun, _C, _n)                            \
    ({                                                           \
        for (int _i = 0; _i < (_n); ++_i) (_A)[_i] = _fun((_C)); \
    })

#define ARRAY_1_FUN(_A, _fun, _C, _n)                                \
    ({                                                               \
        for (int _i = 0; _i < (_n); ++_i) (_A)[_i] = _fun((_C)[_i]); \
    })


// specicals for complex array // @noused
#define ARRAY_SIN(_A, _C, _n)  (ARRAY_1_FUN(_A, std::sin, _C, _n))
#define ARRAY_COS(_A, _C, _n)  (ARRAY_1_FUN(_A, std::cos, _C, _n))
#define ARRAY_EXP(_A, _C, _n)  (ARRAY_1_FUN(_A, std::exp, _C, _n))
#define ARRAY_LOG(_A, _C, _n)  (ARRAY_1_FUN(_A, std::log, _C, _n))
#define ARRAY_CONJ(_A, _C, _n) (ARRAY_1_FUN(_A, CONJ_OF, _C, _n))
#define ARRAY_REAL(_A, _C, _n) (ARRAY_1_FUN(_A, REAL_OF, _C, _n))
#define ARRAY_IMAG(_A, _C, _n) (ARRAY_1_FUN(_A, IMAG_OF, _C, _n))
#define ARRAY_NORM(_A, _C, _n) (ARRAY_1_FUN(_A, NORM_OF, _C, _n))
#define ARRAY_PHAS(_A, _C, _n) (ARRAY_1_FUN(_A, PHAS_OF, _C, _n))
#define ARRAY_ABS(_A, _C, _n)  (ARRAY_1_FUN(_A, ABS_OF, _C, _n))

#define ARRAY_EYE(_A, _n)                                                                          \
    ({                                                                                             \
        for (int _i = 0, _idx = 0; _i < (_n); ++_i) {                                              \
            for (int _j = 0; _j < (_n); ++_j, ++_idx) { (_A)[_idx] = ((_i == _j) ? 1.0f : 0.0f); } \
        }                                                                                          \
    })

// @noused
#define ARRAY_ONES(_A, _n1, _n2)                                              \
    ({                                                                        \
        for (int _i = 0, _idx = 0; _i < (_n1); ++_i) {                        \
            for (int _j = 0; _j < (_n2); ++_j, ++_idx) { (_A)[_idx] = 1.0f; } \
        }                                                                     \
    })

// @noused
#define ARRAY_ZEROS(_A, _n1, _n2)                                             \
    ({                                                                        \
        for (int _i = 0, _idx = 0; _i < (_n1); ++_i) {                        \
            for (int _j = 0; _j < (_n2); ++_j, ++_idx) { (_A)[_idx] = 0.0f; } \
        }                                                                     \
    })

// @noused
// array A is set to diagonal of matrix B
#define ARRAY_VDIAGM(_A, _B, _n)            \
    ({                                      \
        int _idx = 0, _add = (_n) + 1;      \
        for (int _i = 0; _i < (_n); ++_i) { \
            (_A)[_i] = (_B)[_idx];          \
            _idx += _add;                   \
        }                                   \
    })

// @noused
// diagonal of array A is set to array B
#define ARRAY_MDIAGV(_A, _B, _n)            \
    ({                                      \
        int _idx = 0, _add = (_n) + 1;      \
        for (int _i = 0; _i < (_n); ++_i) { \
            (_A)[_idx] = (_B)[i];           \
            _idx += _add;                   \
        }                                   \
    })

// return Tr[B*C]:  B[n1,n2], C[n2, n1]
#define ARRAY_TRACE2(_B, _C, _n1, _n2)                               \
    ({                                                               \
        auto _res = (std::remove_reference<decltype(*_B)>::type)(0); \
        int _idxB = 0;                                               \
        for (int _i = 0; _i < (_n1); ++_i) {                         \
            int _idxC = _i;                                          \
            for (int _j = 0; _j < (_n2); ++_j) {                     \
                _res += (_B)[_idxB++] * (_C)[_idxC];                 \
                _idxC += (_n1);                                      \
            }                                                        \
        }                                                            \
        _res;                                                        \
    })

// A = B*C^, B[n1], c[n2]
#define ARRAY_OUTER_CONJ2(_A, _B, _C, _n1, _n2)                                                  \
    ({                                                                                           \
        int _idxA = 0;                                                                           \
        for (int _i = 0; _i < (_n1); ++_i) {                                                     \
            for (int _j = 0; _j < (_n2); ++_j) { (_A)[_idxA++] = (_B)[_i] * CONJ_OF((_C)[_j]); } \
        }                                                                                        \
    })

// @NOTE: it is well optimized.
// @brief: A[n1,n3] = B[n1,n2] * C[n2,n3]
#define ARRAY_MATMUL(_A, _B, _C, _n1, _n2, _n3)                                               \
    ({                                                                                        \
        int _size_A = (_n1) * (_n3);                                                          \
        int _idxB   = 0;                                                                      \
        ARRAY_CLEAR(_A, _size_A);                                                             \
        for (int _i = 0; _i < (_n1); ++_i) {                                                  \
            int _idxC = 0;                                                                    \
            for (int _j = 0; _j < (_n2); ++_j) {                                              \
                auto _tmp = (_B)[_idxB++];                                                    \
                int _idxA = _i * (_n3);                                                       \
                for (int _k = 0; _k < (_n3); ++_k) { (_A)[_idxA++] += _tmp * (_C)[_idxC++]; } \
            }                                                                                 \
        }                                                                                     \
    })

// @brief: A[n1,n3] = B^[n1,n2] * C[n2,n3], so B is B[n2,n1]
#define ARRAY_MATMUL_TRANS1(_A, _B, _C, _n1, _n2, _n3)                                        \
    ({                                                                                        \
        int _size_A = (_n1) * (_n3);                                                          \
        ARRAY_CLEAR(_A, _size_A);                                                             \
        for (int _i = 0; _i < (_n1); ++_i) {                                                  \
            int _idxB = _i;                                                                   \
            int _idxC = 0;                                                                    \
            for (int _j = 0; _j < (_n2); ++_j) {                                              \
                auto _tmp = CONJ_OF((_B)[_idxB]);                                             \
                _idxB += (_n1);                                                               \
                int _idxA = _i * (_n3);                                                       \
                for (int _k = 0; _k < (_n3); ++_k) { (_A)[_idxA++] += _tmp * (_C)[_idxC++]; } \
            }                                                                                 \
        }                                                                                     \
    })

// @brief: A[n1,n3] = B[n1,n2] * C^[n2,n3], so C is C[n3,n2]
#define ARRAY_MATMUL_TRANS2(_A, _B, _C, _n1, _n2, _n3)            \
    ({                                                            \
        int _size_A = (_n1) * (_n3);                              \
        int _idxB   = 0;                                          \
        ARRAY_CLEAR(_A, _size_A);                                 \
        for (int _i = 0; _i < (_n1); ++_i) {                      \
            for (int _j = 0; _j < (_n2); ++_j) {                  \
                auto _tmp = (_B)[_idxB++];                        \
                int _idxC = _j;                                   \
                int _idxA = _i * (_n3);                           \
                for (int _k = 0; _k < (_n3); ++_k) {              \
                    (_A)[_idxA++] += _tmp * CONJ_OF((_C)[_idxC]); \
                    _idxC += (_n2);                               \
                }                                                 \
            }                                                     \
        }                                                         \
    })

// @brief: A = B^*C*D, here C is diagonal, B[n2,n1], C[n2,n2], D[n2,n3]
// if C[n2,n2], please set n0 = n2, if C[n2] please set n0 = 0.
#define ARRAY_MATMUL3_TRANS1(_A, _B, _C, _D, _n1, _n2, _n0, _n3)                              \
    ({                                                                                        \
        int _size_A = (_n1) * (_n3);                                                          \
        ARRAY_CLEAR(_A, _size_A);                                                             \
        int _iaddC = (_n0) + 1;                                                               \
        for (int _i = 0; _i < (_n1); ++_i) {                                                  \
            int _idxB = _i;                                                                   \
            int _idxC = 0, _idxD = 0;                                                         \
            for (int _j = 0; _j < (_n2); ++_j) {                                              \
                auto _tmp = CONJ_OF((_B)[_idxB]) * (_C)[_idxC];                               \
                _idxB += (_n1);                                                               \
                _idxC += _iaddC;                                                              \
                int _idxA = _i * (_n3);                                                       \
                for (int _k = 0; _k < (_n3); ++_k) { (_A)[_idxA++] += _tmp * (_D)[_idxD++]; } \
            }                                                                                 \
        }                                                                                     \
    })

// @brief: A = B*C*D^, here C is diagonal, B[n1,n2], C[n2, n2], D[n3,n2]
// if C[n2,n2], please set n0 = n2, if C[n2] please set n0 = 0.
#define ARRAY_MATMUL3_TRANS2(_A, _B, _C, _D, _n1, _n2, _n0, _n3)  \
    ({                                                            \
        int _size_A = (_n1) * (_n3);                              \
        ARRAY_CLEAR(_A, _size_A);                                 \
        int _idxB  = 0;                                           \
        int _iaddC = (_n0) + 1;                                   \
        for (int _i = 0; _i < (_n1); ++_i) {                      \
            int _idxC = 0;                                        \
            for (int _j = 0; _j < (_n2); ++_j) {                  \
                auto _tmp = (_B)[_idxB++] * (_C)[_idxC];          \
                _idxC += _iaddC;                                  \
                int _idxD = _j;                                   \
                int _idxA = _i * (_n3);                           \
                for (int _k = 0; _k < (_n3); ++_k) {              \
                    (_A)[_idxA++] += _tmp * CONJ_OF((_D)[_idxD]); \
                    _idxD += (_n2);                               \
                }                                                 \
            }                                                     \
        }                                                         \
    })

#endif  // Array_Macro_H
