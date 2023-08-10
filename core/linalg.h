/**
 * @file linalg.h
 * @author xshinhe
 * @date 2019-01-01
 * @brief linear algebra for array
 * @details
 *      here we use eigen for linear algebra (with EIGEN_USE_MKL)
 *
 * @see
 *  eigen: https://eigen.tuxfamily.org/dox/
 *  mkl  : https://www.intel.com/content/www/us/en/develop/documentation/onemkl-lapack-examples/top.html
 */

#ifndef LINALG_H
#define LINALG_H

#include <complex>

#include "../thirdpart/Eigen/Dense"
#include "../thirdpart/Eigen/QR"
#include "types.h"

using namespace PROJECT_NS;

#define EigMajor Eigen::RowMajor

template <class T>
using EigVX = Eigen::Matrix<T, Eigen::Dynamic, 1, EigMajor>;

template <class T>
using EigMX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, EigMajor>;

template <class T>
using EigAX = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, EigMajor>;

typedef EigVX<num_real> EigVXr;
typedef EigVX<num_complex> EigVXc;
typedef EigMX<num_real> EigMXr;
typedef EigMX<num_complex> EigMXc;
typedef EigMX<num_real> EigAXr;
typedef EigMX<num_complex> EigAXc;
typedef Eigen::Map<EigVXr> MapVXr;
typedef Eigen::Map<EigVXc> MapVXc;
typedef Eigen::Map<EigMXr> MapMXr;
typedef Eigen::Map<EigMXc> MapMXc;
typedef Eigen::Map<EigAXr> MapAXr;
typedef Eigen::Map<EigAXc> MapAXc;

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
        Eigen::Map<EigMX<T>> MapB(B, N2, N1);
        Eigen::Map<EigMX<TC>> MapC(C, N2, 1);
        Eigen::Map<EigMX<T>> MapD(D, N2, N3);
        MapA = MapB.adjoint() * (MapC.asDiagonal() * MapD);
    } else {  // N0 == N2
        Eigen::Map<EigMX<TA>> MapA(A, N1, N3);
        Eigen::Map<EigMX<T>> MapB(B, N2, N1);
        Eigen::Map<EigMX<TC>> MapC(C, N2, N2);
        Eigen::Map<EigMX<T>> MapD(D, N2, N3);
        MapA = MapB.adjoint() * (MapC * MapD);
    }
}

template <class TA, class T, class TC>
void ARRAY_MATMUL3_TRANS2(TA* A, T* B, TC* C, T* D, size_t N1, size_t N2, size_t N0, size_t N3) {
    if (N0 == 0) {
        Eigen::Map<EigMX<TA>> MapA(A, N1, N3);
        Eigen::Map<EigMX<T>> MapB(B, N1, N2);
        Eigen::Map<EigMX<TC>> MapC(C, N2, 1);
        Eigen::Map<EigMX<T>> MapD(D, N3, N2);
        MapA = MapB * (MapC.asDiagonal() * MapD.adjoint());
    } else {  // N0 == N2
        Eigen::Map<EigMX<TA>> MapA(A, N1, N3);
        Eigen::Map<EigMX<T>> MapB(B, N1, N2);
        Eigen::Map<EigMX<TC>> MapC(C, N2, N2);
        Eigen::Map<EigMX<T>> MapD(D, N3, N2);
        MapA = MapB * (MapC * MapD.adjoint());
    }
}

template <class TB, class TC>
TB ARRAY_TRACE2(TB* B, TC* C, size_t N1, size_t N2) {
    Eigen::Map<EigMX<TB>> MapB(B, N1, N2);
    Eigen::Map<EigMX<TC>> MapC(C, N2, N1);
    TB res = (MapB.array() * (MapC.transpose()).array()).sum();
    return res;
}

template <class T>
void ARRAY_EYE(T* A, size_t n) {
    Eigen::Map<EigMX<T>> MapA(A, n, n);
    MapA = EigMX<T>::Identity(n, n);
}

/*===========================================================
=            realize interface of linear algebra            =
===========================================================*/

inline void LinearSolve(num_real* x, num_real* A, num_real* b, const int& N) {
    MapMXr Mapx(x, N, 1);
    MapMXr MapA(A, N, N);
    MapMXr Mapb(b, N, 1);
    Mapx = MapA.householderQr().solve(Mapb);
}

inline void EigenSolve(num_real* E, num_real* T, num_real* A, const int& N) {
    MapMXr MapE(E, N, 1);
    MapMXr MapT(T, N, N);
    MapMXr MapA(A, N, N);
    Eigen::SelfAdjointEigenSolver<EigMXr> eig(MapA);
    MapE = eig.eigenvalues().real();
    MapT = eig.eigenvectors().real();
}

inline void EigenSolve(num_real* E, num_complex* T, num_complex* A, const int& N) {
    MapMXr MapE(E, N, 1);
    MapMXc MapT(T, N, N);
    MapMXc MapA(A, N, N);
    Eigen::SelfAdjointEigenSolver<EigMXc> eig(MapA);
    MapE = eig.eigenvalues().real();
    MapT = eig.eigenvectors();
}

inline void EigenSolve(num_complex* E, num_complex* T, num_complex* A, const int& N) {
    MapMXc MapE(E, N, 1);
    MapMXc MapT(T, N, N);
    MapMXc MapA(A, N, N);
    Eigen::ComplexEigenSolver<EigMXc> eig(MapA);
    MapE = eig.eigenvalues();
    MapT = eig.eigenvectors();
}

inline void PseudoInverse(num_real* A, num_real* invA, const int& N, num_real e = 1E-5) {
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

/*****************************************
    An brief introduction to Eigen library
------------------------------------------

// Eigen        // Matlab           // comments
x.size()        // length(x)        // vector size
C.rows()        // size(C,1)        // number of rows
C.cols()        // size(C,2)        // number of columns
x(i)            // x(i+1)           // Matlab is 1-based
C(i, j)         // C(i+1,j+1)       //

A.resize(4, 4);   // Runtime error if assertions are on.
B.resize(4, 9);   // Runtime error if assertions are on.
A.resize(3, 3);   // Ok; size didn't change.
B.resize(3, 9);   // Ok; only dynamic cols changed.

A << 1, 2, 3,     // Initialize A. The elements can also be
     4, 5, 6,     // matrices, which are stacked along cols
     7, 8, 9;     // and then the rows are stacked.
B << A, A, A;     // B is three horizontally stacked A's.
A.fill(10);       // Fill A with all 10's.


// Eigen                            // Matlab
MatrixXd::Identity(rows,cols)       // eye(rows,cols)
C.setIdentity(rows,cols)            // C = eye(rows,cols)
MatrixXd::Zero(rows,cols)           // zeros(rows,cols)
C.setZero(rows,cols)                // C = ones(rows,cols)
MatrixXd::Ones(rows,cols)           // ones(rows,cols)
C.setOnes(rows,cols)                // C = ones(rows,cols)
MatrixXd::Random(rows,cols)         // rand(rows,cols)*2-1
C.setRandom(rows,cols)              // C = rand(rows,cols)*2-1
VectorXd::LinSpaced(size,low,high)  // linspace(low,high,size)'
v.setLinSpaced(size,low,high)       // v = linspace(low,high,size)'


// Matrix slicing and blocks. All expressions listed here are read/write.
// Templated size versions are faster. Note that Matlab is 1-based (a size N
// vector is x(1)...x(N)).
// Eigen                           // Matlab
x.head(n)                          // x(1:n)
x.head<n>()                        // x(1:n)
x.tail(n)                          // x(end - n + 1: end)
x.tail<n>()                        // x(end - n + 1: end)
x.segment(i, n)                    // x(i+1 : i+n)
x.segment<n>(i)                    // x(i+1 : i+n)
P.block(i, j, rows, cols)          // P(i+1 : i+rows, j+1 : j+cols)
P.block<rows, cols>(i, j)          // P(i+1 : i+rows, j+1 : j+cols)
P.row(i)                           // P(i+1, :)
P.col(j)                           // P(:, j+1)
P.leftCols<cols>()                 // P(:, 1:cols)
P.leftCols(cols)                   // P(:, 1:cols)
P.middleCols<cols>(j)              // P(:, j+1:j+cols)
P.middleCols(j, cols)              // P(:, j+1:j+cols)
P.rightCols<cols>()                // P(:, end-cols+1:end)
P.rightCols(cols)                  // P(:, end-cols+1:end)
P.topRows<rows>()                  // P(1:rows, :)
P.topRows(rows)                    // P(1:rows, :)
P.middleRows<rows>(i)              // P(i+1:i+rows, :)
P.middleRows(i, rows)              // P(i+1:i+rows, :)
P.bottomRows<rows>()               // P(end-rows+1:end, :)
P.bottomRows(rows)                 // P(end-rows+1:end, :)
P.topLeftCorner(rows, cols)        // P(1:rows, 1:cols)
P.topRightCorner(rows, cols)       // P(1:rows, end-cols+1:end)
P.bottomLeftCorner(rows, cols)     // P(end-rows+1:end, 1:cols)
P.bottomRightCorner(rows, cols)    // P(end-rows+1:end, end-cols+1:end)
P.topLeftCorner<rows,cols>()       // P(1:rows, 1:cols)
P.topRightCorner<rows,cols>()      // P(1:rows, end-cols+1:end)
P.bottomLeftCorner<rows,cols>()    // P(end-rows+1:end, 1:cols)
P.bottomRightCorner<rows,cols>()   // P(end-rows+1:end, end-cols+1:end)




// Of particular note is Eigen's swap function which is highly optimized.
// Eigen                           // Matlab
R.row(i) = P.col(j);               // R(i, :) = P(:, i)
R.col(j1).swap(mat1.col(j2));      // R(:, [j1 j2]) = R(:, [j2, j1])

// Views, transpose, etc; all read-write except for .adjoint().
// Eigen                           // Matlab
R.adjoint()                        // R'
R.transpose()                      // R.' or conj(R')
R.diagonal()                       // diag(R)
x.asDiagonal()                     // diag(x)
R.transpose().colwise().reverse(); // rot90(R)
R.conjugate()                      // conj(R)

// All the same as Matlab, but matlab doesn't have *= style operators.
// Matrix-vector.  Matrix-matrix.   Matrix-scalar.
y  = M*x;          R  = P*Q;        R  = P*s;
a  = b*M;          R  = P - Q;      R  = s*P;
a *= M;            R  = P + Q;      R  = P/s;
                   R *= Q;          R  = s*P;
                   R += Q;          R *= s;
                   R -= Q;          R /= s;


// Vectorized operations on each element independently
// Eigen                  // Matlab
R = P.cwiseProduct(Q);    // R = P .* Q
R = P.array() * s.array();// R = P .* s
R = P.cwiseQuotient(Q);   // R = P ./ Q
R = P.array() / Q.array();// R = P ./ Q
R = P.array() + s.array();// R = P + s
R = P.array() - s.array();// R = P - s
R.array() += s;           // R = R + s
R.array() -= s;           // R = R - s
R.array() < Q.array();    // R < Q
R.array() <= Q.array();   // R <= Q
R.cwiseInverse();         // 1 ./ P
R.array().inverse();      // 1 ./ P
R.array().sin()           // sin(P)
R.array().cos()           // cos(P)
R.array().pow(s)          // P .^ s
R.array().square()        // P .^ 2
R.array().cube()          // P .^ 3
R.cwiseSqrt()             // sqrt(P)
R.array().sqrt()          // sqrt(P)
R.array().exp()           // exp(P)
R.array().log()           // log(P)
R.cwiseMax(P)             // max(R, P)
R.array().max(P.array())  // max(R, P)
R.cwiseMin(P)             // min(R, P)
R.array().min(P.array())  // min(R, P)
R.cwiseAbs()              // abs(P)
R.array().abs()           // abs(P)
R.cwiseAbs2()             // abs(P.^2)
R.array().abs2()          // abs(P.^2)
(R.array() < s).select(P,Q);  // (R < s ? P : Q)

// Reductions.
int r, c;
// Eigen                  // Matlab
R.minCoeff()              // min(R(:))
R.maxCoeff()              // max(R(:))
s = R.minCoeff(&r, &c)    // [s, i] = min(R(:)); [r, c] = ind2sub(size(R), i);
s = R.maxCoeff(&r, &c)    // [s, i] = max(R(:)); [r, c] = ind2sub(size(R), i);
R.sum()                   // sum(R(:))
R.colwise().sum()         // sum(R)
R.rowwise().sum()         // sum(R, 2) or sum(R')'
R.prod()                  // prod(R(:))
R.colwise().prod()        // prod(R)
R.rowwise().prod()        // prod(R, 2) or prod(R')'
R.trace()                 // trace(R)
R.all()                   // all(R(:))
R.colwise().all()         // all(R)
R.rowwise().all()         // all(R, 2)
R.any()                   // any(R(:))
R.colwise().any()         // any(R)
R.rowwise().any()         // any(R, 2)




// Dot products, norms, etc.
// Eigen                  // Matlab
x.norm()                  // norm(x).    Note that norm(R) doesn't work in Eigen.
x.squaredNorm()           // dot(x, x)   Note the equivalence is not true for complex
x.dot(y)                  // dot(x, y)
x.cross(y)                // cross(x, y) Requires #include <Eigen/Geometry>

//// Type conversion
// Eigen                           // Matlab
A.cast<num_real>();                  // num_real(A)
A.cast<float>();                   // single(A)
A.cast<int>();                     // int32(A)
A.real();                          // real(A)
A.imag();                          // imag(A)
// if the original type equals destination type, no work is done


// Solve Ax = b. Result stored in x. Matlab: x = A \ b.
x = A.ldlt().solve(b));  // A sym. p.s.d.    #include <Eigen/Cholesky>
x = A.llt() .solve(b));  // A sym. p.d.      #include <Eigen/Cholesky>
x = A.lu()  .solve(b));  // Stable and fast. #include <Eigen/LU>
x = A.qr()  .solve(b));  // No pivoting.     #include <Eigen/QR>
x = A.svd() .solve(b));  // Stable, slowest. #include <Eigen/SVD>
// .ldlt() -> .matrixL() and .matrixD()
// .llt()  -> .matrixL()
// .lu()   -> .matrixL() and .matrixU()
// .qr()   -> .matrixQ() and .matrixR()
// .svd()  -> .matrixU(), .singularValues(), and .matrixV()

// Eigenvalue problems
// Eigen                          // Matlab
A.eigenvalues();                  // eig(A);
EigenSolver<Matrix3d> eig(A);     // [vec val] = eig(A)
eig.eigenvalues();                // diag(val)
eig.eigenvectors();               // vec
// For self-adjoint matrices use SelfAdjointEigenSolver<>


// Use Map of Eigen, associate Map class with C++ pointer

MapMd W(pW,1,N/F);
std::cout << std::setiosflags(std::ios::scientific)
          << std::setprecision(8) << W << std::endl;

*****************************************/


#endif  // LINALG_H
