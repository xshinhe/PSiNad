/**@file        linalg.h
 * @brief       Provide linalg APIs
 * @details
 *  unified APIs, and realization can be based on either Eigen or MKL
 * @see
 *  eigen: https://eigen.tuxfamily.org/dox/
 *  mkl  : https://www.intel.com/content/www/us/en/develop/documentation/onemkl-lapack-examples/top.html
 *
 * @author      Xin He
 * @date        2024-03
 * @version     1.0
 * @copyright   GNU Lesser General Public License (LGPL)
 *
 *              Copyright (c) 2024 Xin He, Liu-Group
 *
 *  This software is a product of Xin's PhD research conducted by Professor Liu's
 *  Group at the College of Chemistry and Molecular Engineering, Peking University.
 *  All rights are reserved by Peking University.
 *  You should have received a copy of the GNU Lesser General Public License along
 *  with this software. If not, see <https://www.gnu.org/licenses/lgpl-3.0.en.html>
 **********************************************************************************
 * @par revision:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-04-25  <td> split from tpl APIs, which is slow in compilation
 * </table>
 **********************************************************************************
 */

#ifndef KIDS_LINALG_H
#define KIDS_LINALG_H

#include "kids/Types.h"

#define KIDS_LINALG_BIND_EIGEN_NOT_USE_TEMPLATE

#ifndef KIDS_LINALG_BIND_EIGEN_NOT_USE_TEMPLATE
#include "kids/linalg_tpl.h"
#else  // KIDS_LINALG_BIND_EIGEN_NOT_USE_TEMPLATE

namespace PROJECT_NS {

/**
 * Check if all elements of a real array are finite.
 *
 * @param A Pointer to the real array.
 * @param n Size of the array.
 * @return True if all elements are finite, false otherwise.
 */
bool ARRAY_ISFINITE(kids_real* A, size_t n);

/**
 * Check if all elements of a complex array are finite.
 *
 * @param A Pointer to the complex array.
 * @param n Size of the array.
 * @return True if all elements are finite, false otherwise.
 */
bool ARRAY_ISFINITE(kids_complex* A, size_t n);

/**
 * Set all elements of an integer array to zero.
 *
 * @param A Pointer to the integer array.
 * @param N Size of the array.
 */
void ARRAY_CLEAR(kids_int* A, size_t N);

/**
 * Set all elements of a real array to zero.
 *
 * @param A Pointer to the real array.
 * @param N Size of the array.
 */
void ARRAY_CLEAR(kids_real* A, size_t N);

/**
 * Set all elements of a complex array to zero.
 *
 * @param A Pointer to the complex array.
 * @param N Size of the array.
 */
void ARRAY_CLEAR(kids_complex* A, size_t N);

/**
 * Perform matrix multiplication where A = B * C for real matrices.
 *
 * @param A Pointer to the result real matrix.
 * @param B Pointer to the first real matrix.
 * @param C Pointer to the second real matrix.
 * @param N1 Number of rows in B and A.
 * @param N2 Number of columns in B and rows in C.
 * @param N3 Number of columns in C and A.
 */
void ARRAY_MATMUL(kids_real* A, kids_real* B, kids_real* C, size_t N1, size_t N2, size_t N3);

/**
 * Perform matrix multiplication where A = B * C for complex matrices.
 *
 * @param A Pointer to the result complex matrix.
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second complex matrix.
 * @param N1 Number of rows in B and A.
 * @param N2 Number of columns in B and rows in C.
 * @param N3 Number of columns in C and A.
 */
void ARRAY_MATMUL(kids_complex* A, kids_complex* B, kids_complex* C, size_t N1, size_t N2, size_t N3);

/**
 * Perform matrix multiplication where A = B * C for complex and real matrices.
 *
 * @param A Pointer to the result complex matrix.
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second real matrix.
 * @param N1 Number of rows in B and A.
 * @param N2 Number of columns in B and rows in C.
 * @param N3 Number of columns in C and A.
 */
void ARRAY_MATMUL(kids_complex* A, kids_real* B, kids_complex* C, size_t N1, size_t N2, size_t N3);

/**
 * Perform matrix multiplication where A = B * C for complex and real matrices.
 *
 * @param A Pointer to the result complex matrix.
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second real matrix.
 * @param N1 Number of rows in B and A.
 * @param N2 Number of columns in B and rows in C.
 * @param N3 Number of columns in C and A.
 */
void ARRAY_MATMUL(kids_complex* A, kids_complex* B, kids_real* C, size_t N1, size_t N2, size_t N3);

/**
 * Perform matrix multiplication where A = B^T * C for real matrices.
 *
 * @param A Pointer to the result real matrix.
 * @param B Pointer to the first real matrix.
 * @param C Pointer to the second real matrix.
 * @param N1 Number of rows in B^T and A.
 * @param N2 Number of columns in B^T and rows in C.
 * @param N3 Number of columns in C and A.
 */
void ARRAY_MATMUL_TRANS1(kids_real* A, kids_real* B, kids_real* C, size_t N1, size_t N2, size_t N3);

/**
 * Perform matrix multiplication where A = B^T * C for complex matrices.
 *
 * @param A Pointer to the result complex matrix.
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second complex matrix.
 * @param N1 Number of rows in B^T and A.
 * @param N2 Number of columns in B^T and rows in C.
 * @param N3 Number of columns in C and A.
 */
void ARRAY_MATMUL_TRANS1(kids_complex* A, kids_complex* B, kids_complex* C, size_t N1, size_t N2, size_t N3);

/**
 * Perform matrix multiplication where A = B^T * C for complex and real matrices.
 *
 * @param A Pointer to the result complex matrix.
 * @param B Pointer to the first real matrix.
 * @param C Pointer to the second complex matrix.
 * @param N1 Number of rows in B^T and A.
 * @param N2 Number of columns in B^T and rows in C.
 * @param N3 Number of columns in C and A.
 */
void ARRAY_MATMUL_TRANS1(kids_complex* A, kids_real* B, kids_complex* C, size_t N1, size_t N2, size_t N3);

/**
 * Perform matrix multiplication where A = B^T * C for complex and real matrices.
 *
 * @param A Pointer to the result complex matrix.
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second real matrix.
 * @param N1 Number of rows in B^T and A.
 * @param N2 Number of columns in B^T and rows in C.
 * @param N3 Number of columns in C and A.
 */
void ARRAY_MATMUL_TRANS1(kids_complex* A, kids_complex* B, kids_real* C, size_t N1, size_t N2, size_t N3);

/**
 * Perform matrix multiplication where A = B * C^T for real matrices.
 *
 * @param A Pointer to the result real matrix.
 * @param B Pointer to the first real matrix.
 * @param C Pointer to the second real matrix.
 * @param N1 Number of rows in B and A.
 * @param N2 Number of columns in B and A^T.
 * @param N3 Number of columns in C^T and A.
 */
void ARRAY_MATMUL_TRANS2(kids_real* A, kids_real* B, kids_real* C, size_t N1, size_t N2, size_t N3);

/**
 * Perform matrix multiplication where A = B * C^T for complex matrices.
 *
 * @param A Pointer to the result complex matrix.
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second complex matrix.
 * @param N1 Number of rows in B and A.
 * @param N2 Number of columns in B and A^T.
 * @param N3 Number of columns in C^T and A.
 */
void ARRAY_MATMUL_TRANS2(kids_complex* A, kids_complex* B, kids_complex* C, size_t N1, size_t N2, size_t N3);

/**
 * Perform matrix multiplication where A = B * C^T for complex and real matrices.
 *
 * @param A Pointer to the result complex matrix.
 * @param B Pointer to the first real matrix.
 * @param C Pointer to the second complex matrix.
 * @param N1 Number of rows in B and A.
 * @param N2 Number of columns in B and A^T.
 * @param N3 Number of columns in C^T and A.
 */
void ARRAY_MATMUL_TRANS2(kids_complex* A, kids_real* B, kids_complex* C, size_t N1, size_t N2, size_t N3);

/**
 * Perform matrix multiplication where A = B * C^T for complex and real matrices.
 *
 * @param A Pointer to the result complex matrix.
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second real matrix.
 * @param N1 Number of rows in B and A.
 * @param N2 Number of columns in B and A^T.
 * @param N3 Number of columns in C^T and A.
 */
void ARRAY_MATMUL_TRANS2(kids_complex* A, kids_complex* B, kids_real* C, size_t N1, size_t N2, size_t N3);

/**
 * Perform outer product with transpose where A = B^T * C for real matrices.
 *
 * @param A Pointer to the result real matrix.
 * @param B Pointer to the first real vector.
 * @param C Pointer to the second real vector.
 * @param N1 Number of rows in A.
 * @param N2 Number of columns in A.
 */
void ARRAY_OUTER_TRANS2(kids_real* A, kids_real* B, kids_real* C, size_t N1, size_t N2);

/**
 * Perform outer product with transpose where A = B^T * C for complex matrices.
 *
 * @param A Pointer to the result complex matrix.
 * @param B Pointer to the first complex vector.
 * @param C Pointer to the second complex vector.
 * @param N1 Number of rows in A.
 * @param N2 Number of columns in A.
 */
void ARRAY_OUTER_TRANS2(kids_complex* A, kids_complex* B, kids_complex* C, size_t N1, size_t N2);

/**
 * Perform matrix multiplication where A = B^T * C * D for real matrices.
 *
 * @param A Pointer to the result real matrix.
 * @param B Pointer to the first real matrix.
 * @param C Pointer to the second real matrix.
 * @param D Pointer to the third real matrix.
 * @param N1 Number of rows in A.
 * @param N2 Number of columns in B^T and rows in C.
 * @param N0 Number of columns in B^T and rows in D; N0==0 when C is a diagonal vector and N2 is for rows of D.
 * @param N3 Number of columns in D and A.
 */
void ARRAY_MATMUL3_TRANS1(kids_real* A, kids_real* B, kids_real* C, kids_real* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3);

/**
 * Perform matrix multiplication where A = B^T * C * D for complex matrices.
 *
 * @param A Pointer to the result complex matrix.
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second complex matrix.
 * @param D Pointer to the third complex matrix.
 * @param N1 Number of rows in A.
 * @param N2 Number of columns in B^T and rows in C.
 * @param N0 Number of columns in B^T and rows in D; N0==0 when C is a diagonal vector and N2 is for rows of D.
 * @param N3 Number of columns in D and A.
 */
void ARRAY_MATMUL3_TRANS1(kids_complex* A, kids_complex* B, kids_complex* C, kids_complex* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3);

/**
 * Perform matrix multiplication where A = B^T * C * D for complex and real matrices.
 *
 * @param A Pointer to the result complex matrix.
 * @param B Pointer to the first real matrix.
 * @param C Pointer to the second complex matrix.
 * @param D Pointer to the third real matrix.
 * @param N1 Number of rows in A.
 * @param N2 Number of columns in B^T and rows in C.
 * @param N0 Number of columns in B^T and rows in D; N0==0 when C is a diagonal vector and N2 is for rows of D.
 * @param N3 Number of columns in D and A.
 */
void ARRAY_MATMUL3_TRANS1(kids_complex* A, kids_real* B, kids_complex* C, kids_real* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3);

/**
 * Perform matrix multiplication where A = B^T * C * D for complex and real matrices.
 *
 * @param A Pointer to the result complex matrix.
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second real matrix.
 * @param D Pointer to the third complex matrix.
 * @param N1 Number of rows in A.
 * @param N2 Number of columns in B^T and rows in C.
 * @param N0 Number of columns in B^T and rows in D; N0==0 when C is a diagonal vector and N2 is for rows of D.
 * @param N3 Number of columns in D and A.
 */
void ARRAY_MATMUL3_TRANS1(kids_complex* A, kids_complex* B, kids_real* C, kids_complex* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3);

/**
 * Perform matrix multiplication where A = B * C * D^T for real matrices.
 *
 * @param A Pointer to the result real matrix.
 * @param B Pointer to the first real matrix.
 * @param C Pointer to the second real matrix.
 * @param D Pointer to the third real matrix.
 * @param N1 Number of rows in A.
 * @param N2 Number of columns in B and rows in C.
 * @param N0 Number of columns in C and rows in D^T; N0==0 when C is a diagonal vector and N2 is for rows of D^T.
 * @param N3 Number of columns in D^T and A.
 */
void ARRAY_MATMUL3_TRANS2(kids_real* A, kids_real* B, kids_real* C, kids_real* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3);

/**
 * Perform matrix multiplication where A = B * C * D^T for complex matrices.
 *
 * @param A Pointer to the result complex matrix.
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second complex matrix.
 * @param D Pointer to the third complex matrix.
 * @param N1 Number of rows in A.
 * @param N2 Number of columns in B and rows in C.
 * @param N0 Number of columns in C and rows in D^T; N0==0 when C is a diagonal vector and N2 is for rows of D^T.
 * @param N3 Number of columns in D^T and A.
 */
void ARRAY_MATMUL3_TRANS2(kids_complex* A, kids_complex* B, kids_complex* C, kids_complex* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3);

/**
 * Perform matrix multiplication where A = B * C * D^T for complex and real matrices.
 *
 * @param A Pointer to the result complex matrix.
 * @param B Pointer to the first real matrix.
 * @param C Pointer to the second complex matrix.
 * @param D Pointer to the third real matrix.
 * @param N1 Number of rows in A.
 * @param N2 Number of columns in B and rows in C.
 * @param N0 Number of columns in C and rows in D^T; N0==0 when C is a diagonal vector and N2 is for rows of D^T.
 * @param N3 Number of columns in D^T and A.
 */
void ARRAY_MATMUL3_TRANS2(kids_complex* A, kids_real* B, kids_complex* C, kids_real* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3);

/**
 * Perform matrix multiplication where A = B * C * D^T for complex and real matrices.
 *
 * @param A Pointer to the result complex matrix.
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second real matrix.
 * @param D Pointer to the third complex matrix.
 * @param N1 Number of rows in A.
 * @param N2 Number of columns in B and rows in C.
 * @param N0 Number of columns in C and rows in D^T; N0==0 when C is a diagonal vector and N2 is for rows of D^T.
 * @param N3 Number of columns in D^T and A.
 */
void ARRAY_MATMUL3_TRANS2(kids_complex* A, kids_complex* B, kids_real* C, kids_complex* D,  //
                          size_t N1, size_t N2, size_t N0, size_t N3);

/**
 * Compute the trace of the matrix product B * C for real matrices.
 *
 * @param B Pointer to the first real matrix.
 * @param C Pointer to the second real matrix.
 * @param N1 Number of rows in B and columns in C.
 * @param N2 Number of columns in B and rows in C.
 * @return The trace of the matrix product B * C.
 */
kids_real ARRAY_TRACE2(kids_real* B, kids_real* C, size_t N1, size_t N2);

/**
 * Compute the trace of the matrix product B * C for complex matrices.
 *
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second complex matrix.
 * @param N1 Number of rows in B and columns in C.
 * @param N2 Number of columns in B and rows in C.
 * @return The trace of the matrix product B * C.
 */
kids_complex ARRAY_TRACE2(kids_complex* B, kids_complex* C, size_t N1, size_t N2);

/**
 * Compute the trace of the matrix product B * C for complex and real matrices.
 *
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second real matrix.
 * @param N1 Number of rows in B and columns in C.
 * @param N2 Number of columns in B and rows in C.
 * @return The trace of the matrix product B * C.
 */
kids_complex ARRAY_TRACE2(kids_complex* B, kids_real* C, size_t N1, size_t N2);

/**
 * Compute the trace of the matrix product B * C for real and complex matrices.
 *
 * @param B Pointer to the first real matrix.
 * @param C Pointer to the second complex matrix.
 * @param N1 Number of rows in B and columns in C.
 * @param N2 Number of columns in B and rows in C.
 * @return The trace of the matrix product B * C.
 */
kids_complex ARRAY_TRACE2(kids_real* B, kids_complex* C, size_t N1, size_t N2);

/**
 * Compute the trace of the diagonal elements of the matrix product B * C for real matrices.
 *
 * @param B Pointer to the first real matrix.
 * @param C Pointer to the second real matrix.
 * @param N1 Number of rows in B and columns in C.
 * @param N2 Number of columns in B and rows in C.
 * @return The trace of the diagonal elements of the matrix product B * C.
 */
kids_real ARRAY_TRACE2_DIAG(kids_real* B, kids_real* C, size_t N1, size_t N2);

/**
 * Compute the trace of the diagonal elements of the matrix product B * C for complex matrices.
 *
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second complex matrix.
 * @param N1 Number of rows in B and columns in C.
 * @param N2 Number of columns in B and rows in C.
 * @return The trace of the diagonal elements of the matrix product B * C.
 */
kids_complex ARRAY_TRACE2_DIAG(kids_complex* B, kids_complex* C, size_t N1, size_t N2);

/**
 * Compute the trace of the diagonal elements of the matrix product B * C for complex and real matrices.
 *
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second real matrix.
 * @param N1 Number of rows in B and columns in C.
 * @param N2 Number of columns in B and rows in C.
 * @return The trace of the diagonal elements of the matrix product B * C.
 */
kids_complex ARRAY_TRACE2_DIAG(kids_complex* B, kids_real* C, size_t N1, size_t N2);

/**
 * Compute the trace of the diagonal elements of the matrix product B * C for real and complex matrices.
 *
 * @param B Pointer to the first real matrix.
 * @param C Pointer to the second complex matrix.
 * @param N1 Number of rows in B and columns in C.
 * @param N2 Number of columns in B and rows in C.
 * @return The trace of the diagonal elements of the matrix product B * C.
 */
kids_complex ARRAY_TRACE2_DIAG(kids_real* B, kids_complex* C, size_t N1, size_t N2);

/**
 * Compute the trace of the off-diagonal elements of the matrix product B * C for real matrices.
 *
 * @param B Pointer to the first real matrix.
 * @param C Pointer to the second real matrix.
 * @param N1 Number of rows in B and columns in C.
 * @param N2 Number of columns in B and rows in C.
 * @return The trace of the off-diagonal elements of the matrix product B * C.
 */
kids_real ARRAY_TRACE2_OFFD(kids_real* B, kids_real* C, size_t N1, size_t N2);

/**
 * Compute the trace of the off-diagonal elements of the matrix product B * C for complex matrices.
 *
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second complex matrix.
 * @param N1 Number of rows in B and columns in C.
 * @param N2 Number of columns in B and rows in C.
 * @return The trace of the off-diagonal elements of the matrix product B * C.
 */
kids_complex ARRAY_TRACE2_OFFD(kids_complex* B, kids_complex* C, size_t N1, size_t N2);

/**
 * Compute the trace of the off-diagonal elements of the matrix product B * C for complex and real matrices.
 *
 * @param B Pointer to the first complex matrix.
 * @param C Pointer to the second real matrix.
 * @param N1 Number of rows in B and columns in C.
 * @param N2 Number of columns in B and rows in C.
 * @return The trace of the off-diagonal elements of the matrix product B * C.
 */
kids_complex ARRAY_TRACE2_OFFD(kids_complex* B, kids_real* C, size_t N1, size_t N2);

/**
 * Compute the trace of the off-diagonal elements of the matrix product B * C for real and complex matrices.
 *
 * @param B Pointer to the first real matrix.
 * @param C Pointer to the second complex matrix.
 * @param N1 Number of rows in B and columns in C.
 * @param N2 Number of columns in B and rows in C.
 * @return The trace of the off-diagonal elements of the matrix product B * C.
 */
kids_complex ARRAY_TRACE2_OFFD(kids_real* B, kids_complex* C, size_t N1, size_t N2);

/**
 * Compute the inner product of the transpose of B with C for real matrices.
 *
 * @param B Pointer to the real matrix B.
 * @param C Pointer to the real matrix C.
 * @param N1 Number of elements in B and C (length of vectors).
 * @return The inner product of transpose(B) and C.
 */
kids_real ARRAY_INNER_TRANS1(kids_real* B, kids_real* C, size_t N1);

/**
 * Compute the inner product of the conjugate transpose of B with C for complex matrices.
 *
 * @param B Pointer to the complex matrix B.
 * @param C Pointer to the complex matrix C.
 * @param N1 Number of elements in B and C (length of vectors).
 * @return The inner product of conjugate transpose(B) and C.
 */
kids_complex ARRAY_INNER_TRANS1(kids_complex* B, kids_complex* C, size_t N1);

/**
 * Compute the inner product of the conjugate transpose of B with C for complex and real matrices.
 *
 * @param B Pointer to the complex matrix B.
 * @param C Pointer to the real matrix C.
 * @param N1 Number of elements in B and C (length of vectors).
 * @return The inner product of conjugate transpose(B) and C.
 */
kids_complex ARRAY_INNER_TRANS1(kids_complex* B, kids_real* C, size_t N1);

/**
 * Compute the inner product of the transpose of B with C for real and complex matrices.
 *
 * @param B Pointer to the real matrix B.
 * @param C Pointer to the complex matrix C.
 * @param N1 Number of elements in B and C (length of vectors).
 * @return The inner product of transpose(B) and C.
 */
kids_complex ARRAY_INNER_TRANS1(kids_real* B, kids_complex* C, size_t N1);

/**
 * Compute the inner product of the transpose of the real vector B with the matrix C,
 * element-wise multiplication with the real vector D.
 *
 * @param B Pointer to the real vector B.
 * @param C Pointer to the real matrix C.
 * @param D Pointer to the real vector D.
 * @param N1 Number of elements in B and rows in C and D.
 * @param N2 Number of columns in C and elements in D.
 * @return The inner product of transpose(B) * C * D.
 */
kids_real ARRAY_INNER_VMV_TRANS1(kids_real* B, kids_real* C, kids_real* D, size_t N1, size_t N2);

/**
 * Compute the inner product of the conjugate transpose of the complex vector B with
 * the complex matrix C, element-wise multiplication with the complex vector D.
 *
 * @param B Pointer to the complex vector B.
 * @param C Pointer to the complex matrix C.
 * @param D Pointer to the complex vector D.
 * @param N1 Number of elements in B and rows in C and D.
 * @param N2 Number of columns in C and elements in D.
 * @return The inner product of conjugate transpose(B) * C * D.
 */
kids_complex ARRAY_INNER_VMV_TRANS1(kids_complex* B, kids_complex* C, kids_complex* D, size_t N1, size_t N2);

/**
 * Compute the inner product of the conjugate transpose of the complex vector B with the real matrix C,
 * element-wise multiplication with the complex vector D.
 *
 * @param B Pointer to the complex vector B.
 * @param C Pointer to the real matrix C.
 * @param D Pointer to the complex vector D.
 * @param N1 Number of elements in B and rows in C and D.
 * @param N2 Number of columns in C and elements in D.
 * @return The inner product of conjugate transpose(B) * C * D.
 */
kids_complex ARRAY_INNER_VMV_TRANS1(kids_complex* B, kids_real* C, kids_complex* D, size_t N1, size_t N2);

/**
 * Compute the inner product of the transpose of the real vector B with the complex matrix C,
 * element-wise multiplication with the real vector D.
 *
 * @param B Pointer to the real vector B.
 * @param C Pointer to the complex matrix C.
 * @param D Pointer to the real vector D.
 * @param N1 Number of elements in B and rows in C and D.
 * @param N2 Number of columns in C and elements in D.
 * @return The inner product of transpose(B) * C * D.
 */
kids_complex ARRAY_INNER_VMV_TRANS1(kids_real* B, kids_complex* C, kids_real* D, size_t N1, size_t N2);

/**
 * Generate an identity matrix of size n for real numbers.
 *
 * @param A Pointer to the real matrix A to store the result.
 * @param n Size of the identity matrix.
 */
void ARRAY_EYE(kids_real* A, size_t n);

/**
 * Generate an identity matrix of size n for complex numbers.
 *
 * @param A Pointer to the complex matrix A to store the result.
 * @param n Size of the identity matrix.
 */
void ARRAY_EYE(kids_complex* A, size_t n);

/**
 * Copy the diagonal elements from matrix B to matrix A for real matrices.
 *
 * @param A Pointer to the real matrix A.
 * @param B Pointer to the real matrix B.
 * @param N1 Number of rows and columns in A and B.
 */
void ARRAY_MAT_DIAG(kids_real* A, kids_real* B, size_t N1);

/**
 * Copy the diagonal elements from matrix B to matrix A for complex matrices.
 *
 * @param A Pointer to the complex matrix A.
 * @param B Pointer to the complex matrix B.
 * @param N1 Number of rows and columns in A and B.
 */
void ARRAY_MAT_DIAG(kids_complex* A, kids_complex* B, size_t N1);

/**
 * Copy the diagonal elements from matrix B to matrix A for complex and real matrices.
 *
 * @param A Pointer to the complex matrix A.
 * @param B Pointer to the real matrix B.
 * @param N1 Number of rows and columns in A and B.
 */
void ARRAY_MAT_DIAG(kids_complex* A, kids_real* B, size_t N1);

/**
 * Copy the off-diagonal elements from matrix B to matrix A for real matrices.
 *
 * @param A Pointer to the real matrix A.
 * @param B Pointer to the real matrix B.
 * @param N1 Number of rows and columns in A and B.
 */
void ARRAY_MAT_OFFD(kids_real* A, kids_real* B, size_t N1);

/**
 * Copy the off-diagonal elements from matrix B to matrix A for complex matrices.
 *
 * @param A Pointer to the complex matrix A.
 * @param B Pointer to the complex matrix B.
 * @param N1 Number of rows and columns in A and B.
 */
void ARRAY_MAT_OFFD(kids_complex* A, kids_complex* B, size_t N1);

/**
 * Copy the off-diagonal elements from matrix B to matrix A for complex and real matrices.
 *
 * @param A Pointer to the complex matrix A.
 * @param B Pointer to the real matrix B.
 * @param N1 Number of rows and columns in A and B.
 */
void ARRAY_MAT_OFFD(kids_complex* A, kids_real* B, size_t N1);

/**
 * Solve a linear system Ax = b for real matrices.
 *
 * @param x Pointer to the solution vector x.
 * @param A Pointer to the coefficient matrix A.
 * @param b Pointer to the right-hand side vector b.
 * @param N Size of the system (number of rows or columns in A).
 */
void LinearSolve(kids_real* x, kids_real* A, kids_real* b, size_t N);

/**
 * Solve the eigenvalue problem for real matrices.
 *
 * @param E Pointer to the array to store the eigenvalues.
 * @param T Pointer to the temporary matrix.
 * @param A Pointer to the matrix for eigenvalue computation.
 * @param N Size of the matrix (number of rows or columns in A).
 */
void EigenSolve(kids_real* E, kids_real* T, kids_real* A, size_t N);

/**
 * Solve the eigenvalue problem for (hermite) complex matrices.
 *
 * @param E Pointer to the array to store the eigenvalues.
 * @param T Pointer to the temporary complex matrix.
 * @param A Pointer to the complex matrix for eigenvalue computation.
 * @param N Size of the matrix (number of rows or columns in A).
 */
void EigenSolve(kids_real* E, kids_complex* T, kids_complex* A, size_t N);

/**
 * Solve the eigenvalue problem for general complex matrices.
 *
 * @param E Pointer to the array to store the eigenvalues.
 * @param T Pointer to the temporary complex matrix.
 * @param A Pointer to the complex matrix for eigenvalue computation.
 * @param N Size of the matrix (number of rows or columns in A).
 */
void EigenSolve(kids_complex* E, kids_complex* T, kids_complex* A, size_t N);

/**
 * Compute the pseudo-inverse of a matrix for real numbers.
 *
 * @param A Pointer to the matrix A.
 * @param invA Pointer to the matrix to store the pseudo-inverse of A.
 * @param N Size of the matrix (number of rows or columns in A).
 * @param e Tolerance value for singular values.
 */
void PseudoInverse(kids_real* A, kids_real* invA, size_t N, kids_real e = 1E-5);

/**
 * Compute the inverse of a matrix for real numbers.
 *
 * @param invA Pointer to the matrix to store the inverse of A.
 * @param A Pointer to the matrix A.
 * @param N Size of the matrix (number of rows or columns in A).
 */
void MatrixInverse(kids_real* invA, kids_real* A, size_t N);

/**
 * Compute the inverse of a matrix for complex numbers.
 *
 * @param invA Pointer to the matrix to store the inverse of A.
 * @param A Pointer to the matrix A.
 * @param N Size of the matrix (number of rows or columns in A).
 */
void ARRAY_INV_MAT(kids_complex* invA, kids_complex* A, size_t N);

void ARRAY_EXP_MAT_GENERAL(kids_complex* expkA, kids_complex* A, kids_complex k, size_t N);

void ARRAY_CORRECT_U(kids_complex* U, size_t N);

void ARRAY_TRANSPOSE(kids_real* A, size_t N1, size_t N2);

void ARRAY_TRANSPOSE(kids_complex* A, size_t N1, size_t N2);


};  // namespace PROJECT_NS

#endif  // KIDS_LINALG_BIND_EIGEN_NOT_USE_TEMPLATE

#endif  // KIDS_LINALG_H
