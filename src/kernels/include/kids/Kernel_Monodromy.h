/**@file        Kerner_Monodromy.h
 * @brief       this file provides Kerner_Monodromy class for update p
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
 * <tr><td> 2024-04-02  <td> Initial version. Added detailed commentary by ChatGPT.
 * </table>
 **********************************************************************************
 */

#ifndef Kernel_Monodromy_H
#define Kernel_Monodromy_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

class Kernel_Monodromy  // : public Kernel
{
   public:
    static bool enable;

    Kernel_Monodromy(){};

    /*
        void MQC_prefactor() {
            for (int iP = 0; iP < Dimension::P; ++iP) {
                Mpp = ;  //
                Mpq = ;  //
                Mqp = ;  //
                Mqq = ;  //
            }

            M1pp = 0;  //??

            for (int i = 0, ik = 0; i < Dimension::N + Dimension::F; ++i)
                for (int k = 0; k < Dimension::N + Dimension::F; ++k, ++ik)
                    Mat1[ik] = mono1[(N2 + i) * N4 + (N2 + k)] - im * gt * mono1[(N0 + i) * N4 + (N2 + k)];
            for (int i = 0, ik = 0; i < Dimension::N + Dimension::F; ++i)
                for (int k = 0; k < Dimension::N + Dimension::F; ++k, ++ik)
                    Mat2[ik] = (i == k) ? 0.5e0 + 0.5e0 / ((cq[i] + g0) * cp[i] + cq[i] * (1.0e0 / g0 + cp[i])) : 0.0e0;
            for (int i = 0, ik = 0; i < Dimension::N + Dimension::F; ++i)
                for (int k = 0; k < Dimension::N + Dimension::F; ++k, ++ik)
                    Mat3[ik] = mono2[(N2 + i) * N4 + (N2 + k)] * gt + im * mono2[(N2 + i) * N4 + (N0 + k)];
            ARRAY_MATMUL3(Mat4, Mat1, Mat2, Mat3, Nx, Nx, Nx, Nx);

            for (int i = 0, ik = 0; i < Dimension::N + Dimension::F; ++i)
                for (int k = 0; k < Dimension::N + Dimension::F; ++k, ++ik)
                    Mat1[ik] = mono1[(N0 + i) * N4 + (N0 + k)] * gt + im * mono1[(N2 + i) * N4 + (N0 + k)];
            for (int i = 0, ik = 0; i < Dimension::N + Dimension::F; ++i)
                for (int k = 0; k < Dimension::N + Dimension::F; ++k, ++ik)
                    Mat2[ik] =
                        (i == k) ? (0.5e0 / g0 + cp[i]) / ((cq[i] + g0) * cp[i] + cq[i] * (1.0e0 / g0 + cp[i])) : 0.0e0;
            for (int i = 0, ik = 0; i < Dimension::N + Dimension::F; ++i)
                for (int k = 0; k < Dimension::N + Dimension::F; ++k, ++ik)
                    Mat3[ik] = mono2[(N2 + i) * N4 + (N2 + k)] * gt + im * mono2[(N2 + i) * N4 + (N0 + k)];
            ARRAY_MATMUL3(Mat5, Mat1, Mat2, Mat3, Nx, Nx, Nx, Nx);

            for (int i = 0, ik = 0; i < Dimension::N + Dimension::F; ++i)
                for (int k = 0; k < Dimension::N + Dimension::F; ++k, ++ik)
                    Mat1[ik] = mono1[(N0 + i) * N4 + (N0 + k)] * gt + im * mono1[(N2 + i) * N4 + (N0 + k)];
            for (int i = 0, ik = 0; i < Dimension::N + Dimension::F; ++i)
                for (int k = 0; k < Dimension::N + Dimension::F; ++k, ++ik)
                    Mat2[ik] = (i == k) ? 0.5e0 + 0.5e0 / ((cq[i] + g0) * cp[i] + cq[i] * (1.0e0 / g0 + cp[i])) : 0.0e0;
            for (int i = 0, ik = 0; i < Dimension::N + Dimension::F; ++i)
                for (int k = 0; k < Dimension::N + Dimension::F; ++k, ++ik)
                    Mat3[ik] = mono2[(N0 + i) * N4 + (N0 + k)] - im * gt * mono2[(N0 + i) * N4 + (N2 + k)];
            ARRAY_MATMUL3(Mat6, Mat1, Mat2, Mat3, Nx, Nx, Nx, Nx);

            for (int i = 0, ik = 0; i < Dimension::N + Dimension::F; ++i)
                for (int k = 0; k < Dimension::N + Dimension::F; ++k, ++ik)
                    Mat1[ik] = mono1[(N2 + i) * N4 + (N2 + k)] - im * gt * mono1[(N0 + i) * N4 + (N2 + k)];
            for (int i = 0, ik = 0; i < Dimension::N + Dimension::F; ++i)
                for (int k = 0; k < Dimension::N + Dimension::F; ++k, ++ik)
                    Mat2[ik] =
                        (i == k) ? (0.5e0 * g0 + cq[i]) / ((cq[i] + g0) * cp[i] + cq[i] * (1.0e0 / g0 + cp[i])) : 0.0e0;
            for (int i = 0, ik = 0; i < Dimension::N + Dimension::F; ++i)
                for (int k = 0; k < Dimension::N + Dimension::F; ++k, ++ik)
                    Mat3[ik] = mono2[(N2 + i) * N4 + (N2 + k)] * gt + im * mono2[(N2 + i) * N4 + (N0 + k)];
            ARRAY_MATMUL3(Mat7, Mat1, Mat2, Mat3, Nx, Nx, Nx, Nx);
        }
    */
};


};  // namespace PROJECT_NS


#endif  // Kernel_Monodromy_H
