/**@file        Kernel_Elec_Functions.h
 * @brief       this file provides Kernel_Elec_Functions class for public electonic data
 *
 * @author      Xin He
 * @date        2024-03
 * @version     1.0
 * @copyright   GNU Lesser General Public License (LGPL)
 *
 *              Copyright (c) 2024 Xin He, Liu-Group
 *
 *  This software is a product of academic research conducted by Professor Liu's
 *  Group at the College of Chemistry and Molecular Engineering, Peking University.
 *  All rights are reserved by Peking University.
 *  You should have received a copy of the GNU Lesser General Public License along
 *  with this software. If not, see <https://www.gnu.org/licenses/lgpl-3.0.en.html>
 **********************************************************************************
 * @par revision:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-06     <td> updated.
 * </table>
 **********************************************************************************
 */

#ifndef Kernel_Elec_Functions_H
#define Kernel_Elec_Functions_H

#include "psnd/Kernel.h"

namespace PROJECT_NS {

/**
 * this class implements calculation/utils for electronic DOFs:
 */
class Kernel_Elec_Functions final : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    psnd_real gamma1, gamma2, xi1, xi2;
    psnd_bool use_fall;
    psnd_bool use_fssh;
    psnd_bool use_sqc;
    psnd_bool use_cv;
    psnd_bool check_cxs;
    psnd_bool use_strange_win;
    psnd_int  sqc_init;

    int                occ0;
    span<psnd_int>     occ_nuc;
    span<psnd_complex> rho_ele, rho_ele_init;  ///< electronic density
    span<psnd_real>    T, T_init;

    span<psnd_complex> w;  ///< initial measurement of the phase point
    span<psnd_complex> wz_A;
    span<psnd_complex> wz_D;
    span<psnd_complex> ww_A;
    span<psnd_complex> ww_D;
    span<psnd_complex> ww_A_init;
    span<psnd_complex> ww_D_init;
    span<psnd_complex> w_AA, w_AD, w_DD, w_CC, w_CP, w_PP;
    span<psnd_complex> K0;    ///< partial version of K0
    span<psnd_complex> K1;    ///< partial version of K1
    span<psnd_complex> K2;    ///< partial version of K2
    span<psnd_complex> K1QA;  ///< Simplex Quantization
    span<psnd_complex> K2QA;  ///< Heaviside Quantization
    span<psnd_complex> K1DA;
    span<psnd_complex> K2DA;
    span<psnd_complex> K1QD;  ///< Simplex Quantization
    span<psnd_complex> K2QD;  ///< Heaviside Quantization
    span<psnd_complex> K1DD;
    span<psnd_complex> K2DD;

    span<psnd_complex> KSHA;
    span<psnd_complex> KTWA;
    span<psnd_complex> KTWD;

    span<psnd_real> sqcw, trKTWA, trKTWD;

    span<psnd_complex> OpA, OpB;
    span<psnd_complex> TrK1A, TrK2B;

   private:
    void setInputParam_impl(std::shared_ptr<Param> PM);

    void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    Status &initializeKernel_impl(Status &stat);

    Status &executeKernel_impl(Status &stat);
};

};  // namespace PROJECT_NS

#endif  // Kernel_Elec_Functions_H
