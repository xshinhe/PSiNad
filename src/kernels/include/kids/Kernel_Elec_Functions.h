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
 *  This software is a product of Xin's PhD research conducted by Professor Liu's
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

#include "kids/Kernel.h"

namespace PROJECT_NS {

/**
 * this class implements calculation/utils for electronic DOFs:
 */
class Kernel_Elec_Functions final : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    kids_real gamma1, gamma2, xi1, xi2;
    kids_bool use_fall;
    kids_bool use_fssh;
    kids_bool use_sqc;
    kids_bool use_cv;
    kids_bool check_cxs;
    kids_bool use_strange_win;
    kids_int  sqc_init;

    int                occ0;
    span<kids_int>     occ_nuc;
    span<kids_complex> rho_ele, rho_ele_init;  ///< electronic density
    span<kids_real>    T, T_init;

    span<kids_complex> w;  ///< initial measurement of the phase point
    span<kids_complex> wz_A;
    span<kids_complex> wz_D;
    span<kids_complex> ww_A;
    span<kids_complex> ww_D;
    span<kids_complex> ww_A_init;
    span<kids_complex> ww_D_init;
    span<kids_complex> w_AA, w_AD, w_DD, w_CC, w_CP, w_PP;
    span<kids_complex> K0;    ///< partial version of K0
    span<kids_complex> K1;    ///< partial version of K1
    span<kids_complex> K2;    ///< partial version of K2
    span<kids_complex> K1QA;  ///< Simplex Quantization
    span<kids_complex> K2QA;  ///< Heaviside Quantization
    span<kids_complex> K1DA;
    span<kids_complex> K2DA;
    span<kids_complex> K1QD;  ///< Simplex Quantization
    span<kids_complex> K2QD;  ///< Heaviside Quantization
    span<kids_complex> K1DD;
    span<kids_complex> K2DD;

    span<kids_complex> KSHA;
    span<kids_complex> KTWA;
    span<kids_complex> KTWD;

    span<kids_real> sqcw, trKTWA, trKTWD;

    span<kids_complex> OpA, OpB;
    span<kids_complex> TrK1A, TrK2B;

   private:
    void setInputParam_impl(std::shared_ptr<Param> PM);

    void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    Status &initializeKernel_impl(Status &stat);

    Status &executeKernel_impl(Status &stat);
};

};  // namespace PROJECT_NS

#endif  // Kernel_Elec_Functions_H
