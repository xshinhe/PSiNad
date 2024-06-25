/**@file        Kernel_Elec_NAD.h
 * @brief       this file provides Kernel_Elec_NAD class for electronic dynamics
 *              and properties in nonadiabatic trajectory dynamics.
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
 * <tr><td> 2024-04-02  <td> Initial version.
 * </table>
 **********************************************************************************
 */

#ifndef Kernel_Elec_NAD_H
#define Kernel_Elec_NAD_H

#include "kids/Kernel.h"
#include "kids/Kernel_Elec.h"
#include "kids/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(NADPolicy,
              EHR,    // Ehrenfest Dynamics
              BOSH,   // BO dynamics & hopping according to W(\rho)
              CVSH,   // CV dynamics & hopping acoording to W(\rho)
              BOSD,   // BO dynamics & smoothing acoording to W(\rho)
              CVSD);  // CV dynamics & smoothing acoording to W(\rho)

/**
 * @brief initialization kernel for electonic DOFs in NAD
 */
class Kernel_Elec_NAD final : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Kernel_Elec_NAD(double scale = 1.0e0) : Kernel(), scale{scale} {
        appendChild(std::shared_ptr<Kernel_Elec>(new Kernel_Elec()));  //
    }

   private:
    NADPolicy::_type cmsh_type;

    kids_real gamma1, gamma2, xi1, xi2;
    bool      use_focus = false;
    bool      use_cv    = true;   // adapt cv in rho_nuc
    bool      use_wmm   = false;  // in this case, gamma1 will be used as delta in wMM
    bool      use_fall  = false;
    bool      use_gdtwa = false;
    bool      use_sum   = false;

    bool cread_from_ds        = false;
    bool disable_inner_switch = false;

    bool reflect = true;  // treatment in hopping
    int  hopping_type1;
    int  hopping_type2;

    bool dynamic_alpha;

    double        scale;
    double        dt;
    kids_real     alpha0;
    kids_real*    alpha;
    kids_real*    p;
    kids_real*    m;
    kids_real*    fadd;
    kids_real*    ftmp;
    kids_real*    direction;
    kids_real *   vpes, *V, *E, *dE, *T;
    kids_real*    Epot;
    kids_complex* H;
    kids_complex* wrho;
    kids_real*    sqcw;
    kids_real *   sqcIA, *sqcID;

    kids_bint* at_samplingstep_finally_ptr;

    //
    bool use_sqc;
    int  sqc_init;
    bool only_adjust;
    bool check_cxs;

    bool use_fssh;
    bool use_strange_win;

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // Kernel_Elec_NAD_H
