/**@file        Kernel_NADForce.h
 * @brief       this file provides Kernel_NADForce class enabling force weighting
 *              from electronic properties.
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


#ifndef Kernel_NADForce_H
#define Kernel_NADForce_H

#include "../core/Kernel.h"
#include "../core/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(NADForcePolicy,  //
              EHR,             //
              BO,              //
              CV,              //
              BOSD,            //
              CVSD,            //
              ELSE);           //

namespace FORCE_OPT {
extern bool BATH_FORCE_BILINEAR;
extern int  nbath;
extern int  Nb;
};  // namespace FORCE_OPT

/**
 * this class implements process of calculation of nuclear force for nonadiabatic dynamics
 */
class Kernel_NADForce : public Kernel {
   public:
    static NADForcePolicy::_type     NADForce_type;
    inline virtual const std::string name() { return "Kernel_NADForce"; }

   private:
    bool offd_projected;

    kids_real *f, *grad, *dV, *dE, *Force, *T;
    kids_real *p, *m;
    kids_real *fadd, *fproj;

    bool *succ_ptr;

    virtual void read_param_impl(Param *PM);

    virtual void init_data_impl(DataSet *DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS

#endif  // Kernel_NADForce_H
