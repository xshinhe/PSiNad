/**@file        Kernel_NAForce.h
 * @brief       this file provides Kernel_NAForce class enabling force weighting
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
 * <tr><td> 2024-06-15  <td> Updated version.
 * </table>
 **********************************************************************************
 */


#ifndef Kernel_NAForce_H
#define Kernel_NAForce_H

#include "kids/Kernel.h"
#include "kids/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(NAForcePolicy,  //
              EHR,            //
              BO,             //
              NAF,            //
              BOSD,           //
              NAFSD,          //
              ELSE);          //

namespace FORCE_OPT {
extern bool BATH_FORCE_BILINEAR;
extern int  nbath;
extern int  Nb;
};  // namespace FORCE_OPT

/**
 * this class implements process of calculation of nuclear force for nonadiabatic dynamics
 */
class Kernel_NAForce : public Kernel {
   public:
    static NAForcePolicy::_type NAForce_type;

    virtual const std::string getName();

    virtual int getType() const;

   private:
    bool offd_projected;

    kids_real *   f, *grad, *dV, *dE, *ForceMat, *EMat, *T, *V;
    kids_real *   p, *m;
    kids_real *   fadd, *ftmp, *fproj;
    kids_real*    alpha;
    kids_complex* wrho;

    kids_real *Epot, *vpes;
    kids_real* dt_ptr;

    kids_int*     occ_nuc;
    kids_complex *rho_ele, *rho_nuc;

    kids_bint* succ_ptr;

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // Kernel_NAForce_H
