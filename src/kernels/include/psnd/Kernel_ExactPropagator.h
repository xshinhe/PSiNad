/**@file        Kernel_ExactPropagator.h
 * @brief       this file provides Kernel_NAForce class enabling force weighting
 *              from electronic properties.
 *
 * @author      Baihua Wu, Haocheng Lu
 * @date        2025-08
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
 * <tr><td> 2025-08-01  <td> Initial version.
 * <tr><td> 2025-10-29  <td> Updated version.
 * </table>
 **********************************************************************************
 */


#ifndef Kernel_ExactPropagator_H
#define Kernel_ExactPropagator_H

#include "psnd/Kernel.h"
#include "psnd/Kernel_NAForce.h"

namespace PROJECT_NS {

class Kernel_ExactPropagator : public Kernel {
   public:
    Kernel_ExactPropagator() : Kernel(), scale{1.0} {};  // 默认构造函数
    Kernel_ExactPropagator(double scale) : Kernel(), scale{scale} {};
    virtual ~Kernel_ExactPropagator() = default;  // 虚析构函数

    virtual const std::string getName();

    virtual int getType() const;

   private:
    double             scale;
    span<psnd_real>    p, f, fadd, minv, ve;
    span<psnd_real>    T, eig, dE, dV, ddV, nac, grad, hess;
    span<psnd_complex> mask, dmask;
    span<psnd_complex> c, rho_nuc, rho_ele;
    span<psnd_real>    Ekin;
    span<psnd_real>    dt;

    span<psnd_int>     occ_nuc;
    span<psnd_real>    ForceMat, EMat, V;
    span<psnd_real>    m;
    span<psnd_real>    ftmp, fproj;
    span<psnd_real>    alpha;
    span<psnd_complex> wrho;

    NAForcePolicy::_type NAForce_type;


    span<psnd_real> Epot, vpes;
    span<psnd_real> dt_ptr;

    span<psnd_complex> MFFtmp1, MFFtmp2, MFFtmp3, MFFtmp4, MFFtmp5, MFFtmp6;

    // span<psnd_real> B_vec;

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status &initializeKernel_impl(Status &stat);

    virtual Status &executeKernel_impl(Status &stat);
};
}  // namespace PROJECT_NS

#endif  // Kernel_ExactPropagator_H
