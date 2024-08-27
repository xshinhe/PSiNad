/**@file        Kernel_Iterative.h
 * @brief       this file provides Kernel_Iterative class enabling iteration
 *              and conservation.
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

#ifndef Kernel_Iterative_H
#define Kernel_Iterative_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

/**
 * Minimal iterator for integration of the equations of motion
 */
class Kernel_Iterative final : public Kernel {
   public:
    Kernel_Iterative() { enable_call_child = false; }

    virtual const std::string getName();

    virtual int getType() const;

   private:
    int             sstep, nstep, nsamp;
    double          t0, tend, dt0;
    span<kids_real> t, dt;
    span<kids_bint> at_condition;
    span<kids_int>  istep, isamp;

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Iterative_H
