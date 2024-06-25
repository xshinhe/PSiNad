/**@file        Kernel_Prioritization.h
 * @brief       this file provides Kernel_Prioritization class enabling reordering
 *              different kernels in parsing Param, connecting DataSet and initial
 *              calculation.
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


#ifndef Kernel_Prioritization_H
#define Kernel_Prioritization_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

/**
 * This class specifies the kernels that prioritizing other kernels
 */
class Kernel_Prioritization final : public Kernel {
   public:
    Kernel_Prioritization(std::vector<std::shared_ptr<Kernel>> kers, int ptype_in);

    virtual const std::string getName();

    virtual int getType() const;

   private:
    int                                  ptype;
    std::vector<std::shared_ptr<Kernel>> _ref_kernels;

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& initializeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // Kernel_Prioritization_H