/**@file        Kernel_Read_Dimensions.h
 * @brief       this file provides Kernel_Read_Dimensions class enabling reading
 *              and initializing the sizes and shapes of various variables
 *
 *              please refer to vars_list.h for more details.
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


#ifndef Kernel_Read_Dimensions_H
#define Kernel_Read_Dimensions_H

#include "kids/Kernel.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

/**
 * This class specifies the kernels that initialize the dimensions of data
 */
class Kernel_Read_Dimensions final : public Kernel {
   public:
    Kernel_Read_Dimensions() : Kernel(){};

    inline virtual const std::string name() { return "Kernel_Read_Dimensions"; }

   private:
    virtual void setInputParam_impl(std::shared_ptr<Param>& PM) {
        Dimension::M = PM->get_int("M", LOC(), 1);
        Dimension::P = PM->get_int("P", LOC(), 1);
        Dimension::N = PM->get_int("N", LOC(), 1);
        Dimension::F = PM->get_int("F", LOC(), 1);
        Dimension::static_build_shapes();
    };
};

};  // namespace PROJECT_NS

#endif  // Kernel_Read_Dimensions_H