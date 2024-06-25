/**@file        Model.h
 * @brief       wrapper/interface for model kernel
 * @details
 *
 * @author      Xin He
 * @date        2024-05
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
 * @par [logs]:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-05     <td> initial version.
 * </table>
 *
 **********************************************************************************
 */

#ifndef KIDS_Model_H
#define KIDS_Model_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

class Model : public Kernel, std::enable_shared_from_this<Model> {
   public:
    inline virtual const std::string getName() { return utils::concat("Kernel__", kernel_name); }

    Model(const std::string& customized_name = "");

    virtual ~Model(){};

   protected:
    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual Status& executeKernel_impl(Status& stat);

    virtual Status& finalizeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // KIDS_Model_H
