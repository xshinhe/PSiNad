/**@file        Kernel_Load_DataSet.h
 * @brief       this file provides Kernel_Load_DataSet class enabling load dataset.
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

#ifndef Kernel_Load_DataSet_H
#define Kernel_Load_DataSet_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

/**
 * this class implements a process for loading data in state stucture
 */
class Kernel_Load_DataSet : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

   private:
    std::string load_fn;

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Load_DataSet_H
