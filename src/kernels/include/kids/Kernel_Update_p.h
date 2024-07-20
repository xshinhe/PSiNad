/**@file        Kernel_Update_p.h
 * @brief       this file provides Kernel_Update_p class for update p
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

#ifndef Kernel_Update_p_H
#define Kernel_Update_p_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

class Kernel_Update_p : public Kernel {
   public:
    Kernel_Update_p(double scale) : Kernel(), scale{scale} {};

    virtual const std::string getName();

    virtual int getType() const;

   private:
    double *p, *f, *minv, *ve;
    double* Ekin;
    double  scale, *dt;

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual Status& executeKernel_impl(Status& stat);
};


};  // namespace PROJECT_NS


#endif  // Kernel_Update_p_H
