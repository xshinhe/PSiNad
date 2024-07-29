/**@file        Kernel_Update_x.h
 * @brief       this file provides Kernel_Update_x class for update x
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

#ifndef Kernel_Update_x_H
#define Kernel_Update_x_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

class Kernel_Update_x final : public Kernel {
   public:
    Kernel_Update_x(double scale) : Kernel(), scale{scale} {};

    virtual const std::string getName();

    virtual int getType() const;

   private:
    double *x, *p, *m, *minv;
    double  scale, *dt;

    kids_real *mono, *monodt;

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& executeKernel_impl(Status& stat);

    void update_monodromy();
};

};  // namespace PROJECT_NS


#endif  // Kernel_Update_x_H
