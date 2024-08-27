/**@file        Kernel_Iterative_Adapt.h
 * @brief       this file provides Kernel_Iterative_Adapt class enabling iteration
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

#ifndef Kernel_Iterative_Adapt_H
#define Kernel_Iterative_Adapt_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

class Kernel_Iterative_Adapt final : public Kernel {
   public:
    Kernel_Iterative_Adapt() { enable_call_child = false; }

    virtual const std::string getName();

    virtual int getType() const;

   private:
    int             nstep, sstep, nsamp;
    double          t0, tend, dt0;
    span<kids_real> t, dt;
    span<kids_bint> at_condition;

    int            msize;
    span<kids_int> tsize, dtsize;
    span<kids_int> istep, isamp;
    int            nbackup;

    double time_unit;

    const std::vector<std::string> backup_fields = {"x", "p", "U", "occ_nuc", "f", "Ekin", "Epot"};

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Iterative_Adapt_H
