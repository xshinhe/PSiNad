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
    virtual const std::string getName();

    virtual int getType() const;

   private:
    double     t0, *t0_ptr;
    double     t, *t_ptr;
    double     dt, *dt_ptr;
    double     tend, *tend_ptr;
    kids_bint* at_samplingstep_initially_ptr;
    kids_bint* at_samplingstep_finally_ptr;

    kids_bint* succ_ptr;
    kids_bint* frez_ptr;
    kids_bint* last_attempt_ptr;
    int*       fail_type_ptr;  // record the failure information (longtime keeped)

    int msize, *tsize_ptr, *dtsize_ptr;
    int sstep, *sstep_ptr;
    int istep, *istep_ptr, nstep, *nstep_ptr;
    int isamp, *isamp_ptr, nsamp, *nsamp_ptr;
    int nbackup;

    double time_unit;

    const std::vector<std::string> backup_fields = {"x", "p", "U", "occ_nuc", "f", "Ekin", "Epot"};

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Iterative_Adapt_H
