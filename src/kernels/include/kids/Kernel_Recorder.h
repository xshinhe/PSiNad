/**@file        Kernel_Recorder.h
 * @brief       this file provides Kernel_Recorder class for trace data in dataset
 *              during the dynamics.
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

#ifndef Kernel_Recorder_H
#define Kernel_Recorder_H

#include "kids/Kernel.h"
#include "kids/RuleEvaluator.h"

namespace PROJECT_NS {

class Kernel_Recorder final : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Kernel_Recorder();

    virtual ~Kernel_Recorder();

   private:
    // friend class Kernel_Report;

    int*                     istep_ptr;
    int*                     sstep_ptr;
    int*                     isamp_ptr;
    int*                     nsamp_ptr;
    kids_bint*               at_samplingstep_initially_ptr;
    double                   t0, dt, time_unit;
    std::vector<std::string> opened_files;

    virtual void token(Param::JSON& j);

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual Status& executeKernel_impl(Status& stat);

    virtual Status& finalizeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Recorder_H
