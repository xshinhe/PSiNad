/**@file        Sampling_Elec.h
 * @brief       this file provides Sampling_Elec class for electronic sampling
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

#ifndef Sampling_Elec_H
#define Sampling_Elec_H

#include "kids/Kernel.h"
#include "kids/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(ElectronicSamplingPolicy,
              Focus,         // Focus sampling / Langer Modification
              GDTWA,         // GDTWA sampling
              SQCsqr,        // square SQC sampling
              SQCtri,        // triangle SQC sampling (TWF)
              SQCspx,        // simplex SQC sampling
              SQCtest01,     // used for sqc test 1
              SQCtest02,     // used for sqc test 2
              SQCtest03,     // used for sqc test 3
              Gaussian,      // Gaussian sampling
              Constraint,    // Constraint sampling
              ReadDataSet);  // Read from dataset


class Sampling_Elec final : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Sampling_Elec(){};

   private:
    ElectronicSamplingPolicy::_type sampling_type;
    std::string                     sampling_file;


    kids_int  occ0;
    kids_real gamma1, xi1;
    bool      use_cv   = true;   // adapt cv in rho_nuc
    bool      use_wmm  = false;  // in this case, gamma1 will be used as delta in wMM
    bool      use_sum  = false;
    bool      use_fssh = false;

    kids_real *alpha, *V;

    kids_int*     occ_nuc;
    kids_real*    T;
    kids_complex *c, *rho_ele, *rho_nuc;
    kids_complex* w;


    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // Sampling_Elec_H
