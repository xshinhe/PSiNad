/**@file        Sampling_Nucl.h
 * @brief       this file provides Sampling_Nucl class for nuclear sampling
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

#ifndef Sampling_Nucl_H
#define Sampling_Nucl_H

#include "kids/Kernel.h"
#include "kids/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(NuclearSamplingPolicy,
              Fix,           // from fixed coordinate and zero momentum
              Fix2,          // from fixed coordinate and fixed momentum
              ClassicalHO,   // from classical Harmonic Oscillator
              WignerHO,      // from wigner Harmonic Oscillator
              QcHO,          // from quasi-classical Harmonic Oscillator
              ClassicalNMA,  // from classical normal-mode analysis
              WignerNMA,     // from wigner normal-mode analysis
              QcNMA,         // from quasi-classical normal-mode analysis
              Gaussian,      // from gaussian distribution
              ReadDataSet,   // read from dataset file
              ReadAmberRST,   // read from Amber's RST file
              ReadXYZ);      // read from XYZ format file

class Sampling_Nucl final : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Sampling_Nucl(){};

   private:
    NuclearSamplingPolicy::_type sampling_type;
    std::string                  sampling_file;
    std::string                  ignore_nma;
    kids_real                    beta;
    kids_int                     screen_hfreq_type;

    span<kids_real> x, p;
    span<kids_real> x0, p0, x_sigma, p_sigma, w, mass, Tmod;

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // Sampling_Nucl_H
