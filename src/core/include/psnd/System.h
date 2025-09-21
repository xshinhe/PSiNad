/**@file        System.h
 * @brief       provide System class
 * @details     System is a wrapper over Model Kernel with some private data
 *
 * @author      Xin He
 * @date        2024-04
 * @version     1.0
 * @copyright   GNU Lesser General Public License (LGPL)
 *
 *              Copyright (c) 2024 Xin He, Liu-Group
 *
 *  This software is a product of academic research conducted by Professor Liu's
 *  Group at the College of Chemistry and Molecular Engineering, Peking University.
 *  All rights are reserved by Peking University.
 *  You should have received a copy of the GNU Lesser General Public License along
 *  with this software. If not, see <https://www.gnu.org/licenses/lgpl-3.0.en.html>
 **********************************************************************************
 * @warning    Do not include this file to any header. You'd better include it only
 *  in source files!
 * @par [logs]:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-05-14  <td> moved here for initial version.
 * </table>
 *
 **********************************************************************************
 */


#ifndef PSND_SYSTEM_H
#define PSND_SYSTEM_H

#include "psnd/DataSet.h"
#include "psnd/Kernel.h"
#include "psnd/Model.h"
#include "psnd/Param.h"

namespace PROJECT_NS {

class System {
   public:
    System(std::shared_ptr<Model> model, std::shared_ptr<Param> PM, std::shared_ptr<DataSet> DS);

    inline std::shared_ptr<Model> getModel() { return _model; }

    inline std::shared_ptr<Param> getParam() { return _param; }

    inline std::shared_ptr<DataSet> getDataSet() { return _dataset; }

   private:
    std::shared_ptr<Model>              _model;
    std::shared_ptr<Param>              _param;
    std::shared_ptr<DataSet>            _dataset;
    int                                 _dim   = 1;
    int                                 _natom = 1;
    int                                 _nmole = 1;
    std::vector<std::vector<psnd_real>> _boxvectors;
    std::vector<psnd_real>              _mass;
    // std::vector<std::shared_ptr<ConstraintInfo>> constraints;
    // std::vector<std::shared_ptr<Force>>          forces;
    // std::vector<std::shared_ptr<VirtualSite>>    virtualSites;
};
};  // namespace PROJECT_NS

#endif  // PSND_SYSTEM_H