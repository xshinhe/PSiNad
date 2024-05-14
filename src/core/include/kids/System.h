/**@file        System.h
 * @brief       provide System class
 * @details     System
 *
 * @author      Xin He
 * @date        2024-04
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

#ifndef KIDS_System_H
#define KIDS_System_H

namespace PROJECT_NS {

class System {
    std::shared_ptr<Platform> getPlatfrom() { return _platform; }
    std::shared_ptr<Param>    getParam() { return _param; }
    std::shared_ptr<DataSet>  getDataSet() { return _dataset; }
    std::shared_ptr<System>   getSystem() { return _system; }
    std::shared_ptr<Solver>   getSolver() { return _solver; }

   private:
    double*                                      boxvectors;
    double*                                      masses;
    std::vector<std::shared_ptr<ConstraintInfo>> constraints;
    std::vector<std::shared_ptr<Force>>          forces;
    std::vector<std::shared_ptr<VirtualSite>>    virtualSites;
};
};  // namespace PROJECT_NS

#endif  // KIDS_System_H