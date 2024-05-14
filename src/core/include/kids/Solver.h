/**@file        Solver.h
 * @brief       provide Solver class
 * @details     Solver
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

#ifndef KIDS_Solver_H
#define KIDS_Solver_H

namespace PROJECT_NS {

class Solver {
    std::shared_ptr<Platform> getPlatfrom() { return _platform; }
    std::shared_ptr<Param>    getParam() { return _param; }
    std::shared_ptr<DataSet>  getDataSet() { return _dataset; }
    std::shared_ptr<System>   getSystem() { return _system; }
    std::shared_ptr<Solver>   getSolver() { return _solver; }

   private:
    std::shared_ptr<Platform> _platform;
    std::shared_ptr<Param>    _param;
    std::shared_ptr<DataSet>  _dataset;
    std::shared_ptr<System>   _system;
    std::shared_ptr<Solver>   _solver;
};
};  // namespace PROJECT_NS

#endif  // KIDS_Solver_H