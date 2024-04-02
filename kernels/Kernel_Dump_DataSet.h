/**@file        Kernel_Dump_DataSet.h
 * @brief       this file provides Kernel_Dump_DataSet class enabling dump dataset.
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

#ifndef Kernel_Dump_DataSet_H
#define Kernel_Dump_DataSet_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

/**
 * this class implements a process for dumping current data with stuctured format
 */
class Kernel_Dump_DataSet : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Dump_DataSet"; }

   private:
    std::string directory;  ///< path for dumping
    std::string fn;         ///< filename (stamp) for dumping
    std::string hdlr_str;   ///< handler type specifier

    virtual void read_param_impl(Param* PM);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Dump_DataSet_H
