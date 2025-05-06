import os, sys

def generate(fin, fout):
    lines = open(fin).readlines()

    f2 = open(fout, 'w')
    f2.write('''/**@file        vars_list.h
 * @brief       declaration of variables used in the program.
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
 * <tr><td> 2024-03-29  <td> Initial version.
 * </table>
 **********************************************************************************
 */

#ifndef VARS_LIST_H
#define VARS_LIST_H

#include "kids/Shape.h"
#include "kids/Types.h"
#include "kids/Variable.h"

namespace PROJECT_NS {\n
''')
    for l in lines:
        terms = l.split()
        if len(terms) == 0: continue

        if len(terms[0]) >= 11 and terms[0][0:11] == 'std::size_t':
            name = terms[1].split(';')[0]
            f2.write('namespace Dimension {\n')
            f2.write(f'extern {l}')
            f2.write('}; // namespace Dimension\n')

        if len(terms[0]) >= 5 and terms[0][0:5] == 'Shape':
            name = terms[1].split('(')[0]
            f2.write('namespace Dimension {\n')
            f2.write(f'extern Shape {name};\n')
            f2.write('}; // namespace Dimension\n')

        if len(terms[0]) >= 8 and terms[0][0:8] == 'VARIABLE':
            names = terms[1].split('(')[1].split(',')[0].split('::')
            f2.write('namespace DATA {\n')
            for k in range(len(names)-1):
                f2.write('namespace %s {\n'%names[k])
            f2.write(f'extern {terms[0]} {names[-1]};')
            for k in range(len(names)-1):
                f2.write('}; // namespace %s \n'%names[k])
            f2.write('}; // namespace DATA\n')

    f2.write('namespace Dimension {\n')
    f2.write('extern void static_build_shapes();')
    f2.write('}; // namespace Dimension\n')

    f2.write('}; // namespace PROJECT_NS\n')
    f2.write('#endif  // VARS_LIST_H\n')
    f2.close()


if __name__ == '__main__':
    filein = sys.argv[1]
    fileout = sys.argv[2]
    generate(filein, fileout)
