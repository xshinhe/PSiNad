#!/usr/bin/env python3
#   Coding=utf-8

#   KIDS SCRIPTS
#   Author: xshinhe
#   
#   Copyright (c) 2024 PeKing Univ. - GNUv3 License

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#####################################################################################################

from argparse import ArgumentParser, RawTextHelpFormatter
from enum import Enum, unique
# from pprint import pprint
# from traceback import format_exc

@unique
class CalculationType(Enum):
    NAD = 'nad'
    ENERGY = 'energy'
    GRADIENT = 'gradient'
    HESSIAN = 'hessian'
    HESSLOG2DS = 'hesslog2ds'
    level0 = '0'
    level1 = '1'
    level2 = '2'
    level3 = '3'
    level4 = '4'
    level5 = '5'

@unique
class QMSolverType(Enum):
    ADF = 'adf'
    BAGEL = 'bagel'
    BDF = 'bdf'
    COLUMBUS = 'columbus'
    GAMESS = 'gamess'
    GAUSSIAN = 'gaussian'
    MNDO = 'mndo'
    MOLCAS = 'molcas'
    MOLPRO = 'molpro'
    ORCA = 'orca'
    PYSCF = 'pyscf'
    QCHEM = 'qchem'
    SHARC = 'sharc'
    TURBOMOLE = 'turbomole'

    @property
    def description(self):
        descriptions = {
        'adf': 
'''Amsterdam Density Functional (ADF) program is a powerful computational 
chemistry software for modeling molecules, surfaces, and solids. It uses 
density functional theory (DFT) for electronic structure calculations. 
For more information, visit: https://www.scm.com/amsterdam-modeling-suite/ADF/''',
        'bagel': 
'''BAGEL is an open-source quantum chemistry program that specializes in 
highly accurate ab initio methods. It is designed for parallel computing 
environments. For more information, visit: https://github.com/nubakery/bagel/''',
        'bdf': 
'''Beijing Density Functional (BDF) program is a quantum chemistry software 
package developed in China, which includes a variety of methods for electronic 
structure calculations. 
For more information, visit: https://bdf-manual.readthedocs.io/''',
        'columbus': 
'''COLUMBUS is a quantum chemistry software package for highly accurate multi-reference 
configuration interaction (MRCI) calculations. It is widely used for studying excited 
states and potential energy surfaces. 
For more information, visit: https://columbus-program-system.gitlab.io/columbus/''',
        'gamess': 
'''General Atomic and Molecular Electronic Structure System (GAMESS) is a general-purpose 
quantum chemistry software package. It offers a wide range of methods for electronic 
structure calculations. For more information, visit: http://www.msg.chem.iastate.edu/gamess/''',
        'gaussian': 
'''Gaussian is a widely used computational chemistry software package which provides 
state-of-the-art capabilities for electronic structure modeling. 
For more information, visit: http://gaussian.com/''',
        'mndo': 
'''Modified Neglect of Differential Overlap (MNDO) is a semiempirical quantum chemistry 
program written by W. Thiel, with contributions from M. Beck, S. Billeter, R. Kevorkiants, 
M. Kolb, A. Koslowski, S. Patchkovskii, A. Turner, E.-U. Wallenborn, and W. Weber et al. 
For more information, visit: https://en.wikipedia.org/wiki/MNDO''',
        'molcas': 
'''MOLCAS is a quantum chemistry software package that provides tools for multi-configurational 
SCF calculations and multi-reference perturbation theory. 
For more information, visit: http://www.molcas.org/''',
        'molpro': 
'''MOLPRO is a comprehensive system of ab initio programs for advanced molecular electronic 
structure calculations. It is particularly well-suited for highly accurate calculations on 
small to medium-sized molecules. For more information, visit: https://www.molpro.net/''',
        'orca': 
'''ORCA is an ab initio quantum chemistry program package that contains a wide variety of methods 
for electronic structure calculations. It is known for its efficiency and user-friendly interface.
For more information, visit: https://orcaforum.kofo.mpg.de/''',
        'pyscf': 
'''PySCF is an open-source Python library for quantum chemistry, which supports a range of 
methods including Hartree-Fock, DFT, and post-Hartree-Fock methods. 
For more information, visit: https://pyscf.org/''',
        'qchem': 
'''Q-Chem is a comprehensive ab initio quantum chemistry package. It offers a wide range of 
methods for electronic structure calculations and is widely used in both academic and 
industrial research. For more information, visit: http://www.q-chem.com/''',
        'sharc': 
'''SHARC (Surface Hopping including ARbitrary Couplings) is a program for nonadiabatic dynamics 
simulations. It is designed for studying the dynamics of molecules in excited electronic states. 
For more information, visit: http://sharc-md.org/''',
        'turbomole': 
'''TURBOMOLE is a quantum chemistry program package, particularly well-suited for large-scale 
ground state and excited state calculations. 
For more information, visit: https://www.turbomole.org/''',
        }
        return descriptions[self.value]
    
class MMSolverType(Enum):
    AMBER = 'amber'
    GROMACS = 'gromacs'
    OPENMM = 'openmm'

    @property
    def description(self):
        descriptions = {
        'amber': 
'''Amber is a suite of biomolecular simulation programs. it is used to simulate the physical 
movements of atoms and molecules, allowing researchers to analyze the structures and dynamics 
of complex biomolecular systems. 
For more information, visit: https://ambermd.org/''',
        'gromacs': 
'''GROMACS is a versatile package to perform molecular dynamics, i.e. simulate the Newtonian 
equations of motion for systems with hundreds to millions of particles. it is primarily designed 
for biochemical molecules like proteins, lipids, and nucleic acids. 
For more information, visit: http://www.gromacs.org/''',
        'openmm': 
'''OpenMM is a high-performance toolkit for molecular simulation. it can be used as a library, 
or as an application by itself. it provides a combination of extreme flexibility and high 
performance, and is used in a variety of applications, including drug discovery and materials 
science. 
For more information, visit: https://openmm.org/''',
        }
        return descriptions[self.value]
    

parser = ArgumentParser(description='Execute <KIDS SCRIPTS> Calculation',
    formatter_class=RawTextHelpFormatter)

parser.add_argument('-d', '--directory', dest='directory', nargs='?', 
    default='.', 
    type=str,
    help='specify working directory')
parser.add_argument('-i', '--input', dest='input', nargs='?', 
    default='QM.in',
    type=str,
    help='''input file in TOML format, detailing molecular structure and computational parameters. 
example:
```
    [GEOM] # comment: input molecule geometry
    xyz = """
    2

    H    0.0    0.0    0.0
    H    0.0    0.0    1.0
    """

    [QM]
    exec = "MNDO" # point to parameters in [QM.MNDO]

    [QM.MNDO]
    path = "/usr/bin/mndo99" # exetutable file
    N = 6         # nuclear DoF. (=3*Natom)
    F = 2         # eletronic DoF.

    # default setups for interfaced calculation
    level0 = """JOP=-2 IOP=-6 IGEOM=1 IFORM=1 ICUTS=-1 ICUTG=-1 +
    ISCF=9 IPLSCF=9 DPREC=1D-8 DSTEP=1D-5 IPRINT=1 +
    NCIGRD=2 IEF=1 IPREC=100 +
    IMULT=1 IMOMAP=1 IUVCD=2 +
    KCI=5 IOUTCI=1 MPRINT=1 ICROSS=7 +
    MOVO=0 ICI1=10 ICI2=8 NCIREF=3 MCIREF=0 LEVEXC=2 IROOT=3 KITSCF=2000


    $COORD_XYZ
    0 0.00000000 0 0.00000000 0 0.00000000 0
    1 2
    """

    # next level setups for interfaced calculation
    level1 = "..."
```
''')
parser.add_argument('-qm', '--qmsolver', dest='qmsolver', nargs='?', 
    default=QMSolverType.MNDO.value, type=str,
    choices=[s.value for s in QMSolverType],
    help='ab initio QM solver. choices are:\n' + '\n'.join(
        ['[' + s.value + ']:\n' + s.description for s in QMSolverType]))
parser.add_argument('-mm', '--mmsolver', dest='mmsolver', nargs='?', 
    default=MMSolverType.AMBER.value, type=str,
    choices=[s.value for s in MMSolverType],
    help='molecule mechanics solver. choices are:\n' + '\n'.join(
        ['[' + s.value + ']:\n' + s.description for s in MMSolverType]))
parser.add_argument('-top', '--topology', dest='topology', nargs='?', 
    default='real.top,model-H.top',     # support format: crd, xyz, etc.
    type=str,
    help='topology files for MM (seperated by comma). (default: real.top,model-H.top)')
parser.add_argument('-c', '--coord', dest='coord', nargs='?', 
    default='',
    type=str,
    help='coordinate file. if not provided, coordinate will be initialized from TOML file. (default null)')
parser.add_argument('-l', '--layer', dest='layer', nargs='?', 
    default='layer.info',
    type=str,
    help='layer information file. if layer is null/blank, all particles are treated in QM (default: layer.info)')
parser.add_argument('-t', '--type', dest='type', nargs='?', 
    default=CalculationType.level0.value, 
    type=str,
    choices=[t.value for t in CalculationType],
    help='calculation type. Choices are: ' + ', '.join([t.value for t in CalculationType]))
parser.add_argument('-o', '--output', dest='output', nargs='?', default='QM.out', type=str,
    help='specify the name of the output file (default: QM.out)')
parser.add_argument('-vb', '--verbosity', dest='verbosity', nargs='?', default=1, type=int,
    help='set the verbosity level of the output (default: 1)')
parser.add_argument('-ds', '--datasetfn', dest='datasetfn', nargs='?', default='interface.ds', type=str,
    help='specify the dataset file to dump. (default: interface.ds)')

if __name__ == '__main__':
    args = parser.parse_args()
