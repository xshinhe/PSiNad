#!/usr/bin/env python3
# coding=utf-8

#    COBRAMM
#    Copyright (c) 2019 ALMA MATER STUDIORUM - Università di Bologna

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

# import statments of module from python standard library

import argparse  # command line arguments parser
import sys  # System-specific parameters and functions
import os  # filesystem utilities

# hijack $PYTHONPATH (sys.path) to find cobramm local modules
try:
    sys.path.append(os.path.join(os.environ['COBRAM_PATH'], 'cobramm'))
except KeyError:
    raise RuntimeError('cannot import cobramm module.\n'
                       'Please check if $COBRAM_PATH is defined in your environment, if not run:\n'
                       'export COBRAM_PATH=/path/to/COBRAMM')

# imports of local modules

import logwrt  # write messages and output to log (std out or cobramm.log)
import cobrammenv  # management of the enviroment for running COBRAMM (and utils)

# imports of user-defined classes

from amberCalculator import AmberCalculator  # Amber wrapper


#####################################################################################################

def main():

    # this new class defines a style for the parser message that join the two predefined styles
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    # Define command line arguments with argparse
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="COBRAMM utility script: prepare AMBER input files for a solvated molecule",
        epilog="the script launches a sequence of amberTools programs to construct \n"
               "an MM model of a molecule (described with GAFF2 FF) within a periodic\n"
               "box of solvent (described with specific ad-hoc FF)\n\n",
        formatter_class=CustomFormatter)

    parser.add_argument("--solvent", choices=["water", "methanol", "chloroform", "DMSO", "acetonitrile"], default="water",
                        help="name of the molecule for solvation (default: water)")
    parser.add_argument("-xyz", "--xyzfile", required=True, metavar='XYZ_FILENAME', dest="xyzfile",
                        help="name of the input xyz file defining the molecular structure of the chromophore")
    parser.add_argument("-sz", "--solvationsize", required=True, metavar='SOLVATION_SIZE', dest="boxsize",
                        help="size of the solvation layer in Angstrom")
    parser.add_argument("-q", "--netcharge", metavar="NET CHARGE", dest="netcharge", default=0,
                        type=int, help="net charge of the chromophore")
    parser.add_argument('--version', action='version',
                        version='%(prog)s is part of COBRAMM {}'.format(cobrammenv.getVersion()))
    parser.add_argument("-ch", "--charges", metavar='QM_CALCULATION_FILENAME', dest="chargesfile",
                        help="path to output of a QM single point calculation containing ESP charges to extract for the chromophore")
    args = parser.parse_args()

    # source COBRAMM configuration file (if available)
    confFile = cobrammenv.setCobrammProfile()
    # Print to log information on the environment configuration
    logwrt.writelog(cobrammenv.checkConfigFile(confFile))

    # initialize Amber wrapper class
    AmberCalculator()

    # ===================================================================

    logwrt.writelog("\nReading geometry from file {0} and preparing input file for solvated molecule".format(
        args.xyzfile))

    with open(args.xyzfile) as f:
        f.readline(), f.readline()
        atomlabels, coords = [], []
        for ln in f:
            if ln.split():
                atomlabels.append(ln.split()[0])
                coords.append([float(f) for f in ln.split()[1:]])

    # create solvated molecule AMBER files with the method of ClassCalculator

    if args.chargesfile:
        logwrt.writelog("\nReading external QM charges from calculation {0}".format(args.chargesfile))
        charges = os.path.abspath(args.chargesfile)
        snapshot, topology = AmberCalculator.createSolvatedMolecule(atomlabels, coords, args.boxsize, solvent=args.solvent, gau_charges=True,
                                                                    gau_file=charges, net_charge=args.netcharge)
    else:
        snapshot, topology = AmberCalculator.createSolvatedMolecule(atomlabels, coords, args.boxsize, solvent=args.solvent, net_charge=args.netcharge)

    # define names for the topology and coordinate output files
    basename = os.path.splitext(args.xyzfile)[0]

    # write coordinates and topology
    with open(basename + ".top", "w") as f:
        f.write(topology.topotext)
    with open(basename + ".crd", "w") as f:
        f.write(snapshot.crdtext)

    logwrt.writelog("\nThe AMBER input files for the solvated molecule have been created with success!\n"
                    "Topology written to {0} and coordinates "
                    "with octahedral BPC written to {1}\n\n".format(basename + ".top", basename + ".crd"))


#####################################################################################################

if __name__ == '__main__':
    main()
