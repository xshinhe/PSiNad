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

import argparse    # command line arguments parser
import sys         # System-specific parameters and functions
import os          # filesystem utilities
import shutil      # filesystem utilities

# hijack $PYTHONPATH (sys.path) to find cobramm local modules
try:
    sys.path.insert(0, os.path.join(os.getenv('COBRAM_PATH'),'cobramm'))
except KeyError:
    raise RuntimeError('cannot import cobramm module.\n'
                       'Please check if $COBRAM_PATH is defined in your environment, if not run:\n'
                       'export COBRAM_PATH=/path/to/COBRAMM')

# imports of local modules

import logwrt  # write messages and output to log (std out or cobramm.log)
import cobrammenv  # management of the enviroment for running COBRAMM (and utils)

# imports of user-defined classes

from amberCalculator import AmberCalculator, AmberSnapshot, AmberTopology  # Amber wrapper
from cobrammCalculator import CobrammCalculator  # Cobramm wrapper

#####################################################################################################

# SCRIPT PARAMETERS

_TARGET_DIR = "cobramm_input_files"  # name of the directories where final COBRAMM input files are saved

_DROPLET_DEFAULT_RAD = 15.0  # default radius of the solvent droplet
_MLAYER_DEFAULT_RAD = 7.0  # default radius of the solvent droplet


#####################################################################################################

def main():

    # this new class defines a style for the parser message that join the two predefined styles
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    # Define command line arguments with argparse
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="COBRAMM utility script: create a COBRAMM input for a molecule in a droplet of solvent",
        epilog="the script extract the results from an equilibrated AMBER snapshot \n"
               "and creates input files for a COBRAMM QM/MM calculation in which the\n"
               "molecule is in the QM layer and a variable solvation shell is in the \n"
               "mobile MM layer. \n\n",
        formatter_class=CustomFormatter)

    parser.add_argument("-p", "--topology", required=True, metavar='PRMTOP', dest="topologyfile",
                        help="name of the topology file (AMBER prmtop file)")
    parser.add_argument("-c", "--crdfile", required=True, metavar='INPCRD', dest="coordfile",
                        help="name of the coordinates file (AMBER inpcrd file)")

    parser.add_argument("-sr", "--solv-radius", default=_DROPLET_DEFAULT_RAD, metavar='RADIUS1', dest="solvradius",
                        type=float, help="radius of the final droplet of solvent, in Angstrom")
    parser.add_argument("-mr", "--mlayer-radius", default=_MLAYER_DEFAULT_RAD, metavar='RADIUS2', dest="mradius",
                        type=float, help="radius of the mobile solvent sphere, in Angstrom")

    parser.add_argument('--version', action='version',
                        version='%(prog)s is part of COBRAMM {}'.format(cobrammenv.getVersion()))

    args = parser.parse_args()

    # source COBRAMM configuration file (if available)
    confFile = cobrammenv.setCobrammProfile()
    # Print to log information on the environment configuration
    logwrt.writelog(cobrammenv.checkConfigFile(confFile))

    # initialize Amber and COBRAMM wrapper classes
    AmberCalculator(), CobrammCalculator()

    # ===================================================================

    # read coordinates and topology files
    logwrt.writelog("\nReading topology from AMBER prmtop file {}\n".format(args.topologyfile))
    with open(args.topologyfile) as ftopo:
        topology = AmberTopology(ftopo.read())
    logwrt.writelog("Reading coordinates from AMBER coordinate file {}\n\n".format(args.coordfile))
    with open(args.coordfile) as fcrd:
        snapshotInit = AmberSnapshot.readcrd(fcrd.read())

    # ================================================================================================================

    logwrt.writelog("Constructing COBRAMM input with a droplet of solvent of radius"
                    " {} Ang around the fist residue\n\n".format(args.solvradius))
    logwrt.writelog(" * high layer includes the first residue only\n")
    logwrt.writelog(" * medium layer includes the other residues within "
                    "{} Ang from the QM part\n".format(args.mradius))
    logwrt.writelog(" * lower layer includes the rest of the residues\n\n")

    # use AmberCalculator to create the Amber snapshot and topology objects that represent
    # a molecule surrounded by a solvent droplet of radius give by the argument spherRad
    dropletSnapshot, dropletTopology = AmberCalculator.cutdroplet(topology, snapshotInit, spherRad=args.solvradius)

    # now we can create input for COBRAMM with the qmmmLayers method of CobrammCalculator
    geometry, modelTopology, realTopology = CobrammCalculator.qmmmLayers(
        dropletSnapshot, dropletTopology, mobileThreshold=args.mradius, reorderRes=True)

    logwrt.writelog("Writing input files for COBRAMM in the directory {0}\n\n".format(_TARGET_DIR))

    # store starting dir and move to directory of the calculation
    startDir = os.getcwd()
    # then move to the work dir where final files will be saved, if the dir exists, remove first
    if os.path.isdir(_TARGET_DIR):
        logwrt.writewarning("overwriting previous " + _TARGET_DIR + " directory\n\n")
        shutil.rmtree(_TARGET_DIR)
    os.mkdir(_TARGET_DIR), os.chdir(_TARGET_DIR)

    # write input files
    # write model and full topologies to files
    with open("real.top", "w") as f:
        f.write(realTopology.topotext)
    with open("model-H.top", "w") as f:
        f.write(modelTopology.topotext)
    # prepare real_layers.xyz file
    with open("real_layers.xyz", "w") as f:
        f.write(geometry.reallayertext)
    # write the real.crd file
    geometry.makerealcrd()

    # write template for cobram.command
    with open("cobram.command", "w") as f:
        # keyword group
        f.write("!keyword\ntype=optxg\nnsteps=100\nqm-type=gauss\nqmem=500MB\ngeomem=500MB\n?keyword\n\n")
        # write MM part for AMBER
        f.write("!sander \ncomment line \n&cntrl \nimin   = 1, \nmaxcyc = 0, \nntb    = 0, \nigb    = 0,\n"
                "ntr    = 0, \nibelly = 1, \ncut    = 10 \n/\n?sander\n\n")
        # write QM part for GAUSSIAN
        f.write("!gaussian \n#p STO-3G HF nosym \n\n")
        f.write("gaussian input generated by COBRAMM \n\n0 1 \n?gaussian\n\n")

    # TODO: insert a line to add the RATTLE section to the dummy input file

    # move back to starting directory
    os.chdir(startDir)


#####################################################################################################

if __name__ == '__main__':
    main()
