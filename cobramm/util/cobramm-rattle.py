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

from cobrammCalculator import CobrammCalculator  # Cobramm wrapper
from layers import Layers  # object to store information on the geometry and the layers definition

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
        description="COBRAMM utility script: create the RATTLE section for the COBRAMM input",
        epilog="the script analyze the geometry and QM/MM setup of some pre-existing\n"
               "COBRAMM input files and produces the RATTLE section that is needed\n"
               "to constrain the X-H bonds in a molecular dynamics run.\n\n",
        formatter_class=CustomFormatter)

    parser.add_argument("-r", "--real-layers", default="real_layers.xyz", metavar='REALLAYERS', dest="reallayers",
                        help="name of the real_layers.xyz file (COBRAMM layers file)")
    parser.add_argument("-c", "--crdfile", default="real.crd", metavar='INPCRD', dest="coordfile",
                        help="name of the coordinates file (AMBER inpcrd file)")
    parser.add_argument('--version', action='version',
                        version='%(prog)s is part of COBRAMM {}'.format(cobrammenv.getVersion()))
    args = parser.parse_args()

    # source COBRAMM configuration file (if available)
    confFile = cobrammenv.setCobrammProfile()
    # Print to log information on the environment configuration
    logwrt.writelog(cobrammenv.checkConfigFile(confFile))

    # initialize COBRAMM wrapper classes
    CobrammCalculator()

    # ===================================================================

    # first define a Layers instance with the data contained in the real_layers.xyz file
    logwrt.writelog("\nReading layers from COBRAMM file '{}'\n".format(args.reallayers))
    with open(args.reallayers, "r") as f:
        geometry = Layers.from_real_layers_xyz(f.read())
    if os.path.isfile(args.coordfile):
        geometry.updatereal(filename=args.coordfile)
        logwrt.writelog("Reading coordinates from AMBER inpcrd file '{}'\n\n".format(args.coordfile))
    else:
        logwrt.writelog("Reading coordinates from COBRAMM file '{}'\n\n".format(args.reallayers))

    # now find the coordinates that need to be constrained and write the rattle section
    rattleSect = CobrammCalculator.rattle(geometry)

    if rattleSect is None:
        logwrt.writewarning("\nNo bond to constrain was found in the given geometry")
    else:
        logwrt.writelog("\nWriting RATTLE group for cobramm input in file '{0}'\n".format("rattle"))
        with open("rattle", "w") as f:
            f.write(rattleSect)


#####################################################################################################

if __name__ == '__main__':
    main()
