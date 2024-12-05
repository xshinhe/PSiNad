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
import numpy as np

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
import constants
from cobrammCalculator import CobrammCalculator, CobrammInput, CobrammOutput
from transientCalculator import PumpProbeCalculator, MultiwfnCalculator

#####################################################################################################
def main():

    # this new class defines a style for the parser message that join the two predefined styles
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    def dir_path(string):
        if os.path.isdir(string):
            return string
        else:
            raise NotADirectoryError(f"directory {string} doesn't exist")


    # Define command line arguments with argparser
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="COBRAMM utility script: setup and run time-resolved transient absorption",
        epilog="The script runs vertical excitations on top of geometries along a TSH trajectories\n\n",
        formatter_class=CustomFormatter)

    parser.add_argument("-dir", "--directories", nargs="+", required=True, metavar='DIR', dest="directories_list",
                        type=dir_path, help="path of the directories containing the TSH trajectories")
    parser.add_argument("-m", "--method", required=True, metavar="METHOD", dest="method",
                        type=str, help="QM method of choice (TD-DFT, CASPT2)")
    parser.add_argument("-bs", "--bset", required=True, default="", metavar="BASISSET", dest="bset",
                        type=str, help="basis set of choice")
    parser.add_argument("-f", "--functional", default="", metavar="FUNCTIONAL", dest="functional",
                        type=str, help="functional of choice in case of TD-DFT calculation")
    parser.add_argument("-templ", "--template", default="command_molcas_template", metavar="TEMPLATE", dest="template",
                        type=str, help="name of molcas template file in case of CAS/RAS calculation")
    parser.add_argument("-q", "--charge", default=0, metavar="CHARGE", dest="charge",
                        type=int, help="net charge of the QM region")
    parser.add_argument("-ns", "--nstates", required=True, metavar='NSTATE', dest="nstates",
                        type=int, help="number of electronic states to consider")
    parser.add_argument("-st", "--simtime", required=True, metavar='SIMTIME', dest="simulation_time",
                        type=int, help="simulation time of the dynamics (fs)")
    parser.add_argument("-ts", "--timestep", required=True, metavar='TIMESTEP', dest="time_step",
                        help="time step performed (fs)")
    parser.add_argument("-Dt", "--deltat", required=True, metavar="DELTATIME", dest="delta_t",
                        type=int, help="delta of time to perform vertical excitations (fs)")
    parser.add_argument("--cluster", action="store_true",
                        help="submission of vertical excitation on cluster")
    parser.add_argument("-sub", "--submiss", metavar="SUBSTR", dest="submission_string",
                        type=str, help="submission string for cluster")
    parser.add_argument("-nc", "--ncores", default=1, metavar='NCORES', dest="ncores", type=int,
                        help="number of cpu's used")
    parser.add_argument('--version', action='version',
                        version='%(prog)s is part of COBRAMM {}'.format(cobrammenv.getVersion()))
    # TODO ADD PARSER OUTPUT

    args = parser.parse_args()

    #initialize class
    pp = PumpProbeCalculator(trjdir_list=args.directories_list, method=args.method, simulation_time=args.simulation_time,
                             delta_t=args.delta_t, nstates=args.nstates)

    upper_dir = os.getcwd()
    print("\nSetting up vertical excitations with a time step of {} fs\n".format(pp.delta_t))

    # setup and run single point calculations
    for folder in args.directories_list:

        print("\nEntering in the trajectory folder : {}\n".format(folder))

        pp.check_folder(folder)
        shutil.copy(args.template, folder)
        os.chdir(folder)

        if args.cluster:

            pp.setup_vertical_excitations(chrg=args.charge, functional=args.functional,
                                          template=args.template, ncores=args.ncores, basis_set=args.bset, cluster=True,
                                          submission_string=args.submission_string)
        else:
            pp.setup_vertical_excitations(chrg=args.charge, functional=args.functional,
                                          template=args.template, ncores=args.ncores, basis_set=args.bset)

        os.chdir(upper_dir)

    if args.cluster:
        print("\nVertical excitations submitted!\n")
    else:
        print("\nVertical excitations completed!\n")

#####################################################################################################

if __name__ == '__main__':
    main()
