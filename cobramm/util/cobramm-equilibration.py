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
import optparse #command line option parser
import sys  # System-specific parameters and functions
import os  # filesystem utilities

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
from amberCalculator import AmberInput, AmberCalculator, AmberSnapshot, AmberTopology  # Amber wrapper

#####################################################################################################

# PARAMETERS

# directories where to run AMBER calculations
_CALC01_DIR = "01_opt"
_CALC02_DIR = "02_heat"
_CALC03_DIR = "03_equil"

# output files
_CALC01_PLOT = "optimization.pdf"
_CALC02_PLOT = "thermalization.pdf"
_CALC03_PLOT = "equilibration.pdf"


#####################################################################################################

def main():

    # this new class defines a style for the parser message that join the two predefined styles
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    # Define command line arguments with argparse
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="COBRAMM utility script: equilibrate initial conditions with AMBER molecular mechanics",
        epilog="the script launches a sequence of three AMBER calculations:\n"
               "1) an optimization run to minimize the energy with respect the system coordinates\n"
               "2) a thermalization run to reach a given temperature at constant volume\n"
               "3) an equilibration run to reach given temperature and pressures\n \n",
        formatter_class=CustomFormatter)

    parser.add_argument("-p", "--topology", required=True, metavar='PRMTOP', dest="topologyfile",
                        help="name of the topology file (AMBER prmtop file)")
    parser.add_argument("-c", "--crdfile", required=True, metavar='INPCRD', dest="coordfile",
                        help="name of the coordinates file (AMBER inpcrd file)")

    parser.add_argument("-T", "--temperature", default=300.0, metavar='TEMP', dest="temperature", type=float,
                        help="target temperature of the equilibration in K")
    parser.add_argument("-P", "--pressure", default=1.0, metavar='PRESS', dest="pressure", type=float,
                        help="target pressure of the equilibration in bar")

    parser.add_argument("-nc", "--ncores", default=1, metavar='NCORES', dest="ncores", type=int,
                        help="number of cpu's used, when > 1 run AMBER with MPI")

    parser.add_argument("-ctf", "--cutoff", default=12.0, metavar='CUTOFFVALUE', dest="cutoff", type=float,
                        help="cutoff value for MM potential evaluation in ang")

    parser.add_argument("-opt", "--optimiz-steps", default=1000, metavar='NSTEPS', dest="optsteps", type=int,
                        help="number of optimization steps, before heating and equilbration")
    parser.add_argument("-ht", "--heating-time", default=20.0, metavar='HEATINGTIME', dest="heatingtime", type=float,
                        help="time of the heating dynamics in ps")
    parser.add_argument("-et", "--equilibr-time", default=100.0, metavar='EQUILTIME', dest="equiltime", type=float,
                        help="time of the equilibration dynamics in ps")
    parser.add_argument("-dt", "--time-step", default=0.002, metavar='TSTEP', dest="timestep", type=float,
                        help="time step time of the heating+equilibration dynamics in ps")

    parser.add_argument('--version', action='version',
                        version='%(prog)s is part of COBRAMM {}'.format(cobrammenv.getVersion()))
    parser.add_argument('--frozenQM', action='store_true')
    args = parser.parse_args()

    # source COBRAMM configuration file (if available)
    confFile = cobrammenv.setCobrammProfile()
    # Print to log information on the environment configuration
    logwrt.writelog(cobrammenv.checkConfigFile(confFile))

    # initialize Amber wrapper class
    AmberCalculator()

    # ================================================================================================================

    # read coordinates and topology files
    logwrt.writelog("\nReading topology from AMBER prmtop file {}\n".format(args.topologyfile))
    with open(args.topologyfile) as ftopo:
        topology = AmberTopology(ftopo.read())
    logwrt.writelog("Reading coordinates from AMBER coordinate file {}\n".format(args.coordfile))
    with open(args.coordfile) as fcrd:
        snapshotInit = AmberSnapshot.readcrd(fcrd.read())

    if args.frozenQM:
        logwrt.writelog("Keeping H layer frozen\n")

    # define whether PBC are defined
    CoordhasPBC = snapshotInit.unitcell is not None

    # ================================================================================================================

    # create an instance of the Amber interface
    AmberCalculator()

    logwrt.writelog("\n\nGeometry optimization at the MM level, {} steps\n".format(args.optsteps))
    logwrt.writelog("files are stored in the directory {}\n".format(_CALC01_DIR))

    # prepare the input for a minimization run
    if args.frozenQM:
        input01 = AmberInput(minimize=True, nrSteps=args.optsteps, nprntsteps=1, cutoff=args.cutoff, usePBC=CoordhasPBC, freezesolute=True)
    else:
        input01 = AmberInput(minimize=True, nrSteps=args.optsteps, nprntsteps=1, cutoff=args.cutoff, usePBC=CoordhasPBC)

    # run the optimization, and do not remove directory at the end
    if args.frozenQM:
        output01 = AmberCalculator.run(input01, topology, snapshotInit, calcDir=_CALC01_DIR, nCores=args.ncores, store=True, ref=True)
    else:
        output01 = AmberCalculator.run(input01, topology, snapshotInit, calcDir=_CALC01_DIR, nCores=args.ncores, store=True)

    # analysis of the optimization
    output01.optimizationPlot(_CALC01_PLOT)
    logwrt.writelog("Plot of the energy optimization written to file {0}\n".format(_CALC01_PLOT))
    logwrt.writelog("\n")

    # extract optimized geometry
    optsnapshot = output01.snapshot(-1)

    logwrt.writelog("Heating the system at {0} K for {1} ps\n".format(args.temperature, args.heatingtime))
    logwrt.writelog("files are stored in the directory {0}\n".format(_CALC02_DIR))

    # prepare the input for the heating
    nrSteps = int(args.heatingtime / args.timestep)
    prtSteps = int(0.1 / args.timestep)  # print snapshot every 100 fs
    if args.frozenQM:
        input02 = AmberInput(minimize=False, nrSteps=nrSteps, deltaT=args.timestep, temperature=args.temperature,
                         pressure=None, cutoff=args.cutoff, freezeH=True, nprntsteps=prtSteps, usePBC=CoordhasPBC, freezesolute=True)
    else:
        input02 = AmberInput(minimize=False, nrSteps=nrSteps, deltaT=args.timestep, temperature=args.temperature,
                         pressure=None, cutoff=args.cutoff, freezeH=True, nprntsteps=prtSteps, usePBC=CoordhasPBC)

    # run the optimization, and do not remove directory at the end
    if args.frozenQM:
        output02 = AmberCalculator.run(input02, topology, optsnapshot, calcDir=_CALC02_DIR, nCores=args.ncores, store=True, ref=True)
    else:
        output02 = AmberCalculator.run(input02, topology, optsnapshot, calcDir=_CALC02_DIR, nCores=args.ncores, store=True)

    # analysis of the optimization
    output02.dynamicsPlot(_CALC02_PLOT)
    logwrt.writelog("Plot of the thermalization dynamics written to file {0}\n".format(_CALC02_PLOT))
    logwrt.writelog("\n")

    # extract thermalized geometry
    thermalsnapshot = output02.snapshot(-1)

    logwrt.writelog("Equilibrating the system at {0} K and {1} bar for {2} ps\n".format(
        args.temperature, args.pressure, args.equiltime))
    logwrt.writelog("files are stored in the directory {0}\n".format(_CALC03_DIR))

    # prepare the input for the heating
    nrSteps = int(args.equiltime / args.timestep)
    prtSteps = int(0.1 / args.timestep)  # print snapshot every 100 fs
    if args.frozenQM:
        input03 = AmberInput(minimize=False, nrSteps=nrSteps, deltaT=args.timestep, temperature=args.temperature,
                     pressure=args.pressure, cutoff=args.cutoff, freezeH=True, nprntsteps=prtSteps,
                     usevelocity=True, usePBC=CoordhasPBC, freezesolute=True)
    else:
        input03 = AmberInput(minimize=False, nrSteps=nrSteps, deltaT=args.timestep, temperature=args.temperature,
                     pressure=args.pressure, cutoff=args.cutoff, freezeH=True, nprntsteps=prtSteps,
                     usevelocity=True, usePBC=CoordhasPBC)

    # run the equilibration MD with AMBER
    if args.frozenQM:
        output03 = AmberCalculator.run(input03, topology, thermalsnapshot, calcDir=_CALC03_DIR,
                               nCores=args.ncores, store=True, ref=True)
    else:
        output03 = AmberCalculator.run(input03, topology, thermalsnapshot, calcDir=_CALC03_DIR,
                               nCores=args.ncores, store=True)

    # analysis of the optimization
    output03.dynamicsPlot(_CALC03_PLOT)
    logwrt.writelog("Plot of the equilibration dynamics written to file {0}\n".format(_CALC03_PLOT))
    logwrt.writelog("\n")

    # extract final equilibration geometry and write it to file
    finalsnapshot = output03.snapshot(-1)
    with open("finalsnapshot.crd", "w") as finalf:
        finalf.write(finalsnapshot.crdtext)


#####################################################################################################

if __name__ == '__main__':
    main()
