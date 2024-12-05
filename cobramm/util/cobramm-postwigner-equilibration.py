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
import shutil as sh
import numpy as np

# hijack $PYTHONPATH (sys.path) to find cobramm local modules
try:
    sys.path.insert(0,os.path.join(os.environ['COBRAM_PATH'], 'cobramm'))
except KeyError:
    raise RuntimeError('cannot import cobramm module.\n'
                       'Please check if $COBRAM_PATH is defined in your environment, if not run:\n'
                       'export COBRAM_PATH=/path/to/COBRAMM')

# imports of local modules

import logwrt  # write messages and output to log (std out or cobramm.log)
import cobrammenv  # management of the enviroment for running COBRAMM (and utils)

# imports of user-defined classes
from amberCalculator import AmberInput, AmberCalculator, AmberSnapshot, AmberTopology  # Amber wrapper
from cobrammCalculator import CobrammCalculator

#####################################################################################################

# PARAMETERS

# directories where to run AMBER calculations
_CALC_DIR = "MedLay_equil"

# output files
_CALC_PLOT = "MediumLayerEquilibration.pdf"

_TARGET_DIR = "cobramm_input_files"

#####################################################################################################

def main():

    # this new class defines a style for the parser message that join the two predefined styles
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    # Define command line arguments with argparse
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="COBRAMM utility script: equilibrate initial conditions with AMBER molecular mechanics",
        epilog="the script launches an AMBER calculations:\n"
               "an equilibration run at a given temperature to equilibrate Mobile Layers atoms to Wigner geometries of the H layer\n \n",
        formatter_class=CustomFormatter)

    parser.add_argument("-p", "--topology", required=True, metavar='PRMTOP', dest="topologyfile",
                        help="name of the topology file (AMBER prmtop file)")
    parser.add_argument("-c", "--crdfile", required=True, metavar='INPCRD', dest="coordfile",
                        help="name of the coordinates file (AMBER inpcrd file)")

    parser.add_argument("-mr", "--mlayer-radius", default=4, metavar='RADIUS2', dest="mradius",
                        type=float, help="radius of the mobile solvent sphere, in Angstrom")

    parser.add_argument("-T", "--temperature", default=300.0, metavar='TEMP', dest="temperature", type=float,
                        help="target temperature of the equilibration in K")

    parser.add_argument("-nc", "--ncores", default=1, metavar='NCORES', dest="ncores", type=int,
                        help="number of cpu's used, when > 1 run AMBER with MPI")

    parser.add_argument("-ctf", "--cutoff", default=12.0, metavar='CUTOFFVALUE', dest="cutoff", type=float,
                        help="cutoff value for MM potential evaluation in ang")

    parser.add_argument("-et", "--equilibr-time", default=100.0, metavar='EQUILTIME', dest="equiltime", type=float,
                        help="time of the equilibration dynamics in ps")

    parser.add_argument("-dt", "--time-step", default=0.002, metavar='TSTEP', dest="timestep", type=float,
                        help="time step time of the heating+equilibration dynamics in ps")

    parser.add_argument("-Hvel", "--HighLayer-velocty", required=True, default="Hvel.dat", metavar='HVEL', dest="hvel", type=str,
                        help="name of the file contaning Wigner velocities of H layer")
    
    parser.add_argument('--version', action='version',
                        version='%(prog)s is part of COBRAMM {}'.format(cobrammenv.getVersion()))
    args = parser.parse_args()

    # source COBRAMM configuration file (if available)
    confFile = cobrammenv.setCobrammProfile()
    # Print to log information on the environment configuration
    logwrt.writelog(cobrammenv.checkConfigFile(confFile))

    # initialize Amber wrapper class
    AmberCalculator()

    # ================================================================================================================

    with open(args.hvel, "r") as f:
        hvel=f.readlines()

    # read coordinates and topology files
    logwrt.writelog("\nReading topology from AMBER prmtop file {}\n".format(args.topologyfile))
    with open(args.topologyfile) as ftopo:
        topology = AmberTopology(ftopo.read())
    logwrt.writelog("Reading coordinates from AMBER coordinate file {}\n".format(args.coordfile))
    with open(args.coordfile) as fcrd:
        snapshotInit = AmberSnapshot.readcrd(fcrd.read())

    # define whether PBC are defined
    CoordhasPBC = snapshotInit.unitcell is not None

    # ================================================================================================================

    # create an instance of the Amber interface
    AmberCalculator()

    # prepare the input for the equilibration
    nrSteps = int(args.equiltime / args.timestep)
    prtSteps = int(0.1 / args.timestep)  # print snapshot every 100 fs
        
    input_eq = AmberInput(minimize=False, nrSteps=nrSteps, deltaT=args.timestep, temperature=args.temperature,
                     cutoff=args.cutoff, freezeH=True, freezesolute=True, nprntsteps=prtSteps,
                     usevelocity=False, usePBC=CoordhasPBC, freezeLowLayer=True, mradius=args.mradius)

    # run the equilibration MD with AMBER
    output_eq = AmberCalculator.run(input_eq, topology, snapshotInit, calcDir=_CALC_DIR,
                               nCores=args.ncores, store=True, ref=True)

    # analysis of the optimization
    output_eq.dynamicsPlot(_CALC_PLOT)
    logwrt.writelog("Plot of the equilibration dynamics written to file {0}\n".format(_CALC_PLOT))
    logwrt.writelog("\n")

    # extract final equilibration geometry and write it to file
    finalsnapshot = output_eq.snapshot(-1)
    with open("finalsnapshot.crd", "w") as finalf:
        finalf.write(finalsnapshot.crdtext)

    startDir = os.getcwd()
    # then move to the work dir where final files will be saved, if the dir exists, remove first
    if os.path.isdir(_TARGET_DIR):
        logwrt.writewarning("overwriting previous " + _TARGET_DIR + " directory\n\n")
        sh.rmtree(_TARGET_DIR)
    os.mkdir(_TARGET_DIR), os.chdir(_TARGET_DIR)

    geometry, modelTopology, realTopology = CobrammCalculator.qmmmLayers(
        finalsnapshot, topology, mobileThreshold=args.mradius, reorderRes=True)

    # write input files
    # write model and full topologies to files
    with open("real.top", "w") as f:
        f.write(realTopology.topotext)
    with open("model-H.top", "w") as f:
        f.write(modelTopology.topotext)
    # prepare real_layers.xyz file
    with open("real_layers.xyz", "w") as f:
        f.write(geometry.reallayertext)

    # read number of M atoms
    layertext = open("real_layers.xyz", "r")
    text = layertext.read()
    nMatoms = text.count("M")

    # write the real.crd file
    geometry.makerealcrd()

    # write the velocity.dat file combining H velocities from Wigner and M velocities from the equilibration
    firstMatom = len(hvel)
    lastMatom = nMatoms+firstMatom
    velocities = np.array(finalsnapshot.velocity).transpose().reshape((-1))
    velM = ""
    for i, c in enumerate(velocities):
        if i < firstMatom*3:
            pass
        elif i >= firstMatom*3 and i < lastMatom*3:
            velM += "{0:12.7f}".format(c)
            if (i + 1) % 3 == 0: velM += "\n"
        elif i >= nMatoms*3:
            break
    if len(velocities) % 3 != 0: velM += "\n"

    with open("velocity.dat", "w") as vel:
        for Hatom in range(len(hvel)):
            vel.write(hvel[Hatom])
        vel.write(velM)

    # write template for cobram.command for TSH
    with open("cobram.command", "w") as f:
        # keyword group
        f.write("!keyword\ntype=mdv\nnproc=1\nnumproc=1\nnsteps=100\nqm-type=gauss\nqmem=2000MB\ntstep=0.5\n"
                "surhop=persico\nnacs=tdc\ntdctype=0\nhoptogs=0\nbackhop=0\nvelafterhop=1\n?keyword\n\n")
        # write MM part for AMBER
        f.write("!sander \ncomment line \n&cntrl \nimin   = 1, \nmaxcyc = 0, \nntb    = 0, \nigb    = 0,\n"
                "ntr    = 0, \nibelly = 1, \ncut    = 10 \n/\n?sander\n\n")
        # write QM part for GAUSSIAN
        f.write("!gaussian \n#p def2svp cam-b3lyp nosym tda=(nstates=5) \n\n")
        f.write("gaussian input generated by COBRAMM \n\n0 1 \n?gaussian\n\n")

    # TODO: insert a line to add the RATTLE section to the dummy input file

    # move back to starting directory
    os.chdir(startDir)

#####################################################################################################

if __name__ == '__main__':
    main()
