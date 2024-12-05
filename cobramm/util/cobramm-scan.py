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
import subprocess  # run external program as child process
import re          # process output with regular expressions

# other modules

try:
    from openbabel import openbabel
except:
    import openbabel  # openbabel API
import numpy as np  # numpy numerical library

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
import constants  # physical constants and conversion factors

# imports of user-defined classes

from cobrammCalculator import CobrammCalculator  # Cobramm wrapper classes
from layers import Layers  # object to store information on the geometry and the layers definition
from amberCalculator import AmberTopology  # AMBER wrapper classes

#####################################################################################################

# SCRIPT PARAMETERS

_REALLAYERS_FNAME = "real_layers.xyz"
_REALCRD_FNAME = "real.crd"
_REALTOP_FNAME = "real.top"
_MODELHTOP_FNAME = "model-H.top"
_COMMAND_FNAME = "cobram.command"

_SCANVALUE_FILE = "SCAN_VALUE"


######################################################################################################################

def main():

    # this new class defines a style for the parser message that join the two predefined styles
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    # Define command line arguments with argparse
    # create the top-level parser
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="COBRAMM utility script: coordinate scan with COBRAMM",
        epilog="this script realize a sequence of COBRAMM calculation starting from a single set of \n"
               "input files, by changing the value of a chosen coordinate - atom distance, angle or \n"
               "torsion - in a uniform grid. The full workflow is divided in a sequence of three steps \n"
               "that include the preparation of the input files (action 'prepare'), the execution of \n"
               "COBRAMM calculations (action 'run') and the plot of the results (action 'plot')\n\n",
        formatter_class=CustomFormatter)

    parser.add_argument('--version', action='version',
                        version='%(prog)s is part of COBRAMM {}'.format(cobrammenv.getVersion()))

    # define subparsers depending for each of the script actions
    subparsers = parser.add_subparsers(title='valid actions',
                                       description='choose one of the steps of the scanning procedure')

    # create the parser for the "prepare" command
    parser_prepare = subparsers.add_parser('prepare', help='prepare COBRAMM input files for the scan snapshots',
                                           formatter_class=CustomFormatter)
    parser_prepare.add_argument("-c", "--coordinate", required=True, metavar="n1,n2,[n3,n4]", dest="coordinate",
                                help="coordinate to scan (2 integers = bond, 3 int = angle, 4 int = torsion")
    parser_prepare.add_argument("-s", "--scan-limit", required=True, metavar="vi,vf,dv", dest="scanlimit",
                                help="scan values: starting value, (excluded) ending value, step")
    parser_prepare.add_argument("-d", "--input-dir", default="template", metavar='PATH', dest="inputdir",
                                help="path where the template input files for COBRAMM are stored")
    parser_prepare.add_argument("--force", action='store_true', dest="overwrite",
                                help="overwrite previous scan directories")
    parser_prepare.set_defaults(func=prepare)

    # create the parser for the "run" command
    parser_run = subparsers.add_parser('run', help='launch COBRAMM calculations', formatter_class=CustomFormatter)
    parser_run.add_argument('-cmd', "--command", default="cobram.py > cobramm.log", metavar="CMD", dest="command",
                            help='shell command for running COBRAMM')
    parser_run.add_argument('-bg', "--background", action='store_true', dest="runinbackground",
                            help='launch COBRAMM commands in background (use with caution!)')
    parser_run.set_defaults(func=run)

    # create the parser for the "plot" command
    parser_plot = subparsers.add_parser('plot', help='process results, plot QMMM energy vs the scan coordinate',
                                        formatter_class=CustomFormatter)
    parser_plot.add_argument('-p', '--plot', help='automatic plot to screen using matplotlib', action='store_true')
    parser_plot.add_argument('-pdf', '--save-pdf', metavar='FILENAME', dest="pdffile",
                             help='save matplotlib plot to pdf file')
    parser_plot.add_argument('-png', '--save-png', metavar='FILENAME', dest="pngfile",
                             help='save matplotlib plot to png file')
    parser_plot.set_defaults(func=plot)

    args = parser.parse_args()

    # call the specific action requested by command argument
    try:
        args.func(args)
    except AttributeError:
        parser.print_help()


######################################################################################################################

def prepare(args):

    # source COBRAMM configuration file (if available)
    confFile = cobrammenv.setCobrammProfile()
    # Print to log information on the environment configuration
    logwrt.writelog(cobrammenv.checkConfigFile(confFile))

    # initialize COBRAMM wrapper classes
    CobrammCalculator()

    # PROCESS COMMAND LINE ARGUMENTS

    coordinate = [int(s) for s in args.coordinate.split(",")]
    scanlimit = [float(s) for s in args.scanlimit.split(",")]
    scanvalue = np.arange(scanlimit[0], scanlimit[1], scanlimit[2])

    # ================================================================================================================

    # CHECK THAT ALL THE NECESSARY INPUT FILES ARE PRESENT IN THE inputdir DIRECTORY

    reallayers = os.path.join(args.inputdir, _REALLAYERS_FNAME)
    coordfile = os.path.join(args.inputdir, _REALCRD_FNAME)
    modelHtop = os.path.join(args.inputdir, _MODELHTOP_FNAME)
    realtop = os.path.join(args.inputdir, _REALTOP_FNAME)
    commandfile = os.path.join(args.inputdir, _COMMAND_FNAME)

    if os.path.isfile(commandfile):
        logwrt.writelog("\nUsing COBRAMM command from input file '{}'\n".format(commandfile))
    else:
        logwrt.fatalerror("cannot find required file {0}".format(commandfile))

    # first define a Layers instance with the data contained in the real_layers.xyz file
    logwrt.writelog("Reading layers from COBRAMM file '{}'\n".format(reallayers))
    with open(reallayers, "r") as f:
        geometry = Layers.from_real_layers_xyz(f.read())

    # depending on the availability of the real.crd file, update coordinates
    if os.path.isfile(coordfile):
        geometry.updatereal(filename=coordfile)
        logwrt.writelog("Reading coordinates from AMBER inpcrd file '{}'\n".format(coordfile))
    else:
        logwrt.writelog("Reading coordinates from COBRAMM file '{}'\n".format(reallayers))

    # when the calculation is not only QM, a real topology is needed
    if geometry.calculationType != "H":
        if os.path.isfile(realtop):
            logwrt.writelog("Using topology for the full system from AMBER partop file '{}'\n".format(realtop))
            with open(realtop, "r") as f:
                _ = AmberTopology(f.read())
        else:
            logwrt.fatalerror("cannot find required file {0}".format(realtop))

        # when the calculation is QM/MM, also a modelH topology is needed
        if "H" in geometry.calculationType:
            if os.path.isfile(modelHtop):
                logwrt.writelog("Using topology for the QM system from AMBER partop file '{}'\n".format(modelHtop))
                with open(modelHtop, "r") as f:
                    _ = AmberTopology(f.read())
            else:
                logwrt.fatalerror("cannot find required file {0}".format(modelHtop))

    # ================================================================================================================

    # PREPARE MOLECULE IN OB CLASS

    xyzinp = str(geometry.NatomHM) + "\n\n" + geometry.to_string(geom_model='MEDIUM_HIGH', labels=True)

    mol = openbabel.OBMol()  # OpenBabel class to store molecular information
    obConversion = openbabel.OBConversion()  # OpenBabel class for input/output
    obConversion.SetInAndOutFormats("xyz", "xyz")

    check = obConversion.ReadString(mol, xyzinp)
    if not check:
        logwrt.fatalerror("Error reading the molecular structure with OpenBabel.")
    mol.ConnectTheDots()

    # ================================================================================================================

    # PRINT INFO ON SCANNING COORDINATE (+ preliminary actions)

    at = [mol.GetAtom(ind) for ind in coordinate]

    if len(at) == 2:
        logwrt.writelog("\nScanning along the distance between atoms {0}-{1} and {2}-{3}\n".format(
            at[0].GetType(), at[0].GetIdx(), at[1].GetType(), at[1].GetIdx()))
        check = mol.AddBond(coordinate[0], coordinate[1], 1)
        if check:
            logwrt.writelog("Bond {0}-{1} added!\n".format(coordinate[0], coordinate[1]))
        else:
            logwrt.writelog("Bond {0}-{1} was already present!\n".format(coordinate[0], coordinate[1]))

    elif len(at) == 3:
        logwrt.writelog("\nScanning along the angle between atoms {0}-{1}, {2}-{3} and {4}-{5}\n".format(
            at[0].GetType(), at[0].GetIdx(), at[1].GetType(), at[1].GetIdx(), at[2].GetType(), at[2].GetIdx()))

    elif len(at) == 4:
        logwrt.writelog("\nScanning along the torsion between atoms {0}-{1}, {2}-{3}, {4}-{5} and {6}-{7}\n".format(
            at[0].GetType(), at[0].GetIdx(), at[1].GetType(), at[1].GetIdx(),
            at[2].GetType(), at[2].GetIdx(), at[3].GetType(), at[3].GetIdx()))

    else:
        logwrt.fatalerror("invalid coordinate: please give 2 (distance), 3 (angle) or 4 (torsions) integer numbers")

    # make a string to accumulate the snapshots xyz, that will be written in a unique file at the end
    snapmovie = ""

    # ================================================================================================================

    # LOOP OVER THE SCAN VALUES

    for isnap, value in enumerate(scanvalue):

        # define the new geometry with OB and store it in the string newxyz in xyz format
        newxyz = ""

        if len(at) == 2:

            logwrt.writelog("Setting the bond length to {0:.3f} Ang in snapshot {1}\n".format(value, isnap))
            bond = mol.GetBond(*at)  # identify the bond
            bond.SetLength(value)  # set the length of the bond to the given value
            newxyz = obConversion.WriteString(mol)

        elif len(at) == 3:

            # in the case of the angle, the procedure is a bit complicated, and does not always work:
            # we get the Z-matrix from the cartesian coordinates, and in case the angle that we want
            # to change is defined within the Z-matrix, we modified, otherwise there is not much
            # to do and we exit with error

            logwrt.writelog("Setting the angle to {0:.1f} degrees in snapshot {1}\n".format(value, isnap))

            coords = mol.GetInternalCoord()  # transform to Z matrix

            # loop over the lines of the Z matrix, searching for the angle
            angleFound = False
            for i, item in enumerate(coords):
                if item is not None:
                    try:
                        setIdx = {i, item._a.GetIdx(), item._b.GetIdx()}
                    except AttributeError:
                        setIdx = set()
                    if set(coordinate) == setIdx:  # use set comparison to disregard the order of the elements
                        item._ang = value  # when the angle is found, change its value and break the loop
                        angleFound = True
                        break

            # if the angle was not found, exit with error
            if not angleFound:
                logwrt.fatalerror("the angle is not defined in OB Z-matrix... try to reordering the xyz")

            # now transform back to cartesian coordinate, and write the structure in xyz format
            openbabel.InternalToCartesian(openbabel.vectorpOBInternalCoord(coords[1:]), mol)
            newxyz = obConversion.WriteString(mol)

        elif len(at) == 4:

            logwrt.writelog("Setting the torsion angle to {0:.1f} degrees in snapshot {1}\n".format(value, isnap))
            mol.SetTorsion(at[0], at[1], at[2], at[3], value*constants.Deg2Rad)
            newxyz = obConversion.WriteString(mol)

        # accumulate the xyz snapshot to print the whole movie at the end
        snapmovie += newxyz

        # read coordinates from xyz format and update geometry
        coords = [[], [], []]
        for line in newxyz.splitlines()[2:]:
            coords[0].append(float(line.split()[1]))
            coords[1].append(float(line.split()[2]))
            coords[2].append(float(line.split()[3]))
        geometry.updateHMlayers(coords)

        # this is the name of the directory where the files will be stored
        dirname = "scan_{0:03}".format(isnap)
        # if a directory with the same name is present,
        # stop unless the overwrite option is given... in such case remove the directory and go on
        if os.path.isdir(dirname):
            if args.overwrite:
                shutil.rmtree(dirname)
            else:
                logwrt.writelog("\nA directory of a previous scan is present... stopping the script! \n"
                                "If you want to overwrite it, please consider using the --force option\n")
                sys.exit()

        # create the directory with the input file for cobramm
        shutil.copytree(args.inputdir, dirname)

        # create a SCAN_VALUE file
        with open(os.path.join(dirname, _SCANVALUE_FILE), "w") as f:
            f.write(str(value))

        # write a new real.crd there to update the coordinates of the snapshot
        geometry.makerealcrd(filename=os.path.join(dirname, _REALCRD_FNAME))

    with open("scan_movie.xyz", "w") as f:
        f.write(snapmovie)


######################################################################################################################

def run(args):

    # source COBRAMM configuration file (if available)
    confFile = cobrammenv.setCobrammProfile()
    # Print to log information on the environment configuration
    logwrt.writelog(cobrammenv.checkConfigFile(confFile))

    # initialize COBRAMM wrapper classes
    CobrammCalculator()

    logwrt.writelog("\nRunning COBRAMM calculations from sub-directories scan_* \n")

    # define the list of scan_* directories
    # list all the subdirectories
    listDir = [name for name in os.listdir(".") if os.path.isdir(os.path.join(".", name))]
    # define those directories that match the "scan_*" pattern
    matchingDir = [name for name in listDir if re.match(r'.*scan_[0-9]+.*', name)]
    matchingDir = sorted(matchingDir)

    for isnap, dirname in enumerate(matchingDir):

        # move to the snapshot dir saving the name of the starting directory
        currentDir = os.getcwd()
        os.chdir(dirname)

        # run cobramm
        if args.runinbackground:  # launch process and go on, do not wait for the single calculation to finish
            logwrt.writelog("Launching COBRAMM calculation from directory {0} in background\n".format(dirname))
            subprocess.Popen(args.command, shell=True, stdout=subprocess.PIPE)
        else:  # run one command after another
            logwrt.writelog("Running COBRAMM calculation from directory {0}\n".format(dirname))
            subprocess.run(args.command, shell=True, stdout=subprocess.PIPE)

        # move back to the starting directory
        os.chdir(currentDir)


######################################################################################################################

def plot(args):

    # source COBRAMM configuration file (if available)
    confFile = cobrammenv.setCobrammProfile()
    # Print to log information on the environment configuration
    logwrt.writelog(cobrammenv.checkConfigFile(confFile))

    # initialize COBRAMM wrapper classes
    CobrammCalculator()

    # define the list of scan_* directories
    # list all the subdirectories
    listDir = [name for name in os.listdir(".") if os.path.isdir(os.path.join(".", name))]
    # define those directories that match the "scan_*" pattern
    matchingDir = [name for name in listDir if re.match(r'.*scan_[0-9]+.*', name)]
    matchingDir = sorted(matchingDir)

    # name of the output dat file
    energydatfile = "scan_qmmm_energy.dat"
    logwrt.writelog("\nWriting QM/MM energies vs the scan coordinate to file {}\n\n".format(energydatfile))

    # remove file from previous run
    if os.path.isfile(energydatfile):
        os.remove(energydatfile)

    plt_data = []

    for isnap, dirname in enumerate(matchingDir):

        # create a SCAN_VALUE file
        with open(os.path.join(dirname, _SCANVALUE_FILE), "r") as f:
            value = float(f.read())

        # read cobramm.log file
        with open(os.path.join(dirname, "cobramm.log"), "r") as f:
            cobrammlog = f.read()

        # extract energy values with regular expressions
        regex_QMMMEnergy = 'STATE *[0-9]*.*?E\(tot\)\=E\(Model\-H QM\)\+\(Real MM\)\-' \
                           '\(Model\-H MM\)\-\(Emb\-emb\)\= *(-?[0-9]*.[0-9]*)'
        E_QMMM = re.findall(regex_QMMMEnergy, cobrammlog, flags=re.DOTALL)

        # accumulate data to plot
        plt_data.append([value] + [float(e) for e in E_QMMM])

        # update the scan.dat file with the QM/MM energies depending on the scanning coordinate
        with open(energydatfile, "a") as f:
            f.write("{0:10.3f} ".format(value) + " ".join(["{0}".format(e) for e in E_QMMM]) + "\n")

    if args.plot or args.pdffile or args.pngfile:
        # import matplotlib and create a graph
        import matplotlib.pyplot as plt
        fig, graph = plt.subplots(1, 1)

        # transpose plt_data
        plt_data = list(zip(*plt_data))
        # plot data
        legend_list = []
        for i in range(1, len(plt_data)):
            g = graph.plot(plt_data[0], plt_data[i], label='el. state {0}'.format(i))
            legend_list.append(g[0])

        # legend, axis labels, title
        graph.legend(legend_list, [p.get_label() for p in legend_list])
        graph.set_xlabel('scan coordinate')
        graph.set_ylabel('QM/MM energy / Hartree')

        # plot to screen
        if args.plot:
            plt.show()
        if args.pngfile:
            plt.savefig(args.pngfile, dpi=None, facecolor='w', edgecolor='w', orientation='landscape', format="png")
        if args.pdffile:
            plt.savefig(args.pdffile, dpi=None, facecolor='w', edgecolor='w', orientation='landscape', format="pdf")


######################################################################################################################

if __name__ == '__main__':
    main()
