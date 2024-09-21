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

# import statements of module from python standard library

import os  # filesystem utilities
import shutil  # filesystem utilities
import subprocess  # run external program as child process
import fileinput  # fileinput is used for inplace editing of files
import shlex  # simple lexical analysis for unix syntax
import math  # import mathematical functions
import re  # process output with regular expressions
import copy  # shallow and deep copy operations
import time  # provides various time-related functions
import urllib.request  # Extensible library for opening URLs

# imports of local modules

import logwrt  # write messages and output to log (std out or cobramm.log)
import constants  # physical and mathematical constants
import qmmmenv  # function to check if AMBER executables are available

# math libraries

import numpy as np  # numpy: arrays and math utilities
from scipy.io import netcdf  # read NetCDF file using scipy


#####################################################################################################

class AmberError(Exception):
    """A custom exception used to report errors in use of Amber"""


#####################################################################################################

class AmberDriver:
    """ wrapper for the execution of AMBER """

    # class variables defined by the constructor
    amberhome = ""  # value of the env variables that points to the AMBER root directory
    amberVersion = None  # version of AMBER installed in $AMBERHOME
    mpirun = ""
    # class variable that controls when the check on the environment has been done
    _amberCheck = False

    # constant class variables with the names of the Amber I/O files
    INPNAME = "inp"
    TOPNAME = "top"
    CRDNAME = "crd"
    OUTNAME = "out"
    MDCRDNAME = "mdcrd"
    RESTRTNAME = "restrt"

    # Data to be used for force-field definitions in each version of Amber
    _TIP3PBOX_18 = ["", "TIP3PBOX", {}, "Standard AMBER force field"]
    _OPC3_WAT_18 = ["loadAmberParams frcmod.opc3", "OPC3BOX", {}, "standard AMBER force field"]
    _SPCF_WAT_18 = ["WAT = SPF\nloadAmberParams frcmod.spcfw\nset default FlexibleWater on", "SPCFWBOX", {}, "standard AMBER force field"]
    _SPCE_WAT_18 = ["loadAmberParams frcmod.spce", "SPCBOX", {}, "standard AMBER force field"]
    _METHANOL_18 = ["loadAmberParams frcmod.meoh", "MEOHBOX", {}, "standard AMBER force field"]
    _CHLOROFO_18 = ["loadAmberParams frcmod.chcl3", "CHCL3BOX", {}, "standard AMBER force field"]
    _DMSO_ONLINE = ["loadOff dmsobox.off\nloadAmberParams frcmod.dmso", "d",
                    {"dmsobox.off": "http://personalpages.manchester.ac.uk/staff/Richard.Bryce/amber/box/dmsobox.off",
                     "frcmod.dmso": "http://personalpages.manchester.ac.uk/staff/Richard.Bryce/amber/box/frcmod.dmso"},
                    "Fox T, Kollman PA, J. Phys. Chem. B 1998, 102, 8070"]
    _ACNT_ONLINE = ["loadAmberprep acn.prep\nloadAmberParams frcmod.acn\nACNBOX = loadpdb acn.pdb", "ACNBOX",
                    {"acn.pdb": "http://amber.manchester.ac.uk/box/ch3cn_210.pdb",
                     "acn.prep": "http://amber.manchester.ac.uk/box/prep.ch3cn",
                     "frcmod.acn": "http://amber.manchester.ac.uk/box/frcmod.ch3cn"},
                    "Grabuleda, X, Jaime, C and Kollman, PA, 2000, J. Comput. Chem., 21: 901-908"]

    _AMBERSOLVENTS = {18: {"water": _TIP3PBOX_18, "methanol": _METHANOL_18, "chloroform": _CHLOROFO_18,
                           "DMSO": _DMSO_ONLINE, "acetonitrile": _ACNT_ONLINE, "opc3": _OPC3_WAT_18, "spcfw":_SPCF_WAT_18, "spce":_SPCE_WAT_18},
                      19: {"water": _TIP3PBOX_18, "methanol": _METHANOL_18, "chloroform": _CHLOROFO_18,
                           "DMSO": _DMSO_ONLINE, "acetonitrile": _ACNT_ONLINE, "opc3": _OPC3_WAT_18, "spcfw":_SPCF_WAT_18, "spce":_SPCE_WAT_18},
                      20: {"water": _TIP3PBOX_18, "methanol": _METHANOL_18, "chloroform": _CHLOROFO_18,
                           "DMSO": _DMSO_ONLINE, "acetonitrile": _ACNT_ONLINE, "opc3": _OPC3_WAT_18, "spcfw":_SPCF_WAT_18, "spce":_SPCE_WAT_18},
                      21: {"water": _TIP3PBOX_18, "methanol": _METHANOL_18, "chloroform": _CHLOROFO_18,
                           "DMSO": _DMSO_ONLINE, "acetonitrile": _ACNT_ONLINE, "opc3": _OPC3_WAT_18, "spcfw":_SPCF_WAT_18, "spce":_SPCE_WAT_18},
                      22: {"water": _TIP3PBOX_18, "methanol": _METHANOL_18, "chloroform": _CHLOROFO_18,
                           "DMSO": _DMSO_ONLINE, "acetonitrile": _ACNT_ONLINE, "opc3": _OPC3_WAT_18, "spcfw":_SPCF_WAT_18, "spce":_SPCE_WAT_18},
                      }

    _LEAP_PREAMBLE = {18: "source leaprc.protein.ff14SB\nsource leaprc.DNA.OL15\nsource leaprc.lipid17"
                          "\nsource leaprc.water.tip3p\nsource leaprc.gaff2",
                      19: "source leaprc.protein.ff14SB\nsource leaprc.DNA.OL15\nsource leaprc.lipid17"
                          "\nsource leaprc.water.tip3p\nsource leaprc.gaff2",
                      20: "source leaprc.protein.ff14SB\nsource leaprc.DNA.OL15\nsource leaprc.lipid17"
                          "\nsource leaprc.water.tip3p\nsource leaprc.gaff2",
                      21: "source leaprc.protein.ff14SB\nsource leaprc.DNA.OL15\nsource leaprc.lipid17"
                          "\nsource leaprc.water.tip3p\nsource leaprc.gaff2",
                      22: "source leaprc.protein.ff14SB\nsource leaprc.DNA.OL15\nsource leaprc.lipid17"
                          "\nsource leaprc.water.tip3p\nsource leaprc.gaff2",
                      }

    # ================================================================================

    def __init__(self):
        """ the constructor of the AmberDriver class checks if the environment for AMBER execution
        is properly defined and if all the necessary executables are available """

        if not AmberDriver._amberCheck:

            # environmental variable with the AMBER path
            try:
                AmberDriver.amberhome = os.environ['AMBERHOME']
            except KeyError:
                raise AmberError("environment variable $AMBERHOME is not defined")

            # environmental variable with the MPI command
            try:
                AmberDriver.mpirun = os.environ['PARA_EXE']
            except KeyError:
                AmberDriver.mpirun = 'mpirun -np'

            # call update_amber to check the current version of amber
            updateExe = os.path.join(AmberDriver.amberhome, "update_amber")
            if os.path.isfile(updateExe) and os.access(updateExe, os.X_OK):
                command = shlex.split("python3 " + updateExe + " -v")
                proc = subprocess.run(command, stdout=subprocess.PIPE)
                # process the std out line-by-line
                for line in proc.stdout.splitlines():
                    if b"AmberTools version" in line:
                        AmberDriver.amberVersion = line.decode().split()[2].split(".")
            else:
                AmberDriver.amberVersion = ["20"]

            # define the list of executables depending on the version of amber
            if int(AmberDriver.amberVersion[0]) <= 12:
                exeList = ["antechamber", "tleap", "sander", "ambpdb", "parmchk", "cpptraj", "ambmask"]
            else:
                exeList = ["antechamber", "tleap", "sander", "ambpdb", "parmchk2", "cpptraj", "ambmask"]

            # check if amber executables are available
            for exe in exeList:
                if not qmmmenv.which(exe):
                    raise AmberError(exe + " executable is not available")

            # check that openbabel is available
            if not qmmmenv.which("obabel"):
                raise AmberError("obabel executable is not available")

            # write message about the AMBER version currently in use
            logwrt.writelog("Running AMBER {0} from {1}\n".format(
                AmberDriver.amberVersion[0], AmberDriver.amberhome))

            # set class variable to true, to skip this check at next initialization
            AmberDriver._amberCheck = True

    # ================================================================================

    @staticmethod
    def run(inputfile, topologyfile, snapshot, calcDir="amberCalc", overwrite=False, nCores=1, GPU=False, store=False, ref=False):
        """ run AMBER with the initial conditions defined in inputData and store the results in outputData.

        :rtype: AmberOutput
        :param inputfile: AmberInput class instance with the input file for AMBER
        :param topologyfile: AmberTopology class instance with the topology file
        :param snapshot: AmberSnapshot class instance with the coordinates and velocities of the system
        :param calcDir: string that specify where the calculation is run
        :param overwrite: if True overwrite previous calculations in calcDir, otherwise attempt to read existing files
        :param nCores: number of cores used in a parallel run
        :param GPU: use the GPU compiled version of the AMBER code
        :param store: do not remove output files when execution ends
        :return: AmberOutput class instance with the results of the calculation
        """

        # check if inputfile is of AmberInput type
        if not isinstance(inputfile, AmberInput):
            raise AmberError("amberDriver.run: first argument should be of amberInput type")
        # check if topologyfile is of AmberTopology type
        if not isinstance(topologyfile, AmberTopology):
            raise AmberError("amberDriver.run: second argument should be of AmberTopology type")
        # check if snapshot is of AmberSnapshot type
        if not isinstance(snapshot, AmberSnapshot):
            raise AmberError("amberDriver.run: third argument should be of AmberSnapshot type")

        # store starting dir
        startDir = os.getcwd()

        # initialize a logical flag readExisting (when true, read results from the dir and do not re-run existing calc)
        readExisting = False

        # check if the existing directory contains an identical calculation
        if os.path.isdir(calcDir) and not overwrite:
            os.chdir(calcDir)

            # check whether an identical calculation setup has been used, by comparing input files
            identicalSetup = False
            try:
                with open(AmberDriver.INPNAME, "r") as oldinp:
                    if oldinp.read().strip() == inputfile.inputtext.strip(): identicalSetup = True
            except IOError:
                pass

            # check whether an identical topology has been used, by comparing topology files
            identicalTopology = False
            try:
                with open(AmberDriver.TOPNAME, "r") as oldtop:
                    if oldtop.read().strip() == topologyfile.topotext.strip(): identicalTopology = True
            except IOError:
                pass

            # check whether identical initial conditions have been used
            identicalInitCond = False
            try:
                with open(AmberDriver.CRDNAME, "r") as f:
                    oldcrd = AmberSnapshot.readcrd(f.read())
                if oldcrd == snapshot: identicalInitCond = True
            except IOError:
                pass

            # read only when identical setup AND identical initial conditions
            readExisting = identicalSetup and identicalInitCond and identicalTopology
            os.chdir(startDir)

        # when the existing working directory is not re-used, move there and remove it
        if os.path.isdir(calcDir) and not readExisting:
            logwrt.writewarning("overwriting previous " + calcDir + " directory ")
            shutil.rmtree(calcDir)

        # if needed create directory and move there
        if not os.path.isdir(calcDir): os.mkdir(calcDir)

        # prepare input file and run calculation only when previous files are not reused
        if not readExisting:

            # write topology to file
            with open(os.path.join(calcDir, AmberDriver.TOPNAME), "w") as f:
                f.write(topologyfile.topotext)

            # write coordinates and unit cell constants to file
            with open(os.path.join(calcDir, AmberDriver.CRDNAME), "w") as f:
                f.write(snapshot.crdtext)

            # write sander input
            with open(os.path.join(calcDir, AmberDriver.INPNAME), "w") as f:
                f.write(inputfile.inputtext)

            # defines command for running SANDER
            if GPU:  # run calculation with cuda support
                amberCommand = "pmemd.cuda "
            elif nCores > 1:  # run calculation with MPI (cuda has precedence)
                amberCommand = "{0} {1} sander.MPI ".format(AmberDriver.mpirun, nCores)
            else:  # run serially
                amberCommand = "sander "
            # add options for input and output files
            if ref:
                amberCommand += "-O -i {0} -o {1} -p {2} -c {3} -ref {3}".format(
                        AmberDriver.INPNAME, AmberDriver.OUTNAME, AmberDriver.TOPNAME, AmberDriver.CRDNAME)
            else:
                amberCommand += "-O -i {0} -o {1} -p {2} -c {3}".format(
                        AmberDriver.INPNAME, AmberDriver.OUTNAME, AmberDriver.TOPNAME, AmberDriver.CRDNAME)
            # split command in pieces with shlex
            amberCommand = shlex.split(amberCommand)

            # run the amber calculation
            os.chdir(calcDir)
            subprocess.run(amberCommand, stdout=None, stderr=None)
            os.chdir(startDir)

            # now create the outputData instance reading from the output of the AMBER calculation
            # AMBER files are stored only if requested with the argument store
            outputData = AmberOutput(calcDir, storeFiles=store)

        else:

            logwrt.writewarning("reading results from pre-existing {0} directory".format(calcDir))

            # now create the outputData instance reading from the output of the AMBER calculation
            # it is assumed that one wants to keep the original AMBER files
            outputData = AmberOutput(calcDir, storeFiles=True)

        return outputData

    # =============================================================================================================

    @staticmethod
    def createSolvatedMolecule(atomlabels: list, coords: list, solvsize: float, calcDir: str = "solvatedMolecule",
                               residueName: str = "CHR", solvent: str = "water", charges: bool = False,
                               newcharges: list = [], dry: bool = False, net_charge=0, gau_charges = False, gau_file=""):
        """ Use AMBER tools to generate the topology and an initial snapshot for a molecule in solvent.

        :rtype: (AmberSnapshot, AmberTopology)
        :param atomlabels:      list of atomic labels
        :param coords:          list of 3d coordinates for each atom [[x0,y0,z0], ... [x1,y1,z1]]
        :param solvsize:        size of the solvation layer added to the molecule (in Angstrom)
        :param calcDir:         name of the directory where the processing will be done
        :param residueName:     three-letter string with the name of the residue
        :param solvent:         string with the common name of the solvent
        :param charges:         True if QM charges have to be replaced with other charges read from external file
        :param newcharges:      list of QM charges provided externally to substitute the ones automatically defined
        :param dry:             True if needed to run a COBRAMM calculation in gas phase
        :param net_charge:      value of chromophore net charge
        :return:                snapshot and topology of the box with the molecule and the solvent
        """

        # check whether the solvent is available
        amberVersion = int(AmberDriver.amberVersion[0])
        if amberVersion not in AmberDriver._AMBERSOLVENTS:
            raise AmberError("amberDriver.createSolvatedMolecule: AMBER v. {} not supported".format(amberVersion))
        if solvent not in AmberDriver._AMBERSOLVENTS[amberVersion]:
            raise AmberError("amberDriver.createSolvatedMolecule: solvent {} not available".format(solvent))

        # store starting dir
        startDir = os.getcwd()

        # if the dir already exists, remove first, then create it from scratch
        if os.path.isdir(calcDir):
            logwrt.writewarning("overwriting previous " + calcDir + " directory ")
            shutil.rmtree(calcDir)
        os.mkdir(calcDir)

        # then move to the working directory
        os.chdir(calcDir)

        # write atoms and coords to xyz format
        with open("molecule.xyz", "w") as fxyz:
            fxyz.write("{}\n{}\n".format(len(atomlabels), "geometry of the molecule for AMBER input preparation"))
            for lab, c in zip(atomlabels, coords):
                fxyz.write("{0}{1:15.6f}{2:15.6f}{3:15.6f}\n".format(lab, *c))

        # convert with openbabel to define the connectivity of the molecule
        with open("babel.log", "w") as fstdout:
            subprocess.run(["obabel", "molecule.xyz", "-omol2", "-Omolecule.mol2"], stdout=fstdout,stderr=subprocess.STDOUT)

        # TODO: check existence of molecule.pdb
        # change the name of the residue
        f = fileinput.FileInput("molecule.mol2", inplace=True)
        for line in f:
            print(line.replace("UNL1", residueName).rstrip())
        f.close()

        with open("antechamber.log", "w") as fstdout:
            if gau_charges:
                command = shlex.split("antechamber -i {0} -fi gout -o molecule.mol2 -fo mol2 -c resp -at gaff2 -nc {1} -dr no".format(gau_file, net_charge))
            else:
                command = shlex.split("antechamber -i molecule.mol2 -fi mol2 -o molecule.mol2 -fo mol2 -c bcc -at gaff2 -nc {} -dr no".format(net_charge))
            subprocess.run(command, stdout=fstdout, stderr=subprocess.STDOUT)

        # TODO: check existence of molecule.mol2

        # check force field and define missing parameters with parmchk
        with open("parmchk.log", "w") as fstdout:
            if int(AmberDriver.amberVersion[0]) <= 12:
                command = shlex.split("parmchk -a Y -i molecule.mol2 -f mol2 -o molecule.frcmod")
            else:
                command = shlex.split("parmchk2 -a Y -i molecule.mol2 -f mol2 -o molecule.frcmod")
            subprocess.run(command, stdout=fstdout, stderr=subprocess.STDOUT)

        # TODO: check existence of molecule.frcmod

        # define the command to source solvent parameters and the name of the solvent box
        paramCmd = AmberDriver._AMBERSOLVENTS[amberVersion][solvent][0]
        solvName = AmberDriver._AMBERSOLVENTS[amberVersion][solvent][1]
        # set a name for the solvated system
        sysName = residueName + "_SOLV"

        # print a message on the force field that will be used
        logwrt.writelog("\nUsing force field and solvent box information obtained \nfrom {0}\n".format(
            AmberDriver._AMBERSOLVENTS[amberVersion][solvent][3]))

        # get the files that are needed with the force field and solvent box definitions
        for filename, address in AmberDriver._AMBERSOLVENTS[amberVersion][solvent][2].items():
            urllib.request.urlretrieve(address, filename)

        # create input file for leap and run leap
        with open("leap.script", "w") as f:
            inputText = AmberDriver._LEAP_PREAMBLE[amberVersion] + "\n" + \
                        "{2}\n" + \
                        "loadamberparams molecule.frcmod\n" + \
                        "{0} = loadmol2 molecule.mol2\n" + \
                        "check {0}\n"
            if dry:
                inputText += "saveAmberParm {0} molecule_dry.top molecule_dry.crd\n" + \
                "quit\n"
            else:
                if net_charge > 0:
                    inputText += "addIons {0} Cl- 0\n"
                elif net_charge < 0:
                    inputText += "addIons {0} Na+ 0\n"
                inputText += "solvateOct {0} {1} {3}\n" + \
                             "saveAmberParm {0} molecule_solv.top molecule_solv.crd\n" + \
                             "quit\n"
            f.write(inputText.format(sysName, solvName, paramCmd, solvsize))

        with open("leap.log", "w") as fstdout:
            subprocess.run(["tleap", "-f", "leap.script"], stdout=fstdout, stderr=subprocess.STDOUT)

        # TODO: check existence of topology and coordinate file, if they exist give success message

        # read coordinates and topology
        if dry:
            with open("molecule_dry.top") as f:
                topology = AmberTopology(f.read())
        else:
            with open("molecule_solv.top") as f:
                topology = AmberTopology(f.read())
                if charges:
                    topology = topology.substituteQMcharges(newcharges)

        if dry:
            with open("molecule_dry.crd") as f:
                snapshot = AmberSnapshot.readcrd(f.read())
        else:
            with open("molecule_solv.crd") as f:
                snapshot = AmberSnapshot.readcrd(f.read())

        # move back to starting directory
        os.chdir(startDir)

        # return topology and snapshot file
        return snapshot, topology


    # =============================================================================================================

    @staticmethod
    def cutdroplet(topologyfile, snapshot, spherRad: float = None, nrResidue: int = 1, calcDir: str = "tmpDroplet"):
        """ run AMBER with the initial conditions defined in inputData and store the results in outputData.

        :rtype:                 (AmberSnapshot, AmberTopology)
        :param topologyfile:    topology of the AMBER data to process in this method
        :param snapshot:        snapshot of the AMBER data to process in this method
        :param spherRad:        radius of the final droplet
        :param nrResidue:       ordinal number (1,2, ...) of the residue at the center of the droplet
        :param calcDir:         name of the directory where the processing will be done
        :return:                snapshot and topology of the droplet, cut from the original structure
        """

        # store starting dir and move to directory of the calculation
        startDir = os.getcwd()
        # then move to the work d where preprocessing files are stored, if the dir exists, remove first
        if os.path.isdir(calcDir):
            logwrt.writewarning("overwriting previous " + calcDir + " directory ")
            shutil.rmtree(calcDir)
        os.mkdir(calcDir), os.chdir(calcDir)

        # write topology to file
        with open("input00.top", "w") as f:
            f.write(topologyfile.topotext)
        # write coordinates and unit cell constants to file
        with open("input00.crd", "w") as f:
            f.write(snapshot.crdtext)

        # write first input file for cpptraj (extract snapshot, center and remap coordinates in the first unit cell)
        with open("cpptrajscript.01", "w") as f:
            f.write(" parm {0}\n trajin {1}\n".format("input00.top", "input00.crd"))
            f.write(" center :{0} \n".format(nrResidue))
            f.write(" image familiar\n")
            f.write(" trajout {0} restart\n go\n".format("input01.crd"))

        # run cpptraj
        with open("cpptrajscript.stdout", "w") as fstdout:
            subprocess.run(["cpptraj", "-i", "cpptrajscript.01"], stdout=fstdout, stderr=subprocess.STDOUT)

        # write second input file for cpptraj (create droplet with spherical radius)
        with open("cpptrajscript.02", "w") as f:
            f.write(" parm {0}\n trajin {1}\n reference {1}\n".format("input00.top", "input01.crd"))
            if spherRad: f.write(" strip :1>:{0} outprefix droplet\n".format(spherRad))
            f.write(" trajout {0} restart\n go\n".format("input02.crd"))

        # run cpptraj
        with open("cpptrajscript.stdout", "a") as fstdout:
            subprocess.run(["cpptraj", "-i", "cpptrajscript.02"], stdout=fstdout, stderr=subprocess.STDOUT)
        # handling errors, when the processed topology is not present
        if not os.path.isfile("droplet.input00.top"):
            with open("cpptrajscript.stdout") as fstdout:
                cpptrajout = fstdout.read()
            if "No atoms to strip. Skipping 'strip' for topology 'input00.top'" in cpptrajout:
                # in this case, it means that the new topology is identical to the old one
                os.rename("input00.top", "input02.top")
            else:
                raise AmberError("cpptraj unknown failure in topology stripping")
        else:
            os.rename("droplet.input00.top", "input02.top")

        # now read the processed topology and coordinate file
        with open("input02.crd", "r") as f:
            dropletSnapshot = AmberSnapshot.readcrd(f.read())
        with open("input02.top", "r") as f:
            dropletTopology = AmberTopology(f.read())

        # move back to starting directory
        os.chdir(startDir)
        # clean up processing directory
        #shutil.rmtree(calcDir)

        # return snapshot and topology for the droplet
        return dropletSnapshot, dropletTopology

    # =============================================================================================================

    @staticmethod
    def extractAtomsAndResidues(topology, snapshot):
        """Processing the PDB output obtained from the input topology and the input snapshot, this method
        extracts complete information about the atomic structure: the atomic labels, the atomic coordinates
        and the assignment of each atom to the residues defined in the topology

        :rtype:             (list, list, list)
        :param topology:    topology of the AMBER data to process in this method
        :param snapshot:    snapshot of the AMBER data to process in this method
        :return:            lists with: 1) the atomic symbols, 2) the atomic coords, 3) the residue the atom belongs to
        """

        # write topology to file
        with open("tmp.top", "w") as f:
            f.write(topology.topotext)
        # write coordinates and unit cell constants to file
        with open("tmp.crd", "w") as f:
            f.write(snapshot.crdtext)

        # create PDB file with ambpdb
        ambtext = subprocess.run(["ambpdb", "-p", "tmp.top", "-c", "tmp.crd"], check=True, stdout=subprocess.PIPE)

        # move back to starting directory and remove temporary one
        os.remove("tmp.top"), os.remove("tmp.crd")

        # initialize lists to store the data
        symbolList = []  # list containing the atomic labels
        coordsList = []  # list with the atomic coordinates
        residueList = []  # list with the number of residue to which each atom belongs

        # loop over the lines of the PDB file
        for line in ambtext.stdout.splitlines():
            if b"END" in line:  # loop until the "END" code is found
                break
            elif b"ATOM" in line:  # store the data corresponding to current atom
                symbolList.append(line[76:78].strip().decode())
                coordsList.append([float(x) for x in line[30:54].split()])
                residueList.append(int(line[22:26]))
            else:  # otherwise, there is nothing to do
                pass

        return symbolList, coordsList, residueList

    # =============================================================================================================

    @staticmethod
    def atomlist2amberformat(atomlist):
        """Convert a list of atoms to an Amber-style string,
        eg. [1,2,3,5,9,10,11,12] -> '1-3,5,9-12' """

        # first order atomlist
        atomlist = sorted(atomlist)

        # initialization of lists for start and end of the sequences
        start, end = [], []
        # the fist serie of consecutive numbers start with the first element
        start.append(atomlist[0])
        # then at each discontinuity in the numbers of atomlist, close the current sequence and start another one
        for i in range(len(atomlist) - 1):
            if atomlist[i] != atomlist[i + 1] - 1:
                start.append(atomlist[i + 1])
                end.append(atomlist[i])
        # close the last sequence with the last number
        end.append(atomlist[-1])

        # now join together the intervals and return the comma-separeted interval string
        listpieces = [str(s) if s == e else "{0}-{1}".format(s, e) for s, e in zip(start, end)]
        return "':" + ",".join(listpieces) + "'"

    # =============================================================================================================

    @staticmethod
    def amberformat2atomlist(amberformat):
        """Convert an Amber-style string that defines a list of atoms to the list of atoms with integer indices,
        eg. '1-3,5,9-12' -> [1,2,3,5,9,10,11,12] """

        # remove starting and ending "/'
        amberformat = amberformat.strip(" \'\":")
        # split the string to a list of intervals using comma as separator
        intervallist = amberformat.split(",")

        # initialize the list
        atomlist = []
        # now loop over the intervals
        for interval in intervallist:
            if "-" in interval:  # this means that the interval contains more than one element
                begin, end = [int(i) for i in interval.split("-")]
                atomlist += list(range(begin, end + 1))
            else:
                atomlist.append(int(interval))

        # now return the list
        return intervallist

    # =============================================================================================================

    @staticmethod
    def ambermask2atomlist(topologyfile, snapshot, mask, calcDir="amberMask"):
        """Convert an Amber-style string that defines a list of atoms to the list of atoms with integer indices,
        eg. '1-3,5,9-12' -> [1,2,3,5,9,10,11,12]. Use the ambmask Amber utility for the conversion. """

        # check if topologyfile is of AmberTopology type
        if not isinstance(topologyfile, AmberTopology):
            raise AmberError("amberDriver.ambermask2atomlist: first argument should be of AmberTopology type")
        # check if snapshot is of AmberSnapshot type
        if not isinstance(snapshot, AmberSnapshot):
            raise AmberError("amberDriver.ambermask2atomlist: second argument should be of AmberSnapshot type")

        # when the existing working directory is not re-used, move there and remove it
        if os.path.isdir(calcDir):
            logwrt.writewarning("overwriting previous " + calcDir + " directory ")
            shutil.rmtree(calcDir)
        os.mkdir(calcDir)

        # write topology to file
        with open(os.path.join(calcDir, AmberDriver.TOPNAME), "w") as f:
            f.write(topologyfile.topotext)
        # write coordinates and unit cell constants to file
        with open(os.path.join(calcDir, AmberDriver.CRDNAME), "w") as f:
            f.write(snapshot.crdtext)

        # define command for ambmaks
        amberCommand = shlex.split("ambmask -p {0} -c {1} -out pdb -find '{2}'".format(
            os.path.join(calcDir, AmberDriver.TOPNAME), os.path.join(calcDir, AmberDriver.CRDNAME), mask))
        # run ambmaks
        result = subprocess.run(amberCommand, check=True, stdout=subprocess.PIPE)

        # check if the run terminated with error
        if b"ERROR" in result.stdout:
            raise AmberError("amberDriver.ambermask2atomlist: ambmask error termination")
        # remove temporary directory
        shutil.rmtree(calcDir)

        # process standard output, and return  list of atoms with their integer index (starting from 1)
        return [int(n) for n in re.findall(b"ATOM *([0-9]*)", result.stdout)]


#####################################################################################################

class AmberInput:
    """This class defines the input keywords for an AMBER run. It can be instantiated in two ways,
    either with the standard constructor that sets the values of the keywords based on a series of
    input arguments, or with the readinput method that set the keywords by reading an already
    existing Amber input file. """

    def __init__(self, minimize=True, nrSteps=1, deltaT=0.001, temperature=None, pressure=None,
                 cutoff=999.0, freezeH=False, freezesolute=False, freezeLowLayer=False, mradius=999, nprntsteps=50, usevelocity=False, usePBC=False):
        """ Prepare input file by setting reasonable options depending on the input arguments of this function

        :param minimize: True if this is a minimization run, False for a MD run
        :param nrSteps: total number of steps
        :param deltaT: time step for MD (in picoseconds)
        :param temperature: temperature of the MD run (in K)
        :param pressure: pressure of the MD run (in bar)
        :param cutoff: cutoff value for evaluating the MM interaction potentials
        :param freezeH: freeze H atoms
        :param freezesolute: freeze QM atoms
        :param freezeLowLayer: freeze all the atoms above a certain distance (:param mradius) from the solute (Low-layer-like)
        :param nprntsteps: print average values every nprntsteps steps
        :param usevelocity: set initial velocity for the calculations (if False, initial velocity are zero)
        :param usePBC: use periodic boundary conditions
        """

        self._input = {}
        self._freezeoptions = {}

        # OPTIONS FOR OPTIMIZATION
        if minimize:
            self._input["imin"] = 1  # this is a minimization calculation
            self._input["ntx"] = 1  # use only coordinate, no velocity restart
            self._input["irest"] = 0
            self._input["maxcyc"] = nrSteps  # maxcyc = max nr of opt cycles
            self._input["ncyc"] = 500  # initially use steepest descent, after ncyc switch on conjugate gradient
            if not usePBC:  # if there is no unit cell, do not use periodic boundary conditions
                self._input["ntb"] = 0
            else:
                self._input["ntb"] = 1

        # OPTIONS FOR MOLECULAR DYNAMICS
        else:
            self._input["imin"] = 0  # this is a molecular dynamics run
            if usevelocity:
                self._input["ntx"] = 5  # use initial velocity
                self._input["irest"] = 1
            else:
                self._input["ntx"] = 1  # start from coordinates only
                self._input["irest"] = 0
            self._input["nstlim"] = nrSteps  # total nr of steps and ...
            self._input["dt"] = deltaT  # ... timestep
            if temperature is not None and temperature > 0.0:  # put thermostat
                self._input["ntt"] = 3
                self._input["gamma_ln"] = 2.
                self._input["temp0"] = temperature
            else:
                self._input["ntt"] = 0
            if not usePBC:  # if there are no periodic boundary conditions, pressure makes no sense
                self._input["ntb"] = 0
                self._input["ntp"] = 0
            else:
                if pressure is not None and pressure > 0.0:  # put barostat
                    self._input["ntb"] = 2
                    self._input["ntp"] = 1
                    self._input["pres0"] = pressure
                    self._input["taup"] = 2.
                else:  # simulation with fixed periodic box
                    self._input["ntb"] = 1
                    self._input["ntp"] = 0

        # OTHER COMMON OPTIONS
        self._input["cut"] = cutoff  # cutoff in potential terms evaluation
        if freezeH:
            self._input["ntf"] = 2  # freeze bonds involving H atoms
            self._input["ntc"] = 2
        if freezesolute:
            self._input["ntr"] = 1  # set keyword to freeze residue
            self._input["restraint_wt"] = 500.0
            if freezeLowLayer:
                self._input["restraintmask"] = "':1 | :1>:{}'".format(mradius)
            else:
                self._input["restraintmask"] = "':1'"


        # COMMON PRINT OPTIONS
        self._input["ig"] = -1  # random seed generator for stochastic runs
        self._input["ntxo"] = 2  # final restart file is in NETCDF format
        self._input["ntpr"] = nprntsteps  # print energy info every nprntsteps steps
        self._input["ioutfm"] = 1  # write coords, vel and forces every nprntsteps in NETCDF mdcrd
        self._input["ntwx"] = nprntsteps
        self._input["ntwv"] = -1
        self._input["ntwf"] = -1

    # ================================================================================

    def __del__(self):
        """Destroy the AmberInput data"""
        del self._input

    # ================================================================================

    def __getitem__(self, item):
        """Returns variables stored in the self.input dictionary with the
        "item" key, by using the syntax self[item]. """
        return self._input[item]

    def __setitem__(self, key, value):
        """Change the value of a directive of the sander input file: if the key already
        exists, its value will be changed, otherwise, the new key will be added"""
        self._input[key] = value

    # ================================================================================

    @property
    def inputtext(self):
        """ Return a string with sander input file constructed with
        the directives stored in the self instance of the AmberInput class

        :return: string with the input file text
        """

        inputtext = "AMBER input generated by COBRAMM wrapper\n &cntrl\n"
        for key, value in self._input.items():
            inputtext += "   {0} = {1},\n".format(key, value)
        inputtext += "/\n"
        for key, value in self._freezeoptions.items():
            inputtext += "{0}\n".format(value)

        return inputtext

    # ================================================================================

    @classmethod
    def readinput(cls, inputtext):
        """Parse the text of an amber input file, extracting the keywords and their value
         to set up a new instance of the AmberInput class

         :param inputtext: Amber input file text
         :return: new AmberInput instance with the input file directives """

        # first breaks the text in three parts: the title, the &cntrl section and what follows the &cntrl section
        regexpr = '(.*)&cntrl(.*)/(.*)'
        results = re.findall(regexpr, inputtext, re.MULTILINE | re.DOTALL)
        try:
            title = results[0][0].strip()
            cntrlsection = results[0][1].strip()
            rest = results[0][2].strip()
        except IndexError:
            raise AmberError("error parsing sander input file ...\n" +
                             "======================================= \n" +
                             "{0}".format(inputtext) +
                             "======================================= \n")

        # in the current parsing scheme, the part that follows the cntrl section is ignored,
        # in the case in which this is not null, print a warning to screen to let the user know
        if rest:
            logwrt.writewarning("ignoring sander keywords outside the cntrl section in file '{0}'".format(title))

        # initialize dictionary to save results of the parsing
        sanderdirectives = {}

        # the following regular expression looks for
        # 1) expressions of the type keyword = "whatever there is inside the double quotes"
        # 2) expressions of the type keyword = 'whatever there is inside the double quotes'
        # 3) expressions of the type keyword = value
        # making sure that what is recognized as the first (or second pattern) is not processed
        # as one of the following expressions
        regexpr = '[\w-]+ *= *".+?"|[\w-]+ *= *\'.+?\'|[\w-]+ *= *[\w-]+'
        # save the keys and their values to the dictionary
        for r in re.findall(regexpr, cntrlsection):
            sanderdirectives[r.split("=")[0].strip()] = r.split("=")[1].strip()

        # create a new instance of AmberInput and overwrite the dictionary with the keywords
        newinput = cls()
        newinput._input = sanderdirectives
        # return the AmberInput instance
        return newinput


# ####################################################################################################7

class AmberTopology:
    """This small class defines the storage for the topology file (as a string) and implements
    methods to extract information or change the topology."""

    def __init__(self, topotext):
        """Initialize the instance by storing the text of the topology file."""

        self.topotext = topotext

    # =============================================================================================================

    def __del__(self):
        """Destroy the AmberTopology data"""
        del self.topotext

    # ================================================================================

    @property
    def charges(self):
        """Extract the list of atomic charges from an Amber topology file

        :return: list of floating point numbers with the charge values for the atoms
        """

        topcharges = []  # initialize container for charges
        reading = False
        for line in self.topotext.splitlines():  # process text line by line
            if "%FLAG" in line and reading:  # *) we are reading, and there is a "%FLAG" label, indicating
                break  # _ _ _ _ _ _ _ _ _ _ _ _ _that the CHARGE section is finished
            if "CHARGE" in line:  # _ _ _ _ _ _*) we find "CHARGE" in the line, the charge section is
                reading = True  # _ _ _ _ _ _ _ _ starting here, so we can start reading
            if reading:  # _ _ _ _ _ _ _ _ _ _ *) we are reading, so try to cast the line to float numbers
                try:
                    topcharges += [float(c) * constants.MDcharge2au for c in line.split()]
                except ValueError:
                    pass

        return topcharges

    # =============================================================================================================

    def substitutecharges(self, newcharges):
        """Returns a new Amber topology file in which the charges have been substituted with an input list

        :param newcharges: list of the charges to include in the new topology file
        :return: new AmberTopology instance in which the charges have been substituted
        """

        newtopotext = ""  # initialize string for the text of the new topology file
        chargesect = False
        for line in self.topotext.splitlines():  # process text line by line
            if "%FLAG" in line and chargesect:  # *) we are reading, and there is a "%FLAG" label, indicating
                chargesect = False  # _ _ _ _ _ _ _ _that the CHARGE section is finished
            if "CHARGE" in line:  # _ _ _ _ _ _ _ *) we find "CHARGE" in the line, the charge section is
                chargesect = True  # _ _ _ _ _ _ _ _ starting here, so we can put here the new charge section
                newtopotext += line + "\n"
                newtopotext += "%FORMAT(5E16.8)" + "\n"
                for i, ch in enumerate(newcharges):
                    newtopotext += "{0:16.8E}".format(ch)
                    if (i + 1) % 5 == 0:
                        newtopotext += "\n"
                if len(newcharges) % 5 != 0: newtopotext += "\n"
            if not chargesect:  # _ _ _ _ _ _ _ _*) we are outside the charge section, so save line to newtopotext
                newtopotext += line + "\n"

        return AmberTopology(newtopotext)

    # =============================================================================================================

    def substituteQMcharges(self, newcharges):
        """Returns a new Amber topology file in which the QM charges have been substituted with an input list

        :param newcharges: list of the QM charges to include in the new topology file
        :return: new AmberTopology instance in which the QM charges have been substituted
        """

        newtopotext = ""  # initialize string for the text of the new topology file
        chargesect = False
        for line in self.topotext.splitlines():  # process text line by line
            if "%FLAG" in line and chargesect:  # *) we are reading, and there is a "%FLAG" label, indicating
                chargesect = False  # _ _ _ _ _ _ _ _that the CHARGE section is finished
                if newtopotext[-1] != "\n": #if last line of charges ended with less than 5 fields, we need to add a "\n"
                    newtopotext += "\n"
            if "CHARGE" in line:  # _ _ _ _ _ _ _ *) we find "CHARGE" in the line, the charge section is
                chargesect = True  # _ _ _ _ _ _ _ _ starting here, so we can put here the new charge section
                newtopotext += line + "\n"
                newtopotext += "%FORMAT(5E16.8)" + "\n"
            elif chargesect :
                if "%FORMAT" in line: #first line of charge section does not contain numbers but FORMAT specifications: just initialize counter to 0
                    counter = 0
                else: #we are in chargesect but not in the line containing "%FORMAT"
                    l = line.split() #works on line fields
                    for i in range(len(l)):
                        if counter < len(newcharges): #substitute the field until counter is below the number of QM charges
                            newtopotext +=  "{0:16.8E}".format(newcharges[counter]) #replace field with new charge
                        else: 
                            newtopotext += "{0:16.8E}".format(float(l[i])) #QM charges ended, so keep original charges from top file
                        counter += 1
                        if counter % 5 == 0 and counter != 0:
                            newtopotext += "\n"
            if not chargesect:  # _ _ _ _ _ _ _ _*) we are outside the charge section, so save line to newtopotext
                newtopotext += line + "\n"

        return AmberTopology(newtopotext)

    # =============================================================================================================

    @property
    def atoms(self):
        """Extract the list of Amber atomic labels from an Amber topology file

        :return: list of atomic labels
        """

        topatoms = []  # initialize container for charges
        reading = False
        for line in self.topotext.splitlines():  # process text line by line
            if "%FLAG" in line and reading:  # *) we are reading, and there is a "%FLAG" label, indicating
                break  # _ _ _ _ _ _ _ _ _ _ _ _ _that the AMBER_ATOM_TYPE section is finished
            if "%FORMAT" in line and reading:  # *) we are reading, and there is a "%FORMAT" label, this
                continue  # _ _ _ _ _ _ _ _ _ _ _ _ line can be skipped
            if "AMBER_ATOM_TYPE" in line:  # _ *) we find "AMBER_ATOM_TYPE" in the line, the section to read is
                reading = True  # _ _ _ _ _ _ _ _ starting here, so we can start reading
                continue
            if reading:  # _ _ _ _ _ _ _ _ _ _ *) we are reading, so store the line
                topatoms += line.split()

        return topatoms


# ####################################################################################################7

class AmberSnapshot:
    """This simple class defines a convenient way to store a single snapshot of
    an Amber calculation, including the coordinates of the atoms, and possibly
    their velocities and the unit cell in a calculation with PBCs.
    The object can be initialized with the standard initializer, using
    the three arguments explicitely, or the class function 'readfromcrd'
    can be used, that instantiate AmberSnapshot by reading the data from
    an input crd file. """

    def __init__(self, coords, unitcell=None, velocity=None):
        """Define an instance of AmberSnapshot from arguments given explicitely:
        the coordinates, the unit cell and the velocities. The unit-cell and the
        velocities can be omitted, in this case it is assumed that the calculation
        is not periodic (for the unit cell) and that velocity information is
        not present (for the velocities).

        :param coords: input coordinates of the atoms [ [x_0, y_0, z_0] ... [x_N, y_N, z_N] ]
        :param unitcell: unit cell definition, list of six floats, None if PBC are not present
        :param velocity: input velocities (same format as coords), None if not velocity is present
        """

        # store the input variables, coordinates, velocity and topology file
        self.coords = coords
        self.unitcell = unitcell
        self.velocity = velocity

    # =============================================================================================================

    def __del__(self):
        """Destroy the AmberSnapshot data"""
        del self.coords
        del self.unitcell
        del self.velocity

    # =============================================================================================================

    def __eq__(self, other):
        """Equality of two instance of the AmberSnapshot class"""

        if isinstance(other, AmberSnapshot):
            return self._identicalData(self.coords, other.coords) and \
                   self._identicalData(self.unitcell, other.unitcell) and \
                   self._identicalData(self.velocity, other.velocity)
        return NotImplemented

    # ================================================================================

    @staticmethod
    def _identicalData(var1, var2):
        """Extend the comparison between two variables, as defined in the numpy function np.allclose(),
         by considering the case in which one of the two variables or both can assume a None value.
         When the two variables are both None, the variables are considered equal."""

        if var1 is None and var2 is None:  # when both variables are None, return True
            return True
        elif var1 is None or var2 is None:  # when only one of the two variables is None, return False
            return False
        else:  # otherwise, make the comparison as two lists of numbers
            try:
                return np.allclose(var1, var2, atol=1.e-6)
            except ValueError:
                return False

    # =============================================================================================================

    @classmethod
    def readcrd(cls, crdText):
        """Read CRD file and return an instance of AmberSnapshot with the data read from
        the crd file: the coordiantes, + unit cell and velocities when available.

        :param crdText: text of the CRD file
        :return: AmberSnapshot with the data read from the crd file
        """

        crdLines = crdText.splitlines()  # process text line by line

        try:
            crdLines.pop(0)  # skip first line
            nAtoms = int(crdLines.pop(0).split()[0])  # second line has the number of atoms

            # Decide whether the file contains unitcell and velocity information, depending on the lines
            if len(crdLines) == int(math.ceil(nAtoms * 3.0 / 6.0)):
                readVelocity, readUnitCell = False, False
            elif len(crdLines) == int(math.ceil(nAtoms * 3.0 / 6.0)) + 1:
                readVelocity, readUnitCell = False, True
            elif len(crdLines) == 2 * int(math.ceil(nAtoms * 3.0 / 6.0)):
                readVelocity, readUnitCell = True, False
            elif len(crdLines) == 2 * int(math.ceil(nAtoms * 3.0 / 6.0)) + 1:
                readVelocity, readUnitCell = True, True
            else:
                raise AmberError("Unexpected number of lines reading an Amber coordinate file")

            # read the coordinates (keep them in Angstrom, as in the crd files)
            coords = []
            for i in range(int(math.ceil(nAtoms * 3.0 / 6.0))):
                coords += [float(c) for c in crdLines.pop(0).split()]
            coords = np.array(coords).reshape((-1, 3)).transpose()

            # when requested, file should contain velocities as well (convert them to au from the strange Amber units)
            if readVelocity:
                vels = []
                for i in range(int(math.ceil(nAtoms * 3.0 / 6.0))):
                    vels += [float(c) / constants.Bohr2Ang / constants.MDtime2au for c in crdLines.pop(0).split()]
                vels = np.array(vels).reshape((-1, 3)).transpose()
            else:
                vels = None

            # read the unit cell sizes and angles
            if readUnitCell:
                unitcellline = crdLines.pop(0).split()
                unitcell = np.array([float(c) for c in unitcellline[0:3]] +
                                    [float(c) * constants.Deg2Rad for c in unitcellline[3:]])
            else:
                unitcell = None

        except IndexError:
            coords, vels, unitcell = None, None, None

        # return coords, vels and unitcell (vels can be None if no velocity was found)
        return cls(coords, unitcell, vels)

    # =============================================================================================================

    @property
    def crdtext(self):
        """ Write to a CRD file the lists of atom coordinates, velocities and unit cell constants
        defined in the instance of AmberSnapshot

        :return: string with the text of the CRD file with the data from AmberSnapshot
        """

        crdText = ""  # initialize string of the crd text

        # first line is a comment
        crdText += " crd generated by COBRAMM - {0}\n".format(
            time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
        crdText += "  {0:d}\n".format(len(self.coords[0]))  # write number of atoms

        # write coordinates
        coords = np.array(self.coords).transpose().reshape((-1))
        for i, c in enumerate(coords):
            crdText += "{0:12.7f}".format(c)
            if (i + 1) % 6 == 0: crdText += "\n"
        if len(coords) % 6 != 0: crdText += "\n"

        # write velocities
        if self.velocity is not None:
            velocities = np.array(self.velocity).transpose().reshape((-1))
            for i, c in enumerate(velocities):
                crdText += "{0:12.7f}".format(c * constants.Bohr2Ang * constants.MDtime2au)
                if (i + 1) % 6 == 0: crdText += "\n"
            if len(velocities) % 6 != 0: crdText += "\n"

        # write unit cell constants
        if self.unitcell is not None:
            for c in self.unitcell[0:3]:
                crdText += "{0:12.7f}".format(c)
            for c in self.unitcell[3:]:
                crdText += "{0:12.7f}".format(c / constants.Deg2Rad)
            crdText += "\n"

        return crdText


#####################################################################################################

class AmberOutput:
    """The AmberOutput class defines the object that stores output information collected from an Amber output file."""

    def __init__(self, outDir, storeFiles=False):
        """ Constructor of the AmberOutput instance. Needs as argument the path of the directory where the
         calculation is stored. The logical variables storeFiles defines whether the Amber files are removed or not
         when cleaning the AmberOutput instance. It should be set to true when the user is interested in keeping
          the original Amber files for future analysis of the results.

         :param outDir:  path of the directory where the calculation is stored
         :param storeFiles: Amber files are removed or not when cleaning the AmberOutput instance
         """

        # initialize storeFiles to True, in this way when the constructor dyes unexpectedly we can save the files
        self.storeFiles = True

        # store the path of the directory where the Amber files are read
        self.outDir = outDir

        # setup a dictionary that will be used to store the data read from Amber output
        self.dataDict = {
            "type": None,  # the calculation is an optimization or a MD run
            "topology": None,  # topology file of the system
            "logtext": None,  # log file of the amber calculation
            "mdcrdstart": None,  # nr of the first step that should be found in the mdcrd file
            "mdcrdfile": None,  # path of the mdcrd file
            "restrtfile": None,  # path of the restrt file
            "energy": None,  # dictionary of the values of total energy {step(nr/time):Energy}
            "temperature": None,  # dictionary of the values of temperature (only MD calculation)
            "pressure": None,  # dictionary of the values of pressure (only MD calculation)
            "volume": None,  # dictionary of the values of volume (only MD calculation)
        }

        # now populate the dictionary, starting from the topology
        with open(os.path.join(outDir, AmberDriver.TOPNAME), "r") as top:
            self.dataDict["topology"] = top.read()

        # now read the output file
        with open(os.path.join(outDir, AmberDriver.OUTNAME), "r") as out:
            self.dataDict["logtext"] = out.read()

        # check whether the calculation finished with success
        if "Maximum number of minimization cycles reached." not in self.dataDict["logtext"] and \
                "5.  TIMINGS" not in self.dataDict["logtext"]:
            raise AmberError("sander calculation in {0} ended with error".format(outDir))

        # check whether the calculation is an optimization or a molecular dynamics run
        regExpr = "imin *= *([0-9]*),"
        results = re.findall(regExpr, self.dataDict["logtext"])
        if results:

            if int(results[0]) == 1:  # this calculation is an optimization
                self.dataDict["type"] = "opt"
                # extract the values from the out file
                steps, energies = AmberOutput._grepAmberOptimization(self.dataDict["logtext"])
                # store the values in the dictionary
                self.dataDict["energy"] = dict(zip(steps, energies))

            else:  # this calculation is a molecular dynamics simulation
                self.dataDict["type"] = "MD"
                # extract the values from the out file
                times, energies, temperatures, pressures, volumes = \
                    AmberOutput._grepAmberDynamics(self.dataDict["logtext"])
                # store the values in the dictionary
                self.dataDict["energy"] = dict(zip(times, energies))
                self.dataDict["temperature"] = dict(zip(times, temperatures))
                if pressures is not None: self.dataDict["pressure"] = dict(zip(times, pressures))
                if volumes is not None: self.dataDict["volume"] = dict(zip(times, volumes))

        else:
            raise AmberError("error reading Amber file {0}: cannot determine 'imin' value".format(
                AmberDriver.OUTNAME))

        # define the list of steps that are expected in the mdcrd file
        regExpr = "ntwx *= *([0-9]*),"
        results = re.findall(regExpr, self.dataDict["logtext"])
        self.dataDict["mdcrdstart"] = int(results[0])

        # store the path of the mdcrd file
        self.dataDict["mdcrdfile"] = os.path.join(self.outDir, AmberDriver.MDCRDNAME)
        # store the path of the restrtfile file
        self.dataDict["restrtfile"] = os.path.join(self.outDir, AmberDriver.RESTRTNAME)

        # store the value of storeFiles
        self.storeFiles = storeFiles

    # ================================================================================

    @staticmethod
    def _grepAmberOptimization(fileText):
        """ use regexpr to extract the list of total energies from AMBER output file

        :param fileText: text of the AMBER output file
        :return: two lists, for each step the first contains the step number and the second the MM energy
        """

        reString = "NSTEP.*?\n *([0-9]* *\-?[0-9]*.[0-9]*E[\+\-][0-9]*).*?\n"
        results = re.findall(reString, fileText)

        enList, stepList = [], []
        for strg in results:
            stepList.append(int(strg.split()[0]))
            enList.append(float(strg.split()[1]) / constants.Hatree2kcalmol)

        if stepList[-1] == stepList[-2]:
            stepList.pop()
            enList.pop()

        return stepList, enList

    # ================================================================================

    @staticmethod
    def _grepAmberDynamics(fileText):
        """ use regexpr to extract information from AMBER output for a Molecular Dynamics run

        :param fileText: text of the AMBER output file
        :return: lists of time, energy, temperature, pressure, volume for each step of the output file
        """

        # select the first part of the text, up to the string "A V E R A G E S   O V E R"
        nmax = fileText.find("A V E R A G E S   O V E R")

        regExpr = "TIME\(PS\) = *([0-9]*\.[0-9]*) *TEMP\(K\) = *(\-?[0-9]*\.[0-9]*) *" \
                  "PRESS = *(\-?[0-9]*\.[0-9]*).*\n.*Etot   = *(\-?[0-9]*\.[0-9]*)"
        results = re.findall(regExpr, fileText[0:nmax])

        tList, enList, tempList, presList = [], [], [], []
        for x in results:
            tList.append(float(x[0]) * constants.ps2au)
            tempList.append(float(x[1]))
            presList.append(float(x[2]))
            enList.append(float(x[3]) / constants.Hatree2kcalmol)

        # when the pressure is always == 0, this is a constant cell / non periodic calculatio, pressure means nothing
        if np.all(np.array(presList) == 0): presList = None

        regExpr = "VOLUME     = *(\-?[0-9]*\.[0-9]*)"
        results = re.findall(regExpr, fileText[0:nmax])

        volList = []
        if results:
            for x in results:
                volList.append(float(x) / (constants.Bohr2Ang ** 3))
        else:
            try:
                regExpr = "Box X = *(\-?[0-9]*\.[0-9]*) *Box Y = *(\-?[0-9]*\.[0-9]*) *Box Z = *(\-?[0-9]*\.[0-9]*) *" \
                          "\n *Alpha = *(\-?[0-9]*\.[0-9]*) *Beta  = *(\-?[0-9]*\.[0-9]*) *Gamma = *(\-?[0-9]*\.[0-9]*)"
                unitCellStrings = re.findall(regExpr, fileText[0:nmax])[0]
                unitCell = [float(x) / constants.Bohr2Ang for x in unitCellStrings[0:3]] + \
                           [float(x) * constants.Deg2Rad for x in unitCellStrings[3:6]]
                volume = unitCell[0] * unitCell[1] * unitCell[2] * math.sqrt(
                    1.0 - math.cos(unitCell[3]) ** 2 - math.cos(unitCell[4]) ** 2 - math.cos(unitCell[5]) ** 2
                    + 2.0 * math.cos(unitCell[3]) * math.cos(unitCell[4]) * math.cos(unitCell[5]))
                volList = [volume] * len(tList)
            except IndexError:
                volList = None

        return tList, enList, tempList, presList, volList

    # ================================================================================

    def __del__(self):
        """ destroy the data stored in the class instance, and possibly
        remove the directory of the calculation unless it is requested otherwise """

        if not self.storeFiles:
            shutil.rmtree(self.outDir)

        del self.dataDict
        del self.outDir
        del self.storeFiles

    # ================================================================================

    def __getattr__(self, item):
        """Catches non-standard class attributes, by using the "item" key to access the
        datadictionary self.dataDict that store the results of the calculation."""
        return self.dataDict[item]

    def __getitem__(self, item):
        """Returns variables stored in the dataDict dictionary, by using
        the syntax self[item]. When item is an integer, return one of the
        snapshots of the calculation, otherwise try to access the
        dataDict dictionary using item as key. """
        if isinstance(item, int):
            return self.snapshot(item)
        else:
            return self.dataDict[item]

    # ================================================================================

    def snapshot(self, nrSnap=-1):
        """Extract one of the snapshots of the Amber calculation, using the nrSnap integer
         index. Returns the three usual arguments that define a snapshot: the coordinates (Ang),
         the unit cell constants (Ang + rad) and the velocities (au), stored in an instance
         of the object AmberSnapshot. If any of these variables is not present, it is set to None.
         When nrSnap is equal to -1 or to the last snapshot, the restart file is read instead
         of the mdcrd file, and the final results are returned.

         :param nrSnap: requested snapshot to extract
         :return: AmberSnapshot object with the requested snapshot"""

        # this calculation is an optimization
        if self.type == "opt":
            # first check if the requested snapshot is present
            laststep = sorted(self.energy.items())[-1][0]
            if nrSnap == -1 or nrSnap == laststep:
                lastSnapshot = True  # when last snapshot is requested, results will be read from restart file
                snapKey = None
            else:
                lastSnapshot = False  # otherwise, we need to ensure that the step is present in the mdcrd file
                snapList = [n for n in self.energy.keys() if n % self.mdcrdstart == 0]
                if nrSnap not in snapList:
                    logwrt.writewarning("requested snapshot is not present in Amber results")
                    return None
                snapKey = snapList.index(nrSnap)
        else:
            # in the case of an md run, do thing in an immediate way:
            # -1 means read snapshot from the restart file
            # n > 0 means read n-th snapshot stored in the mdcrd
            if nrSnap == -1:
                lastSnapshot = True  # results will be read from restart file
                snapKey = None  # there is no need of an integer index
            else:
                lastSnapshot = False  # results will be read from mdcrd file
                snapKey = nrSnap  # the key is the integer argument of the function

        #  ############### READ COORDINATES AND VELOCITIES ###############
        if lastSnapshot:  # the final step is requested, use the restart file
            with netcdf.netcdf_file(self.restrtfile) as restrt:
                # read coordinates
                coords = copy.deepcopy(restrt.variables["coordinates"].data)
                coords = self._toAtomicUnits(coords, restrt.variables["coordinates"].units)
                # read velocities
                try:
                    vels = copy.deepcopy(restrt.variables["velocities"].data)
                    if np.all(vels == 0.0):
                        vels = None
                    else:
                        vels = self._toAtomicUnits(vels, restrt.variables["velocities"].units)
                except KeyError:
                    vels = None
        else:  # an intermediate snapshot is requested, read coordinates and velocities from mdcrd
            with netcdf.netcdf_file(self.mdcrdfile) as mdcrd:
                # read coordinates
                coords = copy.deepcopy(mdcrd.variables["coordinates"].data[snapKey])
                coords = self._toAtomicUnits(coords, mdcrd.variables["coordinates"].units)
                # read velocities
                try:
                    vels = copy.deepcopy(mdcrd.variables["velocities"].data[snapKey])
                    if np.all(vels == 0.0):
                        vels = None
                    else:
                        vels = self._toAtomicUnits(vels, mdcrd.variables["velocities"].units)
                except KeyError:
                    vels = None

        #  ############### READ UNIT CELL CONSTANTS ###############
        # for the unit cell is a bit more complicated:
        # - for an optimization run, data is always in the restrt
        # - for a dynamical run, either use restrt or mdcrd depending on the step requested
        if self.type == "opt" or lastSnapshot:
            with netcdf.netcdf_file(self.restrtfile) as restrt:
                try:
                    # unit cell information is always read from the restrt file, regardless the step requested
                    cell_lengths = copy.deepcopy(restrt.variables["cell_lengths"].data)
                    cell_lengths = self._toAtomicUnits(cell_lengths, restrt.variables["cell_lengths"].units)
                    cell_angles = copy.deepcopy(restrt.variables["cell_angles"].data)
                    cell_angles = self._toAtomicUnits(cell_angles, restrt.variables["cell_angles"].units)
                    unitcell = np.array(list(cell_lengths) + list(cell_angles))
                except KeyError:  # when keyerror is raised, it means that this is not a periodic calculation
                    unitcell = None
        else:
            with netcdf.netcdf_file(self.mdcrdfile) as mdcrd:
                try:
                    # unit cell information is always read from the restrt file, regardless the step requested
                    cell_lengths = copy.deepcopy(mdcrd.variables["cell_lengths"].data[snapKey])
                    cell_lengths = self._toAtomicUnits(cell_lengths, mdcrd.variables["cell_lengths"].units)
                    cell_angles = copy.deepcopy(mdcrd.variables["cell_angles"].data[snapKey])
                    cell_angles = self._toAtomicUnits(cell_angles, mdcrd.variables["cell_angles"].units)
                    unitcell = np.array(list(cell_lengths) + list(cell_angles))
                except KeyError:  # when keyerror is raised, it means that this is not a periodic calculation
                    unitcell = None

        # now transpose the data to store them in the COBRAMM typical format [Xlist, Ylist, Zlist]
        if coords is not None:
            coords = np.transpose(coords)
        if vels is not None:
            vels = np.transpose(vels)

        return AmberSnapshot(coords, unitcell, vels)

    # ================================================================================

    @property
    def gradient(self):
        """Extract the forces stored in the mdcrd output file of an Amber MD calculation,
        using the nrSnap integer index. If the force is not present in the calculation,
        return None. """

        if self.type == "opt":  # this calculation is an optimization, the forces are not stored
            return None

        else:  # the calculation is an MD run, we can get the forces from the mdcrd file
            with netcdf.netcdf_file(self.mdcrdfile) as mdcrd:
                forces = copy.deepcopy(mdcrd.variables["forces"].data)
                return - self._toAtomicUnits(forces, mdcrd.variables["forces"].units)

    # ================================================================================

    @staticmethod
    def _toAtomicUnits(data, units):
        """ Convert data read from NetCDF file to the units used internally by COBRAMM.
            Internal units are: Angstrom for length, a.u. for velocity, radiants for angles, au for forces

        :param data: numerical data that is converted
        :param units: byte object with the units of the input data
        :return:
        """
        if units.decode() == 'angstrom':
            return data
        elif units.decode() == 'angstrom/picosecond':
            # WARNING!!!!! despite the units that are written to the NETCDF file, the numbers for the
            # atom velocities are stored in the usual AKMA system of units.
            # this means that the conversion factor is  1 / constants.Bohr2Ang / constants.MDtime2au
            # and not 1 / constants.Bohr2Ang / constants.ps2au
            return data / constants.Bohr2Ang / constants.MDtime2au
        elif units.decode() == 'degree':
            return data * constants.Deg2Rad
        elif units.decode() == "kilocalorie/mole/angstrom":
            return data / constants.Hatree2kcalmol * constants.Bohr2Ang
        else:
            raise AmberError("{0} unit not yet implemented".format(units.decode()))

    # ================================================================================

    def optimizationPlot(self, pdfFileName=None):
        """ print graph that shows the the energy over the optimization steps

        :param pdfFileName: name of the output pdf file
        """

        if self.type == "opt":

            # plot data with matplotlib
            import matplotlib.pyplot as plt

            # extract results in a format that can be plot
            lists = sorted(self.energy.items())
            x, y = zip(*lists)

            # plot data with matplotlib
            plt.plot(x, y, label="optimization energy", c="#e6194B", ls="-", lw=2, marker="o")
            plt.ylabel("energy / Hartree"), plt.xlabel("optimization step")

            # write to PDF file or show to screen
            if pdfFileName:
                plt.savefig(pdfFileName, bbox_inches='tight')
            else:
                plt.show()

        else:
            logwrt.writewarning("the calculation is not an optimization run, skipping analysis... ")

    # ================================================================================

    def dynamicsPlot(self, pdfFileName=None):
        """ print graph that shows some variables over the times of the dynamics

        :param pdfFileName: name of the output pdf file
        """

        if self.type == "MD":

            # plot data with matplotlib
            import matplotlib.pyplot as plt
            # create a figure
            fig = plt.figure(num=None, figsize=(8, 10), dpi=100, facecolor='w', edgecolor='k')

            # define the plots that are needed
            if self.volume is None and self.pressure is None:
                nrgraphs = 2
            else:
                nrgraphs = 3
            graphA = fig.add_subplot(nrgraphs, 1, nrgraphs)  # this will be the graph with the energy
            graphB = fig.add_subplot(nrgraphs, 1, nrgraphs - 1,
                                     sharex=graphA)  # this will be the temperature
            graphC1, graphC2 = None, None
            if self.volume is not None or self.pressure is not None:
                graphC1 = fig.add_subplot(3, 1, 1, sharex=graphA)
                graphC2 = graphC1.twinx()

            # plot energy
            x, y = zip(*sorted(self.energy.items()))
            ln1 = graphA.plot(np.array(x) / constants.ps2au, y, label="total energy", c="#4363d8", ls="-", lw=2)

            # plot temperature and average
            x, y = zip(*sorted(self.temperature.items()))
            y2 = [sum(y[0:i + 1]) / float(i + 1) for i in range(len(y))]
            ln2 = graphB.plot(np.array(x) / constants.ps2au, y, label="inst. temperature", c="#e6194B", ls="--", lw=1)
            ln3 = graphB.plot(np.array(x) / constants.ps2au, y2, label="avg.  temperature", c="#e6194B", ls="-", lw=2)

            # plot pressure and average
            if self.pressure is not None:
                x, y = zip(*sorted(self.pressure.items()))
                y2 = [sum(y[0:i + 1]) / float(i + 1) for i in range(len(y))]
                ln4 = graphC1.plot(np.array(x) / constants.ps2au, y, label="inst. pressure", c="#3cb44b", ls="--", lw=1)
                ln5 = graphC1.plot(np.array(x) / constants.ps2au, y2, label="avg.  pressure", c="#3cb44b", ls="-", lw=2)
            else:
                ln4, ln5 = [], []
            # plot volume
            if self.volume is not None:
                x, y = zip(*sorted(self.volume.items()))
                ln6 = graphC2.plot(np.array(x) / constants.ps2au, np.array(y) * constants.Bohr2Ang ** 3,
                                   label="volume", c="#ffe119", ls="-", lw=2)
            else:
                ln6 = []

            # define labels
            graphA.set_ylabel("E / Hartree"), graphA.set_xlabel("t / ps")
            graphB.set_ylabel("T / K")
            if graphC1: graphC1.set_ylabel("p / bar")
            if graphC2: graphC2.set_ylabel("V / $\AA^3$")

            # create legends
            graphA.legend(ln1, [line.get_label() for line in ln1])
            graphB.legend(ln2 + ln3, [line.get_label() for line in ln2 + ln3])
            if graphC1: graphC1.legend(ln4 + ln5 + ln6, [line.get_label() for line in ln4 + ln5 + ln6])

            # write to PDF file or show to screen
            if pdfFileName:
                plt.savefig(pdfFileName, bbox_inches='tight')
            else:
                plt.show()

        else:
            logwrt.writewarning("the calculation is not a dynamics simulation, skipping analysis... ")


# #####################################################################################################

if __name__ == '__main__':
    # test run of the classes contained in this amberDriver module

    # use files from example calculations
    topofile = "../examples/TURBOMOLE/QMMM_OPTXG_DSCF/real.top"
    crdfile = "../examples/TURBOMOLE/QMMM_OPTXG_DSCF/real.crd"

    # create an instance of the Amber interface
    AmberDriver()

    # read coordinates and topology
    with open(topofile) as ftopo:
        topologyTest = AmberTopology(ftopo.read())
    with open(crdfile) as fcrd:
        snapshotTest = AmberSnapshot.readcrd(fcrd.read())

    # prepare the input for a minimization run
    optinput = AmberInput(minimize=True, nrSteps=10, nprntsteps=1, cutoff=9.)
    # run the optimization
    optoutput = AmberDriver.run(optinput, topologyTest, snapshotTest, calcDir="amberOpt", overwrite=False)
    optoutput.storeFiles = False  # remove directories at the end

    # plot a graph of the optimization energy to file
    optoutput.optimizationPlot("testopt.pdf")

    # define input coordinates and unit cell for restart with the MD
    optsnapshot = optoutput.snapshot(-1)

    # prepare the input for a molecular dynamics
    mdinput = AmberInput(minimize=False, freezeH=True, nrSteps=30, deltaT=0.001, nprntsteps=2,
                         temperature=300., pressure=1.0, cutoff=9.)
    # run the calculation
    mdoutput = AmberDriver.run(mdinput, topologyTest, optsnapshot, calcDir="amberMD", overwrite=False)
    mdoutput.storeFiles = False  # remove directories at the end

    # plot an analysis of the molecular dynamics trajectory
    mdoutput.dynamicsPlot("testmd.pdf")

    # extract final results
    mdfinalsnapshot = mdoutput.snapshot(-1)
