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

######################################################################################################################

# import statements of module from python standard library

import os  # filesystem utilities

# imports of local modules

import logwrt  # manages log file output + start/end procedures
import constants  # values of physical constants and conversion factors

# import of local classes
from QMOutput import QMOutput


#####################################################################################################

class SharcInterfaceInput:
    """Store the data for the execution of sharc QM interface """

    def __init__(self, qmcode, coords, symbols, otheropt):
        """Constructor of the SharcInterfaceInput class, requires these arguments:
        - qmcode: is a string with the name of the QM code interfaced through the SHARC package
        - coods: list of the actual coordinates of the atoms of the QM part
        - symbols: list of the atomic symbols of the atoms of the QM part
        """

        # initialize attributes of the SharcInterfaceInput class instance
        # name of the QM program
        self.qmcode = None
        # sharc-qm interface name
        self.interface = None
        # coordinates and symbols of the molecule
        self.coords, self.symbols = None, None
        # dictionary with other options
        self.otheropt = {}

        # store the input variable
        self.qmcode = qmcode
        # construct the name of the sharc interface python executable
        self.interface = "SHARC_{0}.py".format(qmcode.upper())
        # store the information about the molecular geometry
        self.coords = coords
        self.symbols = symbols

        # process the dictionary with the other options, to check and complete it
        for key in ["charges", "field"]:
            if key in otheropt:
                self.otheropt[key] = otheropt[key]
            else:
                self.otheropt[key] = None

        # read the QM.in file which should be present in the current directory where COBRAMM was launched
        # process the QM.in file to get the sections that need to be reproduced in the sharc interface
        rdsymbols, rdcoords = [], []
        with open("QMMM.in") as fin:
            # first line is the number of atoms (geometry written here in QM.in is the REAL one)
            natom = int(fin.readline().strip())
            # second line is a title
            self.title = fin.readline().strip()
            # then there are: atom symbol, atom coords, atom velocities
            for i in range(natom):
                line = fin.readline().strip()
                rdsymbols.append(line.split()[0])
                rdcoords.append(line.split()[1:4])
                # rdvel.append(line.split()[4:])
            # then the rest is a specification of the type of calculation
            self.qmininfo = fin.read()

    # =============================================================================================================

    def fileText(self):
        """Prepare the text of QM.in file for the modelH structure, adapted from the QM.in file
        for the real structure that has been used to set up the calculation"""

        # initialize the text of the input file
        inputText = ""

        # first line of QM.in file is the number of atoms
        inputText += "{0}\n".format(len(self.coords[0]))
        # second line is the title
        inputText += "{0}\n".format(self.title)
        # following lines are atom symbols and geometry
        for s, x, y, z in zip(self.symbols, *self.coords):
            inputText += "{0:2s}{1:15.6f}{2:15.6f}{3:15.6f}\n".format(s, x, y, z)

        # the rest is the specification of hte QM options for the SHARC interface
        # process the text, to ensure that the line "unit angstrom" is present
        unitLine = False
        for line in self.qmininfo.splitlines():
            if line.split() and line.split()[0] == "unit":
                inputText += "unit angstrom\n"
                unitLine = True
            else:
                inputText += line + "\n"
        if not unitLine:
            inputText += "unit angstrom\n"

        # return the complete text of the input file
        return inputText

    # =============================================================================================================

    def chargeText(self):
        """Prepare a dat file with information on the charges and on the points where to
        compute electric field. """

        # initialize the text of the input file
        chargeText = ""

        if self.qmcode.upper() == "MOLCAS":
            if self.otheropt["charges"]:
                chargeText += "XField\n{0}\n".format(len(self.otheropt["charges"][1]))
                for ch, x, y, z in zip(self.otheropt["charges"][1], *self.otheropt["charges"][0]):
                    chargeText += "{1:14.8f}{2:14.8f}{3:14.8f}{0:16.8f} 0.0 0.0 0.0\n".format(
                        ch, x / constants.Bohr2Ang, y / constants.Bohr2Ang, z / constants.Bohr2Ang)
            if self.otheropt["field"]:
                chargeText += "EFLD\n{0}\n".format(len(self.otheropt["field"][1]))
                for x, y, z in zip(*self.otheropt["field"][0]):
                    chargeText += "{0:14.8f}{1:14.8f}{2:14.8f}\n".format(
                        x / constants.Bohr2Ang, y / constants.Bohr2Ang, z / constants.Bohr2Ang)

        elif self.qmcode.upper() == "ORCA":
            if self.otheropt["charges"]:
                chargeText += "{0}\n".format(len(self.otheropt["charges"][1]))
                for ch, x, y, z in zip(self.otheropt["charges"][1], *self.otheropt["charges"][0]):
                    chargeText += "{1:14.8f}{2:14.8f}{3:14.8f}{0:16.8f}\n".format(
                    ch, x/constants.Bohr2Ang, y/constants.Bohr2Ang, z/constants.Bohr2Ang)

        else:
           logwrt.fatalerror('Cannot write charge.dat for SHARC_{0}.py SHARC interface'.format(self.qmcode))

        # return the complete text of the dat file
        return chargeText


###################################################################################################################

class SharcInterfaceOutput(QMOutput):

    def __init__(self, name, calcdir):

        # define options of the base class
        QMOutput.__init__(self)

        # Add additional SHARC-specific attributes to the data dictionary
        # FILE NAMES, DIRECTORIES, INPUT/OUTPUT
        self.dataDict["outfile"] = None  # string with the text of the QM.out file, initialized to None
        self.dataDict["name"] = name  # name of the calc, it is the base name of the I/O files (.com, .chk, .log)
        self.dataDict["dir"] = calcdir  # name of the directory where the output files are stored
        self.dataDict["efieldfile"] = None  # string with the text of the efield.dat, when present
        self.dataDict["grad_chargesfile"] = None  # string with the text of the grad_charges file, when present
        # INFO ON QM EXECUTION
        self.dataDict["termination"] = None   # final termination: 0 = Normal termination, 1 = Error termination
        self.dataDict["errormsg"] = None   # in case of Error termination, error message
        # ADDITIONAL PHYSICAL INFORMATION
        self.dataDict["statelabels"] = {}  # dictionary to store the definitions of the electronic states
        self.dataDict["elfield"] = None   # electric field at the position of hte point charges
        self.dataDict["natoms"] = None   # number of atoms in the QM calculation

        logwrt.writelog("\n", 1)
        logwrt.writelog("Reading output files from SHARC interface\n", 1)
        logwrt.writelog("Currently working in directory {0}\n".format(os.getcwd()), 1)
        logwrt.writelog("Reading files from directory {0}\n".format(calcdir), 1)
        logwrt.writelog("Content of the calculation dir:\n{0}\n".format(os.listdir(calcdir)), 1)
        logwrt.writelog("\n", 1)

        # read the QM.out file text
        try:
            with open(os.path.join(calcdir, name+".out")) as fout:
                self.dataDict["outfile"] = fout.read()
        except IOError:
            logwrt.fatalerror('Cannot find the output file from SHARC QM interface {0}...\n'
                              'terminating COBRAMM execution'.format(name+".out"))

        # loop over the lines of the QM.out file
        linebyline = self.dataDict["outfile"].splitlines()
        for i, line in enumerate(linebyline):

            # the line is the starting line of a new section
            if line.split() and line.split()[0] == "!":

                if line.split()[1] == "1":  # hamiltonian matrix block
                    nstates = int(linebyline[i+1].split()[0])
                    # THE FOLLOWING IS A DUMMY DEFINITION OF THE OPTIMIZED STATE,
                    # IT DOES NOT REALLY MATTER BUT IS NEEDED BY COBRAMM
                    self.dataDict["optstate"] = nstates-1
                    self.log += 'Number of QM electronic states: {0:3d}\n'.format(nstates)
                    for istate in range(nstates):
                        QMenergy = float(linebyline[i+2+istate].split()[2*istate])
                        self.dataDict["energy"][istate] = QMenergy
                        self.log += 'QM energy for state {0:3d} is {1:12.8f} Hartree\n'.format(istate+1, QMenergy)

                # elif line.split()[1] == "2":  # dipole moment
                #     continue

                elif line.split()[1] == "3":  # gradient vectors
                    nstates = int(line.split("(")[1].split("x")[0])
                    j = i+1
                    # loop to read all the gradients available in the output
                    for istate in range(nstates):
                        # first line is the size of the gradient and the state label
                        size, statelabel = linebyline[j].split("!")
                        self.dataDict["natoms"] = int(size.split()[0])
                        x, y, z = [], [], []
                        # read the gradient
                        for jatom in range(self.dataDict["natoms"]):
                            linesplit = linebyline[j + 1 + jatom].split()
                            x.append(float(linesplit[0])), y.append(float(linesplit[1])), z.append(float(linesplit[2]))
                        # save results that have been extracted
                        self.dataDict["statelabels"][istate] = statelabel
                        self.dataDict["gradient"][istate] = x, y, z
                        # move to next gradient section
                        j += self.dataDict["natoms"] + 1

        self.log += "\n"

        # read the grad_charges file when present (gradient of the point charges)
        if os.path.exists(os.path.join(calcdir, "grad_charges")):
            with open(os.path.join(calcdir, "grad_charges")) as fout:
                self.dataDict["grad_chargesfile"] = fout.read()

        # read the charge.dat file when present
        if os.path.exists(os.path.join(calcdir, "efield.dat")):
            with open(os.path.join(calcdir, "efield.dat")) as fout:
                self.dataDict["efieldfile"] = fout.read()

        # process the file with the gradient of the point charges (only when gradient is available)
        if self.dataDict["grad_chargesfile"] is not None and self.dataDict["gradient"] is not None:
            # initialize dictionary to store electric field
            self.dataDict["fullgradcharge"] = {}
            # split lines of the grad_charges file
            linebyline, j = self.dataDict["grad_chargesfile"].splitlines(), 0
            # loop over the electronic states (the number of states should be the same of the atom gradient)
            for istate in range(len(self.dataDict["gradient"])):
                size, statelabel = linebyline[j].split("!")  # first line is the state description
                natomMM = int(size.split()[0])  # extract number of M layer atoms
                j += 1  # increment line counter
                if natomMM == 0:
                    self.dataDict["fullgradcharge"][istate] = None
                else:
                    # get electric field computed for the state (the number of points is the number of M atoms)
                    x, y, z = [], [], []
                    for jatom in range(natomMM):
                        linesplit = linebyline[j].split()
                        x.append(float(linesplit[0])), y.append(float(linesplit[1])), z.append(float(linesplit[2]))
                        j += 1  # increment line counter
                    self.dataDict["fullgradcharge"][istate] = x, y, z

        if self.dataDict["efieldfile"] is not None and self.dataDict["gradient"] is not None:
            # initialize dictionary to store electric field
            self.dataDict["elfield"] = {}
            # split lines of the efield.dat file
            linebyline, j = self.dataDict["efieldfile"].splitlines(), 0
            # loop over the electronic states (the number of states should be the same of the atom gradient)
            for istate in range(len(self.dataDict["gradient"])):
                size, statelabel = linebyline[j].split("!")  # first line is the state description
                natomMM = int(size.split()[0])  # extract number of M layer atoms
                j += 1  # increment line counter
                if natomMM == 0:
                    self.dataDict["elfield"][istate] = None
                else:
                    # get electric field computed for the state (the number of points is the number of M atoms)
                    x, y, z = [], [], []
                    for jatom in range(natomMM):
                        linesplit = linebyline[j].split()
                        x.append(float(linesplit[0])), y.append(float(linesplit[1])), z.append(float(linesplit[2]))
                        j += 1  # increment line counter
                    self.dataDict["elfield"][istate] = x, y, z

        # finish collecting data from .log file
        # now collect some info on the job execution

        # check if QM termination is normal or error, and store excerpt of log file with error message
        self.dataDict["termination"] = 1
        try:
            with open(os.path.join(calcdir, name + ".log")) as flog:
                for line in flog.readlines():
                    if line.strip() == "===> Writing output to file QM.out in SHARC Format":
                        self.dataDict["termination"] = 0
        except IOError:
            self.dataDict["termination"] = 1
        if self.dataDict["termination"] == 1:
            with open(os.path.join(calcdir, name + ".err")) as ferr:
                self.dataDict["errormsg"] = ferr.read()

    # =============================================================================================================

    def __del__(self):
        """Destructor for the gaussianOutput class: not only the memory allocated for the oject attributes
        (the dictionary self.dataDict) needs to be released, also the Gaussian I/O files stored on disk
        can be safely removed... """

        # permanently remove the QM.out file
        if not logwrt.DEBUG_COBRAMM_RUN: os.remove(os.path.join(self.dataDict["dir"], self.dataDict["name"]+".out"))
        if not logwrt.DEBUG_COBRAMM_RUN: os.remove("charge.dat")
        # clean up the log
        del self.log
        # destroy the dictionary that contains the output data
        del self.dataDict

    # =============================================================================================================

    def orbfiles(self):
        """Check their existance and then return a list with the names of the files to be
         saved to store the orbitals to disk. In this case, no file needs to be saved"""

        return []

    # =============================================================================================================

    def restartfile(self):
        """Check their existance and then return a the names of the files to be
         saved to store the orbitals to disk for restart purpose. In this case, there is no QM restart."""

        return None


###################################################################################################################
