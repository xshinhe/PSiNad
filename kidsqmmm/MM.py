#!/usr/bin/env python3
# coding=utf-8

#    COBRAMM 这个文件还没有被使用!!!!
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
import copy  # shallow and deep copy operations

# imports of local objects

from Timer import Timer  # keep timings of the different sections of the code
from amberDriver import AmberInput, AmberDriver, AmberSnapshot, AmberTopology  # Amber wrapper

# imports of local modules

import logwrt  # manages log file output + start/end procedures
import CBF


###################################################################################################################

class MM:
    """Object that generate, store and process data for the MM part of the QM/MM calculation"""

    # variable with the name of the directory where to store the MM output files (logs)
    _MM_DATA_STORAGE = 'MM_data'

    # CLASS variables that controls the namings for the files of a Amber calculation
    _MODELH_TOP = "model-H.top"
    _REAL_TOP = "real.top"

    @Timer("MM setup")
    def __init__(self, cobcom, command, geometry, auxfiles=None):
        """Constructor of the MM class. The class store all the data that is necessary to setup a
        calculation, but the actual input files of the MM calculation are written when
        the class instance is created. """

        # auxfiles should be a dictionary, if None make it an empty dictionary
        if auxfiles is None: auxfiles = {}

        # Function to initialize the _MM_DATA_STORAGE directory: clean up old storage dir and remove
        # files that might remain from previous run of COBRAMM

        # remove existing QM_DATA_STORAGE dir, which has been left by another cobramm execution
        if MM._MM_DATA_STORAGE in os.listdir('.'):
            shutil.rmtree(MM._MM_DATA_STORAGE)
        # create a new directory
        os.mkdir(MM._MM_DATA_STORAGE)

        # depending on the type of MM chosen with the command[71] keyword, create the
        # instances of the appropriate classes that store the information on the input
        # and topology files.

        if command[71] == '0':  # MM calculation by Amber

            # define options for the amber run
            if command[20] == '1':
                self.ncores, self.GPU = 1, True
            elif command[20] == '2':
                self.ncores, self.GPU = command[7], False
            else:
                self.ncores, self.GPU = 1, False

            # the real top file is alway needed
            if MM._REAL_TOP not in auxfiles:
                logwrt.fatalerror('mandatory input file {0} not available!'.format(MM._REAL_TOP))
            self.realTop = AmberTopology(auxfiles[MM._REAL_TOP])

            # extract the real model charge from the topology file
            self.realCharge = self.realTop.charges

            # check whether the modelH top is present when required (when the calculation has an H layer)
            if "H" in geometry.calculationType:
                if MM._MODELH_TOP not in auxfiles:
                    logwrt.fatalerror('mandatory input file {0} not available!'.format(MM._MODELH_TOP))
                self.modelHTop = AmberTopology(auxfiles[MM._MODELH_TOP])
            else:
                self.modelHTop = None

            # extract list of atoms from the topology file
            realAtoms = self.realTop.atoms
            # check that the number of atoms corresponds
            if len(realAtoms) != geometry.atomNum:
                logwrt.fatalerror("wrong number of atoms in the {0} topology".format(MM._REAL_TOP))
            for labelTop, labelGeom in zip(realAtoms, geometry.atomLabel):
                if labelTop[:len(labelGeom)] != labelGeom and labelTop != "IP":
                    logwrt.writewarning("atom type mismatch in the {0} topology ({1} /= {2})".format(
                        MM._REAL_TOP, labelTop, labelGeom))

            # when the modelH topology is available, extract the list of atoms and make additional checks
            if self.modelHTop is not None:
                modelHAtoms = self.modelHTop.atoms
                # check that the number of atoms corresponds
                if len(modelHAtoms) != geometry.NatomQM + geometry.NsubH:
                    logwrt.fatalerror("wrong number of atoms in the {0} topology".format(MM._MODELH_TOP))
                for i in range(geometry.NatomQM):
                    if modelHAtoms[i] != realAtoms[geometry.list_HIGH[i] - 1]:
                        logwrt.writewarning("atom type mismatch in the {0} and {1} topologies ({3} /= {2})".format(
                            MM._REAL_TOP, MM._MODELH_TOP, modelHAtoms[i], realAtoms[geometry.list_HIGH[i] - 1]))

            # extract sander input for the "first" calculation from the command file
            commandInput = CBF.ReadCobramCommand(cobcom, 'sander', 'real-sander.inp')

            # the missing "sander" section is not a problem, since this section controls a calc that can be skipped
            if len(commandInput) == 0:
                self.inputFromCMND = None

            else:  # when the input is present, however, we need to make sure that the relevant options are present
                self.inputFromCMND = AmberInput.readinput("\n".join(commandInput))

                # 1. bellymaks must contain all the low layer atoms, because the low layer is allowed to move
                #    in the microiterative scheme
                try:
                    # extract the current geometry and store it as an AmberSnapshot instance
                    bellymaskstring = self.inputFromCMND["bellymask"]
                    freeatoms = AmberDriver.ambermask2atomlist(self.realTop, AmberSnapshot(geometry.cartesian),
                                                                   bellymaskstring) + list(geometry.list_LOW)
                    freeatoms = sorted(list(set(freeatoms)))  # get unique entries in increasing order
                except KeyError:
                    freeatoms = list(geometry.list_LOW)
                # now get the amber-style string with the list of free atoms
                freestring = AmberDriver.atomlist2amberformat(freeatoms)
                # put bellymask statement in self.sanderInput
                self.inputFromCMND["bellymask"], self.inputFromCMND["ibelly"] = freestring, 1
                # 2. the input should be of the type:  "ntx = 1, irest = 0" to read TXT coordinates from input crd
                self.inputFromCMND["ntx"], self.inputFromCMND["irest"] = 1, 0
                # 3. the restart should be of the type: "ntxo = 2" to get the final coords from the NETCDF restrt file
                self.inputFromCMND["ntxo"] = 2
                # 4. there should be no unit cell nor periodic boundary conditions, "ntb = 0"
                self.inputFromCMND["ntb"] = 0

            # now if a MM calc for a QM/MM is being setup, create topo with 0 charge for the QM part of the model
            if "H" in geometry.calculationType:

                # create the list of MM charges with zeroed QM charges
                realZeroedCharges = copy.deepcopy(self.realCharge)
                for i in geometry.list_QM:
                    realZeroedCharges[i-1] = 0.0
                # substitute the list in the real topology, and save the new topology as instance attribute
                self.realTopZeroQMCh = self.realTop.substitutecharges(realZeroedCharges)

                # create the list of modelH zero charge
                modelHZeroedCharges = [0.0 for _ in range(geometry.NatomQM + geometry.NsubH)]
                # substitute the list in modelH topology, and save the new topology as instance attribute
                self.modelHTopZeroQMCh = self.modelHTop.substitutecharges(modelHZeroedCharges)

            # create generic input files for a single point calculation
            self.SPInput = AmberInput(minimize=False, nrSteps=1, nprntsteps=1, deltaT=0.0, cutoff=float(command[41]),
                                      usevelocity=False, usePBC=False)

            # define whether the "first" AMBER calculation should be done or not
            if self.inputFromCMND is None or command[1] == 'freqxg':
                self.doFirstStep = False
            else:
                self.doFirstStep = True

        else:
            logwrt.fatalerror("MM software (command 71) '{0}' is not defined. ".format(command[71]))

        # define instance attributes that will be populated later when running the actual calculation
        self.firstOutput = None  # output results of the preliminary full MM calculation
        self.realOutput = None  # output results of the SP MM calculation for the full system
        self.modelHOutput = None  # output results of the SP MM calculation for the modelH system

    # =============================================================================================================

    @Timer("MM run")
    def run(self, geometry):

        # the layer scheme does not contain any MM atom, so nothing needs to be done
        if "M" not in geometry.calculationType and "L" not in geometry.calculationType:

            # leave the output instance variables initialiazed to none and go on
            self.firstOutput, self.realOutput, self.modelHOutput = None, None, None

        # the scheme contains MM atoms, we need to do a sequence of MM calculations
        else:

            # depending on the type of MM chosen with the command[71] keyword,
            # run the MM calculations that are needed to construct the QM/MM gradient.

            if isinstance(self.SPInput, AmberInput):  # MM calculation by Amber
                logwrt.writelog('Starting MM calculation using AMBER\n')

                # extract the current geometry and store it as an AmberSnapshot instance
                realCoords = AmberSnapshot(geometry.cartesian)

                # Optionally, a first calculation is done prior to the computation
                # of the gradient to perform some low layer displacements (e.g. in the microiterative scheme)
                if self.doFirstStep:
                    logwrt.writelog("\nPerforming a preliminary full calculation of the entire system "
                                    "to move the low layer ...")
                    self.firstOutput = AmberDriver.run(self.inputFromCMND, self.realTop, realCoords, GPU=self.GPU,
                                                           nCores=self.ncores, calcDir=MM._availablePathForMM())
                    logwrt.writelog(' done!\n')
                    # extract new coordinates obtained at the end of the calculation,
                    realCoords = self.firstOutput.snapshot(-1)
                    # and update geometry with the new coordinates
                    geometry.cartesian = realCoords.coords

                # the layer scheme does not contain any QM atom, so only one MM for the full system needs to be done
                if "H" not in geometry.calculationType:

                    # amber full MM single point on real
                    logwrt.writelog("\nPerforming a single point calculation and collecting energy and gradient ...")
                    self.realOutput = AmberDriver.run(self.SPInput, self.realTop, realCoords, GPU=self.GPU,
                                                          nCores=self.ncores, calcDir=MM._availablePathForMM())
                    logwrt.writelog(' done!\n')
                    logwrt.writelog('MM energy is {0} Hartree\n'.format(self.realOutput.energy[0]))

                # the layer scheme contains both QM and MM atoms, so we need to compute the gradient and
                # the energy with the subtractive scheme
                else:

                    # amber MM single point on full system with ZERO charges on model-H
                    logwrt.writelog("\nPerforming a single point calc. of the entire system with ZERO charge on the \n"
                                    "High layer, and collecting energy and gradient ...")
                    self.realOutput = AmberDriver.run(self.SPInput, self.realTopZeroQMCh, realCoords, GPU=self.GPU,
                                                          nCores=self.ncores, calcDir=MM._availablePathForMM())
                    logwrt.writelog(' done!\n')
                    logwrt.writelog('MM energy is {0} Hartree\n'.format(self.realOutput.energy[0]))

                    # extract coordinates of modelH and store them in an instance of AmberSnapshot
                    modelHCoords = AmberSnapshot(geometry.modelH)
                    # amber MM single point on model-H system with ZERO charges
                    logwrt.writelog("\nPerforming a single point calc. of the High layer, with H-saturated bonds and \n"
                                    "ZERO charge, and collecting energy and gradient ...")
                    self.modelHOutput = AmberDriver.run(self.SPInput, self.modelHTopZeroQMCh, modelHCoords,
                                                            calcDir=MM._availablePathForMM(), GPU=self.GPU,
                                                            nCores=self.ncores)
                    logwrt.writelog(' done!\n')
                    logwrt.writelog('MM energy is {0} Hartree\n'.format(self.modelHOutput.energy[0]))

        logwrt.writelog('\n')

    # =============================================================================================================

    def __getattr__(self, item):
        """Catches non-standard class attributes, by using the "item" key to access the
        some selected results that are necessary to construct the QMMM energy
        and gradient: the energy and the gradient of the single point calculations
        for the modelH system and the real system (possibly with zeroed QM charges)."""

        if item in ["E_real", "E_modelH", "Grad_real", "Grad_modelH"]:
            try:
                if item == "E_real":
                    return self.realOutput.energy[0]
                elif item == "E_modelH":
                    return self.modelHOutput.energy[0]
                elif item == "Grad_real":
                    return self.realOutput.gradient[0]
                elif item == "Grad_modelH":
                    return self.modelHOutput.gradient[0]
            except AttributeError:
                return None

    # =============================================================================================================

    def __del__(self):
        """Destructor for the QM class: destroys the input data (defined on instance construction)
        and the output data, when defined"""

        # destroy all attributes
        del self.GPU, self.ncores, self.doFirstStep
        del self.SPInput, self.inputFromCMND
        del self.firstOutput
        del self.realCharge
        del self.realOutput, self.realTop, self.realTopZeroQMCh
        del self.modelHOutput, self.modelHTop, self.modelHTopZeroQMCh

    # =============================================================================================================

    def archiveStep(self, step):
        """ save output defined in the instance of MM. Step is an integer number that is used to label
        the information that is saved, and tipically is the step number. """

        # file where to store QM results
        firstAll = os.path.join(MM._MM_DATA_STORAGE, 'first_ALL.log')
        realAll = os.path.join(MM._MM_DATA_STORAGE, 'real_ALL.log')
        modelHAll = os.path.join(MM._MM_DATA_STORAGE, 'modelH_ALL.log')

        # decorator of the log file to distinguish the steps
        stepseparation = '=' * 80 + '\n' +\
                         "Start MM {0} calculation output of STEP : {1}\n" +\
                         ' || ' * 20 + '\n' + ' \/ ' * 20 + '\n' +\
                         "{2}" + \
                         ' /\ ' * 20 + '\n' + ' || ' * 20 + '\n' + \
                         "End MM {0} calculation output of STEP : {1}\n" + \
                         '=' * 80 + '\n'

        for output, allfile, name in zip([self.firstOutput, self.modelHOutput, self.realOutput],
                                         [firstAll, modelHAll, realAll],
                                         ["full", "model-H SP", "real SP"]):
            # print only if the output is present
            if output is not None:
                # write the content of the AMBER log file to the ALL file, decorated with the step number
                with open(allfile, 'a') as f:
                    f.write(stepseparation.format(name, step, output.logtext))

    # =============================================================================================================

    @staticmethod
    def _availablePathForMM():
        """Returns a path corresponding to an empty directory where the MM calculation
        can be run by the runMM function. The directory is created by this function."""

        mmID = 0
        while True:
            # increment an ID variable that labels the mm directory
            mmID += 1
            # define the name of the directory
            mmdir = "mmCalc{0:05d}".format(mmID)
            # check whether the directory exists, if it does not returns it to caller function
            if not os.path.isdir(mmdir): break

        # create the directory
        return mmdir


# #####################################################################################################

if __name__ == '__main__':
    # test run of the MM class

    # use files from example calculations
    realtopofile = "../examples/TURBOMOLE/QMMM_OPTXG_DSCF/real.top"
    modelHtopofile = "../examples/TURBOMOLE/QMMM_OPTXG_DSCF/model-H.top"
    crdfile = "../examples/TURBOMOLE/QMMM_OPTXG_DSCF/real.crd"
    reallayersfile = "../examples/TURBOMOLE/QMMM_OPTXG_DSCF/real_layers.xyz"
    commandfile = "../examples/TURBOMOLE/QMMM_OPTXG_DSCF/cobram.command"

    # initialize geometry
    from Layers import Layers
    with open(reallayersfile, "r") as f:
        testgeometry = Layers.from_real_layers_xyz(f.read())

    # initialize the list with the COBRAMM options
    testcommand = CBF.makehard()
    # get the content of the cobramm command file
    testcobcom = CBF.getCobrammCommand(commandfile)
    # merge command soft and hard
    key_soft = CBF.ReadCobramCommand(testcobcom, 'keyword', 'keywords')
    testcommand = CBF.key2hard(key_soft, testcommand)
    command_soft = CBF.ReadCobramCommand(testcobcom, 'command', 'commands')
    testcommand = CBF.soft2hard(command_soft, testcommand)

    # read topology file
    topofiles = {}
    with open(realtopofile) as f:
        topofiles["real.top"] = f.read()
    with open(modelHtopofile) as f:
        topofiles["model-H.top"] = f.read()

    # try to read charges from real.top
    MMtest = MM(testcobcom, testcommand, testgeometry, topofiles)
    MMtest.run(testgeometry)

    # print results for QM/MM energy and gradient
    print("E real = {0} Hartree".format(MMtest.E_real))
    print("grad real = {0} a.u.".format(MMtest.Grad_real))
    print("E modelH = {0} Hartree".format(MMtest.E_modelH))
    print("grad modelH = {0} a.u.".format(MMtest.Grad_modelH))

    # save log files
    MMtest.archiveStep(0)
