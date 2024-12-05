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
# import shelve
import math  # mathematical functions
import multiprocessing as mp  # run child processes in parallel

# math libraries
import numpy as np  # numpy library for scientific computation

# imports of local objects

from timer import Timer  # keep timings of the different sections of the code
from QMOutput import QMOutput  # template class for the output of a QM calculation
from gaussianCalculator import GaussianInput, GaussianOutput  # objects to handle Gaussian Input and Output
from molcasCalculator import MolcasInput, MolcasOutput
from sharcQMCalculator import SharcInterfaceInput, SharcInterfaceOutput  # objects to handle QM with SHARC QM Interfaces
from orbitals import Orbitals  # object to parse and store orbital information from QM output files

# imports of local modules

import logwrt  # manages log file output + start/end procedures
import CBF
import cobrammenv  # functions to handle third-party software environment
import constants  # values of physical constants and conversion factors

###################################################################################################################


class QM:
    """Object that generate, store and process data for the QM part of the QM/MM calculation"""

    # CLASS variables that controls the preparation of the QM_DATA_STORAGE directory
    # variable with the name of the directory where to store the QM output files (orbitals and logs)
    _QM_DATA_STORAGE = 'QM_data'
    # initialization status of the QM data storage
    _QM_DATA_SETUP = False

    @Timer("QM input")
    def __init__(self, cobcom, command, geometry, charges, step, restart=None, displacement=None,
                 setState=None, forceState=None, MS_second_run=False, SS_second_run=False, gradOnly=False, NAC=False): #--> F: added step and tstep attribute
        """Constructor of the QM class. The input of the QM calculation is defined when
        the class instance is created. """

        # define a string attribute to store the log related to QM input and output, to be printed when convenient
        self.log = ""

        # when the dictionary displacement is given, use its arguments to set up a displaced modelH geometry
        # otherwise use the regular modelH defined in geometry
        if displacement is not None:
            modelHgeom = geometry.modelHdisplace(displacement["iAt"], displacement["iCoord"], displacement["iDir"],
                                                 displacement["length"], displacement["linkfollow"])
        else:
            modelHgeom = geometry.modelH

        self.tstep = None
        # following attributes are needed to activate special actions in case of PT2 calculations
        # for QM/MM X/R/MS-PT2 calculations different from SP, a second run is needed on the roatted CASSCF states to retrieve electric field
        self.PT2_QMMM = False
        #in case of SS-CASPT2 calculations requiring gradients, the &CASPT2(GRDT) and &ALASKA routines must be executed in a secondary calculation
        #(after the root order is retrieved from the first run)
        self.is_SSPT2_first_run = False
        # for PT2 nonadiabatic dynamics with THS scheme, the PT2 WF must be saved for overlap calculation at the next step
        self.save_PT2_WF = False

        # depending on the type of QM chosen with the command[51] keyword, create the
        # instance of the appropriate class that stores the information on the input file

        if command[51] == '1':  # QM calculation by Gaussian

            # extract sections of the command file that are related to the QM code
            keys = CBF.ReadCobramCommand(cobcom, 'gaussian', '')  # gaussian keywords
            gen = CBF.ReadCobramCommand(cobcom, 'gen', '')  # gaussian basis set
            gaussadd = CBF.ReadCobramCommand(cobcom, 'gaussadd', '')  # gaussian additional parts
            gaussweights = CBF.ReadCobramCommand(cobcom, 'gaussweights', '')  # gaussian weights of the CASSCF

            # create dictionary with additional options for the gaussian calculation
            optdict = {}
            # background charges, only for a calculation with atoms in the H or M layers
            if geometry.NatomMM > 0: optdict["charges"] = geometry.pod, charges.CRG_emb
            # when this is not a single point calculation, put option to compute forces
            if command[60] != '1': #-->F: is '1' and 'sp' really equivalent everywhere in the code? CHECK!!!
                optdict["forces"] = True
                optdict["SP"] = False
                # when this is not a single step and there are M atoms, electric field computation is required
                if geometry.NatomM > 0: optdict["field"] = geometry.MEDIUM, charges.CRG_MEDIUM
            else:
                optdict["SP"] = True
            # increase verbosity in gaussianALL.log, do not delete AO, MO and CI deffinitions
            if command[2] == '2':
                optdict["verbose"] = True
            else:
                optdict["verbose"] = False
            # when the restart input value contains the path of an available chk file, set restart to true
            if restart and os.path.exists(restart):
                self.log += '{0} file found! \nIt will be used as chk file ' \
                            'for initial WF guess.\n\n'.format(restart)
                optdict["restart"] = restart
            # option to disable nosymm control on input file
            if command[192] != '0': optdict["suppress_nosymm_error"] = True

            # when TDcouplings are requested, require tamm-damkov and force printing of all excitations
            if command[14] == '1':
                optdict["tdcouplings"] = True

            if command[103] == '0':
                optdict["suppress_TDDFT"] = True

            # when a specific state is requested, define the "nstate" option in optdic
            if forceState is not None:
                optdict["nstate"] = int(forceState)
                optdict["forcestate"] = True
            # this is an alternative option to set the initial state, that gives error when there is a conflict with
            # what is defined in the QM input text
            elif setState is not None:
                optdict["nstate"] = int(setState)
                optdict["forcestate"] = False

            # when a CI optimization is requested, set an option for the calculation of both gradients (upper and lower state)
            # CI optimization at (TD)DFT level is only available with branching plane based on gradient mean and gradient difference (no NACs)
            if command[1] == 'ci':
                if command[19] == '0':
                    self.inputData = None
                    self.outputData = None
                    logwrt.fatalerror('CI optimization with gaussian is only available with gmean branching plane! Please set command 19=1 (keyword BranchPlane=gmean)\n')
                else:
                    optdict["CIopt"] = True
                    optdict["upper_state"] = int(command[13])
            else:
                optdict["CIopt"] = False

            # save instance of GaussianInput as the input file for the QM calculation
            self.inputData = GaussianInput(keys, gen, gaussadd, gaussweights, modelHgeom,
                                           geometry.getAtomLabels("modelH"), optdict)

        elif command[51] == '6':  # QM calulation by Molcas
            # extract sections of the command file that are related to the QM code
            keys = CBF.ReadCobramCommand(cobcom, 'molcas', '')  # molcas keywords
            basisset = CBF.ReadCobramCommand(cobcom, 'basisset', '')  # user-defined molcas basis set (optional)
            sewardkey = CBF.ReadCobramCommand(cobcom, 'seward', '')  # additional seward keywords added through !sewars section of cobram.command

            # create dictionary with additional options for the gaussian calculation
            optdict = {}
            # background charges, only for a calculation with atoms in the H or M layers
            if geometry.NatomMM > 0: optdict["charges"] = geometry.pod, charges.CRG_emb
            # when this is not a single point calculation, put option to compute forces
            if command[60] != '1': #-->F: is '1' and 'sp' really equivalent everywhere in the code? CHECK!!!
                optdict["forces"] = True
                optdict["SP"] = False
                # when this is not a single step and there are M atoms, electric field computation is required
                if geometry.NatomM > 0: optdict["field"] = geometry.MEDIUM, charges.CRG_MEDIUM
            else:
                optdict["SP"] = True

            # when the restart input value contains the path of an available restart file, set restart to true
            if restart and os.path.exists(restart):
                self.log += '{0} file found! \nIt will be used ' \
                            'for initial WF guess.\n\n'.format(os.path.basename(restart))
                optdict["restart"] = restart

            #basis set from !command (or !keyword) section, to use for molcas (unless non standard basis set is specified through !bassset section)
            optdict["basis"] = command[197]

            # type of basis functions
            if command[194] == '1':
                optdict["cartesian"] = True
            else:
                optdict["cartesian"] = False

            # use Cholesky decomposition
            if command[196] == '1':
                optdict["ricd"] = True
                optdict["cdth"] = command[195]

            # pass command[1] (i.e. type of calculatio ) to MolcasInput to build appropriate input file
            optdict['command1'] = command[1]
            #in case of MDV, pass hopping scheme to MolcasInput for appropriate input construction
            if command[1] == 'mdv':
                if int(command[85]) == 1:
                    # Tully with THS time-derivative coupling
                    if int(command[14]) == 1:
                        optdict['mdv_type'] = 1
                    # Tully with NACs 
                    else:
                        optdict['mdv_type'] = 2
                # energy gap-based hop
                elif int(command[85]) == 2:
                    optdict['mdv_type'] = 3
                else:
                    optdict['mdv_type'] = 0


            # for CI optimizations: define type of branching plane (NAC-based or gradients-only)
            if command[19] == '1':
                optdict["newBP"] = True

            # when a specific state is requested, define the "nstate" option in optdic
            if forceState is not None:
                optdict["nstate"] = int(forceState)
                optdict["forcestate"] = True
            # this is an alternative option to set the initial state, that gives error when there is a conflict with
            # what is defined in the QM input text
            elif setState is not None:
                optdict["nstate"] = int(setState)
                optdict["forcestate"] = False
            # True if this is the second QM run for MS-XMS-RMS PT2 QM/MM (to retrieve electric field)
            if MS_second_run or SS_second_run:
                optdict["ROST"] = True
            if command[81] != '0':
                optdict["highest_NAC_state"] = int(command[81])
            if gradOnly:
                optdict["gradOnly"] = True
            if NAC:
                #argument NAC (when present) is a list of two elements, identifying the old and new state of a hopping event
                #these are passed to MolcasInput through optdict so that NAC is calculated after hop and we can rescale velocity
                optdict["NAC_for_vel_rescale"] = NAC

            # save instance of MolcasInput as the input file for the QM calculation
            self.inputData = MolcasInput(keys, basisset, sewardkey, modelHgeom,
                                           geometry.getAtomLabels("modelH"), optdict, step)
            # set an attributes to remember if a second PT2 run is needed
            # (in case of QM/MM calculations different from SP a second run is needed to retrieve electric field)
            if (self.inputData.calctype in ['MSPT2', 'XMSPT2', 'RMSPT2']) and not optdict["SP"] and (geometry.NatomMM > 0) and not MS_second_run:
                self.PT2_QMMM = True
            # (in case of SS-CASPT2 calculations different from SP a second run is needed to retrieve gradient of the correct state, after checking possible root flipping)
            if (self.inputData.calctype == 'SSPT2') and not optdict["SP"] and not SS_second_run:
                self.is_SSPT2_first_run = True
            # in case of PT2 nonadiabatic dynamics (command 1 + command 85) with THS scheme (command 14)
            if (self.inputData.calctype in ['MSPT2', 'XMSPT2', 'RMSPT2']) and (command[1] == 'mdv') and (int(command[85]) > 0) and (command[14] == '1') and not MS_second_run:
                self.save_PT2_WF = True
            elif self.inputData.calctype == 'SSPT2' and (command[1] == 'mdv') and (int(command[85]) > 0) and (command[14] == '1') and not self.is_SSPT2_first_run:
                self.save_PT2_WF = True

        elif command[51] == '7':  # QM calulation by MOLPRO
            # Molpro is not yet supported by the QM class
            self.inputData = None

        elif command[51] == '11':  # QM calulation with SHARC INTERFACES

            qmcode = command[110]

            if forceState is not None:
                logwrt.fatalerror("cannot change electronic state in QM calculation with SHARC interface")

            # create dictionary with additional options for the gaussian calculation
            optdict = {}
            # background charges and electric field computation, only for a calc with atoms in the H or M layers
            if geometry.NatomMM > 0:
                optdict["charges"] = geometry.pod, charges.CRG_emb, geometry.pod_mobileatoms
                if qmcode.upper() in ["MOLCAS"]:
                    optdict["field"] = geometry.MEDIUM, charges.CRG_MEDIUM

            self.inputData = SharcInterfaceInput(qmcode, modelHgeom, geometry.getAtomLabels("modelH"), optdict)

        else:
            # other type of calculations are not supported by this QM class
            self.inputData = None

        # initialize an empty outputData container, that will be used to store the result of the QM calculation
        self.outputData = QMOutput()
        # initialize another attribute to store the orbitals, for the QM for which orbital parsing is available
        self.orbitals = None

    # =============================================================================================================
    
    def __del__(self):
        """Destructor for the QM class: destroys the input data (defined on instance construction)
        and the output data, when defined"""

        # destroy input data
        del self.inputData
        # when available, destroy output data
        if self.outputData: del self.outputData

    # =============================================================================================================

    def setOutputData(self, energy, gradient, state, charges, dipole, selfenergy, gradcharges, NAC):
        """Workaround to define the QM instance from QM results that have been already obtained
        elsewhere, in order to use the QM class as a simple data storage object in cases where the QM
        code is not fully implemented"""

        # store list of the electronic energies and the corresponding state ID (GS = 1, ES = 2,3,4... )
        # data is stored as a dictionary {state_number : energy_of_the_state}
        self.outputData.set("energy", {i: en for i, en in enumerate(energy)})
        # store the gradient of the electronic state of interest
        self.outputData.set("gradient", {i: grad for i, grad in zip(state, gradient)})
        # store information on the state that has been optimized in the QM calculation
        self.outputData.set("optstate", state[0])
        # store the electrostatic output information
        self.outputData.set("charges", charges)
        self.outputData.set("dipole", dipole)
        # store the energy/gradient components of the electrostatic embedding
        self.outputData.set("selfenergy", selfenergy)
        self.outputData.set("gradcharges", {i: grad for i, grad in zip(state, gradcharges)})
        self.outputData.set("nac", {i: {j: NAC[i][j] for j in range(len(NAC[i]))} for i in range(len(NAC))})

    # =============================================================================================================

    def setTstep(self, tstep):
        self.tstep = tstep

    # =============================================================================================================

    def setTDC(self, tdc):
        self.outputData.set("tdc", tdc)

    # =============================================================================================================

    # TODO: extend the behaviour of the __getattribute__ to handle in more generic and
    #  rational way the code to get access to the output data stored in self.outputData...
    #  the following methods energy(), gradient(), gradcharges(), energydict(), ... dipole()
    #  should be revisited.
    #  - For most of the attributes, one can simply use the name of the
    #  attribute as argument of a self.outputData.get() call
    #  - In the case of energy/gradient/gradcharges, it should be decided what is
    #  the default behaviour of the attribute (full dictionary vs optimized element)

    def energy(self, nstate=None):
        """Return a list containing the electronic state energy (or energies) computed with the QM calculation """

        # get the dictionary that contains the QM energy of the states
        energydict = self.outputData.get("energy")

        try:  # exception handler, return None if that state requested from input is not available
            # when nstate is None, no state is requested... return the full list in order of increasing energy
            if nstate is None:
                return [energydict[n] for n in sorted(energydict.keys())]
            # if nstate = "opt", then return the energy of the state that has been optimized in QM calculation
            elif nstate == "opt":
                return energydict[self.outputData.get("optstate")]
            # otherwise return the energy of the state requested as input
            else:
                return energydict[nstate]
        except (KeyError, TypeError):
            return None

    # =============================================================================================================

    def gradient(self, nstate=None):
        """Return the gradient on the QM part (or gradients) computed with the QM calculation """

        # get the dictionary that contains the QM gradient of the states computed with QM
        graddict = self.outputData.get("gradient")

        try:  # exception handler, return None if that state requested from input is not available
            # when nstate is None, no state is requested... return the full list in order of increasing energy
            if nstate is None:
                return [graddict[n] for n in sorted(graddict.keys())]
            # if nstate = "opt", then return the gradoemt of the state that has been optimized in QM calculation
            elif nstate == "opt":
                return graddict[self.outputData.get("optstate")]
            # otherwise return the gradient of the state requested as input
            else:
                return graddict[nstate]
        except (KeyError, TypeError):
            return None

    # =============================================================================================================

    def nac(self, nstate1=None, nstate2=None):
        """Return the NAC on the QM part computed with the QM calculation """

        # get the dictionary that contains the QM nac of the states computed with QM
        nacdict = self.outputData.get("nac")

        try:  # exception handler, return None if that state requested from input is not available
            # when nstate is None, no state is requested... return the full list in order of increasing energy
            if nstate1 is None or nstate2 is None:
                logwrt.fatalerror("Looks like this calculation needs a NAC that is not available...please check input\n")
            # otherwise return the nac of the state requested as input
            else:
                return nacdict[nstate1][nstate2]
        except (KeyError, TypeError):
            return None

    # =============================================================================================================

    def gradcharges(self, nstate=None):
        """Return the gradient of the QM-MM electrostatic energy wrt the position of the background charges """

        # get the dictionary that contains the QM gradient of the states computed with QM
        graddict = self.outputData.get("gradcharges")

        try:  # exception handler, return None if that state requested from input is not available
            # when nstate is None, no state is requested... return the full list in order of increasing energy
            if nstate is None:
                return [graddict[n] for n in sorted(graddict.keys())]
            # if nstate = "opt", then return the energy of the state that has been optimized in QM calculation
            elif nstate == "opt":
                return graddict[self.outputData.get("optstate")]
            # otherwise return the energy of the state requested as input
            else:
                return graddict[nstate]
        except (KeyError, TypeError):
            return None

    # =============================================================================================================

    @property
    def energydict(self):
        """Return the dictionary with the energy results of the QM calculation """
        return self.outputData.get("energy")

    # =============================================================================================================

    @property
    def gradientdict(self):
        """Return the dictionary with the gradient results of the QM calculation """
        return self.outputData.get("gradient")

    # =============================================================================================================

    @property
    def charges(self):
        """Return a list containing the atomic charges computed with the QM calculation """
        return self.outputData.get("charges")

    # =============================================================================================================

    @property
    def selfenergy(self):
        """Return a list containing the self-energy of the background charge distribution """
        return self.outputData.get("selfenergy")

    # =============================================================================================================

    @property
    def dipole(self):
        """Return a list containing the x,y,z components + norm of the dipole moment computed with the QM calc """
        return self.outputData.get("dipole")

    # =============================================================================================================

    @staticmethod
    def _prepareQMDataStorage():
        """Function to initialize the QM_DATA_STORAGE directory, to be called at the
        before running any QM calculation. The function cleans up the storage dir and remove
        files that might remain from previous run of COBRAMM. """

        if not QM._QM_DATA_SETUP:  # the data storage is not set up yet...

            # remove existing QM_DATA_STORAGE dir, which has been left by another cobramm execution
            if QM._QM_DATA_STORAGE in os.listdir('.'):
                shutil.rmtree(QM._QM_DATA_STORAGE)

            # create a new directory
            os.mkdir(QM._QM_DATA_STORAGE)

            # set initialization status variable to True
            QM._QM_DATA_SETUP = True

    # =============================================================================================================

    def archiveStep(self, copyLog, copyOrb, step):
        """ save output and orbital files for the calculation defined in the instance of QM:
        the log files is appended to the QMall.log file, and the orbitals files are copied
        depending on the type of QM used. copyLog and copyOrb are logical flags: when True,
        log and orbital files are saved, respectively. Step is an integer number that is used to label
        the information that is saved, and tipically is the step number. """

        # initialize data storage
        QM._prepareQMDataStorage()

        # file where to store QM results
        allName = os.path.join(QM._QM_DATA_STORAGE, 'qmALL.log')

        # check if the QM calculation output exists, otherwise print a warning to screen
        if self.outputData is None:
            logwrt.writewarning("QM output is not available, log and orbital files will not be saved")
            return

        if copyLog:
            # write the content of the QM single point calc in the ALL file, decorated with
            # comments that highligh the step number of the QM single point
            with open(allName, 'a') as qmall:
                qmall.write('=' * 80 + '\n')
                qmall.write("Start QM calculation output of STEP : " + str(step) + "\n")
                qmall.write(' || ' * 20 + '\n' + ' \/ ' * 20 + '\n')
                qmall.write(self.outputData.get("outfile"))
                qmall.write(' /\ ' * 20 + '\n' + ' || ' * 20 + '\n')
                qmall.write("End QM calculation output of STEP : " + str(step) + "\n")
                qmall.write('=' * 80 + '\n' + '\n')

        # copy orbitals  when required
        if copyOrb:
            # loop over the files that one needs to copy to save the orbitals
            for srcf in self.outputData.orbfiles():
                # process name to extract path, name and extension of the file
                path, filename = os.path.split(srcf)
                basename, ext = os.path.splitext(filename)
                # define the name of the destination file
                tarf = os.path.join(QM._QM_DATA_STORAGE, "{0}_{1}".format(basename, step) + ext)
                # save orbital file with a name that contains the step number, and then gzip it
                shutil.copy(srcf, tarf)
                CBF.GZIP(tarf)

    # =============================================================================================================

    def saveRestartFile(self):
        """Save a copy of the file for orbital restart that is appropriate for the type
        of QM calculation that has been done, and return to caller the full path of the
        file."""

        # name of the file for restart
        srcf = self.outputData.restartfile()

        # if the restart file is available
        if srcf is not None:

            # process name to extract path, name and extension of the file
            path, filename = os.path.split(srcf)
            basename, ext = os.path.splitext(filename)

            # define the name of the destination file
            tarf = os.path.join(os.getcwd(), "restart" + ext)
            # save orbital file with the name defined above and then return the path
            shutil.copy(srcf, tarf)
            return tarf

        # if the file is not available, do nothing and return none
        else:
            return None

    # =============================================================================================================

    @staticmethod
    def _availablePathForQM():
        """Returns a path corresponding to an empty directory where the QM calculation
        can be run by the runQM function. The directory is created by this function."""

        qmID = 0
        while True:
            # increment an ID variable that labels the qm directory
            qmID += 1
            # define the name of the directory
            qmdir = "qmCalc{0:05d}".format(qmID)
            # check whether the directory exists, if it does not returns it to caller function
            if not os.path.isdir(qmdir): break

        # create the directory
        os.mkdir(qmdir)
        return qmdir

    # =============================================================================================================

    @staticmethod
    def _lauchMolcasQM(rundir, filename, qmexe):
        """ This simple function executes the command for running Gaussian QM
         It is needed to store the functions to execute for a parallel execution """

        # store starting dir
        startdir = os.getcwd()
        # move to the directory of the calculation
        os.chdir(rundir)

        # environmental variables for MOLCAS execution
        os.putenv('Project', 'molcas')
        os.putenv('WorkDir', os.getcwd())
        os.putenv('MOLCAS_OUTPUT', os.getcwd())

        # run qm calculation
        finp = filename + ".com"
        with open(filename + ".log", "w") as fout:
                with open(filename + ".err", "w") as ferr:
                    subprocess.call(["nice", "-1", qmexe, finp], stdout=fout, stderr=ferr)

        # move to the directory of the calculation
        os.chdir(startdir)

    # =============================================================================================================

    @staticmethod
    def _lauchGaussianQM(rundir, filename, qmexe):
        """ This simple function executes the command for running Gaussian QM
         It is needed to store the functions to execute for a parallel execution """

        # store starting dir
        startdir = os.getcwd()
        # move to the directory of the calculation
        os.chdir(rundir)

        # run qm calculation
        with open(filename + ".com", "r") as finp:
            with open(filename + ".log", "w") as fout:
                with open(filename + ".err", "w") as ferr:
                    subprocess.call(["nice", "-1", qmexe], stdin=finp, stdout=fout, stderr=ferr)

        # move to the directory of the calculation
        os.chdir(startdir)

    # =============================================================================================================

    @staticmethod
    @Timer("QM run")
    def runQM(QMcalculations, memory="500MB", nprocQM=1, nprocParall=1):
        """Given a list of input files, defined as instances of the QM class, run the QM program
        and extract the relevant results from the output files, that are then saved in the QM instances
        themselves and can be later read as attributes of the instances"""

        # make argument a list
        if type(QMcalculations) is not list: QMcalculations = [QMcalculations]

        # for SHARC, the parallel calculation is not possible
        if np.any([isinstance(qm.inputData, SharcInterfaceInput) for qm in QMcalculations]):
            logwrt.writewarning("parallel execution is not available with SHARC, setting nproc = 1")
            nprocQM = 1

        # if the QM calculations are run in parallel, initialize Pool object from the multiprocessing library
        if nprocParall > 1 and len(QMcalculations) > 1:
            pool = mp.Pool(processes=nprocParall)
        else:
            pool = None

        # store starting dir
        startdir = os.getcwd()

        # lists to store the names of the directory and of the file for the calculations
        rundirList, filenameList = [], []

        # in the standard case, run the calculations serially
        for qm in QMcalculations:

            # input should be run with gaussian!
            if isinstance(qm.inputData, GaussianInput):

                # check if there is a working environment for gaussian
                checkResults = cobrammenv.checkGaussianQMEnv()
                if not checkResults[0]:
                    # if not, stop execution with error message
                    logwrt.fatalerror(checkResults[1])

                # name of the executable and of the I/O files
                qmexe, filename = os.environ["GAUSSIAN_EXE_QM"], "gaussian-QM"
                # set the name where the qm calculation will be run, store it, and move there
                rundir = QM._availablePathForQM()
                # store the names of the dir and files for the current calculation
                rundirList.append(rundir), filenameList.append(filename)

                # write the input text to file
                with open(os.path.join(rundir, filename + ".com"), "w") as finp:
                    finp.write(qm.inputData.fileText(filename, memory=memory, nproc=nprocQM, gversion=qmexe))

                # copy restart file to the working directory, to use the file for restart
                if qm.inputData.otheropt["restart"] is not None:
                    shutil.copy(qm.inputData.otheropt["restart"], os.path.join(rundir, filename + ".chk"))

                # when parallel calculation, use apply_async method to run asincronously
                if pool is not None:
                    pool.apply_async(QM._lauchGaussianQM, args=(rundir, filename, qmexe))
                else:  # when serial calculation, run normally
                    QM._lauchGaussianQM(rundir, filename, qmexe)

            # input is of type SHARC INTERFACE
            elif isinstance(qm.inputData, SharcInterfaceInput):

                # check SHARC environment
                # TODO: add here some checks to control whether SHARC has been properly set

                # name of the executable and of the I/O files
                qmexe, filename = os.path.join(os.environ["SHARC"], qm.inputData.interface), "QM"
                # the SHARC INTERFACE should be run from the parent directory
                rundir = '.'
                os.chdir(rundir)
                # store the names of the dir and files for the current calculation
                rundirList.append(rundir), filenameList.append(filename)

                # write the input text to file
                with open(filename + ".in", "w") as finp:
                    finp.write(qm.inputData.fileText())
                # write the input charge info to file
                with open("charge.dat", "w") as fcharge:
                    fcharge.write(qm.inputData.chargeText())

                # run qm calculation
                with open(filename + ".log", "w") as fout:
                    with open(filename + ".err", "w") as ferr:
                        subprocess.call(["nice", "-1", qmexe, filename + ".in"], stdout=fout, stderr=ferr)

                # move back to starting directory
                os.chdir(startdir)

            elif isinstance(qm.inputData, MolcasInput):

                # check if there is a working environment for Molcas
                checkResults = cobrammenv.checkMolcasEnv()
                if not checkResults[0]:
                    # if not, stop execution with error message
                    logwrt.fatalerror(checkResults[1])

                # name of the executable and of the I/O files
                qmexe, filename = os.environ["MOLCAS_SCRIPT"], "molcas"
                # set the name where the qm calculation will be run, store it, and move there
                rundir = QM._availablePathForQM()
                # store the names of the dir and files for the current calculation
                rundirList.append(rundir), filenameList.append(filename)

                # write the input text to file
                with open(os.path.join(rundir, filename + ".com"), "w") as finp:
                    finp.write(qm.inputData.fileText(filename, memory=memory, nproc=nprocQM))

                # copy restart file to the working directory, to use the file for restart
                if qm.inputData.otheropt["restart"] is not None:
                    if qm.inputData.otheropt["restart"] in ("INPORB", "molcas.RasOrb"):
                        shutil.copy(qm.inputData.otheropt["restart"], os.path.join(rundir, "INPORB"))
                    #elif qm.inputData.otheropt["restart"] in ("molcas.JobIph", "restart.JobIph"):
                    else: #restart is not of RasOrb type, so it is assumed to be JobIph type
                        shutil.copy(qm.inputData.otheropt["restart"], os.path.join(rundir, "JOBIPH"))
                # in case of QM/MM X/R/MS-PT2 calculation (second run needed) Do_Rotate is automatically created in main directory
                # so if it is present, we move it to rundir
                if os.path.exists("Do_Rotate.txt"):
                    shutil.move("Do_Rotate.txt", rundir)
                # in case of PT2 dynamics with THS Prev.JobMix is automatically created in main directory
                # so if it is present, we move it to rundir
                if os.path.exists("Prev.JobMix") and qm.save_PT2_WF:
                    shutil.move("Prev.JobMix", rundir)
                if os.path.exists("Prev.JobIph") and qm.save_PT2_WF:
                    shutil.move("Prev.JobIph", rundir)

                # when parallel calculation, use apply_async method to run asincronously
                if pool is not None:
                    pool.apply_async(QM._lauchMolcasQM, args=(rundir, filename, qmexe))
                else:  # when serial calculation, run normally
                    QM._lauchMolcasQM(rundir, filename, qmexe)

        # ---------------------------------------------------------------------------------------------------------

        # parallel execution of the jobs
        if pool is not None:
            pool.close(), pool.join()

        # ---------------------------------------------------------------------------------------------------------

        for qm, filename, rundir in zip(QMcalculations, filenameList, rundirList):

            # extract results from the output files and save them in qm
            if isinstance(qm.inputData, GaussianInput):
                # input is of type GAUSSIAN INTERFACE
                qm.outputData = GaussianOutput(filename, rundir, qm.inputData.otheropt["SP"], qm.inputData.otheropt["verbose"])
                # read orbitals using the ad-hoc class
                try:
                    qm.orbitals = Orbitals.readMOfromGaussianOutput(os.path.join(rundir, filename+".log"))
                except (ValueError, IndexError) as e:
                    qm.orbitals = None

            elif isinstance(qm.inputData, SharcInterfaceInput):
                # input is of type SHARC INTERFACE
                qm.outputData = SharcInterfaceOutput(filename, rundir)

            elif isinstance(qm.inputData, MolcasInput):
                #input is of type MOLCAS INTRFACE
                qm.outputData = MolcasOutput(filename, rundir, qm.inputData.calctype, qm.inputData.otheropt["SP"])
                if qm.PT2_QMMM or qm.is_SSPT2_first_run:
                        with open('Do_Rotate.txt', 'w') as DoRotate:
                            # Do-Rotate.txt file contains the transpose of the eigenvector matrix printed by Molcas
                            for row in qm.outputData.get("eigenvectors"):
                                for col in row:
                                    DoRotate.write("{0:12.8f} ".format(col))
                                DoRotate.write("\n")
                            DoRotate.write(" CMS-PDFT")
                if qm.save_PT2_WF:
                    #in case of SS-CASPT2, save WF only from second PT2 run (identified by is_SSPT2_first_run attribute)
                    #the first PT2 run is indeed only used to retrieve root energy order, the second run is instead the real CASPT2 run
                    if qm.inputData.calctype == 'SSPT2':
                        shutil.copy(os.path.join(rundir, "JOBIPH"), "Prev.JobIph")
                    elif qm.inputData.calctype in ['MSPT2', 'XMSPT2', 'RMSPT2']:
                        shutil.copy(os.path.join(rundir, "molcas.JobMix"), "Prev.JobMix")

            # append log from output file processing to the QM log
            qm.log += qm.outputData.log

            # in case, handle QM execution errors
            if qm.outputData.get("termination") is None:
                logwrt.fatalerror('Cannot determine the final status of the QM calculation ...\n'
                                  'terminating COBRAMM execution')
            else:
                if qm.outputData.get("termination") == 1:
                    errmsg = "QM terminated abnormally. Here is an excerpt of the QM log file:\n" \
                             "-----------------------------------------------------------------------------------\n"
                    for line in qm.outputData.get("errormsg"):
                        errmsg += ">   " + line + "\n"
                    errmsg += "-----------------------------------------------------------------------------------\n" \
                              "For further details please check QM log file. "
                    logwrt.fatalerror(errmsg)

            # when the electric field has been computed at the point charges and the gradient is not available
            if qm.outputData.get("elfield") is not None and qm.outputData.get("gradcharges") is None:
                # process electric field to compute the forces on the point charges dependent on QM-MM interaction
                gradcharges = {}
                for nstate, elfield in qm.outputData.get("elfield").items():
                    # when the electric field is set to None for a specific state, just set the gradient to None
                    if elfield is None:
                        gradcharges[nstate] = None
                    else:
                        g = [], [], []
                        for charge, Ex, Ey, Ez in zip(qm.inputData.otheropt["field"][1], *elfield):
                            g[0].append(-Ex * charge), g[1].append(-Ey * charge), g[2].append(-Ez * charge)
                        gradcharges[nstate] = g
                # store the computed forces
                qm.outputData.set("gradcharges", gradcharges)

            # when the gradient of the whole set of point charges is available, process the data
            # to extract the gradient on the mobile charges alone
            if qm.outputData.get("fullgradcharge") is not None and qm.outputData.get("gradcharges") is None:
                gradcharges = {}
                for nstate, fullgrad in qm.outputData.get("fullgradcharge").items():
                    # when for a state the gradient is None for a specific state, just set the gradient to None
                    if fullgrad is None:
                        gradcharges[nstate] = None
                    else:
                        g = [], [], []
                        for mobile, gx, gy, gz in zip(qm.inputData.otheropt["charges"][2], *fullgrad):
                            if mobile:
                                g[0].append(gx), g[1].append(gy), g[2].append(gz)
                        gradcharges[nstate] = g
                # store the gradient on the mobile charges
                qm.outputData.set("gradcharges", gradcharges)

            # if the force on the point charges has been computed, process the data to acompute when needed the
            # self-energy of the point charges (M-ML interaction) and to remove the self-interaction counter term
            if qm.outputData.get("gradcharges") is not None and qm.outputData.get("gradcharges_extraterm"):
                # initialize the energy and the gradient corresponding to the selfinteraction on the M charges
                FX, FY, FZ = [], [], []
                E = 0.0
                # loop over the charges of the M layer
                for ch, xAng, yAng, zAng in zip(qm.inputData.otheropt["field"][1], *qm.inputData.otheropt["field"][0]):
                    # convert to atomic units
                    x, y, z = xAng / constants.Bohr2Ang, yAng / constants.Bohr2Ang, zAng / constants.Bohr2Ang
                    # initialize component of the force
                    Fx, Fy, Fz = 0.0, 0.0, 0.0
                    # loop over the charges of the M and L layers
                    for ch1, x1Ang, y1Ang, z1Ang in zip(qm.inputData.otheropt["charges"][1],
                                                        *qm.inputData.otheropt["charges"][0]):
                        # convert to atomic units
                        x1, y1, z1 = x1Ang / constants.Bohr2Ang, y1Ang / constants.Bohr2Ang, z1Ang / constants.Bohr2Ang
                        # use this couple of charges only when the two charges are not the same
                        if x1 != x or y1 != y or z1 != z:
                            # distance of the two charges
                            D2 = (x - x1) ** 2 + (y - y1) ** 2 + (z - z1) ** 2
                            # increment energy
                            E = (ch1 * ch) / math.sqrt(D2) + E
                            # increment force components
                            F = (ch1 * ch) / D2
                            N = math.sqrt(D2)
                            Fx = ((x - x1) / N) * F + Fx
                            Fy = ((y - y1) / N) * F + Fy
                            Fz = ((z - z1) / N) * F + Fz
                    # now store the component computed for the M layer atom
                    FX.append(Fx), FY.append(Fy), FZ.append(Fz)
                # store the self energy
                qm.outputData.set("selfenergy", E / 2.0)
                # get the actual force read from outputData
                newgradcharges = {}
                for nstate, Fmol in qm.outputData.get("gradcharges").items():
                    FmolX = list(-np.array(Fmol[0]) + np.array(FX))
                    FmolY = list(-np.array(Fmol[1]) + np.array(FY))
                    FmolZ = list(-np.array(Fmol[2]) + np.array(FZ))
                    newgradcharges[nstate] = [FmolX, FmolY, FmolZ]
                # store the computed forces
                qm.outputData.set("gradcharges", newgradcharges)

            # save all the stuff to sef
            # sef = shelve.open("cobram-sef")
            # sef['gradient'] = qm.outputData.get("gradient")[qm.outputData.get("optstate")]
            # sef['gradch'] = qm.outputData.get("gradcharges")[qm.outputData.get("optstate")]
            # sef['dipole'] = qm.outputData.get("dipole")
            # sef['charges'] = qm.outputData.get("charges")
            # if qm.outputData.get("optstate") is not None:
            #     sef['state'] = qm.outputData.get("optstate")
            #     sef['newstate'] = qm.outputData.get("optstate")
            # sef.close()
