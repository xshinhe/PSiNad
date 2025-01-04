#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# KIDS SCRIPTS (adapted from COMBRAMM)
# Author: xshinhe
#
# Copyright (c) 2024 Peking Univ. - GNUv3 License
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

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

import kids_env  # functions to handle third-party software environment
import constants  # values of physical constants and conversion factors
from kids_log import Timing, Log  # keep timings of the different sections of the code
from QMOutput import QMOutput     # template class for the output of a QM calculation
from drivers.adfDriver import AdfInput, AdfOutput
from drivers.bagelDriver import BagelInput, BagelOutput
from drivers.bdfDriver import BdfInput, BdfOutput
from drivers.columbusDriver import ColumbusInput, ColumbusOutput
from drivers.DFTbabyDriver import DFTbabyInput, DFTbabyOutput
from drivers.gamessDriver import GamessInput, GamessOutput  # objects to handle Gaussian Input and Output
from drivers.gaussianDriver import GaussianInput, GaussianOutput  # objects to handle Gaussian Input and Output
from drivers.mndoDriver import MNDOInput, MNDOOutput
from drivers.molcasDriver import MolcasInput, MolcasOutput
from drivers.molproDriver import MolproInput, MolproOutput
from drivers.orcaDriver import OrcaInput, OrcaOutput
from drivers.pyscfDriver import PyscfInput, PyscfOutput
from drivers.qchemDriver import QchemInput, QchemOutput
from drivers.sharcQMDriver import SharcQMInput, SharcQMOutput  # objects to handle QM with SHARC QM Interfaces
from drivers.turbomoleDriver import TurbomoleInput, TurbomoleOutput
from orbitals import Orbitals  # object to parse and store orbital information from QM output files

# imports of local modules


###################################################################################################################


class QM:
    """Object that generate, store and process data for the QM part of the QM/MM calculation"""

    # CLASS variables that controls the preparation of the QM_DATA_STORAGE directory
    # variable with the name of the directory where to store the QM output files (orbitals and logs)
    _QM_DATA_STORAGE = 'QM_data'
    # initialization status of the QM data storage
    _QM_DATA_SETUP = False

    @Timing("QM input")
    def __init__(self, ks_config, geometry, charges, displacement=None, restart=None,
                 setState=None, forceState=None, MS_second_run=False, SS_second_run=False, 
                 gradOnly=False, NAC=False): 
        """Constructor of the QM class. The input of the QM calculation is defined when
        the class instance is created. """

        self.qmsolver = ks_config.args.qmsolver
        self.qmexe = ks_config.get_nested(f'QM.{self.qmsolver.upper()}.path', '')

        # define a string attribute to store the log related to QM input and output, to be printed when convenient
        self.log = ""

        # allowing multiple scheme controlled by level
        level: int = 0 # default to 0
        if ks_config.args.type.isdigit():
            level = int(ks_config.args.type)

        # when the dictionary displacement is given, use its arguments to set up a displaced modelH geometry
        if displacement is not None:
            modelHgeom = geometry.modelHdisplace(
                displacement["iAt"], displacement["iCoord"], displacement["iDir"],
                displacement["length"], displacement["linkfollow"])
        else:
            modelHgeom = geometry.modelH

        # following attributes are needed to activate special actions in case of PT2 calculations
        # for QM/MM X/R/MS-PT2 calculations different from SP, a second run is needed on the roatted CASSCF states to retrieve electric field
        self.PT2_QMMM = False
        #in case of SS-CASPT2 calculations requiring gradients, the &CASPT2(GRDT) and &ALASKA routines must be executed in a secondary calculation
        #(after the root order is retrieved from the first run)
        self.is_SSPT2_first_run = False
        # for PT2 nonadiabatic dynamics with THS scheme, the PT2 WF must be saved for overlap calculation at the next step
        self.save_PT2_WF = False

        # create dictionary with additional options for the global calculation
        optdict = {}
        # background charges, only for a calculation with atoms in the H or M layers
        if geometry.NatomMM > 0: optdict["charges"] = geometry.pod, charges.CRG_emb, geometry.pod_mobileatoms
        if geometry.NatomM > 0: optdict["field"] = geometry.MEDIUM, charges.CRG_MEDIUM
        if ks_config.args.type != 'sp':
            optdict["forces"] = True
            optdict["SP"] = False            
        else:
            optdict["SP"] = True

        # # increase verbosity in gaussianALL.log, do not delete AO, MO and CI deffinitions
        # if command[2] == '2':
        #     optdict["verbose"] = True
        # else:
        #     optdict["verbose"] = False
        
        # when the restart input value contains the path of an available chk file, set restart to true
        if restart and os.path.exists(restart):
            self.log += '{0} file found! \nIt will be used as restart file ' \
                        'for initial WF guess.\n\n'.format(restart)
            optdict["restart"] = restart

        # option to disable nosymm control on input file
        if ks_config.get_nested('QM.nosymm_error', True): 
            optdict["suppress_nosymm_error"] = True

        # when TDcouplings are requested, require tamm-damkov and force printing of all excitations
        if ks_config.get_nested('QM.tdcoupling', False):
            optdict["tdcouplings"] = True

        if ks_config.get_nested('QM.suppress_excited', False):
            optdict["suppress_excited"] = True

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
        if ks_config.args.type == 'ci':
            if ks_config.get_nested('QM.bp_type', 'gmean') == 'nacs':
                self.inputData = None
                self.outputData = None
                Log.fatalError('CI optimization with gaussian is only available with gmean branching plane! Please set bp_type = "gmean" (keyword BranchPlane=gmean)\n')
            elif ks_config.get_nested('QM.bp_type', 'gmean') == 'gmean':
                optdict["CIopt"] = True
                optdict["newBP"] = True
                optdict["upper_state"] = int(ks_config.get_nested('QM.relax_state', 2))
        else:
            optdict["CIopt"] = False
        
        if False:
            pass

        elif ks_config.args.qmsolver == 'adf':
            Log.fatalError(f"{ks_config.args.qmsolver} is not yet supported")

        elif ks_config.args.qmsolver == 'bagel':
            keys = ks_config.get_nested('QM.BAGEL.level%d'%level, '')
            self.inputData = BagelInput(keys, modelHgeom,
                                           geometry.getAtomLabels("modelH"), optdict)

        elif ks_config.args.qmsolver == 'bdf':
            keys = ks_config.get_nested('QM.BDF.level%d'%level, '')
            self.inputData = BdfInput(keys, modelHgeom,
                                           geometry.getAtomLabels("modelH"), optdict)

        elif ks_config.args.qmsolver == 'columbus':
            Log.fatalError(f"{ks_config.args.qmsolver} is not yet supported")

        elif ks_config.args.qmsolver == 'dftbaby':
            Log.fatalError(f"{ks_config.args.qmsolver} is not yet supported")

        elif ks_config.args.qmsolver == 'gamess':
            Log.fatalError(f"{ks_config.args.qmsolver} is not yet supported")

        elif ks_config.args.qmsolver == 'gaussian':  # QM calculation by Gaussian

            # extract sections of the ks_config file that are related to the QM code
            keys =  ks_config.get_nested('QM.GAUSSIAN.level%d'%level, '')
            gen  =  ks_config.get_nested('QM.GAUSSIAN.gen', '') # gaussian basis set
            gaussadd =  ks_config.get_nested('QM.GAUSSIAN.add', '') # gaussian additional parts
            gaussweights =  ks_config.get_nested('QM.GAUSSIAN.weights', '') # gaussian weights of the CASSCF
            # save instance of GaussianInput as the input file for the QM calculation
            self.inputData = GaussianInput(keys, gen, gaussadd, gaussweights, modelHgeom,
                                           geometry.getAtomLabels("modelH"), optdict)
            # print(self.inputData.fileText)

        elif ks_config.args.qmsolver == 'mndo':
            keys = ks_config.get_nested('QM.MNDO.level%d'%level, '')
            self.inputData = MNDOInput(keys, modelHgeom,
                                           geometry.getAtomLabels("modelH"), optdict)

        
        elif ks_config.args.qmsolver == 'molcas':  # QM calulation by Molcas
            keys = ks_config.get_nested('QM.MOLCAS.level%d'%level, '')
            basisset = ks_config.get_nested('QM.MOLCAS.basisset', '') # user-defined molcas basis set (optional)
            sewardkey = ks_config.get_nested('QM.MOLCAS.seward', '')

            # when the restart input value contains the path of an available restart file, set restart to true
            if restart and os.path.exists(restart):
                self.log += '{0} file found! \nIt will be used ' \
                            'for initial WF guess.\n\n'.format(os.path.basename(restart))
                optdict["restart"] = restart

            optdict["basis"] = ks_config.get_nested('QM.MOLCAS.OWBASIS', '6-31g')
            optdict["cartesian"] = ks_config.get_nested('QM.MOLCAS.OWBTYPE', 'cartesian') == 'cartesian'
            optdict["ricd"] = ks_config.get_nested('QM.MOLCAS.OW_USERI', True)
            optdict["cdth"] = ks_config.get_nested('QM.MOLCAS.OW_THRES', 1.0e-4)

            # pass ks_config[1] (i.e. type of calculatio ) to MolcasInput to build appropriate input file
            optdict['args_type'] = ks_config.args.type
            optdict['mdv_type'] = ks_config.get_nested('QM.mdv_type', 0)

            # True if this is the second QM run for MS-XMS-RMS PT2 QM/MM (to retrieve electric field)
            if MS_second_run or SS_second_run:
                optdict["ROST"] = True
            if ks_config.get_nested('QM.highest_NAC_state', 0) != 0:
                optdict["highest_NAC_state"] = int(ks_config.get_nested('QM.highest_NAC_state', 0))
            if gradOnly:
                optdict["gradOnly"] = True
            if NAC:
                optdict["NAC_for_vel_rescale"] = NAC

            # save instance of MolcasInput as the input file for the QM calculation
            self.inputData = MolcasInput(keys, basisset, sewardkey, modelHgeom,
                            geometry.getAtomLabels("modelH"), optdict)
            # set an attributes to remember if a second PT2 run is needed
            # (in case of QM/MM calculations different from SP a second run is needed to retrieve electric field)
            if (self.inputData.calctype in ['MSPT2', 'XMSPT2', 'RMSPT2']) and not optdict["SP"] and (geometry.NatomMM > 0) and not MS_second_run:
                self.PT2_QMMM = True
            # (in case of SS-CASPT2 calculations different from SP a second run is needed to retrieve gradient of the correct state, after checking possible root flipping)
            if (self.inputData.calctype == 'SSPT2') and not optdict["SP"] and not SS_second_run:
                self.is_SSPT2_first_run = True
            # in case of PT2 nonadiabatic dynamics (ks_config 1 + ks_config 85) with THS scheme (ks_config 14)
            if (self.inputData.calctype in ['MSPT2', 'XMSPT2', 'RMSPT2']) and (ks_config.args.type == 'mdv') and (ks_config.get_nested('QM.sh_type', 'nohop') != 'nohop') and (ks_config[14] == '1') and not MS_second_run:
                self.save_PT2_WF = True
            elif self.inputData.calctype == 'SSPT2' and (ks_config[1] == 'mdv') and (int(ks_config.get_nested('QM.sh_type', 'nohop')) > 0) and (ks_config[14] == '1') and not self.is_SSPT2_first_run:
                self.save_PT2_WF = True

        elif ks_config.args.qmsolver == 'molpro':  # QM calulation by MOLPRO
            Log.fatalError("Molpro is not yet supported")
            self.inputData = None

        elif ks_config.args.qmsolver == 'orca':
            keys = ks_config.get_nested('QM.ORCA.level%d'%level, '')
            self.inputData = OrcaInput(keys, modelHgeom,
                                           geometry.getAtomLabels("modelH"), optdict)
        elif ks_config.args.qmsolver == 'pyscf':
            Log.fatalError(f"{ks_config.args.qmsolver} is not yet supported")

        elif ks_config.args.qmsolver == 'qchem':
            Log.fatalError(f"{ks_config.args.qmsolver} is not yet supported")

        elif ks_config.args.qmsolver == 'sharc':
            qmcode = ks_config.get_nested('QM.SHARC.qm')
            if forceState is not None:
                Log.fatalError("cannot change electronic state in QM calculation with SHARC interface")
            self.inputData = SharcQMInput(qmcode, modelHgeom, 
                geometry.getAtomLabels("modelH"), optdict)

        elif ks_config.args.qmsolver == 'turbomole':
            Log.fatalError(f"{ks_config.args.qmsolver} is not yet supported")

        else:
            Log.fatalError(f"{ks_config.args.qmsolver} is not yet supported")

        self.outputData = QMOutput()
        self.orbitals = None

    # =============================================================================================================
    
    def __del__(self):
        """Destructor for the QM class: destroys the input data (defined on instance construction)
        and the output data, when defined"""

        del self.inputData
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
                Log.fatalError("Looks like this calculation needs a NAC that is not available...please check input\n")
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

    @property
    def nacdict(self):
        return self.outputData.get("nac")

    @property
    def nacCRGdict(self):
        return self.outputData.get("naccharges")

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
            Log.writewarning("QM output is not available, log and orbital files will not be saved")
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
    def _launchQM0(rundir, fninp, fnlog, fnerr, qmexe):
        """ This simple function executes the ks_config for running Gaussian QM
         It is needed to store the functions to execute for a parallel execution """

        # store starting dir
        startdir = os.getcwd()
        # move to the directory of the calculation
        os.chdir(rundir)
        # run qm calculation
        with open(fninp, "r") as finp:
            with open(fnlog, "w") as fout:
                with open(fnerr, "w") as ferr:
                    subprocess.call(["nice", "-1", qmexe], stdin=finp, stdout=fout, stderr=ferr)
        os.chdir(startdir)

    @staticmethod
    def _launchQM1(rundir, fninp, fnlog, fnerr, qmexe):
        """ This simple function executes the ks_config for running Gaussian QM
         It is needed to store the functions to execute for a parallel execution """

        # store starting dir
        startdir = os.getcwd()
        # move to the directory of the calculation
        os.chdir(rundir)
        # run qm calculation
        with open(fnlog, "w") as fout:
            with open(fnerr, "w") as ferr:
                subprocess.call(["nice", "-1", qmexe, fninp], stdout=fout, stderr=ferr)
        os.chdir(startdir)

    # =============================================================================================================
    @staticmethod
    @Timing("QM run")
    def runQM(QMcalculations, memory="500MB", nprocQM=1, nprocParall=1):
        """Given a list of input files, defined as instances of the QM class, run the QM program
        and extract the relevant results from the output files, that are then saved in the QM instances
        themselves and can be later read as attributes of the instances"""

        # make argument a list
        if type(QMcalculations) is not list: QMcalculations = [QMcalculations]

        # for SHARC, the parallel calculation is not possible
        if np.any([isinstance(qm.inputData, SharcQMInput) for qm in QMcalculations]):
            Log.writewarning("parallel execution is not available with SHARC, setting nproc = 1")
            nprocQM = 1

        # if the QM calculations are run in parallel, initialize Pool object from the multiprocessing library
        if nprocParall > 1 and len(QMcalculations) > 1:
            pool = mp.Pool(processes=nprocParall)
        else:
            pool = None

        # store starting dir
        startdir = os.getcwd()

        # lists to store the names of the directory and of the file for the calculations
        rundirList, fninpList, fnlogList = [], [], []

        # in the standard case, run the calculations serially
        for qm in QMcalculations:
            qmsolver = qm.qmsolver

            checkResults = qm.inputData.checkEnv()
            if not checkResults[0]: Log.fatalError(checkResults[1])

            qmexe = qm.qmexe 
            if qmexe == '':
                qmexe = kids_env.checkQMExec(qmsolver)
            
            fninp = f"{qmsolver}-QM.inp"
            fnlog = f"{qmsolver}-QM.log"
            fnerr = f"{qmsolver}-QM.err"
            rundir = QM._availablePathForQM()

            if isinstance(qm.inputData, SharcQMInput):
                fninp = 'QM.in'
                fnlog = 'QM.log'
                fnerr = 'QM.err'
                rundir = '.'
                with open(os.path.join(rundir, 'charge.dat'), "w") as fcharge:
                    fcharge.write(qm.inputData.chargeText())

            elif isinstance(qm.inputData, GaussianInput):
                fninp = f"{qmsolver}-QM.com"
                # copy restart file to the working directory, to use the file for restart
                if qm.inputData.otheropt["restart"] is not None:
                    shutil.copy(qm.inputData.otheropt["restart"], 
                        os.path.join(rundir, fninp.replace('.com', '.chk')))

            elif isinstance(qm.inputData, BagelInput):
                fninp = f"{qmsolver}-QM.json"

            elif isinstance(qm.inputData, MolcasInput):
                # copy restart file to the working directory, to use the file for restart
                if qm.inputData.otheropt["restart"] is not None:
                    if qm.inputData.otheropt["restart"] in ("INPORB", "molcas.RasOrb"):
                        shutil.copy(qm.inputData.otheropt["restart"], os.path.join(rundir, "INPORB"))
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

            else:
                if 'restart' in qm.inputData.otheropt and qm.inputData.otheropt["restart"] is not None:
                    shutil.copy(qm.inputData.otheropt["restart"], rundir)
  
            fninpList.append(fninp)
            fnlogList.append(fnlog)
            rundirList.append(rundir)

            # write the input text to file
            with open(os.path.join(rundir, fninp), "w") as finp:
                finp.write(qm.inputData.fileText(fninp, memory=memory, nproc=nprocQM, version=qmexe))
            # when parallel calculation, use apply_async method to run asincronously
            if pool is not None:
                if qmsolver in ['mndo']:
                    pool.apply_async(QM._launchQM0, args=(rundir, fninp, fnlog, fnerr, qmexe))
                else:
                    pool.apply_async(QM._launchQM1, args=(rundir, fninp, fnlog, fnerr, qmexe))
            else:  # when serial calculation, run normally
                if qmsolver in ['mndo']:
                    QM._launchQM0(rundir, fninp, fnlog, fnerr, qmexe)
                else:
                    QM._launchQM1(rundir, fninp, fnlog, fnerr, qmexe)

        # ---------------------------------------------------------------------------------------------------------

        # parallel execution of the jobs
        if pool is not None:
            pool.close(), pool.join()

        # ---------------------------------------------------------------------------------------------------------

        for qm, fnlog, rundir in zip(QMcalculations, fnlogList, rundirList):

            # extract results from the output files and save them in qm
            if False:
                pass

            elif isinstance(qm.inputData, AdfInput):
                qm.outputData = AdfOutput(fnlog, rundir)

            elif isinstance(qm.inputData, BagelInput):
                qm.outputData = BagelOutput(fnlog, rundir)

            elif isinstance(qm.inputData, BdfInput):
                qm.outputData = BdfOutput(fnlog, rundir)

            elif isinstance(qm.inputData, ColumbusInput):
                qm.outputData = ColumbusOutput(fnlog, rundir)

            elif isinstance(qm.inputData, DFTbabyInput):
                qm.outputData = DFTbabyOutput(fnlog, rundir, qm.inputData.calctype)

            elif isinstance(qm.inputData, GamessInput):
                qm.outputData = GamessOutput(fnlog, rundir, qm.inputData.calctype)
        
            elif isinstance(qm.inputData, GaussianInput):
                qm.outputData = GaussianOutput(fnlog, rundir, 
                    qm.inputData.otheropt["SP"], qm.inputData.otheropt["verbose"])
                try:
                    qm.orbitals = Orbitals.readMOfromGaussianOutput(os.path.join(rundir, fnlog))
                except (ValueError, IndexError) as e:
                    qm.orbitals = None

            elif isinstance(qm.inputData, MNDOInput):
                qm.outputData = MNDOOutput(fnlog, rundir, qm.inputData.calctype)

            elif isinstance(qm.inputData, MolcasInput):
                qm.outputData = MolcasOutput(fnlog, rundir, 
                    qm.inputData.calctype, qm.inputData.otheropt["SP"])
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
            
            elif isinstance(qm.inputData, MolproInput):
                qm.outputData = MolproOutput(fnlog, rundir, qm.inputData.calctype, qm.inputData.otheropt["SP"])

            elif isinstance(qm.inputData, OrcaInput):
                qm.outputData = OrcaOutput(fnlog, rundir, qm.inputData.otheropt["SP"])

            elif isinstance(qm.inputData, PyscfInput):
                qm.outputData = PyscfOutput(fnlog, rundir, qm.inputData.calctype)

            elif isinstance(qm.inputData, QchemInput):
                qm.outputData = QchemOutput(fnlog, rundir, qm.inputData.calctype)

            elif isinstance(qm.inputData, SharcQMInput):
                qm.outputData = SharcQMOutput(filename, rundir)

            elif isinstance(qm.inputData, TurbomoleInput):
                qm.outputData = TurbomoleOutput(fnlog, rundir, qm.inputData.calctype)

            
            # append log from output file processing to the QM log
            qm.log += qm.outputData.log

            # in case, handle QM execution errors
            if qm.outputData.get("termination") is None:
                Log.fatalError('Cannot determine the final status of the QM calculation ...\n'
                                  'terminating COBRAMM execution')
            else:
                if qm.outputData.get("termination") == 1:
                    errmsg = "QM terminated abnormally. Here is an excerpt of the QM log file:\n" \
                             "-----------------------------------------------------------------------------------\n"
                    for line in qm.outputData.get("errormsg"):
                        errmsg += ">   " + line + "\n"
                    errmsg += "-----------------------------------------------------------------------------------\n" \
                              "For further details please check QM log file. "
                    Log.fatalError(errmsg)

            # 如果量化程序自己已经计算了embeded电荷上的受力, 应该优先使用这个 (note: doesn't include MM coulumb interaction)
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
                # store the gradient on the mobile charges (only store mobile force!!!)
                qm.outputData.set("gradcharges", gradcharges)

            # so does the nac (if fullnaccharge is accessible)
            if qm.outputData.get("fullnaccharge") is not None and qm.outputData.get("naccharges") is None:
                naccharges = {}
                for nstate in qm.outputData.get("fullnaccharge"):
                    naccharges[nstate] = {}
                    for kstate, fullnac in qm.outputData.get("fullnaccharge")[nstate].items():
                        # when for a state the gradient is None for a specific state, just set the gradient to None
                        if fullnac is None:
                            naccharges[nstate][kstate] = None
                        else:
                            g = [], [], []
                            for mobile, gx, gy, gz in zip(qm.inputData.otheropt["charges"][2], *fullnac):
                                if mobile:
                                    g[0].append(gx), g[1].append(gy), g[2].append(gz)
                            naccharges[nstate][kstate] = g
                # store the gradient on the mobile charges (only store mobile force!!!)
                qm.outputData.set("naccharges", naccharges)

            # 否则, 如果知道(某个态对应的)计算的电场, 那么embed电荷的受力可以计算 (调换次序)
            # 因为MNDO的不能提供与态相关的电场, 所以请使用上面的计算方式
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
