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

import os, shutil # operating system utilities
import sys  # system-specific parameters and functions
import copy  # shallow and deep copy operations
import shelve  # python object persistence

# check that COBRAMM is running on python 3
if sys.version_info < (3, 0): raise RuntimeError("Python 3 is required for running COBRAMM")

# imports of local modules

import CBF
import amber
import optxg
import parallel_numerics
import Vverlet
import tullyNEW
import logwrt  # manages log file output + start/end procedures
import inpdata  # input data read from files
import cobrammenv  # environmental variable for COBRAMM and 3rd-party software
import constants  # physical constants and conversion factors

# imports of local objects

from layers import Layers  # Layers class to manage geometries
from timer import Timer  # keep timings of the different sections of the code
from charge import Charge  # Charge class that stores real and modelH charges
from output import Output  # Output and Step classes to read/write cobramm.xml file
from QM import QM  # QM class controls the QM calculation, define its input and stores the output
from QMMM import QMMM  # join QM and MM data to construct QMMM results
from QMMM import GradientProjectionOnNormalModes
from tullyNEW import Tully
try:
    from tdcoupling import TDCoupling, MolcasTDCoupling  # compute TD couplings from the wavefunction overlap
except ImportError:
    pass

# math libraries

import numpy as np  # numpy library for scientific computation

#####################################################################################################

# source COBRAMM configuration file (if available)
confFile = cobrammenv.setCobrammProfile()
# start COBRAMM with message to log, and start timer for total calculation time

totalTimer = Timer("total")
totalTimer.start()
logwrt.cobramstart()

# check that COBRAMM environment is properly defined
envDefined, errorMsg = cobrammenv.checkCobrammEnv()
if not envDefined: logwrt.fatalerror(errorMsg)
# Print to log information on the environment configuration
logwrt.writelog(cobrammenv.checkConfigFile(confFile))

# start FILE CONTROL section
logwrt.startSection("FILE CONTROL")

# define starting input file
fileinp_1 = "real_layers.xyz"
fileinp_2 = "cobram.command"
# check  for the presence of the two mandatory input files
for fname in [fileinp_1, fileinp_2]:
    if not os.path.isfile(fname): logwrt.fatalerror('mandatory input file {0} not found!'.format(fname))

# make a cpy of the input files
CBF.saveinputs()

# initialize the list with the COBRAMM options
commandhard = CBF.makehard()
print(command)

commandhard.insert(0, '0')
# get the content of the cobramm command file
cobcom = CBF.getCobrammCommand(fileinp_2)
# merge command soft and hard
key_soft = CBF.ReadCobramCommand(cobcom, 'keyword', 'keywords')
commandhard = CBF.key2hard(key_soft, commandhard)
command_soft = CBF.ReadCobramCommand(cobcom, 'command', 'commands')
commandhard = CBF.soft2hard(command_soft, commandhard)
# list of commands
command = commandhard

print(key_soft)
print('2', command)


# set the level of verbosity of the log file
logwrt.setverbositylevel(command[2])

# read the list of the layers
with open(fileinp_1, "r") as f:
    geometry = Layers.from_real_layers_xyz(f.read())
# update geometry with crd if present, otherwise write crd file
geometryfile = fileinp_1  # by default, geometry is read from real_layers.xyz
if command[51] != '11':  # this run is a normal COBRAMM run
    # update the geometry by reading a real.crd file when present
    if os.path.isfile("real.crd"):
        geometry.updatereal()
        geometryfile = "real.crd"
elif command[51] == '11':  # this run is not coupled with SHARC
    # in this case the QMMM.in file is mandatory and is read to update the geometry
    geometry.updatereal("QMMM.in")
    geometryfile = "QMMM.in"
geometry.makerealcrd()

# depending on the type of calculation, set up number of steps
if command[1] in ['freqxg', 'freqxgp']:
    # for a frequency run, the number of steps is given by the number of finite difference to compute
    command[60] = str((len(geometry.list_MEDIUM_HIGH) * 6) + 2)  # numerical freq requires 6N+1 steps (N=nr of atoms)
    # plus 1 because the gaussian optimizer makes a 0th step in which nothing is done
elif command[1] == "sp":  # new option to define SP calculation, type = "sp"
    command[1] = "optxg"
    command[60] = '1'
elif command[60] == 'sp':  # old option to define SP calculation, nsteps = "sp"
    command[60] = '1'
elif command[1] == 'nac':
    command[1] = 'mdvp' 
    command[60] = 1
    command[211] = 1
else:  # otherwise, use the value of command[60] + 1
    command[60] = str(int(command[60]) + 1)
max_step = int(command[60])

# a global variable "smart_numerics" is initialized, which controls the accuracy of the numerical computations
# when surface hopping is inactive (smart_numerics == 0) GRADs are computed by + displacement only;
# near the seam GRADs and NACs are computed by +/- displacements
smart_numerics = 0
if command[1] == 'mdvp' and command[10] == '2':
    command[10] = '0'
    smart_numerics = 1

# the global variable "par_num" is initialized, which controls the behavior of Cobram in the CBF.QM routine
# 0: sequential computation
# 1: computation at the reference geometry during parallel numerics (behaves like a sequential computation)
# 2: collecting data for computing FREQs, GRADs and NACs during parallel numerics
# 3: perform surface hopping
par_num = 0
if command[1] == 'freqxgp':
    par_num = 1
    command[8] = 1
    command[1] = 'freqxg'
if command[1] == 'optxgp':
    par_num = 1
    command[8] = 1
    command[1] = 'optxg'
if command[1] == 'mdvp':
    par_num = 1
    command[8] = 1
    command[1] = 'mdv'
if command[1] == 'ircp':
    par_num = 1
    command[8] = 1
    command[1] = 'irc'
if command[1] == 'cip':
    par_num = 1
    command[8] = 1
    command[1] = 'ci'
if command[1] == 'tsp':
    par_num = 1
    command[8] = 1
    command[1] = 'ts'

# in case of parallel numerical mdv Molcas at SS-PT2 or CASSCF level NACs are computed numerically through the overlap
# of the CASSCF WFs at r and r+dr (key 14 defaults to 0) unless the user specifies explicitly that the NACs should be
# computed analytically through Molpro (i.e. setting key 14 to 1)

# when the calculation is a Gaussian MD run with Surface Hopping, require computation of TD couplings
if command[51] == '1' and command[1] == "mdv" and int(command[85]) > 0:
    if command[14] != '1':
        logwrt.writewarning("Only time-derivative couplings are available with Gaussian, setting 'nacs' (key 14) to 1")
    command[14] = '1'

if command[1] == "mdv" and int(command[85]) > 0 and command[14] == '1':
    if float(command[84]) != float(command[83]):
        logwrt.writewarning("With time-derivative couplings only the time step defined by 'tstep' (key 83) will be used")
    command[84] = command[83]

# initialise system exchange file cobram-sef
sef = shelve.open("cobram-sef", flag="n")
sef['old_step'] = -1
sef['DEarray'] = []
sef['DE_oldarray'] = []
sef['NAC'] = []
sef['NAC_old'] = []
sef['TDC'] = []
sef['TDC_old'] = []
sef['nroots'] = 1
sef['AM'] = []
#sef['AM1'] = []
sef['CIold'] = []
sef['calc_coupl'] = []
sef['smart_numerics'] = smart_numerics
sef['MP2active'] = 0
sef['stop'] = 0.0
sef['tstep_old'] = 0
sef['state'] = 0
sef['newstate'] = 0
sef.close()

# ####################################################
#          print summary of input options
# ####################################################

# construct a string that describes the layers that are defined
logwrt.writelog("\n{0} is requested, with {1}.\n".format(inpdata.getCalcType(command), inpdata.getLayers(geometry)))
if "H" in geometry.calculationType:
    logwrt.writelog("QM third party software : {0}\n".format(inpdata.getQMCode(command)))
if "M" in geometry.calculationType or "L" in geometry.calculationType:
    logwrt.writelog("MM third party software : {0}\n".format(inpdata.getMMCode(command)))
logwrt.writelog("\n")

logwrt.startSubSection("COMMAND OPTIONS SUMMARY")

logwrt.writelog(inpdata.getOptionSummary(command))
logwrt.writelog("\n\nMax number of stpes is set to {} for this calculation\n\n".format(max_step))

# info on number of steps done for a frequency run
if command[1] == 'freqxg':
    if command[18] == '0':
        logwrt.writelog(
            "A standard frequency run has been requested: approximations are enforced to speed up calulation:\n" +
            "QM calculation is performed only when a QM atom moves, equlibrium WF is used when a MM atom moves\n" +
            "1+{0}*6 = {1} QM wf.s will be computed\n".format(geometry.NatomHM, 1 + geometry.NatomHM * 6) +
            "but QM calc.s will be done for {0} atoms ({1} steps)\n".format(geometry.NatomQM, 1 + geometry.NatomQM * 6) +
            "and will be skipped for {0} atoms ({1} steps)\n\n".format(geometry.NatomM, geometry.NatomM * 6))
    else:
        logwrt.writelog(
            "A standard frequency run has been requested with a full QM calculation at all displacemnets (including those of MM atoms):\n" +
            "1+{0}*6 = {1} QM wf.s will be computed\n".format(geometry.NatomHM, 1 + geometry.NatomHM * 6) +
            "and QM calc.s will be done for all atoms ({1} steps)\n\n".format(geometry.NatomQM, 1 + geometry.NatomHM * 6))
    # print a message that explains the behaviour with link atoms
    if geometry.NsubH > 0:  # atom links are present, decide what to do based on command[16] option
        if command[16] == '0':
            logwrt.writelog("In numerical differentiation, H link atoms will be displaced "
                            "together with the neighboring QM atom.\n\n")
        else:
            logwrt.writelog("In numerical differentiation, H link atoms will be kept fixed in the "
                            "position computed at the equilibrium geometry.\n\n")

logwrt.startSection('INPUT MOLECULAR DESCRIPTION')

# print formatted output on layer and molecular geometry definition
logwrt.writelog("Geometry has been read from the file {0}\n\n".format(geometryfile))
logwrt.printLayers(geometry.list_HIGH, geometry.list_MEDIUM, geometry.list_LOW, geometry.atomLink_BA)
logwrt.writelog('Calculation type is {0}, '.format(geometry.calculationType))
if geometry.calculationType in ['HL', 'HM', 'HML']:
    logwrt.writelog('a QM/MM calculation is then requested\n\n')
if geometry.calculationType in ['ML', 'M']:
    logwrt.writelog('a MM calculation is then requested\n\n')
if geometry.calculationType in ['H']:
    logwrt.writelog('a QM calculation is then requested\n\n')
logwrt.printGeom(geometry)

logwrt.startSection('ATOMIC CHARGES')

# create a vector to store charges for the real system
if geometry.calculationType == "H":
    # when the calculation is QM only, create a vector of zeros
    CRG_real = np.zeros(geometry.NatomQM)
else:
    # when the calculation has MM part, read them from AMBER topology file
    # amber.prepare also creates the input file for sander calculation
    CRG_real = amber.prepare(geometry, cobcom, command)

# create a vector to store charges for the model-H system
if geometry.calculationType in ['HML', 'HM', 'HL']:
    CRG_model_H = amber.read_crgA()
elif geometry.calculationType in ["H"]:
    CRG_model_H = np.zeros(geometry.NatomQM)
else:
    CRG_model_H = np.array([])

# create object Charge to store model and embedding charges
charges = Charge(geometry, CRG_real, CRG_model_H)
charges.checkConsistency()

# ####################################################
#                    write xml file
# ####################################################

# get cobramm version
version = cobrammenv.getVersion()
# read model-H.top (if exists) and real.top
try:
    with open('model-H.top') as f:
        modelH_top = f.read()
except IOError:
    modelH_top = None
try:
    with open('real.top') as f:
        real_top = f.read()
except IOError:
    real_top = None
# init Output instance and write first part of the cobramm.xml file
with Timer("xml output"):
    xmlfile = Output()
    xmlfile.write_init(version, geometry, charges, command, modelH_top, real_top)

############################################################
#          Define QM restart file
############################################################

# set the name of the orbital restart file depending on the type of QM
if command[51] == '1':
    #Gaussian
    restartFileName = ["gaussian-QM.chk", "gaussian.chk"]
elif command[51] == '6':
    #Molcas
    #files are in order of increasing precedence (molcas.JobIph has higher priority)
    restartFileName = ["molcas.RasOrb", "INPORB", "molcas.JobIph"]
else:
    restartFileName = ""

# when the file is present in the main dir, set the variable QMRestart for later use
QMRestart = None
for nm in restartFileName:
    if os.path.exists(nm): QMRestart = nm

############################################################
#          Starting optimization/MD
############################################################

step, substep = 0, 1
if command[1] == 'mdv':
    if xmlfile.step:
        _out = Output(parse=True)
        actualTime = _out.get_step(_out.steps).time
        logwrt.startSection('RESTART MD')
        logwrt.writelog('Actual time set to {0:.2f}fs\n\n'.format(actualTime / constants.fs2au))
        del _out
    else:
        actualTime = 0.0

# initialization of variables
displQM = []
QMPrevious = None
prevx1, prevx2 = None, None

while True:
    if max_step > 1:
        if command[1] == "irc":
            logwrt.startSection('ENTER STEP {0}.{1}'.format(step, substep))
        else:
            logwrt.startSection('ENTER STEP {0}'.format(step))
    else:
        logwrt.startSection('QM/MM SINGLE POINT')

    sef = shelve.open("cobram-sef")
    sef['MDstep'] = step
    if command[203] == '1' and sef['MP2active'] == 1:
        par_num = 0
        sef['MP2active'] = 0
    sef.close()

    ############################################################
    #                  MM calculations
    ############################################################

    with Timer("MM section"):
        Fxyz_second, E_second, E_modelnoc, Fxyz_modelnoc, Fxyz_modelH, E_modelH = CBF.MM(
            step, geometry, command, cobcom)
        MM_Results = [Fxyz_second, E_second, E_modelnoc, Fxyz_modelnoc, Fxyz_modelH, E_modelH]

    sef = shelve.open("cobram-sef")
    sef['Fxyz_modelH'] = Fxyz_modelH
    sef.close()

    ##########################################################
    #                 QM calculations
    ##########################################################

    with Timer("QM section"):

        # define the type of QM calculation: old style interface or new OOP interface
        if command[51] == '1' or command[51] == '11' or command[51] == '6':
            NewQM = True
        else:
            NewQM = False

        # new QM calculation
        if NewQM:

            #for normal run (no DEBUG run), clean QM directories to avoid mess:
            if not logwrt.DEBUG_COBRAMM_RUN:
                QM_DIRECTORIES = [item for item in os.listdir() if os.path.isdir(item) and item.startswith('qmCalc')]
                if QM_DIRECTORIES:
                    for qm_dir in QM_DIRECTORIES:
                        shutil.rmtree(qm_dir)

            # parallel, freqxg calculation, in a step which is not the first
            if command[8] == 1 and command[1] == 'freqxg' and step > 1:
                QMCalc = displQM[step - 1]
                logwrt.writelog(QMCalc.log)

            else:    # in all the other case, really do the calculation

                # for a SH calculation, we need to select the initial state and force-set the elect state in step > 0
                stateargs = {}
                if int(command[85]) > 0:
                    if step == 0:  # if this is the first step, define the initial state using command[13]
                        if command[13] != "1":
                            stateargs["setState"] = int(command[13])-1
                    else:  # if this is not the first step, read the actual state from sef
                        sef = shelve.open('cobram-sef')
                        stateargs["forceState"] = sef['state']
                        sef.close()
                elif command[51] == '6': #Molcas needs to know opt state
                    stateargs["setState"] = int(command[13])-1
                
                #check if QM calculation can be skipped from freq run, and in case this is True take a copy of reference QMCalc
                if command[1] == 'freqxg' and command[18] == '0' and parallel_numerics.freq_QM_skipper(step - 1, geometry, command) and step != 0:
                    logwrt.writelog("QM computation is skipped in this step!\n\n")
                    QMCalc = copy.deepcopy(QMCalc_ref)
                else:
                    # initialize QM instance to define input of the QM calculation
                    QMCalc = QM(cobcom, command, geometry, charges, step, QMRestart, **stateargs) #--> F: added step as attribute to QM class in order to build correct Molcas input
                    # run the calculation
                    QM.runQM(QMCalc, memory=command[53], nprocQM=int(command[7])), logwrt.writelog(QMCalc.log)
                    '''for QM/MM calculations (other then SP) of (X/R)MSPT2 type, we need to run two QM calculations in series
                    first, we run the proper QM PT2 calculation and produce energies and gradient of QM atoms
                    in a subsequent run, we rotate the CASSCF states accoridng to PT2 eigenvalues (ROST keyword) and perform an SS-PT2 calc to obtain
                    the properties of the (rotated) state'''
                    if QMCalc.PT2_QMMM:
                        #use restart files from first run
                        QMRestart = QMCalc.saveRestartFile()
                        QMCalc2 = QM(cobcom, command, geometry, charges, step, QMRestart, **stateargs, MS_second_run=True)
                        # run the calculation
                        QM.runQM(QMCalc2, memory=command[53], nprocQM=int(command[7])), logwrt.writelog(QMCalc2.log)
                        #overwrite properties
                        QMCalc.outputData.dataDict["elfield"] = QMCalc2.outputData.dataDict["elfield"]
                        QMCalc.outputData.dataDict["dipole"] = QMCalc2.outputData.dataDict["dipole"]
                        QMCalc.outputData.dataDict["charges"] = QMCalc2.outputData.dataDict["charges"]
                        #join output of first and second run (for printout in qmALL.log)
                        # DEBUG: I THINK THERE IS AN ERROR HERE...CHECK
                        QMCalc.outputData.dataDict["outfile"] = QMCalc.outputData.dataDict["outfile"] + QMCalc2.outputData.dataDict["outfile"]
                    elif QMCalc.is_SSPT2_first_run:
                        logwrt.writewarning("SS-CASPT2 state order might be different from root order! States and relative properties will be reordered according to CASPT2 energies...")
                        #use restart files from first run
                        QMRestart = QMCalc.saveRestartFile()
                        QMCalc2 = QM(cobcom, command, geometry, charges, step, QMRestart, **stateargs, SS_second_run=True)
                        #run the calculation
                        QM.runQM(QMCalc2, memory=command[53], nprocQM=int(command[7])), logwrt.writelog(QMCalc2.log)
                        #join output of first and second run (for printout in qmALL.log)
                        QMCalc2.outputData.dataDict["outfile"] = QMCalc.outputData.dataDict["outfile"] + QMCalc2.outputData.dataDict["outfile"]
                        #overwrite data
                        #QMCalc.outputData.dataDict = QMCalc2.outputData.dataDict
                        QMCalc = QMCalc2

                #print gradient 
                if command[1] != 'irc' and QMCalc.gradient(QMCalc.outputData.get("optstate")) is not None:
                    logwrt.writelog("\nQM gradient of High Layer (+ link atoms) for state {0:2d}\n".format(QMCalc.outputData.get("optstate")+1), 2)
                    if geometry.NsubH == 0:
                        logwrt.writelog(logwrt.matrix_prettystring(np.array(QMCalc.gradient(QMCalc.outputData.get("optstate"))).T, ".6f", geometry.getAtomLabels("HIGH")), 2)
                    else:
                        logwrt.writelog(logwrt.matrix_prettystring(np.array(QMCalc.gradient(QMCalc.outputData.get("optstate"))).T, ".6f", geometry.getAtomLabels("modelH")), 2)
                    logwrt.writelog("\n", 2)
                if command[1] == 'ci' and (QMCalc.gradient(QMCalc.outputData.get("optstate")-1) is not None):
                    lowerstate = QMCalc.outputData.get("optstate") - 1
                    logwrt.writelog("\nQM gradient of High Layer (+ link atoms) for state {0:2d}\n".format(lowerstate+1), 2)
                    if geometry.NsubH == 0:
                        logwrt.writelog(logwrt.matrix_prettystring(np.array(QMCalc.gradient(lowerstate)).T, ".6f", geometry.getAtomLabels("HIGH")), 2)
                    else:
                        logwrt.writelog(logwrt.matrix_prettystring(np.array(QMCalc.gradient(lowerstate)).T, ".6f", geometry.getAtomLabels("modelH")), 2)
                    logwrt.writelog("\n", 2)

                #print NAC
                #in case of ciopt, only NAC between the desired states are printed
                if (command[1] == 'ci' and command[19] == '0'):
                    logwrt.writelog("\nQM NAC of High Layer (+ link atoms) for states {0:2d} {1:2d}\n".format(QMCalc.outputData.get("optstate"), QMCalc.outputData.get("optstate")+1), 2)
                    if geometry.NsubH == 0:
                        logwrt.writelog(logwrt.matrix_prettystring(np.array(QMCalc.nac(QMCalc.outputData.get("optstate")-1, QMCalc.outputData.get("optstate"))).T, ".6f", geometry.getAtomLabels("HIGH")), 2)
                    else:
                        logwrt.writelog(logwrt.matrix_prettystring(np.array(QMCalc.nac(QMCalc.outputData.get("optstate")-1, QMCalc.outputData.get("optstate"))).T, ".6f", geometry.getAtomLabels("modelH")), 2)
                    logwrt.writelog("\n", 2)
                    #in case of mnd with NACs, all computed NACs are printed
                if (command[1] == 'mdv') and (command[85] == '1') and (command[14] != '1'):
                    NAC = QMCalc.outputData.get("nac")
                    for i in NAC:
                        for j in NAC[i]:
                            if j > i:
                                logwrt.writelog("\nQM NAC of High Layer (+ link atoms) for states {0:2d} {1:2d} (NAC_ij = -NAC_ji)\n".format(i+1, j+1), 2)
                                if geometry.NsubH == 0:
                                    logwrt.writelog(logwrt.matrix_prettystring(np.array(NAC[i][j]).T, ".6f", geometry.getAtomLabels("HIGH")), 2)
                                else:
                                    logwrt.writelog(logwrt.matrix_prettystring(np.array(NAC[i][j]).T, ".6f", geometry.getAtomLabels("modelH")), 2)
                                logwrt.writelog("\n", 2)

                #print overlap matrix (in case of mdv with time-derivatice couplings)
                if (command[1] == 'mdv') and (command[85] == '1') and (command[14] == '1') and command[51] == '6' and step != 0:
                    logwrt.writelog("Overlap matrix with states at previous time step:\n",1)
                    OvMat = QMCalc.outputData.dataDict["psioverlap"]
                    logwrt.writelog(logwrt.matrix_prettystring(OvMat, ".8f"), 1)
                    logwrt.writelog("\n", 1)

                # when this is the first step, initialize the number of the electronic state
                if step == 0:
                    sef = shelve.open('cobram-sef')
                    sef['state'] = QMCalc.outputData.get("optstate")
                    sef['newstate'] = QMCalc.outputData.get("optstate")
                    sef.close()

                # copy and save path of the orbital restart file for the next step
                # in case of frequency calculation, only reference JobIph (step0) is saved to restart all next calculations
                if command[1] == 'freqxg' and step == 0:
                    QMRestart = QMCalc.saveRestartFile()
                    #in case of skipQM, save the reference QMCalc to be used for MM atomic displacemnets
                    if command[18] == '0':
                        QMCalc_ref = copy.deepcopy(QMCalc)
                elif command[1] != 'freqxg':
                    QMRestart = QMCalc.saveRestartFile()

            # in step 0 of a parallel freqxg calculation, compute all the displacements QM in parallel
            if command[8] == 1 and command[1] == 'freqxg' and step == 0:
                logwrt.startSection('ENTER PARALLEL NUMERICS - FREQUENCY CALCULATION')
                displQM = parallel_numerics.newFreqParallelRun(cobcom, command, geometry, charges, QMCalc, **stateargs)

            # at each step of a parallel optimization/MD calculation, compute all the displacements QM in parallel
            if command[8] == 1 and command[1] in ['optxg', 'mdv', 'irc', 'ci', 'ts']:
                logwrt.startSection('ENTER PARALLEL NUMERICS - GRADIENT CALCULATION')
                displQM = parallel_numerics.newGradientParallerRun(cobcom, command, geometry, charges, QMCalc)

            # compute the coupling
            sef = shelve.open('cobram-sef')
            actualstate = sef['state']
            sef.close()
            if command[14] == '1' and not (actualstate == 0 and command [103] == '0') and not (actualstate == 0 and command[90] == '2'):
                if command[51] == '1':
                    tdc = TDCoupling(QMPrevious, QMCalc, float(command[98]), float(command[84])*constants.fs2au, 
                                int(command[7]), float(command[86]),
                                 couplingwithGS=(True if command[90] == "1" else False),
                                 standardHST=(True if command[97] == "1" else False),
                                 transformS=(True if command[97] == "2" else False))
                elif command[51] == '6':
                    tdc = MolcasTDCoupling(QMPrevious, QMCalc, float(command[83])*constants.fs2au)
                
                QMCalc.setTDC(tdc.getTDMatrix())
                #F: I have moved the QMPrevious assignment to the end of QMMM reults
                #QMPrevious = copy.deepcopy(QMCalc)


        # old style QM calculation
        else:

            # initialize QM instance to define input of the QM calculation
            QMCalc = QM(cobcom, command, geometry, charges, QMRestart)

            # for a frequency calculation, decide whether to skip step:
            # step can be skipped if the displacement involves atoms of the M layer
            if command[1] == 'freqxg' and parallel_numerics.freq_QM_skipper(step - 1, geometry, command) and step != 0:
                logwrt.writelog("QM computation is skipped in this step!\n\n")
                QM_Results = copy.deepcopy(QM_result_0th_step)

            # in the standard case, perform QM calculation
            else:
                # run QM calculation
                QM_Results = CBF.QM(command, cobcom, charges, geometry, step, par_num)

            # save results of the first step
            if step == 0: QM_result_0th_step = copy.deepcopy(QM_Results)

            if par_num == 1:
                logwrt.startSection('ENTER PARALLEL NUMERICS')
                QM_Results, par_num = parallel_numerics.run(QM_Results, step, geometry, command, cobcom, charges, par_num)

            sef = shelve.open('cobram-sef')
            if command[1] == 'ci':  # in this case the gradient is available for state and state-1
                _state = [sef['state'], sef['state'] - 1]
                _gradient = [-np.array(QM_Results[1]), -np.array(sef['gradient2'])]
                _gradcharge = [QM_Results[5], sef['gradch2']]
            else:  # in this case only the gradient of the active state is available
                if step > 0:
                    _state = [sef['newstate']]
                else:
                    _state = [sef['state']]
                if len(QM_Results[1][0]) > 0:  # the gradient is not always required (for sp it is not used)
                    _gradient = [-np.array(QM_Results[1])]
                else:
                    _gradient = [None]
                _gradcharge = [QM_Results[5]]
            sef.close()

            #print gradient in case it is not a freqxg calculation and the displaced atom is MM
            if command[1] == 'freqxg' and parallel_numerics.freq_QM_skipper(step - 1, geometry, command) and step != 0:
                pass
            else:
                if int(command[60]) > 1:
                    sef = shelve.open('cobram-sef')
                    logwrt.writelog("\nQM gradient of High Layer (+ link atoms) for state {0:2d}\n".format(sef['state']+1), 2)
                    sef.close()
                    if geometry.NsubH == 0:
                        logwrt.writelog(logwrt.matrix_prettystring(-np.array(QM_Results[1]).T, ".6f", geometry.getAtomLabels("HIGH")), 2)
                    else:
                        logwrt.writelog(logwrt.matrix_prettystring(-np.array(QM_Results[1]).T, ".6f", geometry.getAtomLabels("modelH")), 2)
                    logwrt.writelog("\n", 2)

            # update the QMCalc instance of QM with the results obtained from CBF.QM
            sef = shelve.open('cobram-sef')
            NAC = sef['NAC']
            sef.close()
            QMCalc.setOutputData(QM_Results[0], _gradient, _state, QM_Results[2], QM_Results[4],
                                 QM_Results[3], _gradcharge, NAC)
    if int(command[211]) > 0 :
        break
    ##########################################################
    #         Rebuilding the the charges of the system
    ##########################################################

    # TODO: with the current algorithm, this command update the definition of the charges
    #  only in the instance of the charge object
    #  we should consider a re-implementation of the update of the topology file with QMMM charges,
    #  to work with the microiterative sheme, but the current order of QM and MM should be rethought
    #  since one wants to to the MM dynamics with the values of the charge of the same step, not of the
    #  step before that
    if not list(QMCalc.charges):
        logwrt.writewarning("QM charges are not available: their value will not be updated!\n")
    else:
        charges.rallyCharges(QMCalc.charges)

    ##########################################################
    #         COMPUTING QMMM RESULTS
    ##########################################################

    with Timer("QMMM section"):

        QMMM_Results = QMMM(command, geometry, QMCalc, MM_Results, step, prevx1, prevx2)

        #save QMMM branching plane vectors of previous step in case of ci search with new branching plane
        if command[19] == '1' and command[1] == 'ci':
            prevx1, prevx2 = copy.deepcopy(QMMM_Results.x1), copy.deepcopy(QMMM_Results.x2)

        # getting actual electronic state number from shelve
        sef = shelve.open('cobram-sef')
        actualstate = sef['state']
        if step > 0:
            newstate = sef['newstate']
        else:
            newstate = sef['state']
        sef.close()

        #print gradient
        if command[1] != 'irc' and geometry.calculationType in ['HL', 'HM', 'HML'] and QMCalc.gradient(QMCalc.outputData.get("optstate")) is not None:
            # print High layer gradient
            logwrt.writelog("\nQM/MM gradient of High Layer for state {0:2d}\n".format(newstate+1), 1)
            gradientH=[]
            for iat, ixyz in enumerate(np.array(QMMM_Results.getgradient(newstate)).T):
                if iat+1 in geometry.list_HIGH:
                    gradientH.append(ixyz)
            logwrt.writelog(logwrt.matrix_prettystring(np.array(gradientH), ".6f", geometry.getAtomLabels("HIGH")), 1)
            logwrt.writelog("\n", 1)

            # print High+Medium layer gradient
            logwrt.writelog("\nQM/MM gradient of High+Medium Layer for state {0:2d}\n".format(newstate+1), 2)
            logwrt.writelog(logwrt.matrix_prettystring(np.array(QMMM_Results.getgradient(newstate)).T, ".6f", geometry.getAtomLabels("MEDIUM_HIGH")), 2)
            logwrt.writelog("\n", 2)
        if command[1] == 'ci' and (geometry.calculationType in ['HL', 'HM', 'HML']) and (QMCalc.gradient(QMCalc.outputData.get("optstate") -1 ) is not None):
            lowerstate = QMCalc.outputData.get("optstate") - 1
            # print High layer gradient
            logwrt.writelog("\nQM/MM gradient of High Layer for state {0:2d}\n".format(lowerstate+1), 1)
            gradientH=[]
            for iat, ixyz in enumerate(np.array(QMMM_Results.getgradient(lowerstate)).T):
                if iat+1 in geometry.list_HIGH:
                    gradientH.append(ixyz)
            logwrt.writelog(logwrt.matrix_prettystring(np.array(gradientH), ".6f", geometry.getAtomLabels("HIGH")), 1)
            logwrt.writelog("\n", 1)
    
            # print High+Medium layer gradient
            logwrt.writelog("\nQM/MM gradient of High+Medium Layer for state {0:2d}\n".format(lowerstate+1), 2)
            logwrt.writelog(logwrt.matrix_prettystring(np.array(QMMM_Results.getgradient(lowerstate)).T, ".6f", geometry.getAtomLabels("MEDIUM_HIGH")), 2)
            logwrt.writelog("\n", 2)

        # define gradient for CI optimization or for a regular calculation
        if command[1] != 'ci':
            NewGrad = QMMM_Results.getgradient(newstate)
        else:
            NewGrad = QMMM_Results.cigradient

    # initialize variable to store a new geometry for update (when None, do not update geometry)
    newgeom = None

    if max_step > 1:

        ##########################################################
        #              geometry optimization
        ##########################################################

        if command[1] in ['optxg', 'freqxg', 'irc', 'ci', 'ts']:
            if step == 0: logwrt.startSection('GAUSSIAN EXTERNAL OPTIMIZER SETUP')

            convValues, convTable, newgeom = optxg.GAUSSIAN_optimizator_X(
                    step, geometry, command, cobcom, QMMM_Results.getenergy(newstate), NewGrad, QMCalc.dipole)
            mdItems = [None, None, None, None, None]

        ##########################################################
        #                 molecular dynamics
        ##########################################################

        elif command[1] in ['mdv']:

            # project gradient on the sub-space of normal modes < command[40]
            if float(command[40]) > 0. and ( os.path.isfile("geometry.chk") or os.path.isfile("geometry.chk.gz") ):
                if step == 0:
                    geometry_ref = copy.deepcopy(geometry) 
                NewGrad = GradientProjectionOnNormalModes(geometry_ref, command, NewGrad)

            OldGrad=copy.deepcopy(NewGrad)
            # surface hopping algorithm
            ''' Conditions to enter TSH algorithm:
            - dynamics with hop is requested (command[85] > 0)
            - QM calc using new interfaces (NewQM)
            - NOT (user has requested to switch off hops after hop to GS (command[103]=0) and we are in the GS)
            - NOT (TDDFT dynamics is in GS and GS-ES hop is modeled with deltaE (command[90]=2))'''
            if int(command[85]) > 0 and NewQM and not (actualstate == 0 and command[103] == '0') and not (actualstate == 0 and command[90] == '2'):
                if step == 0:
                    AM = []
                else:
                    AM = tully.AM
                tully = Tully(QMPrevious, QMCalc, actualstate, actualTime, command, cobcom, AM, geometry, step)
                newstate = tully.newstate
                #old version:
                #newstate = tully.tully(command, geometry, actualstate, cobcom, step)

                # propose a hop to the ground state based only on energy difference when command[90] is 2
                if command[90] == "2" and newstate == actualstate and actualstate > 0 and step > 0:
                    sef = shelve.open("cobram-sef", 'r')
                    DEarray = sef['DEarray']
                    sef.close()
                    # when the energy difference of actual state with GS is < command[91], set newstate to 0 (GS)
                    if DEarray[actualstate][0] < float(command[91]):
                        newstate = 0
                        logwrt.writelog('-----------------------------\n')
                        logwrt.writelog('  !!!!Gimme Hop Joanna!!!\n')
                        logwrt.writelog('-----------------------------\n')
                        logwrt.writelog('hopping from state {0} to --> {1}\n'.format(actualstate + 1, newstate + 1))
                        logwrt.writelog('-----------------------------\n\n')
                        # create an empty HOP file
                        with open("HOP", "a") as f:
                            pass
                        # save "newstate" to shelve
                        sef = shelve.open("cobram-sef")
                        sef['newstate'] = newstate
                        sef.close()

                # when a hopping is proposed by the SH algorithm, we need to compute the gradient of the new state
                if newstate != actualstate:
                    logwrt.startSubSection("QM/MM GRADIENT OF THE NEW STATE")
                    # initialize QM instance to define input of the QM calculation (use gradOnly option to avoid useless Molcas routines)
                    # at present, gradOnly has no effect in case of QM software other than Molcas
                    NAC_to_compute = False
                    if command[206] == '0':
                        if command[85] == '2' or command[14] == '1':
                            NAC_to_compute = [actualstate, newstate]
                    if NAC_to_compute:
                        logwrt.writelog("Computing energy gradient for the new state n = {0} as well as <{1}|d/dR|{0}> NAC for velocity rescaling... \n\n".format(newstate + 1, actualstate + 1))
                    else:
                        logwrt.writelog("Computing energy gradient for the new state n = {0}... \n\n".format(newstate + 1))
                    QMCalcNewState = QM(cobcom, command, geometry, charges, step, QMRestart, forceState=newstate, gradOnly=True, NAC=NAC_to_compute)
                    if QMCalcNewState.is_SSPT2_first_run:
                        #in case this is a SS_CASPT2 calculation, actually we need to perform the new QM calc as a SS_second run
                        #so we re-initialize QMCalcNewState in the correct way
                        QMCalcNewState = QM(cobcom, command, geometry, charges, step, QMRestart, forceState=newstate, SS_second_run=True, gradOnly=True, NAC=NAC_to_compute)
                    # run the calculation
                    QM.runQM(QMCalcNewState, memory=command[53], nprocQM=int(command[7])), logwrt.writelog(QMCalc.log)
                    '''just like fpr main QMCalc, in case of QM/MM calculations of (X/R)MSPT2 type, we need to run two QM calculations in series:
                    first run (QMCalcNewState, above) produces energies and gradient of QM atoms
                    second run (QMCalcNewState2, below), rotates the CASSCF states accoridng to PT2 eigenvalues (ROST keyword) and performs SS-PT2 to obtain
                    the properties of the (rotated) new state'''
                    if QMCalcNewState.PT2_QMMM:
                        #use restart files from first run
                        QMRestart = QMCalcNewState.saveRestartFile()
                        QMCalcNewState2 = QM(cobcom, command, geometry, charges, step, QMRestart,forceState=newstate, MS_second_run=True)
                        # run the calculation
                        QM.runQM(QMCalcNewState2, memory=command[53], nprocQM=int(command[7]))
                        #overwrite properties
                        QMCalcNewState.outputData.dataDict["elfield"] = QMCalcNewState2.outputData.dataDict["elfield"]
                        QMCalcNewState.outputData.dataDict["dipole"] = QMCalcNewState2.outputData.dataDict["dipole"]
                        QMCalcNewState.outputData.dataDict["charges"] = QMCalcNewState2.outputData.dataDict["charges"]
                    #print gradient
                    logwrt.writelog("\n QM gradient of High Layer (+ link atoms) for state {0:2d}\n".format(newstate+1), 1)
                    if geometry.NsubH == 0:
                        logwrt.writelog(logwrt.matrix_prettystring(np.array(QMCalcNewState.gradient(newstate)).T, ".6f", geometry.getAtomLabels("HIGH")), 1)
                    else:
                        logwrt.writelog(logwrt.matrix_prettystring(np.array(QMCalcNewState.gradient(newstate)).T, ".6f", geometry.getAtomLabels("modelH")), 1)
                    logwrt.writelog("\n", 1)
                    # compute new QMMM gradient
                    QMMM_ResultsNewState = QMMM(command, geometry, QMCalcNewState, MM_Results, step)
                    NewGrad = QMMM_ResultsNewState.getgradient(newstate)
                    if NAC_to_compute:
                        #save NAC to sef in old format for use in VVerlet (velocity rescaling after hop)
                        NAC = QMCalcNewState.outputData.get("nac")
                        NAC_old_format = []
                        nroots = QMCalcNewState.outputData.get("nroots")
                        for root1 in range(nroots):
                            NAC_old_format.append([])
                            for root2 in range(nroots):
                                try:
                                    NAC_old_format[root1].append([NAC[root1][root2][0], NAC[root1][root2][1], NAC[root1][root2][2]])
                                except:
                                    NAC_old_format[root1].append([[],[],[]])
                        sef = shelve.open("cobram-sef") 
                        sef['NAC'] = NAC_old_format
                        sef.close()

                    # print High layer gradient
                    if command[1] != 'irc' and geometry.calculationType in ['HL', 'HM', 'HML']:
                        logwrt.writelog("\nQM/MM gradient of High Layer for state {0:2d}\n".format(newstate+1), 1)
                        gradientH=[]
                        for iat, ixyz in enumerate(np.array(QMMM_ResultsNewState.getgradient(newstate)).T):
                            if iat+1 in geometry.list_HIGH:
                                gradientH.append(ixyz)
                        logwrt.writelog(logwrt.matrix_prettystring(np.array(gradientH), ".6f", geometry.getAtomLabels("HIGH")), 1)
                        logwrt.writelog("\n", 1)
                        # print High+Medium layer gradient
                        logwrt.writelog("\nQM/MM gradient of High+Medium Layer for state {0:2d}\n".format(newstate+1), 2)
                        logwrt.writelog(logwrt.matrix_prettystring(np.array(QMMM_ResultsNewState.getgradient(newstate)).T, ".6f", geometry.getAtomLabels("MEDIUM_HIGH")), 2)
                        logwrt.writelog("\n", 2)

                    # when using correction of the velocity along the gradient difference
                    if command[206] == '2':
                        logwrt.writelog("Computing gradient difference for velocity scaling... \n")
                        fx_state, fy_state, fz_state = QMMM_Results.getgradient(actualstate)
                        fx_newstate, fy_newstate, fz_newstate = QMMM_ResultsNewState.getgradient(newstate)
                        GD = [[], [], []]
                        # extract from the full gradient the QM part ...
                        for i in range(len(geometry.list_MEDIUM_HIGH)):
                            if geometry.list_MEDIUM_HIGH[i] in geometry.list_HIGH:
                                GD[0].append(fx_state[i] - fx_newstate[i])
                                GD[1].append(fy_state[i] - fy_newstate[i])
                                GD[2].append(fz_state[i] - fz_newstate[i])
                        sef = shelve.open("cobram-sef")
                        sef['GD'] = GD
                        sef.close()

            logwrt.writelog("Integrating Newton's equations of motion with Velocity Verlet... \n")
            if step == 0:
                iniD, CSA = [], []
            MD_results = Vverlet.run(cobcom, geometry, command, NewGrad, OldGrad, QMMM_Results.getenergy(newstate),
                                     step, iniD, CSA)
            mdinfo, mdItems, newgeom, tstep, iniD, CSA = MD_results
            convValues = [None, None, None, None]
            logwrt.writelog("\n")


    ##########################################################
    #           write data in xml file
    ##########################################################

    # preparing arguments for write_step
    # values should NOT be of type string
    energies = [E_second, E_modelnoc, E_modelH, QMCalc.selfenergy, QMCalc.energy(), QMMM_Results.getenergy()]
    if command[1] == 'mdv':
        if os.path.exists('velocity_0.dat'):
            with open('velocity_0.dat', 'r') as f:
                _vel_data = f.readlines()
            vx, vy, vz = [], [], []
            for line in _vel_data:
                line = line.strip().split()
                [l.append(float(line[i])) for i, l in enumerate([vx, vy, vz])]
            _velocity = [np.array(vx), np.array(vy), np.array(vz)]
        else:
            _velocity = None
        _MD_etot, _MD_epot, _MD_ekin = MD_results[1][0], MD_results[1][1], MD_results[1][2]
        _MD_temp = None
        dyn = [actualTime, actualstate, _velocity, _MD_etot, _MD_epot, _MD_ekin, _MD_temp]
    else:
        dyn = None

    if command[1] in ['optxg', 'freqxg', 'ts', 'irc', 'ci'] and max_step > 1:
        _nIRCstep = substep if command[1] == 'irc' else None
        opt = [actualstate, _nIRCstep] + convValues
    else:
        opt = None

    with Timer("xml output"):
        xmlfile.write_step(step, geometry, charges, energies, NewGrad, dyn, opt)

    ##########################################################
    #         Rebuilding the geometry of the system
    ##########################################################

    if newgeom is not None:
        # update the HM coordinates within the geometry object
        geometry.updateHMlayers(newgeom)
        # remove real.crd file and create a new one from current geometry
        if os.path.isfile('real.crd'): os.remove('real.crd')
        geometry.makerealcrd('real.crd')

    ##########################################################
    #                finishing the cycle
    ##########################################################

    # check convergence of IRC point
    if command[1] == 'irc':
        IRCconverged = all(check == "YES" for check in
                           [line.split()[4] for line in convTable.split("\n")[2:] if line.split()])
        if IRCconverged: logwrt.writelog("\n" + ' ' * 27 + '#' * 25 + '\n' + ' ' * 27 +
                                         " IRC POINT {0} CONVERGED!\n".format(step) + ' ' * 27 + '#' * 25 + "\n\n")

    # save orbitals at current step
    if command[1] != 'irc' or IRCconverged or step == 0:

        # copyLog flag to save log file
        copyLog = int(command[99]) > 0 and not step % int(command[99])
        # copyOrb flag to save orbital file if any of the following conditions is true
        # 1. step is multiple of command[100] and it's not a freq calculation
        # 2. step 0 for freq calculation
        # 3. command[100] is negative and it's the last step, still not a frequency calculation
        copyOrb_condition1 = command[1] != 'freqxg' \
                             and int(command[100]) not in [-1, 0] \
                             and not step % int(command[100])
        copyOrb_condition2 = command[1] == 'freqxg' and step == 0
        copyOrb_condition3 = command[1] != 'freqxg' and int(command[100]) < 0 and step == max_step - 1
        copyOrb = np.any([copyOrb_condition1, copyOrb_condition2, copyOrb_condition3])
        
        # the calculation has a QM part, collect files with orbitals
        if "H" in geometry.calculationType:
            if NewQM:
                QMCalc.archiveStep(copyLog, copyOrb, step)
                if newstate != actualstate:
                    QMCalcNewState.archiveStep(copyLog, 0, step)
            else:
                CBF.save_step_temp(step, geometry, command, copyLog, copyOrb)
            
            QMPrevious = copy.deepcopy(QMCalc)

        # the calculation has a MM part, collect AMBER output files
        if "M" in geometry.calculationType or "L" in geometry.calculationType:
            # amber at present is the only option for MM
            amber.save_MM_amber_step(step, command)

    ##############################################
    # PRINT PHYSICAL PROPERTIES AT THIS STEP
    ##############################################

    if command[1] == 'irc' and IRCconverged:
        # print formatted output on molecular geometry
        logwrt.startSubSection('IRC STEP {0} GEOMETRY'.format(step))
        logwrt.printGeom(geometry)

    if (command[1] != 'irc' or IRCconverged) and (command[1] != 'freqxg' or step == 0):
        # print formatted output on QM/MM energy/energies
        logwrt.startSubSection("QM/MM ENERGIES")
        logwrt.printEnergies(QMCalc.energy(), E_modelH, E_modelnoc, E_second, QMCalc.selfenergy,
                             QMMM_Results.getenergy(), geometry.calculationType)

    if ((command[1] != 'irc' or IRCconverged) and (command[1] != 'freqxg' or step == 0) and
            (command[1] != "mdv" or step == 0 or step == max_step - 1)):
        # print formatted output on electrostatic properties (dipole moment and modelH charges)
        if list(QMCalc.charges) and list(QMCalc.dipole):
            logwrt.startSubSection("ELECTROSTATIC PROPERTIES")
            if command[51] == "1" or command[51] == '10':  # for Gaussian, these are ESP charges
                logwrt.writelog("Electrostatic properties - ESP atom charges and dipole moment - are extracted \n" +
                                "from the QM calculation for the H part with H-saturated bonds\n\n")
            elif command[51] == "6" or command[51] == "7":  # for Molcas and Molpro, mulliken charges are extracted
                logwrt.writelog("Electrostatic properties - Mulliken atom charges and dipole moment - are extracted\n" +
                                "from the QM calculation for the H part with H-saturated bonds\n\n")
            logwrt.printModelHCharges(geometry, QMCalc.charges, QMCalc.dipole)

    if command[1] == 'mdv':
        # get the values of the current state and time step
        sef = shelve.open("cobram-sef")
        tStep = sef['tstep_OLD']
        actualstate = sef['state']
        sef.close()
        # at step zero, save the current value of potential that will be used to scale the energy
        # for the printout of MD info (the energy written to output is computed with respect to the initial energy)
        if step == 0: potReference = mdItems[1]
        # print info for molecular dynamics
        logwrt.startSubSection("MOLECULAR DYNAMICS")
        logwrt.printMDinfo(actualTime, actualstate, tStep, mdItems, potReference)

    ##############################################
    # PRINT INFO ON CONVERGENCE, WHEN APPLICABLE
    ##############################################

    if command[1] in ['optxg', 'ci', 'ts', 'irc'] and max_step > 1:  # for optimization print convTable
        logwrt.startSubSection("OPTIMIZATION CONVERGENCE")
        logwrt.writelog(convTable + "\n")

    ##############################################
    # CHECK CONDITIONS FOR LOOP TERMINATION
    ##############################################

    # for a single point calculation, always break loop
    if max_step == 1:
        break

    # for a calculation using GAUSSIAN optimizator, check the geometry.log file
    # and terminate loop when gaussian has ended
    elif command[1] in ['optxg', 'freqxg', 'irc', 'ci', 'ts']:
        # increment step counters (for IRC there are two counters!)
        if command[1] == "irc":
            if step == 0:
                step += 1
            elif IRCconverged:
                step += 1
                substep = 1
            else:
                substep += 1
        else:
            step += 1
        # check termination in geometry.log file, and in case break loop
        gauTerm = optxg.checkGauTermination("geometry.log")
        if gauTerm == 1:
            logwrt.writelog(
                "\n" + ' ' * 27 + '#' * 25 + '\n' + ' ' * 27 + "  CALCULATION COMPLETED\n" + ' ' * 27 + '#' * 25 + "\n\n")
        elif gauTerm == 2:
            logwrt.writelog("\n" + ' ' * 20 + '#' * 38 + '\n' + ' ' * 20 + "  CALCULATION TERMINATED WITH ERRORS\n" +
                            ' ' * 20 + '#' * 38 + "\n\n")
            with open('geometry.log') as geomLog:
                errors = geomLog.readlines()[-11:-3]

            gauFailMsg = 'Optimization terminated abnormally. Here is an excerpt of the geometry.log file:\n' \
                         '-----------------------------------------------------------------------------------\n'
            for line in errors:
                gauFailMsg += '>   ' + line.strip() + '\n'
            gauFailMsg += '-----------------------------------------------------------------------------------\n' \
                          'For further details please check geometry.log file. '
            logwrt.fatalerror(gauFailMsg)
        elif gauTerm == 3:
            logwrt.writelog("\n" + ' ' * 27 + '#' * 25 + '\n' + ' ' * 27 + "  CALCULATION COMPLETED\n" +
                            ' ' * 27 + '#' * 25 + "\n\n")
            convFailMsg = 'Optimization has not converged, number of steps exceeded.\n\n'
            logwrt.writewarning(convFailMsg)
        if gauTerm != 0:
            break

    # for molecular dynamics, run until the total number of steps is reached
    elif command[1] in ['mdv']:
        # increment step counter and actual time of the simulation
        step += 1
        actualTime += tStep
        # when the max number of step is reached, break the loop
        if step == max_step: break
    # flush standard output
    sys.stdout.flush()

##########################################################
# Finish
##########################################################

if command[51] == '11':
    filetext = QMMM_Results.sharcQMMMoutfile(QMCalc.outputData.get("outfile"), geometry)
    with open("QMMM.out", "w") as outf:
        outf.write(filetext)

# save last step data for calculations terminated before max_step (e.g. converged)
if int(command[100]) < 0 and command[1] != 'freqxg' and step != max_step - 1:
    # check if QM-MM data for last step have been already saved
    # replay conditions from "finishing the cycle" section, just in case
    if command[1] != 'irc' or IRCconverged or step == 0:
        # copyOrb conditions have already been defined before breaking the loop
        # if any of the conditions evaluate to True, we have already saved step data
        if not np.any([copyOrb_condition1, copyOrb_condition2, copyOrb_condition3]):
            # switch copyLog boolean value
            copyLog = not copyLog
            if NewQM:
                QMCalc.archiveStep(copyLog, True, step)
            else:
                # set copyMM to False, MM data is already saved after each step
                CBF.save_step_temp(step, geometry, command, copyLog, True)

# print final geometry, when it has changed
if max_step > 1 and command[1] != 'freqxg':
    logwrt.startSection('FINAL GEOMETRY')
    # print formatted output on molecular geometry
    logwrt.printGeom(geometry)

# for a freq calculation print table of normal modes taken from geometry.log
if command[1] == 'freqxg':
    logwrt.startSection('NORMAL MODES')
    logwrt.writelog(optxg.getNormalModes("geometry.log") + "\n")
    logwrt.writelog(
        "For the complete harmonic analysis please check the Gaussian output in the geometry.log file." + "\n")

    # write high precision normal modes in cobramm.xml
    HPNM = optxg.getNormalModes('geometry.log', high_precision=True)

    with Timer("xml output"):
        xmlfile.write_normal_modes(HPNM)

if not logwrt.DEBUG_COBRAMM_RUN: CBF.garbager(geometry, command)

# stop the timer for the main program, and print the report of the timings
totalTimer.stop()
logwrt.cobramend()
