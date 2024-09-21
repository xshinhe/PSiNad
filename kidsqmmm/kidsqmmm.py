#!/usr/bin/env python

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
from pprint import pprint, pformat
from traceback import format_exc
import sys  # system-specific parameters and functions
import copy  # shallow and deep copy operations
import shelve  # python object persistence

# check that COBRAMM is running on python 3
if sys.version_info < (3, 0): raise RuntimeError("Python 3 is required for running COBRAMM")

# imports of local modules

import argparse
import CBF
import amber
import optxg
import parallel_numerics
import logwrt  # manages log file output + start/end procedures
import inpdata  # input data read from files
import qmmmenv  # environmental variable for COBRAMM and 3rd-party software
import constants  # physical constants and conversion factors

# imports of local objects

from Layers import Layers  # Layers class to manage geometries
from Timer import Timer  # keep timings of the different sections of the code
from Charge import Charge  # Charge class that stores real and modelH charges
from Output import Output  # Output and Step classes to read/write cobramm.xml file
from QMCalc import QM  # QM class controls the QM calculation, define its input and stores the output
from QMMMCalc import QMMM  # join QM and MM data to construct QMMM results

# math libraries

import numpy as np  # numpy library for scientific computation

#####################################################################################################

parser = argparse.ArgumentParser(description='Execute QMMM Calculation')
parser.add_argument('-d', '--directory', dest='directory', nargs='?', default='.', type=str,
    help='work directory')
parser.add_argument('-i', '--input', dest='input', nargs='?', default='QMMM.in', type=str,
    help='input file')
parser.add_argument('-c', '--coord', dest='coord', nargs='?', default='real.crd', type=str,
    help='input file')
parser.add_argument('-t', '--task', dest='task', nargs='?', default='0', type=str,
    help='task type')
parser.add_argument('-o', '--output', dest='output', nargs='?', default='QMMM.log', type=str,
    help='output file')

if __name__ == "__main__":
    args = parser.parse_args()
    
    # source COBRAMM configuration file (if available)
    confFile = qmmmenv.setCobrammProfile()

    # store starting dir
    startdir = os.path.abspath(os.getcwd())
    print('####')
    print("startdir is ", startdir)
    print('rundir is ', args.directory)


    rundir = args.directory
    if rundir != '.':
        shutil.copy("real_layers.xyz", rundir)
        shutil.copy("real.top", rundir)
        shutil.copy("model-H.top", rundir)
        shutil.copy(args.input, rundir)
        if os.path.exists('box.info'):
            shutil.copy('box.info', rundir)

    layerfile_rel = os.path.relpath(startdir+'/real_layers.xyz', rundir)
    realtop_rel = os.path.relpath(startdir+'/real.top', rundir)
    modelHtop_rel = os.path.relpath(startdir+'/model-H.top', rundir)
    command_rel = os.path.relpath(args.input, rundir)
    realcrd_rel = os.path.relpath(args.coord, rundir)
    os.chdir(rundir)
    print('switch to rundir and list')
    print('the layer file is: ', layerfile_rel)
    print('the realtop is: ', realtop_rel)
    print('the modelHtop is: ', modelHtop_rel)
    print('the command is: ', command_rel)
    print('the realcrd is: ', realcrd_rel)

    # start COBRAMM with message to log, and start timer for total calculation time
    totalTimer = Timer("total")
    totalTimer.start()
    if os.path.exists('cobramm.xml'): os.remove('cobramm.xml')
    # logwrt.cobramstart()

    # check that COBRAMM environment is properly defined
    envDefined, errorMsg = qmmmenv.checkqmmmenv()
    if not envDefined: logwrt.fatalerror(errorMsg)
    # Print to log information on the environment configuration
    logwrt.writelog(qmmmenv.checkConfigFile(confFile))

    # start FILE CONTROL section
    logwrt.startSection("FILE CONTROL")

    # define starting input file
    fileinp_1 = layerfile_rel
    fileinp_2 = command_rel
    #fileinp_1 = layerfile_rel
    # check  for the presence of the two mandatory input files
    #for fname in [fileinp_1, fileinp_2]:
    #    if not os.path.isfile(fname): logwrt.fatalerror('mandatory input file {0} not found!'.format(fname))

    # make a cpy of the input files
    CBF.saveinputs()

    # initialize the list with the COBRAMM options
    commandhard = CBF.makehard()
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

    # set the level of verbosity of the log file
    logwrt.setverbositylevel(command[2])

    # read the list of the layers
    with open(fileinp_1, "r") as f:
        geometry = Layers.from_real_layers_xyz(f.read())
    geometry.updatereal(realcrd_rel)
    geometryfile = realcrd_rel
    geometry.makerealcrd() #?

    command[60] = 1
    max_step = int(command[60])
    smart_numerics = 0
    par_num = 1
    command[8] = 1
    command[1] = 'nad'

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
    logwrt.writelog("\n\nMax number of steps is set to {} for this calculation\n\n".format(max_step))

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
    # logwrt.printGeom(geometry)

    logwrt.startSection('ATOMIC CHARGES')


    CRG_real = amber.prepare(geometry, cobcom, command)
    CRG_model_H = amber.read_crgA()
    print('inital real charges:', CRG_real)
    print('inital model charges:', CRG_model_H)
    if os.path.exists('laststep.charge'):
        lines = open('laststep.charge', 'r').readlines()
        CRG_model_H = [float(l) for l in lines]
        os.remove('laststep.charge')
    # create object Charge to store model and embedding charges
    charges = Charge(geometry, CRG_real, CRG_model_H)
    charges.checkConsistency()

    # ####################################################
    #                    write xml file
    # ####################################################

    # # get cobramm version
    version = qmmmenv.getVersion()
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
    elif command[51] == '10086':
        # MNDO
        restartFileName = ['fort.7']
    elif command[51] == '10087':
        # BAGEL
        restartFileName = ['laststep.archive']
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

    # initialization of variables
    displQM = []
    QMPrevious = None
    prevx1, prevx2 = None, None

    while True:
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
            # E_second is real system
            # E_modelnoc is real system with no charge on modelH (==? E_second)
            # E_modelH is model H system with no charge
            Fxyz_second, E_second, E_modelnoc, Fxyz_modelnoc, Fxyz_modelH, E_modelH = CBF.MM(
                step, geometry, command, cobcom)
            MM_Results = [Fxyz_second, E_second, E_modelnoc, Fxyz_modelnoc, Fxyz_modelH, E_modelH]
            print("MM Force Norm (real):", np.sum(np.abs(np.array(Fxyz_second))))
            print("MM Force Norm (real noc model):", np.sum(np.abs(np.array(Fxyz_modelnoc))))
            print("MM Force Norm (modelH):", np.sum(np.abs(np.array(Fxyz_modelH))))

        sef = shelve.open("cobram-sef")
        sef['Fxyz_modelH'] = Fxyz_modelH
        sef.close()

        ##########################################################
        #                 QM calculations
        ##########################################################

        with Timer("QM section"):

            # define the type of QM calculation: old style interface or new OOP interface
            if command[51] == '1' or command[51] == '11' or command[51] == '6' or command[51] == '10086' or command[51] == '10087':
                NewQM = True # other should used in NewQM = False, the old interface
        
                #for normal run (no DEBUG run), clean QM directories to avoid mess:
                if not logwrt.DEBUG_COBRAMM_RUN:
                    QM_DIRECTORIES = [item for item in os.listdir() if os.path.isdir(item) and item.startswith('qmCalc')]
                    if QM_DIRECTORIES:
                        for qm_dir in QM_DIRECTORIES:
                            shutil.rmtree(qm_dir)

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
                    
                # initialize QM instance to define input of the QM calculation
                QMCalc = QM(cobcom, command, geometry, charges, step, QMRestart, **stateargs) #--> F: added step as attribute to QM class in order to build correct Molcas input
                # run the calculation
                QM.runQM(QMCalc, memory=command[53], nprocQM=int(command[7])), logwrt.writelog(QMCalc.log)
                
                # #print gradient 
                # if command[1] != 'irc' and command[1] != 'nad' and QMCalc.gradient(QMCalc.outputData.get("optstate")) is not None:
                #     logwrt.writelog("\nQM gradient of High Layer (+ link atoms) for state {0:2d}\n".format(QMCalc.outputData.get("optstate")+1), 2)
                #     if geometry.NsubH == 0:
                #         logwrt.writelog(logwrt.matrix_prettystring(np.array(QMCalc.gradient(QMCalc.outputData.get("optstate"))).T, ".6f", geometry.getAtomLabels("HIGH")), 2)
                #     else:
                #         logwrt.writelog(logwrt.matrix_prettystring(np.array(QMCalc.gradient(QMCalc.outputData.get("optstate"))).T, ".6f", geometry.getAtomLabels("modelH")), 2)
                #     logwrt.writelog("\n", 2)
                # if command[1] == 'nad' and QMCalc.gradient(QMCalc.outputData.get("optstate")) is not None:
                #     logwrt.writelog("\nQM gradient of High Layer (+ link atoms) for all state\n")
                #     for i in range(2):
                #         if geometry.NsubH == 0:
                #             logwrt.writelog(pformat(QMCalc.gradient(i)))
                #         else:
                #             logwrt.writelog(logwrt.matrix_prettystring(np.array(QMCalc.gradient(i)).T, ".6f", geometry.getAtomLabels("modelH")), 2)
                #         logwrt.writelog("\n", 2)
                # if command[1] == 'ci' and (QMCalc.gradient(QMCalc.outputData.get("optstate")-1) is not None):
                #     lowerstate = QMCalc.outputData.get("optstate") - 1
                #     logwrt.writelog("\nQM gradient of High Layer (+ link atoms) for state {0:2d}\n".format(lowerstate+1), 2)
                #     if geometry.NsubH == 0:
                #         logwrt.writelog(logwrt.matrix_prettystring(np.array(QMCalc.gradient(lowerstate)).T, ".6f", geometry.getAtomLabels("HIGH")), 2)
                #     else:
                #         logwrt.writelog(logwrt.matrix_prettystring(np.array(QMCalc.gradient(lowerstate)).T, ".6f", geometry.getAtomLabels("modelH")), 2)
                #     logwrt.writelog("\n", 2)

                #print NAC
                # if (command[1] == 'nad') or ((command[1] == 'mdv') and (command[85] == '1') and (command[14] != '1')):
                #     NAC = QMCalc.outputData.get("nac")
                #     pprint(NAC)
                #     for i in NAC:
                #         for j in NAC[i]:
                #             if j > i:
                #                 logwrt.writelog("\nQM NAC of High Layer (+ link atoms) for states {0:2d} {1:2d} (NAC_ij = -NAC_ji)\n".format(i+1, j+1), 2)
                #                 if geometry.NsubH == 0:
                #                     pprint(NAC[i][j]) # because the old code, nac is only on modelH atoms!

                #                     # logwrt.writelog(logwrt.matrix_prettystring(np.array(NAC[i][j]).T, ".6f", geometry.getAtomLabels("HIGH")), 2)
                #                 else:
                #                     pprint(NAC[i][j]) # because the old code, nac is only on modelH atoms!
                #                     # logwrt.writelog(logwrt.matrix_prettystring(np.array(NAC[i][j]).T, ".6f", geometry.getAtomLabels("modelH")), 2)
                #                 logwrt.writelog("\n", 2)

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
                QMRestart = QMCalc.saveRestartFile()


        ##########################################################
        #         Rebuilding the the charges of the system
        ##########################################################

        # TODO: with the current algorithm, this command update the definition of the charges
        #  only in the instance of the charge object
        #  we should consider a re-implementation of the update of the topology file with QMMM charges,
        #  to work with the microiterative sheme, but the current order of QM and MM should be rethought
        #  since one wants to to the MM dynamics with the values of the charge of the same step, not of the
        #  step before that
        # print(QMCalc.charges)
        if not list(QMCalc.charges):
            logwrt.writewarning("QM charges are not available: their value will not be updated!\n")
        else:
            charges.rallyCharges(QMCalc.charges)
            f = open('laststep.charge', 'w')
            for i in range(len(QMCalc.charges)):
                f.write('{: 12.8e}\n'.format(QMCalc.charges[i]))
            f.flush()
            f.close()

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
            # if command[1] != 'irc' and geometry.calculationType in ['HL', 'HM', 'HML'] and QMCalc.gradient(QMCalc.outputData.get("optstate")) is not None:
            #     # print High layer gradient
            #     logwrt.writelog("\nQM/MM gradient of High Layer for state {0:2d}\n".format(newstate+1), 1)
            #     gradientH=[]
            #     for iat, ixyz in enumerate(np.array(QMMM_Results.getgradient(newstate)).T):
            #         if iat+1 in geometry.list_HIGH:
            #             gradientH.append(ixyz)
            #     logwrt.writelog(logwrt.matrix_prettystring(np.array(gradientH), ".6f", geometry.getAtomLabels("HIGH")), 1)
            #     logwrt.writelog("\n", 1)

            #     # print High+Medium layer gradient
            #     logwrt.writelog("\nQM/MM gradient of High+Medium Layer for state {0:2d}\n".format(newstate+1), 2)
            #     logwrt.writelog(logwrt.matrix_prettystring(np.array(QMMM_Results.getgradient(newstate)).T, ".6f", geometry.getAtomLabels("MEDIUM_HIGH")), 2)
            #     logwrt.writelog("\n", 2)
            # if command[1] == 'ci' and (geometry.calculationType in ['HL', 'HM', 'HML']) and (QMCalc.gradient(QMCalc.outputData.get("optstate") -1 ) is not None):
            #     lowerstate = QMCalc.outputData.get("optstate") - 1
            #     # print High layer gradient
            #     logwrt.writelog("\nQM/MM gradient of High Layer for state {0:2d}\n".format(lowerstate+1), 1)
            #     gradientH=[]
            #     for iat, ixyz in enumerate(np.array(QMMM_Results.getgradient(lowerstate)).T):
            #         if iat+1 in geometry.list_HIGH:
            #             gradientH.append(ixyz)
            #     logwrt.writelog(logwrt.matrix_prettystring(np.array(gradientH), ".6f", geometry.getAtomLabels("HIGH")), 1)
            #     logwrt.writelog("\n", 1)
        
            #     # print High+Medium layer gradient
            #     logwrt.writelog("\nQM/MM gradient of High+Medium Layer for state {0:2d}\n".format(lowerstate+1), 2)
            #     logwrt.writelog(logwrt.matrix_prettystring(np.array(QMMM_Results.getgradient(lowerstate)).T, ".6f", geometry.getAtomLabels("MEDIUM_HIGH")), 2)
            #     logwrt.writelog("\n", 2)

            # define gradient for CI optimization or for a regular calculation
            if command[1] != 'ci':
                NewGrad = QMMM_Results.getgradient(newstate)
            else:
                NewGrad = QMMM_Results.cigradient

        # initialize variable to store a new geometry for update (when None, do not update geometry)
        newgeom = None

        ##########################################################
        #           write data in xml file
        ##########################################################

        stat_number = 0
        with open('interface.ds', 'w') as f:
            # write status 
            f.write('interface.stat\n')
            f.write('kids_int 1\n')
            f.write(f'{stat_number}\n\n')

            # write energy
            f.write('interface.eig\n')
            f.write('kids_real %d\n'%len(QMMM_Results.energies))
            for i in range(len(QMMM_Results.energies)): # sorted order
                f.write('{: 12.8e}\n'.format(QMMM_Results.energies[i]))
            f.write('\n')

            # write energy
            f.write('interface.dE\n')
            f.write('kids_real %d\n'%(len(QMMM_Results.energies) * geometry.atomNum*3))
            jHM = 0 # count for H & M atoms
            for i in range(geometry.atomNum): # sorted order
                if i+1 in geometry.list_MEDIUM_HIGH:
                    for ix in [0,1,2]:
                        for k in range(len(QMMM_Results.energies)):
                            f.write('{: 12.8e} '.format(QMMM_Results.gradient[k][ix][jHM]))
                        f.write('\n')
                    jHM += 1
                if i+1 in geometry.list_LOW:
                    for ix in [0,1,2]:
                        for k in range(len(QMMM_Results.energies)):
                            f.write('{: 12.8e} '.format(0))
                        f.write('\n')
            f.write('\n')

            # write nac
            f.write('interface.nac\n')
            f.write('kids_real %d\n'%(len(QMMM_Results.energies)*len(QMMM_Results.energies)*geometry.atomNum*3) )
            jHM = 0 # count for H & M atoms
            for i in range(geometry.atomNum): # sorted order
                if i+1 in geometry.list_MEDIUM_HIGH:
                    for ix in [0,1,2]:
                        for k1 in range(len(QMMM_Results.energies)):
                            for k2 in range(len(QMMM_Results.energies)):
                                if k2 == k1:
                                    f.write('{: 12.8e} '.format(0))
                                else:
                                    f.write('{: 12.8e} '.format(QMMM_Results.nac[k1][k2][ix][jHM]))
                        f.write('\n')
                    jHM += 1
                if i+1 in geometry.list_LOW:
                    for ix in [0,1,2]:
                        for k1 in range(len(QMMM_Results.energies)):
                            for k2 in range(len(QMMM_Results.energies)):
                                f.write('{: 12.8e} '.format(0))
                        f.write('\n')
            f.write('\n')

            # write ocillation strength
            f.write('interface.strength\n')
            f.write('kids_real %d\n'%len(QMMM_Results.energies))
            for i in range(len(QMMM_Results.energies)): # sorted order
                if i==0:
                    f.write('{: 12.8e}\n'.format(0))
                else:
                    f.write('{: 12.8e}\n'.format(QMCalc.outputData.dataDict["osc_strength"][i]))
            f.write('\n')

            f.flush()
            f.close()

        # pprint(QMCalc.energy())
        # pprint(QMMM_Results.getenergy())
        # pprint(QMMM_Results.getgradient(0))
        # pprint(QMMM_Results.getgradient(1))

        # preparing arguments for write_step
        # values should NOT be of type string
        energies = [E_second, E_modelnoc, E_modelH, QMCalc.selfenergy, QMCalc.energy(), QMMM_Results.getenergy()]
        dyn = None
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

        ##############################################
        # PRINT PHYSICAL PROPERTIES AT THIS STEP
        ##############################################

        # print formatted output on QM/MM energy/energies
        logwrt.startSubSection("QM/MM ENERGIES")
        logwrt.printEnergies(QMCalc.energy(), E_modelH, E_modelnoc, E_second, QMCalc.selfenergy,
                             QMMM_Results.getenergy(), geometry.calculationType)

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

        ##############################################
        # CHECK CONDITIONS FOR LOOP TERMINATION
        ##############################################

        # for a single point calculation, always break loop
        if max_step == 1:
            break

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
                QMCalc.archiveStep(copyLog, True, step)

    if not logwrt.DEBUG_COBRAMM_RUN: CBF.garbager(geometry, command)
    # stop the timer for the main program, and print the report of the timings
    totalTimer.stop()
    logwrt.cobramend()

    os.chdir(startdir)

