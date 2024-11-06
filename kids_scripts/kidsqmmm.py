#!/usr/bin/env python3
#   Coding=utf-8

#   KIDS SCRIPTS
#   Author: xshinhe
#   
#   Copyright (c) 2024 PeKing Univ. - GNUv3 License

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
import sys   # system-specific parameters and functions
import copy  # shallow and deep copy operations
import shelve  # python object persistence
from pprint import pprint, pformat
from traceback import format_exc

if sys.version_info < (3, 0): raise RuntimeError("Python 3 is required")

# imports of local modules (utilities)
import softenv  # environmental variable for COBRAMM and 3rd-party software
import constants  # physical constants and conversion factors
import kids_arg
import kids_log
import kids_io
from kids_log import Timing
from kids_config import Config

# imports of local modules (classes)
from kids_log import Timing, Log    
from Layers import Layers  # Layers class to manage geometries
from Charge import Charge  # Charge class that stores real and modelH charges
from Output import Output  # Output and Step classes to read/write old.xml file
from MMCalc import MM  # QM class controls the QM calculation, define its input and stores the output
import MMCalc  # QM class controls the QM calculation, define its input and stores the output
from QMCalc import QM  # QM class controls the QM calculation, define its input and stores the output
from QMMMCalc import QMMM  # join QM and MM data to construct QMMM results

# math libraries

import numpy as np  # numpy library for scientific computation

#####################################################################################################

if __name__ == "__main__":

    # parse ks_config from command-line arguments & configuration file
    args = kids_arg.parser.parse_args()

    # setup verbosity
    Log.setVerbosityLevel(args.verbosity) # set the level of verbosity of the log file
    Log.start()
    # begin timing
    totalTiming = Timing("Total")

    # check environment
    Log.startSection("ENV CONTROL")
    profile_file = softenv.setKIDSProfile()
    envDefined, errorMsg = softenv.checkKIDSEnv()
    if not envDefined: Log.fatalError(errorMsg)
    Log.writeLog(softenv.checkKIDSProfile(profile_file))

    # store directory information
    startdir = os.path.abspath(os.getcwd())
    rundir = args.directory
    Log.writeLog(f'Start directory is {startdir}\n')
    Log.writeLog(f'Running directory is {rundir}\n')
    if not os.path.isdir(rundir):
        Log.writeLog(f'Making running directory {rundir}\n')
        os.mkdir(rundir)

    # ####################################################
    #          print summary of input options
    # ####################################################
    Log.startSection('INPUT MOLECULAR DESCRIPTION')
    ks_config = Config.load(args.input, args)

    geometry = ks_config['_geom']
    Log.writeLog("\n{0} is requested, with {1}.\n".format(kids_io.getCalcType(args.type), 
        kids_io.getLayers(geometry)))
    if "H" in geometry.calculationType:
        Log.writeLog("QM third party software : {0}\n".format(args.qmsolver))
    if "M" in geometry.calculationType or "L" in geometry.calculationType:
        Log.writeLog("MM third party software : {0}\n".format(args.mmsolver))
    Log.writeLog("\n")

    # Log.writeLog(ks_config.summary())

    # print formatted output on layer and molecular geometry definition
    Log.printLayers(geometry.list_HIGH, geometry.list_MEDIUM, geometry.list_LOW, geometry.atomLink_BA)
    Log.writeLog('Calculation type is {0}, '.format(geometry.calculationType))
    if geometry.calculationType in ['HL', 'HM', 'HML']:
        Log.writeLog('A QM/MM calculation is then requested\n\n')
    if geometry.calculationType in ['ML', 'M']:
        Log.writeLog('A MM calculation is then requested\n\n')
    if geometry.calculationType in ['H']:
        Log.writeLog('A QM calculation is then requested\n\n')
    Log.printGeom(geometry)

    ############################################################
    #          Define QM restart file
    ############################################################

    # set the name of the orbital restart file depending on the type of QM
    if ks_config.args.qmsolver == 'gaussian':
        #Gaussian
        restartFileName = ["gaussian-QM.chk", "gaussian.chk"]
    elif ks_config.args.qmsolver == 'molcas':
        #Molcas
        #files are in order of increasing precedence (molcas.JobIph has higher priority)
        restartFileName = ["molcas.RasOrb", "INPORB", "molcas.JobIph"]
    elif ks_config.args.qmsolver == 'mndo':
        # MNDO
        restartFileName = ['fort.7']
    elif ks_config.args.qmsolver == 'bagel':
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
    charges = []
    mm_results = [None, 0, 0, None, None, 0]
    E_real = 0
    E_modelnoc = 0
    E_modelH = 0
             
    qmcalc     = None
    QMPrevious = None
    prevx1, prevx2 = None, None


    Log.writeLog(f'Changing from {startdir} to {rundir}\n')
    os.chdir(rundir)   
    Log.writeLog(f'Current directory: {os.getcwd()}\n')         
    Log.startSection('QM or QM/MM SINGLE POINT')
    ############################################################
    #                  MM calculations
    ############################################################
    if 'M' in geometry.calculationType or 'L' in geometry.calculationType:        
        with Timing("MM section"):            
            # check files
            if not ks_config.geom_in_toml:
                if not os.path.exists(os.path.join(startdir, ks_config.topo1)):
                    Log.fatalError('Need realtop')
                if not os.path.exists(os.path.join(startdir, ks_config.topo2)):
                    Log.fatalError('Need modelHtop')
                if not os.path.exists(os.path.join(startdir, ks_config.layer)):
                    Log.fatalError('Need layerinfo')  

            if rundir != '.':
                if not ks_config.geom_in_toml:
                    shutil.copy(os.path.join(startdir, ks_config.topo1), '.')
                    shutil.copy(os.path.join(startdir, ks_config.topo2), '.')
                    shutil.copy(os.path.join(startdir, ks_config.layer), '.')
                    shutil.copy(os.path.join(startdir, args.coord), '.')

                shutil.copy(os.path.join(startdir, args.input), '.')

                if os.path.exists(os.path.join(startdir, 'box.info')):
                    shutil.copy(os.path.join(startdir, 'box.info'), '.')

            # Calculate relative paths
            layerfile_rel = os.path.relpath(os.path.join(startdir, ks_config.layer), '.')
            realtop_rel = os.path.relpath(os.path.join(startdir, ks_config.topo1), '.')
            modelHtop_rel = os.path.relpath(os.path.join(startdir, ks_config.topo2), '.')
            ks_config_rel = os.path.relpath(os.path.join(startdir, args.input), '.')
            realcrd_rel = os.path.relpath(os.path.join(startdir, args.coord), '.')

            Log.writeLog('\nInfo files:\n')
            Log.writeLog(f'the layer file is: {layerfile_rel}\n')
            Log.writeLog(f'the realtop is: {realtop_rel}\n')
            Log.writeLog(f'the modelHtop is: {modelHtop_rel}\n')
            Log.writeLog(f'the ks_config is: {ks_config_rel}\n')
            Log.writeLog(f'the realcrd is: {realcrd_rel}\n')
            kids_io.saveinputs() # make a cpy of the input files

            Log.startSection('ATOMIC CHARGES')
            CRG_real, CRG_model_H = MMCalc.prepareCRG(geometry, ks_config)
            if os.path.exists('laststep.charge'):
                lines = open('laststep.charge', 'r').readlines()
                CRG_model_H = [float(l) for l in lines]
                os.remove('laststep.charge')
            # create object Charge to store model and embedding charges
            charges = Charge(geometry, CRG_real, CRG_model_H)
            charges.checkConsistency()

            mm_results = MMCalc.MM(geometry, ks_config)
            Fxyz_real, E_real, E_modelnoc, Fxyz_modelnoc, Fxyz_modelH, E_modelH = mm_results
            # E_real is real system
            # E_modelnoc is real system with no charge on modelH (==? E_real)
            # E_modelH is model H system with no charge

            print("MM Force Norm (real):", np.sum(np.abs(np.array(Fxyz_real))))
            print("MM Force Norm (real noc model):", np.sum(np.abs(np.array(Fxyz_modelnoc))))
            print("MM Force Norm (modelH):", np.sum(np.abs(np.array(Fxyz_modelH))))

    ##########################################################
    #                 QM calculations
    ##########################################################
    if 'H' in geometry.calculationType:
        with Timing("QM section"):
    
            #for normal run (no DEBUG run), clean QM directories to avoid mess:
            # if not Log.DEBUG_COBRAMM_RUN:
            #     QM_DIRECTORIES = [item for item in os.listdir() if os.path.isdir(item) and item.startswith('qmCalc')]
            #     if QM_DIRECTORIES:
            #         for qm_dir in QM_DIRECTORIES:
            #             shutil.rmtree(qm_dir)

            # for a SH calculation, we need to select the initial state and force-set the elect state in step > 0
            stateargs = {}
            if int(ks_config.get_nested('_dyn.actstate', 0)) > 0:
                if step == 0:  # if this is the first step, define the initial state using ks_config[13]
                    if ks_config[13] != "1":
                        stateargs["setState"] = int(ks_config[13])-1
                else:  # if this is not the first step, read the actual state from sef
                    sef = shelve.open('cobram-sef')
                    stateargs["forceState"] = sef['state']
                    sef.close()
            elif ks_config.args.qmsolver == '6': #Molcas needs to know opt state
                stateargs["setState"] = int(ks_config[13])-1
                
            # initialize QM instance to define input of the QM calculation
            qmcalc = QM(ks_config, geometry, charges, **stateargs)
            QM.runQM(qmcalc, 
                memory=ks_config.get_nested('QM.mem', '500MB'), 
                nprocQM=int(ks_config.get_nested('QM.ncore', 1)))
            Log.writeLog(qmcalc.log)

            Log.writeLog('%d\n'%sys._getframe().f_lineno, 2)

            #print overlap matrix (in case of mdv with time-derivatice couplings)
            # if (ks_config[1] == 'mdv') and (ks_config.get_nested('_dyn.actstate', 0) == '1') and (ks_config[14] == '1') and ks_config.args.qmsolver == '6' and step != 0:
            #     Log.writeLog("Overlap matrix with states at previous time step:\n",1)
            #     OvMat = QMCalc.outputData.dataDict["psioverlap"]
            #     Log.writeLog(Log.matrix_prettystring(OvMat, ".8f"), 1)
            #     Log.writeLog("\n", 1)

            # when this is the first step, initialize the number of the electronic state
            # if step == 0:
            #     sef = shelve.open('cobram-sef')
            #     sef['state'] = QMCalc.outputData.get("optstate")
            #     sef['newstate'] = QMCalc.outputData.get("optstate")
            #     sef.close()

            # copy and save path of the orbital restart file for the next step
            QMRestart = qmcalc.saveRestartFile()


        ##########################################################
        #         Rebuilding the the charges of the system
        ##########################################################
        # print(QMCalc.charges)
        if not list(qmcalc.charges):
            Log.writewarning("QM charges are not available: their value will not be updated!\n")
        else:
            charges.rallyCharges(qmcalc.charges)
            f = open('laststep.charge', 'w')
            for i in range(len(qmcalc.charges)):
                f.write('{: 12.8e}\n'.format(qmcalc.charges[i]))
            f.flush()
            f.close()

    ##########################################################
    #         COMPUTING QMMM RESULTS
    ##########################################################
    if True:
        with Timing("QMMM section"):
            global qmmm_results
            qmmm_results = QMMM(ks_config, geometry, qmcalc, mm_results, step, prevx1, prevx2)

            #save QMMM branching plane vectors of previous step in case of ci search with new branching plane
            # if ks_config[19] == '1' and ks_config[1] == 'ci':
                # prevx1, prevx2 = copy.deepcopy(qmmm_results.x1), copy.deepcopy(qmmm_results.x2)

            # getting actual electronic state number from shelve
            # sef = shelve.open('cobram-sef')
            # actualstate = sef['state']
            # if step > 0:
            #     newstate = sef['newstate']
            # else:
            #     newstate = sef['state']
            # sef.close()
            # if ks_config[1] != 'ci':
            #     NewGrad = qmmm_results.getgradient(newstate)
            # else:
            #     NewGrad = qmmm_results.cigradient

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
            f.write('kids_real %d\n'%len(qmmm_results.energies))
            for i in range(len(qmmm_results.energies)): # sorted order
                f.write('{: 12.8e}\n'.format(qmmm_results.energies[i]))
            f.write('\n')

            # write energy
            f.write('interface.dE\n')
            f.write('kids_real %d\n'%(len(qmmm_results.energies) * geometry.atomNum*3))
            jHM = 0 # count for H & M atoms
            for i in range(geometry.atomNum): # sorted order
                if i+1 in geometry.list_MEDIUM_HIGH:
                    for ix in [0,1,2]:
                        for k in range(len(qmmm_results.energies)):
                            f.write('{: 12.8e} '.format(qmmm_results.gradient[k][ix][jHM]))
                        f.write('\n')
                    jHM += 1
                if i+1 in geometry.list_LOW:
                    for ix in [0,1,2]:
                        for k in range(len(qmmm_results.energies)):
                            f.write('{: 12.8e} '.format(0))
                        f.write('\n')
            f.write('\n')

            # write nac
            f.write('interface.nac\n')
            f.write('kids_real %d\n'%(len(qmmm_results.energies)*len(qmmm_results.energies)*geometry.atomNum*3) )
            jHM = 0 # count for H & M atoms
            for i in range(geometry.atomNum): # sorted order
                if i+1 in geometry.list_MEDIUM_HIGH:
                    for ix in [0,1,2]:
                        for k1 in range(len(qmmm_results.energies)):
                            for k2 in range(len(qmmm_results.energies)):
                                if k2 == k1:
                                    f.write('{: 12.8e} '.format(0))
                                else:
                                    f.write('{: 12.8e} '.format(qmmm_results.nac[k1][k2][ix][jHM]))
                        f.write('\n')
                    jHM += 1
                if i+1 in geometry.list_LOW:
                    for ix in [0,1,2]:
                        for k1 in range(len(qmmm_results.energies)):
                            for k2 in range(len(qmmm_results.energies)):
                                f.write('{: 12.8e} '.format(0))
                        f.write('\n')
            f.write('\n')

            # write ocillation strength
            f.write('interface.strength\n')
            f.write('kids_real %d\n'%len(qmmm_results.energies))
            for i in range(len(qmmm_results.energies)): # sorted order
                if i==0:
                    f.write('{: 12.8e}\n'.format(0))
                else:
                    f.write('{: 12.8e}\n'.format(qmcalc.outputData.dataDict["osc_strength"][i]))
            f.write('\n')

            f.flush()
            f.close()

        # pprint(QMCalc.energy())
        # pprint(qmmm_results.getenergy())
        # pprint(qmmm_results.getgradient(0))
        # pprint(qmmm_results.getgradient(1))

        # preparing arguments for write_step
        # values should NOT be of type string
        energies = [E_real, E_modelnoc, E_modelH, qmcalc.selfenergy, qmcalc.energy(), qmmm_results.getenergy()]

        ##########################################################
        #                finishing the cycle
        ##########################################################

        ##############################################
        # PRINT PHYSICAL PROPERTIES AT THIS STEP
        ##############################################

        # print formatted output on QM/MM energy/energies
        Log.startSubSection("QM/MM ENERGIES")
        Log.printEnergies(qmcalc.energy(), E_modelH, E_modelnoc, E_real, qmcalc.selfenergy,
                             qmmm_results.getenergy(), geometry.calculationType)

        # print formatted output on electrostatic properties (dipole moment and modelH charges)
        if list(qmcalc.charges) and list(qmcalc.dipole):
            Log.startSubSection("ELECTROSTATIC PROPERTIES")
            if ks_config.args.qmsolver == "1" or ks_config.args.qmsolver == '10':  # for Gaussian, these are ESP charges
                Log.writeLog("Electrostatic properties - ESP atom charges and dipole moment - are extracted \n" +
                                "from the QM calculation for the H part with H-saturated bonds\n\n")
            elif ks_config.args.qmsolver == "6" or ks_config.args.qmsolver == "7":  # for Molcas and Molpro, mulliken charges are extracted
                Log.writeLog("Electrostatic properties - Mulliken atom charges and dipole moment - are extracted\n" +
                                "from the QM calculation for the H part with H-saturated bonds\n\n")
            Log.printModelHCharges(geometry, qmcalc.charges, qmcalc.dipole)

    ##########################################################
    # Finish
    ##########################################################

    if ks_config.args.qmsolver == '11':
        filetext = qmmm_results.sharcQMMMoutfile(QMCalc.outputData.get("outfile"), geometry)
        with open("QMMM.out", "w") as outf:
            outf.write(filetext)

    if not Log.DEBUG_COBRAMM_RUN: kids_io.garbager(geometry, ks_config)
    # stop the timer for the main program, and print the report of the timings
    Log.end()

    os.chdir(startdir)

