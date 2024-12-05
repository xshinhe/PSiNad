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

# import statements of module from python standard library

import math
import sys  # system-specific parameters and functions
import os  # operating system utilities
from os import system
import shutil
import copy  # shallow and deep copy operations
from fileinput import *
import shelve  # python object persistence
import subprocess

# imports of local modules

import CBF
import molpro
from tullyNEW import Tully
import parallel_numerics
import product
import logwrt  # manages log file output + start/end procedures
import cobrammenv  # environmental variable for COBRAMM and 3rd-party software
import constants  # values of common physical constants
from timer import Timer  # keep timings of the different sections of the code

# math libraries

import numpy as np  # numpy library for scientific computation

#####################################################################################################

MOLCAS_WRKDIR = "MOLCAS"
QM_DATA_STORAGE = 'QM_data'
SINGLEPOINT_DIR = "single_points"


#####################################################################################################


def prepare(command, step):
    # check that the environment to run molcas is properly defined
    envDefined, errorMsg = cobrammenv.checkMolcasEnv()
    if not envDefined: logwrt.fatalerror(errorMsg)

    if step == 0:

        # print message to log, and link molcas script in the working directory
        molcasscript = os.getenv('MOLCAS_SCRIPT')  # path of the MOLCAS launch script
        molcaspath = os.getenv('MOLCAS')  # path of the MOLCAS code
        logwrt.writelog('Using molcas version {} and molcas launch script {}\n'.format(molcaspath, molcasscript))

        # create temporary directories to store MOLCAS calculation files
        if os.path.exists(MOLCAS_WRKDIR): shutil.rmtree(MOLCAS_WRKDIR)
        if os.path.exists(QM_DATA_STORAGE): shutil.rmtree(QM_DATA_STORAGE)
        os.mkdir(MOLCAS_WRKDIR)
        os.mkdir(QM_DATA_STORAGE)

        # environmental variables for MOLCAS execution
        os.putenv('Project', 'molcas')
        os.putenv('MOLCAS_MEM', str(command[53]))
        os.putenv('WorkDir', os.path.join(os.getcwd(), MOLCAS_WRKDIR))
        os.putenv('MOLCAS_OUTPUT', os.path.join(os.getcwd(), MOLCAS_WRKDIR))

        # if a user supply the INPORB or a molcas.RasOrb file, it will be used (INPORB is preferred!)
        if os.path.exists('INPORB'):
            shutil.copy("INPORB", os.path.join(MOLCAS_WRKDIR, "molcas.RasOrb"))
            logwrt.writelog('INPORB found! It will be used to set the initial orbitals.\n')
        elif os.path.exists('molcas.RasOrb'):
            shutil.copy("molcas.RasOrb", os.path.join(MOLCAS_WRKDIR, "molcas.RasOrb"))
            logwrt.writelog('molcas.RasOrb found! It will be used to set the initial orbitals.\n')

    if step != 0 and command[1] in ['freqxg', 'optxg', 'mdv', 'irc', 'ci', 'ts']:

        # start from a clean work directory
        if os.path.isdir("step{}".format(step-1)):
            shutil.copytree(MOLCAS_WRKDIR,"step{}_bis".format(step-1))
        else:
            shutil.copytree(MOLCAS_WRKDIR,"step{}".format(step-1))
        shutil.rmtree(MOLCAS_WRKDIR)
        os.mkdir(MOLCAS_WRKDIR)

        # read available orbitals from INPORB or molcas.RasOrb file
        if os.path.exists('molcas.JobIph'):
            shutil.copy("molcas.JobIph", os.path.join(MOLCAS_WRKDIR, "molcas.JobIph"))
            logwrt.writelog('molcas.JobIph found! It will be used to set the initial wavefunction.\n')
        elif os.path.exists('molcas.RasOrb'):
            shutil.copy("molcas.RasOrb", os.path.join(MOLCAS_WRKDIR, "molcas.RasOrb"))
            logwrt.writelog('molcas.RasOrb found! It will be used to set the initial wavefunction.\n')
        elif os.path.exists('INPORB'):
            shutil.copy("INPORB", os.path.join(MOLCAS_WRKDIR, "molcas.RasOrb"))
            logwrt.writelog('INPORB found! It will be used to set the initial orbitals.\n')

    gen = []
    keys = []
    return [keys, gen]


#####################################################################################################

@Timer("makeinp")
def makeinp(charges, geometry, command, step, cobcom, par_num):
    sef = shelve.open("cobram-sef", 'r')
    if step != 0: nroots = sef['nroots']
    if command[1] == 'mdv' and par_num == 0 and int(command[85]) > 0 and step != 0:
        calc_coupl = sef['calc_coupl']
    DEarray = sef['DEarray']
    sef.close()

    state = None
    ciro_key = 0
    num_key = 0
    caspt2_key = 0

    if step == 0:
        # convention for state labeling in cobram:  0=GS, 1=1st ES, 2=2nd ES"
        # Molcas convention is different: 1=GS, 2=1st ES, 3=2nd ES"
        nroots = 1
        keys = CBF.ReadCobramCommand(cobcom, 'molcas', 'molcas.input')
        for i in range(len(keys)):
            if keys[i].lower()[0:4] == 'ciro':
                try:
                    nroots = int(keys[i + 1].split()[0])
                    ciro_key = 1
                    DEarray = []
                except:
                    nroots = int(keys[i].split('=')[1].split()[0])
                    ciro_key = 1
                    DEarray = []
            if keys[i].lower()[0:4] == 'mult' or keys[i].lower()[0:4] == 'xmul' or keys[i].lower()[0:4] == 'rmul':
                try:
                    nroots = int(keys[i + 1].split()[0])
                    ciro_key = 1
                    DEarray = []
                except:
                    nroots = int(keys[i].split('=')[1].split()[0])
                    ciro_key = 1
                    DEarray = []
            if keys[i].strip().lower()[0:7] == '&alaska':
                for j in range(i + 1, len(keys)):
                    try:
                        if keys[j].strip()[0:1] == '&':
                            break
                        if keys[j].strip().lower()[0:4] == 'nume':
                            num_key = 1
                            break
                    except:
                        pass
            if keys[i].strip().lower()[0:7] == '&caspt2' and (command[1] == 'freqxg' or command[1] == 'freqxgp'):
                caspt2_key = 1

        if ciro_key == 1:
            # generate empty array NAC to take the couplings
            NAC = []
            for i in range(nroots):
                NAC.append([])
            for i in range(nroots):
                for j in range(nroots):
                    NAC[i].append([])
            for i in range(nroots):
                for j in range(nroots):
                    # the multi dimensional matrices must be initialized only with empty lists
                    NAC[i][j].extend([[], [], []])
            NAC_old = copy.deepcopy(NAC)
            # We need an empty array to re-initialise the NAC array every time if a value is given in command 86!
            NAC_empty = copy.deepcopy(NAC)
            # TDNACs will be stored in a different array
            if int(command[85]) > 0 and command[14] == '1':
                TDC = np.zeros((nroots, nroots))
                TDC_old = np.zeros((nroots, nroots))
                # We need an empty array to re-initialize the NAC array every time if a value is given in command 86!
                TDC_empty = np.zeros((nroots, nroots))
                list_of_signs = [1 for i in range(nroots)]
            else:
                TDC, TDC_old, TDC_empty, list_of_signs = None, None, None, None

            sef = shelve.open("cobram-sef")
            sef['NAC_old'] = NAC_old
            sef['NAC'] = NAC
            sef['NAC_empty'] = NAC_empty
            if int(command[85]) > 0 and command[14] == '1':
                sef['TDC_old'] = TDC_old
                sef['TDC'] = TDC
                sef['TDC_empty'] = TDC_empty
                sef['list_of_signs'] = list_of_signs
            sef.close()

            DEarray = np.empty((nroots, nroots))
            DEarray.fill(1000.0)

            # Generate an empty nroots x nroots matrix for usage in tully.py
            zeromat = np.zeros((nroots, nroots))
            sef = shelve.open("cobram-sef")
            sef['zeromat'] = zeromat
            sef.close()
        else:
            # assign dummy values
            DEarray.append([])
            DEarray[0].append(1000.0)
            DEarray[0].append(1000.0)
            NAC = []
            NAC_old = copy.deepcopy(NAC)
            NAC_empty = copy.deepcopy(NAC)
            sef = shelve.open("cobram-sef")
            sef['NAC_old'] = NAC_old
            sef['NAC'] = NAC
            sef['NAC_empty'] = NAC_empty
            sef.close()

        # read additional keys for seward to be added (e.g. PCM etc)
        keys = CBF.ReadCobramCommand(cobcom, 'molcas', 'molcas.input')
        sewardkeys = CBF.ReadCobramCommand(cobcom, 'rfield', 'rfield.input')

        # if ROOT is specified in the ALASKA section of molcas in cobram.command it
        # will take precedence over RLXROOT, STATE and key13
        set_key = 0
        for i in range(len(keys)):
            if keys[i].strip().lower()[0:4] == 'root':
                try:
                    state = int(keys[i + 1]) - 1
                except:
                    state = int(keys[i].split('=')[1]) - 1
                logwrt.writelog("Initial state is {0} and has been taken from ROOT input\n".format(state + 1))
                set_key = 1
        # if ROOT is not specified but RLXROOT is specified the RASSCF section of molcas in cobram.command
        # it will take precedence over STATE and key13
        if set_key == 0:
            for i in range(len(keys)):
                if keys[i].strip().lower()[0:4] == 'rlxr':
                    try:
                        state = int(keys[i + 1]) - 1
                    except:
                        state = int(keys[i].split('=')[1]) - 1
                    logwrt.writelog("Initial state is {0} and has been taken from RLXRoot input.\n".format(state + 1))
                    set_key = 1
        # if not specified in ROOT or RLXRoot, determine the state of interest either through
        # STATE or via key 13 (in that order)
        if set_key == 0:
            if os.path.exists(os.getcwd() + '/STATE'):
                state = int(os.popen('cat STATE').read())
                logwrt.writelog("Initial state is {0} and has been taken from STATE file.\n".format(state + 1))
            else:
                state = int(command[13]) - 1
                logwrt.writelog("Initial state is {0} and has been taken from key 13.\n".format(state + 1))
            # all analytical computations (par_num == 0) will require RLXROOT or NUMERICAL inside &ALASKA
            # numerical frequency computations (par_num == 0 or 1) will require RLXROOT or NUMERICAL inside
            # &ALASKA or a &CASPT2 section
            if ((par_num == 0) or (
                    par_num == 1 and command[1] == 'freqxg' and caspt2_key == 0)) and state != 0 and num_key == 0:
                logwrt.fatalerror('state is different than 0 but there is no RLXROOT section in the RASSCF input.' +
                                  '\nThis is necessary if analytical gradients are requested!')

    else:  # i.e. when step != 0

        sef = shelve.open("cobram-sef")
        state = sef['newstate']
        sef.close()
        cobcom = CBF.getCobrammCommand('cobram.command')

        if command[199] != '0' and float(command[199]) < abs(DEarray[0][1]) and state == 0:
            command[85] = '0'
            nroots = 1
            # read additional keys for seward to be added (e.g. PCM etc)
            keys = CBF.ReadCobramCommand(cobcom, 'molcasSS', 'molcasSS.input')
        else:
            # read additional keys for seward to be added (e.g. PCM etc)
            keys = CBF.ReadCobramCommand(cobcom, 'molcas', 'molcas.input')
        sewardkeys = CBF.ReadCobramCommand(cobcom, 'rfield', 'rfield.input')

    # determine which coupling is to be read from the output or computed numerically during CI optimizations
    if command[1] == 'ci':
        calc_coupl = [[state - 1, state]]
        sef = shelve.open("cobram-sef", 'w')
        sef['calc_coupl'] = calc_coupl
        sef.close()

    # start writing input seward file
    filein = open('seward.input', 'w')
    filein.write(' &SEWARD  &END \n')
    filein.write('Title\n COBRAM\n')
    if command[120] == '1':
        filein.write('EXPERT\n')
    if command[196] == '1':
        filein.write('RICD\n')
        if command[195] != '':
            filein.write('CDTH\n')
            filein.write(str(command[195]) + '\n')
        filein.write('Doana\n')

    # get elements,coordinate and transform in bohr
    elements = [geometry.atomLabel[i - 1] for i in geometry.list_QM] + ["H"] * geometry.NsubH
    X = geometry.modelH[0]
    Y = geometry.modelH[1]
    Z = geometry.modelH[2]

    # displace one of the atoms in case of parallel numerical calculations
    if ((command[1] == 'freqxg' and step != 0 and command[8] == 1) or
            command[1] in ['optxgp', 'mdvp', 'ircp', 'cip', 'tsp']):
        if int(command[210]) > 0:
            iAt, iCoord, iDir, X, Y, Z = geometry.read_displacement(step)
        else: 
            iAt, iCoord, iDir, X, Y, Z = geometry.displace(command, step)

    # use a flexible routine for reading basis set info
    basisset_info = readbassisset(geometry, command, cobcom)

    # construct seward input with basis set definition for each atom (individual definition for each atom is compulsory)
    n = 1
    for i in range(0, len(elements)):
        filein.write('Basis set\n')
        # in the filenames for the Pople's polarized basis set, *'s are substituted with p's
        if basisset_info[i][0] == '6-31G**':
            basisName = '6-31Gpp'
        elif basisset_info[i][0] == '6-31G*':
            basisName = '6-31Gp'
        else:
            basisName = basisset_info[i][0]
        # open file and read content
        try:
            molcaspath = os.getenv('MOLCAS')
            with open(molcaspath + '/basis_library/' + basisName) as basisfile:
                output = basisfile.readlines()
        except IOError:
            logwrt.fatalerror("The given basis set {0} has not been found, please change basis definition".format(
                basisset_info[i][0]))

        # now extract contraction scheme from the basis set file
        for j in range(len(output)):
            # in the files for the Pople's polarized basis set, names with *'s should be used
            if basisset_info[i][0] == '6-31Gpp':
                basisName = '6-31G**'
            elif basisset_info[i][0] == '6-31Gp':
                basisName = '6-31G*'
            else:
                basisName = basisset_info[i][0]
            set2 = output[j].lower().find('/' + elements[i].lower() + '.' + str(basisName.lower()))
            # replace default contraction with the user-defined unless the basis set has only one contraction
            if ((basisset_info[i][0] not in ['6-31G', '6-31G*', '6-31Gp', '6-31G**', '6-31Gpp', 'STO-3G'])
                    and (set2 != -1)):
                tmp_label = output[j][1:].split('.')
                if basisset_info[i][1] != 0:
                    tmp_label[4] = str(basisset_info[i][1])  # set user-defined basis set contraction
                filein.write('.'.join(tmp_label) + '\n')
            elif set2 != -1:
                filein.write(output[j][1:] + '\n')
        filein.write(' ' + elements[i] + str(n) + '   ')  # this format (e.g. N1, C2) is recognized by gv.exe
        filein.write('%14.8f' % (X[i] / constants.Bohr2Ang))
        filein.write('%14.8f' % (Y[i] / constants.Bohr2Ang))
        filein.write('%14.8f' % (Z[i] / constants.Bohr2Ang) + '\n')
        if command[194] == '1':
            filein.write('Cartesian all\n')
        else:
            filein.write('Spherical all\n')
        filein.write('End of Basis\n')
        n = n + 1

    try:
        for i in range(len(sewardkeys)):
            filein.write(str(sewardkeys[i]) + '\n')
    except:
        pass

    # write external charges
    if geometry.calculationType != 'H':

        # ################# USE GHOST ATOMS - COMMAND[120] == "1" ######################
        if command[120] == '1':
            # movable atoms are added as ghost atoms
            if command[1] not in ['freqxgp', 'optxgp', 'mdvp', 'ircp', 'cip', 'tsp']:
                logwrt.writelog('Appending the movable MM atoms as GHOST ATOMS to the seward file\n')
            for n, ch, x, y, z in zip(range(geometry.NatomM), charges.CRG_MEDIUM, *geometry.MEDIUM):
                filein.write('Basis set\nH ..... / inline\n  0.0    0\n0 0\n')
                filein.write(' ' + str(n) + 'H' + '   ')
                filein.write("{0:14.8f}{1:14.8f}{2:14.8f}\n".format(
                    x / constants.Bohr2Ang, y / constants.Bohr2Ang, z / constants.Bohr2Ang))
                filein.write("charge\n  {0:12.8f}\n".format(ch))
                filein.write('End of Basis\n')
            # low layer atoms are added as static charges
            if command[1] not in ['freqxgp', 'optxgp', 'mdvp', 'ircp', 'cip', 'tsp']:
                logwrt.writelog('Appending the atomic point charges as XFIELD to the seward file\n')
            filein.write('XField\n{0}\n'.format(geometry.NatomL))
            for ch, x, y, z in zip(charges.CRG_LOW, *geometry.LOW):
                filein.write("{0:14.8f}{1:14.8f}{2:14.8f}{3:12.8f}   0.0   0.0   0.0\n".format(
                    x / constants.Bohr2Ang, y / constants.Bohr2Ang, z / constants.Bohr2Ang, ch))

        # ################# USE ELECTRIC FIELD (DEFAULT) - COMMAND[120] == "0" ######################
        else:
            if command[1] not in ['freqxgp', 'optxgp', 'mdvp', 'ircp', 'cip', 'tsp']:
                logwrt.writelog('Appending the atomic point charges as XFIELD to the seward file\n')
            filein.write('XField\n{0}\n'.format(geometry.NatomMM))
            for ch, x, y, z in zip(charges.CRG_pod, *geometry.pod):
                filein.write("{0:14.8f}{1:14.8f}{2:14.8f}{3:12.8f}   0.0   0.0   0.0\n".format(
                    x / constants.Bohr2Ang, y / constants.Bohr2Ang, z / constants.Bohr2Ang, ch))
            if command[1] not in ['freqxgp', 'optxgp', 'mdvp', 'ircp', 'cip', 'tsp']:
                logwrt.writelog('Computing elect field at the movable pnt charges to get the coulomb force on the charges\n')
            filein.write('EFLD\n{0}\n'.format(geometry.NatomM))
            for x, y, z in zip(*geometry.MEDIUM):
                filein.write("{0:14.8f}{1:14.8f}{2:14.8f};\n".format(
                    x / constants.Bohr2Ang, y / constants.Bohr2Ang, z / constants.Bohr2Ang))

    filein.write('End of input\n')
    filein.close()

    # automatizing the input generation
    rassi_key = 0
    alaska_key = 0
    num_key = 0
    caspt2_key = 0
    mspt2_key = 0
    #applies to both xms and rms
    rxmspt2_key = 0
    sspt2_key = 0
    link_jobmix_key = 0
    fileout = open('seward.input', 'a')
    for i in range(len(keys)):

        # use JobIph instead of INPORB for faster convergence of the RASSCF
        if step != 0 and keys[i].strip().lower()[0:4] == 'lumo':
            if not os.path.exists('molcas.JobIph'):
                logwrt.writewarning("JOBIPH needed for the restart but it seems to be missing... continue with RASORB")
            elif os.path.exists('molcas.JobIph'):
                # an option for CIONLY during numerical SP computations (if at all should be used with +/- displ.)
                if command[1] in ['mdvp', 'optxgp', 'ircp', 'cip', 'tsp']:
                    if command[11] == '0':
                        keys[i] = 'CIRESTART\nCIONLY'
                    else:
                        keys[i] = 'CIRESTART'
                elif command[1] in ['freqxgp', 'freqxg', 'optxg', 'mdv', 'irc', 'ci', 'ts']:
                    keys[i] = 'CIRESTART'
        # in analytical CI optimizations the RASSCF input must be repeated in order to compute the gradient of the
        # lower state therefore it is first extracted from the user input
        if command[1] == 'ci' and par_num == 0 and keys[i].strip().lower()[0:7] == '&rasscf':
            rasscf_input = []
            rasscf_start = i
            for j in range(i + 1, len(keys)):
                try:
                    if keys[j].strip()[0] == '&':
                        rasscf_end = j - 1
                        break
                    elif keys[j].strip().lower()[0:12] == 'end of input':
                        rasscf_end = j
                        break
                    # in case when the !molcas ?molcas ends with the last line of the rasscf input
                    # i.e., there is neither a following routine nor an 'end of input' we have to make sure that
                    # rasscf_end is assigned, this is why rasscf_end will re-assigend in every iteration
                    else:
                        rasscf_end = j
                except:
                    pass
            # modify the RASSCF input to speed up convergence and to store density of lower state
            skip = -1
            for j in range(rasscf_start, rasscf_end + 1):
                if j == skip:
                    continue
                elif keys[j].strip().lower()[0:4] == 'lumo':
                    rasscf_input.append('CIRESTART\nCIONLY')
                elif keys[j].strip().lower()[0:4] == 'rlxr':
                    try:
                        int(keys[j + 1])
                        rasscf_input.append('RLXROOT\n')
                        rasscf_input.append(str(state))
                        skip = j + 1
                    except:
                        rasscf_input.append('RLXROOT=' + str(state))
                else:
                    rasscf_input.append(keys[j])
            # in the rasscf_input list the counting starts at 0
            rasscf_end = rasscf_end - rasscf_start + 1
            rasscf_start = 0
        # with 'mdvp' and 'cip' there is a special rassi procedure during displacements which is different from the
        # rassi at the reference geometry therefore, delete the rassi input when in 'mdvp' (i.e. when creating the
        # input for the displaced geometries)
        if keys[i].strip().lower()[0:6] == '&rassi' and (command[1] == 'mdvp' or command[1] == 'cip'):
            keys[i] = ''
            for j in range(i + 1, len(keys)):
                try:
                    if keys[j].strip()[0] == '&':
                        break
                    elif keys[j].strip().lower()[0:12] == 'end of input':
                        keys[j] = ''
                        break
                    else:
                        keys[j] = ''
                except:
                    pass
            rassi_key = 0
        if keys[i].strip().lower()[0:7] == '&caspt2' and (command[1] == 'freqxg' or command[1] == 'freqxgp'):
            caspt2_key = 1
        if keys[i].strip().lower()[0:1] == '>' and (command[1] == 'mdvp' or command[1] == 'cip'):
            keys[i] = ''
        # in every other instance the user input has precedence and is not modified
        elif (keys[i].strip().lower()[0:6] == '&rassi') and (command[1] != 'mdvp' and command[1] != 'cip'):
            rassi_key = 1
        if keys[i].strip().lower()[0:7] == '&alaska':
            alaska_key = 1
            """ if Molcas internal numerical gradient algorithm has been called (&ALASKA; nume) but ROOT has not been specified,
            the root is taken from key 13 and added to &ALASKA""" 
            for j in range(i + 1, len(keys)):
                try:
                    if keys[j].strip()[0:1] == '&':
                        break
                    if keys[j].strip().lower()[0:4] == 'nume':
                        num_key = 1
                        for k in range(i + 1, len(keys)):
                            if keys[k].strip()[0:1] == '&':
                                break
                            if keys[k].strip().lower()[0:4] == 'root':
                                num_key = 0
                                break
                        if num_key == 1:
                            if int(command[2]) == 0 and command[1] in ['freqxgp', 'optxgp', 'mdvp', 'ircp', 'cip',
                                                                        'tsp']:
                                pass
                            else:
                                logwrt.writelog('Adding ROOT ' + str(state + 1) + ' to &ALASKA\n')
                            keys[j] = 'NUME\nROOT\n' + str(state + 1)
                except:
                    pass
        if keys[i].strip().lower()[0:6] == '&mbpt2' and command[101] == '0' and int(command[60]) > 1:
            check_grdt = 0
            for j in range(i + 1, len(keys)):
                try:
                    if keys[i].strip().lower()[0:4] == 'grdt':
                        check_grdt = 1
                        break
                    if keys[j].strip()[0:1] == '&' and check_grdt == 0:
                        if int(command[2]) == 0 and command[1] in ['freqxgp', 'optxgp', 'mdvp', 'ircp', 'cip', 'tsp']:
                            pass
                        else:
                            logwrt.writelog("Adding GRDT to MBPT2\n")
                        keys[i] = '&MBPT2\nGRDT'
                        check_grdt = 1
                        break
                except:
                    pass
            if check_grdt == 0:
                if int(command[2]) == 0 and command[1] in ['freqxgp', 'optxgp', 'mdvp', 'ircp', 'cip', 'tsp']:
                    pass
                else:
                    logwrt.writelog("Adding GRDT to MBPT2\n")
                keys[i] = '&MBPT2\nGRDT'
        """ in case of MS-PT2 energies and MediumLayer is present (geometry.NatomMM > 0) only the MS-PT2 properties 
        (dipole, electric field) are taken from an additional CASPT2 section, this is necessary only at the reference geometry 
        (command[1] has no suffix 'p') """ 
        """ if NOMULT or NOMIX is found, then this is a SS-PT2 computations and the properties are available (mspt2_key = 0) """
        if keys[i].strip().lower()[0:7] == '&caspt2' and int(command[60]) > 1 and geometry.NatomMM > 0 and \
                                                                geometry.calculationType not in ['H', 'M', 'ML']:
            mspt2_key = 1
            mspt2_input = ['&CASPT2']
            for j in range(i + 1, len(keys)):
                try:
                    if keys[j].strip().lower()[0:4] == 'nomu' or keys[j].strip().lower()[0:4] == 'nomi':
                        mspt2_key = 0
                        mspt2_input = []
                        break
                    if keys[j].strip()[0] != '>' and keys[j].strip()[0] != '&' and keys[j].split()[0].lower()[0:4] != 'grdt' and keys[j].split()[0].lower()[0:4] != 'sadr' and keys[j].split()[0].lower()[0:4] != 'rlxr':
                        mspt2_input.append(keys[j])
                    elif keys[j].strip()[0] == '>' or keys[j].strip()[0] == '&':
                        mspt2_input.append('end of input')
                except:
                    pass
        """ in case of MS-PT2 energies with MediumLayer (mspt2_key == 1) PROP is set to NOPROP in the first CASPT2 
        iteration to speed up computations; properties will be computed in the second CASPT2 iteration """
        if keys[i].strip().lower()[0:7] == '&caspt2' and mspt2_key == 1:
            for j in range(i + 1, len(keys)):
                try:
                    if (keys[j].strip().lower()[0:12] == 'end of input' or keys[j].strip()[0] == '>'
                            or keys[j].strip()[0] == '&'):
                        break
                    if keys[j].strip().lower()[0:4] == 'prop':
                        keys[j] = 'NOPROP'
                        break
                except:
                    pass
        if keys[i].strip().lower()[0:7] == '&caspt2' and (command[1] == 'mdv' or command[1] == 'ci'):
            link_jobmix_key = 1
            for j in range(i + 1, len(keys)):
                if j == len(keys) - 1:
                    link_jobmix_key = 0
                try:
                    if keys[j].strip().lower()[0:4] == 'nomu' or keys[j].strip().lower()[0:4] == 'nomi':
                        link_jobmix_key = 0
                        break
                    elif keys[j].strip().lower()[0:12] == 'end of input':
                        for k in range(j + 1, len(keys)):
                            if k == len(keys) - 1:
                                link_jobmix_key = 0
                            try:
                                if (keys[k].strip()[0] == '>' and keys[k].split()[-1] == 'JOB001' and
                                        keys[k].split()[-2].split('.')[-1].lower() == 'jobmix'):
                                    link_jobmix_key = 0
                                    logwrt.writelog('Detected link/copy of .JobMix to JOB001 requested' +
                                                    'by user after MS-CASPT2 routine.\n')
                                    break
                                elif (keys[k].strip()[0] == '>' and keys[k].split()[-1] == 'JOB001' and
                                      keys[k].split()[-2].split('.')[-1].lower() == 'jobiph'):
                                    link_jobmix_key = 0
                                    logwrt.writewarning('Detected link/copy of .JobIph to JOB001 requested by user' +
                                                        ' after MS-CASPT2 routine \nalthough NOMULT or NOMIX were ' +
                                                        ' not invoked. Possible source of error!')
                                    break
                                elif keys[k].strip()[0] == '&' and keys[k].strip().lower()[0:6] == '&rassi':
                                    logwrt.writelog('The .JobMix of the MS-CASPT2 computation will be used in RASSI\n')
                                    keys[k - 1] += \
                                        '\n\n>>COPY $WorkDir/molcas.JobMix JOB001\n'
                                    link_jobmix_key = 0
                                    break
                                elif keys[k].strip()[0] == '&' and keys[k].strip().lower()[0:6] != '&rassi':
                                    link_jobmix_key = 0
                                    break
                            except:
                                pass
                        break
                    elif keys[j].strip()[0] == '>' and keys[j].split()[-1] == 'JOB001' and \
                            keys[j].split()[-2].split('.')[-1].lower() == 'jobmix':
                        link_jobmix_key = 0
                        logwrt.writelog(
                            'Detected link/copy of .JobMix to JOB001 requested by user after MS-CASPT2 routine.\n')
                        break
                    elif keys[j].strip()[0] == '>' and keys[j].split()[-1] == 'JOB001' and \
                            keys[j].split()[-2].split('.')[-1].lower() == 'jobiph':
                        link_jobmix_key = 0
                        logwrt.writewarning(
                            'Detected link/copy of .JobIph to JOB001 requested by user after MS-CASPT2 routine ' +
                            '\nalthough NOMULT or NOMIX were not invoked. Possible source of error!')
                        break
                    elif keys[k].strip()[0] == '&' and keys[j].strip().lower()[0:6] == '&rassi':
                        logwrt.writelog('The user is requesting a MULTISTATE CASPT2 computation\n')
                        logwrt.writelog('A JOBMIX file is created and will be used in RASSI\n')
                        keys[j - 1] += '\n\n>>COPY $WorkDir/molcas.JobMix JOB001\n'
                        link_jobmix_key = 0
                        break
                    elif keys[k].strip()[0] == '&' and keys[k].strip().lower()[0:6] != '&rassi':
                        link_jobmix_key = 0
                        break
                except:
                    pass
        if command[1] in ['mdvp', 'optxgp', 'freqxgp', 'freqxg', 'ircp', 'cip', 'tsp'] and keys[i].strip().lower()[
                                                                                           0:7] == '&caspt2':
            for j in range(i + 1, len(keys)):
                try:
                    if keys[j].strip().lower()[0:12] == 'end of input' or \
                            keys[j].strip()[0] == '>' or keys[j].strip()[0] == '&':
                        break
                    if keys[j].strip().lower()[0:4] == 'nomu' or keys[j].strip().lower()[0:4] == 'nomi':
                        sspt2_key = 1
                        break
                except:
                    pass
        """ in case of PT2, properties are skipped at the displaced geometries in numerical compuations """
        if keys[i].strip().lower()[0:7] == '&caspt2' and (
                command[1] in ['mdvp', 'optxgp', 'freqxgp', 'ircp', 'cip', 'tsp'] or (
                command[1] == 'freqxg' and step != 0)):
            for j in range(i + 1, len(keys)):
                try:
                    if keys[j].strip().lower()[0:12] == 'end of input' or \
                            keys[j].strip()[0] == '>' or keys[j].strip()[0] == '&':
                        break
                    if (command[1] != 'freqxgp') and (command[1] != 'freqxg') and keys[j].strip().lower()[0:4] == 'prop':
                        keys[j] = 'noprop'
                    """ in case of SS-CASPT2 all states must be computed at the reference, whereas only one state at 
                    every displaced geometry (i.e. 205 = 1); in order to avoid root swapping between any two states state-1 
                    or/and state+1 will be also computed if gap is < thres 0.004 Hartree (2.5 kcal/mol) for numerical 
                    displacement of 0.001A and < thres 0.04 Hartree (25 kcal/mol) for numerical displacement of 0.01A """
                    if sspt2_key == 1 and keys[j].strip().lower()[0:4] == 'mult' and command[205] == '1' and (
                            (command[1] == 'mdvp' and command[14] == '1') or
                            command[1] in ['optxgp', 'freqxgp', 'freqxg', 'ircp', 'tsp']):
                        sef = shelve.open("cobram-sef")
                        unsort_states = sef['unsort_states']
                        sef.close()
                        try:
                            int(keys[j + 1].split()[0].strip())
                            keys[j + 1] = str(len(unsort_states)) + ' '
                            for istate in unsort_states:
                                keys[j + 1] += ' ' + str(istate + 1) + ' '
                        except:
                            keys[j] = 'multistate=' + str(len(unsort_states)) + ' '
                            for istate in unsort_states:
                                keys[j] += ' ' + str(istate + 1) + ' '
                except:
                    pass
        fileout.write(keys[i] + '\n')
    if command[1] == 'mdv' and par_num == 0 and int(command[85]) > 0 and step != 0 and \
            len(calc_coupl) != 0 and command[14] == '0':
        logwrt.writelog('Adding &ALASKA for non-adiabatic coupling(s)\n')
        for i in range(len(calc_coupl)):
            j = calc_coupl[i][0]
            k = calc_coupl[i][1]
            fileout.write('\n&ALASKA\nNAC\n' + str(j + 1) + ' ' + str(k + 1) + '\n')
    if (command[1] == 'mdv' or command[1] == 'ci') and link_jobmix_key == 1:
        logwrt.writelog('The user is requesting a MULTISTATE CASPT2 computation\n')
        logwrt.writelog('A JOBMIX file is created and will be used in RASSI\n')
        fileout.write('\n>>COPY $WorkDir/molcas.JobMix JOB001\n')
    # rassi input for computing the overlaps between the WF at r (molcas_r.JobIph) and r+dr (molcas.JobIph)
    if (command[1] == 'mdvp' or command[1] == 'cip') and command[14] == '0' and not os.path.exists(
            os.getcwd() + '/MBPT2'):
        for i in range(nroots):
            if abs(DEarray[state][i]) < float(command[86]) or command[1] == 'cip':
                # .JobMix is used in the case of CASPT2
                if os.path.exists(os.getcwd() + '/MOLCAS/molcas_r.JobMix'):
                    fileout.write('>>COPY $WorkDir/molcas_r.JobMix JOB001\n' +
                                  '>>COPY $WorkDir/molcas.JobMix JOB002\n')
                # .JobIph is used in the case of CASSCF
                else:
                    fileout.write('>>COPY $WorkDir/molcas_r.JobIph JOB001\n'
                                  '>>COPY $WorkDir/molcas.JobIph JOB002\n')
                if int(command[2]) == 0:
                    pass
                else:
                    logwrt.writelog(
                        'Adding &RASSI for WF overlaps (numerical computation of non-adiabatic couplings)\n')
                fileout.write('&RASSI\nnrof=2 ' + str(nroots) + ' ' + str(nroots) + '\n')
                for i in range(2):
                    for j in range(1, nroots + 1):
                        fileout.write(str(j) + ' ')
                    fileout.write('\n')
                fileout.write('\nonel\nEnd of input\n')
                break
    """ the &ALASKA section is written to input in case of analytical gradients if it is not a SP (i.e. 60 > 1) and 
    unless already specified by the user """
    if command[1] in ['optxg', 'mdv', 'irc', 'ci', 'ts'] and par_num == 0 and alaska_key == 0 and int(command[60]) > 1:
        logwrt.writelog('Adding &ALASKA for gradient\n')
        fileout.write('\n&ALASKA\nShow\nEnd of input\n')
    """ in case of analytical RASSCF CI optimization add two &RASSCF and &ALASKA sections to the input """
    if command[1] == 'ci' and par_num == 0:
        fileout.write('\n')
        logwrt.writelog('Adding &ALASKA for non-adiabatic coupling\n')
        for i in range(rasscf_start, rasscf_end):
            fileout.write(rasscf_input[i] + '\n')
        fileout.write('\n&ALASKA\nShow\nEnd of input\n')
        fileout.write('\n&ALASKA\nNAC\n' + str(state) + ' ' + str(state + 1) + '\n')
    elif (command[1] == 'freqxg' or command[1] == 'freqxgp') and alaska_key == 0 and caspt2_key == 0:
        if int(command[2]) == 0 and command[1] in ['freqxgp', 'optxgp', 'mdvp', 'ircp', 'cip', 'tsp']:
            pass
        else:
            logwrt.writelog('Adding &ALASKA for gradient\n')
        fileout.write('\n&ALASKA\nShow\nEnd of input\n')
    elif command[1] in ['freqxg', 'freqxgp'] and alaska_key == 0 and caspt2_key == 1:
        if int(command[2]) == 0 and command[1] in ['freqxgp', 'optxgp', 'mdvp', 'ircp', 'cip', 'tsp']:
            pass
        else:
            logwrt.writelog('Adding &ALASKA for numerical gradient\n')
        fileout.write('\n&ALASKA\nNUME\nROOT\n' + str(state + 1) + '\nDELTA\n' + str(
            float(command[12]) / 0.52917720830008323378) + '\nEnd of input\n')
    # the rassi input at the reference compuation is written to seward.input unless already specified by the user
    if rassi_key == 0 and command[1] == 'mdv' and command[85] == '1':
        logwrt.writelog('Adding &RASSI for WF\n')
        fileout.write('&RASSI &END\nNR OF JOBIPHS\n1 ' + str(nroots) + '\n')
        for i in range(1, nroots + 1):
            fileout.write(str(i) + ' ')
        if link_jobmix_key == 1:
            fileout.write('\nEJOB\n')
        fileout.write('\nCIPRint\nTHRS\n0.0\nEnd of Input\n')
    # set a RASSI section to compute TNACS (when repeat_numerics == 0)
    if command[1] == 'mdv' and command[14] == '1' and step > 0 and not os.path.exists(
            os.getcwd() + '/HOP') and not os.path.exists(os.getcwd() + '/MBPT2'):
        # .JobMix is used in the case of CASPT2
        if os.path.exists(os.getcwd() + '/molcas_-dt.JobMix'):
            os.system('cp ' + os.getcwd() + '/molcas_-dt.JobMix ' + os.getcwd() + '/MOLCAS/.')
            fileout.write('>>COPY $WorkDir/molcas_-dt.JobMix JOB001\n'
                          '>>COPY $WorkDir/molcas.JobMix JOB002\n')
        # .JobIph is used in the case of CASSCF
        elif os.path.exists(os.getcwd() + '/molcas_-dt.JobIph'):
            os.system('cp ' + os.getcwd() + '/molcas_-dt.JobIph ' + os.getcwd() + '/MOLCAS/.')
            fileout.write('>>COPY $WorkDir/molcas_-dt.JobIph JOB001\n'
                          '>>COPY $WorkDir/molcas.JobIph JOB002\n')
        else:
            logwrt.fatalerror("COBRAMM expected to find scratch file for time step " + str(step - 1) + " but failed!")
        fileout.write('&RASSI\nnrof=2 ' + str(nroots) + ' ' + str(nroots) + '\n')
        for i in range(2):
            for j in range(1, nroots + 1):
                fileout.write(str(j) + ' ')
            fileout.write('\n')
        fileout.write('\nCIPRint\nonel\nEnd of input\n')

    """ in case of MS-PT2 energies, properties (dipole, electric field) are obtained by computing SS-PT2 on top of the
    PM-CASPT2 WF resulting from the diagonalization of the MS-PT2 Hamiltonian; this is necessary only at reference point 
    (i.e. command[1] has no suffix 'p')"""
    if mspt2_key == 1 and int(command[60]) > 1 and command[1] in ['optxg', 'freqxg', 'irc', 'ts', 'ci', 'mdv']:
        prop_key = 0
        # molcas.JobIph should not be deleted otherwise one cannot make use of the CIRESTART in numerical computations
        fileout.write('\n>>COPY $WorkDir/molcas.JobMix JOBIPH\n')
        for i in range(len(mspt2_input)):
            try:
                if mspt2_input[i] == 'end of input':
                    break
                elif mspt2_input[i] == '':
                    pass
                elif mspt2_input[i].strip().lower()[0:4] == 'prop':
                    prop_key = 1
                elif mspt2_input[i].strip().lower()[0:4] in ['mult', 'xmul', 'rmul']:
                    if mspt2_input[i].strip().lower()[0:4] == 'xmul' or mspt2_input[i].strip().lower()[0:4] == 'rmul':
                        rxmspt2_key = 1
                    """ in case of MS-PT2 energies, in the second &CASPT2 section properties are computed 
                    only for the active state at the reference geometry (i.e. command[1] has no suffix 'p') """
                    try:
                        int(mspt2_input[i + 1].strip()[0])
                        """ in case of CI search we need properties for two states """
                        if command[1] == 'ci':
                            mspt2_input[i + 1] = '2 ' + str(state) + ' ' + str(state + 1)
                        else:
                            mspt2_input[i + 1] = '1 ' + str(state + 1)
                        """ take care that in the 2nd &CASPT2 section one uses mult and not xmul or rmul """
                        mspt2_input[i] = 'multistate'
                    except:
                        mspt2_input[i] = 'multistate=1 ' + str(state + 1)
            except:
                pass
            fileout.write(mspt2_input[i] + '\n')
        fileout.write('NOMULT\nNOMIX')
        if prop_key == 0:
            fileout.write('\nPROP')
        fileout.write('\nend of input\n')
    fileout.close()

    sef = shelve.open("cobram-sef")
    sef['DEarray'] = DEarray
    sef['nroots'] = nroots
    sef['newstate'] = state
    if step == 0:
        sef['state'] = state
    sef.close()

    sys.stdout.flush()
    return


@Timer("QM run")
def launch(command, step):
    # get name of MOLCAS script
    molcasscript = os.getenv('MOLCAS_SCRIPT')
    logwrt.writelog("molcasscript = " + molcasscript + "\n")

    # launch molcas
    fout = open("molcas.log", "w")
    subprocess.call(["nice", "-1", molcasscript, "seward.input"], stdout=fout, stderr=fout)
    fout.close()
    logwrt.writelog("Calculation completed... \n")

    # keeping some file for later use

    # in the case of a frequency calculation, keep only the wavefunction at the reference geometry
    if command[1] == 'freqxg' and step == 0:
        if os.path.exists(os.getcwd() + '/MOLCAS/molcas.JobIph'):
            logwrt.writelog("JOBIPH file saved at the reference geometry \n")
            os.system('cp ' + os.getcwd() + '/MOLCAS/molcas.JobIph ' + os.getcwd())
        elif os.path.exists(os.getcwd() + '/MOLCAS/molcas.RasOrb'):
            logwrt.writelog("RASORB file saved at the reference geometry \n")
            os.system('cp ' + os.getcwd() + '/MOLCAS/molcas.RasOrb ' + os.getcwd())

    # in the other cases, always keep the WF for next step
    elif command[1] in ['optxg', 'mdv', 'irc', 'ci', 'ts']:
        if os.path.exists(os.getcwd() + '/MOLCAS/molcas.JobIph'):
            logwrt.writelog("path = " + os.getcwd() + '/MOLCAS/molcas.JobIph' + "\n")
            logwrt.writelog("JOBIPH file saved for restart in the next step\n")
            os.system('cp ' + os.getcwd() + '/MOLCAS/molcas.JobIph ' + os.getcwd())
        elif os.path.exists(os.getcwd() + '/MOLCAS/molcas.RasOrb'):
            logwrt.writelog("RASORB file saved for restart in the next step\n")
            os.system('cp ' + os.getcwd() + '/MOLCAS/molcas.RasOrb ' + os.getcwd())
        if os.path.exists(os.getcwd() + '/MOLCAS/molcas.JobMix') and (command[1] == 'mdv' or command[1] == 'ci') and \
                command[14] == '0':
            logwrt.writelog("JOBMIX at reference geometry saved for possible numerical NAC computation\n")
            os.system('cp ' + os.getcwd() + '/MOLCAS/molcas.JobMix ' + os.getcwd() + '/MOLCAS/molcas_r.JobMix')
        elif os.path.exists(os.getcwd() + '/MOLCAS/molcas.JobIph') and (command[1] == 'mdv' or command[1] == 'ci') and \
                command[14] == '0':
            logwrt.writelog("JOBIPH at reference geometry saved for possible numerical NAC computation\n")
            os.system('cp ' + os.getcwd() + '/MOLCAS/molcas.JobIph ' + os.getcwd() + '/MOLCAS/molcas_r.JobIph')
        if os.path.exists(os.getcwd() + '/MOLCAS/molcas.JobMix') and not os.path.exists(os.getcwd() + '/HOP') and \
                command[1] == 'mdv' and command[14] == '1':
            logwrt.writelog("JOBMIX saved for possible numerical TDNAC computation\n")
            os.system('cp ' + os.getcwd() + '/MOLCAS/molcas.JobMix ' + os.getcwd() + '/molcas_-dt.JobMix')
        elif os.path.exists(os.getcwd() + '/MOLCAS/molcas.JobIph') and not os.path.exists(os.getcwd() + '/HOP') and \
                command[1] == 'mdv' and command[14] == '1':
            logwrt.writelog("JOBIPH saved for possible numerical TDNAC computation\n")
            os.system('cp ' + os.getcwd() + '/MOLCAS/molcas.JobIph ' + os.getcwd() + '/molcas_-dt.JobIph')
        elif os.path.exists(os.getcwd() + '/SALTO') and command[1] == 'mdv' and command[14] == '1':
            logwrt.writelog("A hop occured, JobIph/JobMix files from previous step are preserved!\n")


@Timer("singlepoint")
def singlepoint(command, step, cobcom):
    """ perform a single step computation after the standard Molcas calculation.
        This option is controlled by command[200] """

    # when command[200] is 0, there is nothing to do
    if command[200] == "0": return

    # at the first step, clean previous calc and create the directory where to run the calculations and to store results
    if step == 0:
        if os.path.isdir(SINGLEPOINT_DIR): shutil.rmtree(SINGLEPOINT_DIR)
        os.mkdir(SINGLEPOINT_DIR)

    # initialize flag that defines whether the point needs to be computed
    doSinglePoint = False

    # IRC case, command[200] = "-1": check convergence of the IRC point
    if command[200] == '-1':
        # always do first step
        if step == 1: doSinglePoint = True
        if step > 1:
            tmp = []
            for line in reversed(open('geometry.log').readlines()):
                tmp.append(line.strip())
            for i in range(len(tmp)):
                try:
                    if tmp[i].find('Converged?') != -1:
                        for j in (range(1, 5)):
                            if tmp[i - j].split()[4] == 'YES':
                                doSinglePoint = True
                            elif tmp[i - j].split()[4] == 'NO':
                                doSinglePoint = False
                                break
                        break
                except:
                    pass
    # when command[200] is positive, do single point every n step, with n = command[200]
    elif int(command[200]) > 0:
        if step % int(command[200]) == 0:
            doSinglePoint = True

    # perform the single point computation when needed
    if doSinglePoint:
        logwrt.writelog("Computing a single point at current step\n")
        logwrt.writelog("Output stored in {0}/molcasALL.log\n".format(SINGLEPOINT_DIR))

        # define env variables to run MOLCAS
        outputDir = os.path.join(os.getcwd(), SINGLEPOINT_DIR)
        os.putenv('MOLCAS_OUTPUT', outputDir)
        tmpDir = os.path.join(outputDir, "tmp")
        os.putenv('WorkDir', tmpDir)
        # create tmpDir
        os.mkdir(tmpDir)

        # prepare input file and restart WF files
        section = CBF.ReadCobramCommand(cobcom, 'singlepoint', 'singlepoint.input')
        for line in section:
            if line.strip().lower()[0:4] == 'cire':
                logwrt.writelog("Restarting CI wavefunction from JobIph file\n")
                shutil.copyfile(os.path.join(MOLCAS_WRKDIR, "molcas.JobIph"), os.path.join(SINGLEPOINT_DIR, "JOBIPH"))
                break
        shutil.copyfile(os.path.join(MOLCAS_WRKDIR, "molcas.RasOrb"), os.path.join(SINGLEPOINT_DIR, "INPORB"))

        # extract the geometry and basis set definitions from main MOLCAS input
        with open("seward.input") as filein:
            newinp = []
            for line in filein:
                newinp.append(line)
                if 'end of input' in line.strip().lower(): break

        # write the geometry and the singlepoint section to the seward file
        with open(os.path.join(SINGLEPOINT_DIR, "seward.input"), "w") as molcasinp:
            for i in range(len(newinp)):
                molcasinp.write(newinp[i])
            for i in range(len(section)):
                molcasinp.write(section[i] + '\n')

        # launch molcas
        os.chdir(SINGLEPOINT_DIR)
        molcasscript = os.getenv('MOLCAS_SCRIPT')
        with open("molcas.log", "w") as fout:
            subprocess.call(["nice", "-1", molcasscript, "seward.input"], stdout=fout, stderr=fout)
        logwrt.writelog("Single point calculation completed... \n")
        os.chdir("..")

        # save molcas log file to the molcasALL.log file
        with open(os.path.join(SINGLEPOINT_DIR, "molcasALL.log"), 'a') as molcastot:
            molcastot.write('=' * 80 + '\n')
            if command[200] == '-1':
                molcastot.write("Start SP Molcas calculation of STEP : " + str(step - 1) + "\n")
            else:
                molcastot.write("Start SP Molcas calculation of STEP : " + str(step) + "\n")
            molcastot.write(' || ' * 20 + '\n')
            molcastot.write(' \/ ' * 20 + '\n')
            with open(os.path.join(SINGLEPOINT_DIR, "molcas.log"), "r") as molcasOut:
                molcastot.write(molcasOut.read())
            molcastot.write(' /\ ' * 20 + '\n')
            molcastot.write(' || ' * 20 + '\n')
            if command[200] == '-1':
                molcastot.write("End SP Molcas calculation of STEP : " + str(step - 1) + "\n")
            else:
                molcastot.write("End SP Molcas calculation of STEP : " + str(step) + "\n")
            molcastot.write('=' * 80 + '\n\n')

        # clean up temporary dir and restore standard MOLCAS environment
        os.putenv('WorkDir', os.path.join(os.getcwd(), MOLCAS_WRKDIR))
        os.putenv('MOLCAS_OUTPUT', os.path.join(os.getcwd(), MOLCAS_WRKDIR))
        shutil.rmtree(tmpDir)


@Timer("numer grad")
def molcasEne(command, geometry):
    """This function computes the numerical gradient from energy calculations for all
    the necessary displacements of the atoms, as computed with the parallel_numerics.run() function"""

    # open shelve data
    sef = shelve.open("cobram-sef")
    nroots = sef['nroots']
    aux_data = sef['aux_data']
    state = sef['newstate']
    DEarray = sef['DEarray']
    MDstep = sef['MDstep']
    sef.close()
    ref_energy = aux_data[19]

    # initialize data
    disp = float(command[12]) / constants.Bohr2Ang
    keep = 0
    gradient = []
    for i in range(nroots):
        gradient.extend([[], [], []])

    # find out the type of calculation
    SCF, MBPT2, CAS, SSPT2, MSPT2 = 'no', 'no', 'no', 'no', 'no'
    append_x, append_y, append_z = 1, 0, 0
    if os.path.exists(os.getcwd() + '/SCF'):
        SCF = 'yes'
    if os.path.exists(os.getcwd() + '/RASSCF'):
        CAS = 'yes'
        SCF = 'no'
    elif os.path.exists(os.getcwd() + '/MBPT2'):
        SCF = 'yes'
        MBPT2 = 'yes'
    elif os.path.exists(os.getcwd() + '/SS-CASPT2'):
        SSPT2 = 'yes'
    elif os.path.exists(os.getcwd() + '/MS-CASPT2'):
        MSPT2 = 'yes'

    # find out whether +/- or only - displacement was performed
    if command[10] == '0':
        tot_disp = 3
        times = 1
    else:
        tot_disp = 6
        times = 2

    # define the number of atoms to displace
    if geometry.NsubH == 0:
        nAtoms, HAtoms = geometry.NatomQM, geometry.NatomQM
    else:
        if command[16] == '0':
            nAtoms, HAtoms = geometry.NatomQM, geometry.NatomQM
        else:
            nAtoms, HAtoms = geometry.NatomQM + geometry.NsubH, geometry.NatomQM

    # cas_order is a (nsteps x nroots) list which contains the order of SS-PT2 states at every displacement
    # it is required for computing the NAC elements properly
    cas_order = []
    for i in range(tot_disp * nAtoms):
        cas_order.append([])
        for j in range(nroots):
            cas_order[i].append(j)

    # loop over the displacements of the atoms
    step = 1
    if int(command[210]) > 0:
        irange = 2*int(command[210])
    else:
        irange = tot_disp * nAtoms
    for i in range(irange):
        sorted_energy, unsorted_energy = [], []
        # read content of output file for geom+str(step) calculation
        logFileName = os.path.join('geom' + str(step), 'molcas.log')
        with open(logFileName) as enefile:
            tmp = [line.strip() for line in enefile]

        # save for debugging
        if command[15] == '1':
            os.system('echo "Displacement step ' + str(i) + ': " >> molcasALLnum.log; ' +
                      'cat ' + logFileName + ' >> molcasALLnum.log')

        # in the case of SS-PT2 computation Cobram has to figure out which CASSCF state (stored in 'cas_order[state]')
        # to read EFIELD and DIPOLES for the variable 'state' (which has to be user-specified via key 13)
        # corresponds to the state of interest at SS-PT2 order
        if SSPT2 == 'yes':
            if SSPT2 == 'yes':
                for j in range(len(tmp)):
                    if tmp[j].find(' CASPT2 Root') != -1:
                        unsorted_energy.append(float(tmp[j].split()[6]))
            if MSPT2 == 'yes':
                tmp_root = 0
                for j in range(len(tmp)):
                    if tmp[j].find(' CASPT2 Root ') != -1 and tmp_root >= nroots:
                        unsorted_energy.append(float(tmp[j].split()[6]))
                    elif tmp[j].find(' CASPT2 Root ') != -1:
                        tmp_root += 1
            if command[205] == '1' and (
                    (command[1] == 'mdv' and command[14] == '1') or command[1] in ['optxg', 'freqxg', 'irc', 'ts']):
                sef = shelve.open("cobram-sef")
                sort_states = sef['sort_states']
                sef.close()
                sorted_energy = list(unsorted_energy)
                sorted_energy.sort()
                the_state = sort_states.index(state)
            else:
                sorted_energy = list(unsorted_energy)
                sorted_energy.sort()
                for j in range(len(sorted_energy)):
                    for k in range(len(unsorted_energy)):
                        if sorted_energy[j] == unsorted_energy[k]:
                            # the state order is stored at each displacement i
                            cas_order[i][j] = k
                            if int(command[2]) > 0 and SSPT2 == 'yes':
                                logwrt.writelog("SSPT2 state " + str(j + 1) + " corresponds to SA-CASSCF state " + str(
                                    cas_order[i][j] + 1) +
                                                " at displaced geometry " + str(step) + "\n")
                            elif int(command[2]) > 0 and MSPT2 == 'yes':
                                logwrt.writelog(
                                    "PM-CASPT2 state " + str(j + 1) + " corresponds to MS-CASPT2 state " + str(
                                        cas_order[i][j] + 1) +
                                    " at displaced geometry " + str(step) + "\n")
                sef = shelve.open("cobram-sef")
                sef['cas_order'] = cas_order
                sef.close()
        elif MSPT2 == 'yes':
            for j in range(len(tmp)):
                if tmp[j].find('MS-CASPT2 Root ') != -1:
                    sorted_energy.append(float(tmp[j].split()[6]))
        elif CAS == 'yes':
            for j in range(len(tmp)):
                if tmp[j].find('root number ') != -1:
                    sorted_energy.append(float(tmp[j].split()[7]))
        elif MBPT2 == 'yes':
            state = 0
            for j in range(len(tmp)):
                if tmp[j].find('Total MBPT2 energy ') != -1:
                    sorted_energy.append(float(tmp[j].split()[4]))
                    break
        elif SCF == 'yes':
            state = 0
            for j in range(len(tmp)):
                if tmp[j].find('Total SCF energy ') != -1:
                    sorted_energy.append(float(tmp[j].split()[4]))
                    break
        if int(command[2]) > 0:
            if command[205] == '1' and command[14] == '1':
                logwrt.writelog("Energy of state {0}: {1} Hartree\n".format(state + 1, sorted_energy[the_state]))
            else:
                for j in range(nroots):
                    logwrt.writelog("Energy of state {0}: {1} Hartree\n".format(j + 1, sorted_energy[j]))

        # a routine to check for WF convergence through inspecting the orbital rotations in case of orbital
        # rotations exceeding the threshold of 0.01 Cobram will print a warming and store the computation log file
        if (SSPT2 == 'yes' or MSPT2 == 'yes' or CAS == 'yes') and command[1] in ['mdv', 'irc', 'ci', 'ts'] and \
                command[11] == '1':
            check_maxrot = 0
            for j in range(len(tmp)):
                if tmp[j].find('max ROT') != -1:
                    k = 1
                    while True:
                        k += 1
                        if tmp[i + k].find('Wave function printout') != -1: break
                        if tmp[j + k].find('No convergence after') != -1:
                            errFileName = os.path.join(QM_DATA_STORAGE, 'molcas_' + str(MDstep) + '.err')
                            with open(errFileName, "a") as errFile:
                                errFile.write("Displacement step " + str(step) + ": \n")
                                with open(logFileName, "r") as logFile:
                                    errFile.write(logFile.read())
                            logwrt.fatalerror(
                                'No convergence in CASSCF numerical computation for displacement {0}.\n'.format(step) +
                                'Check molcas output for step {0} in {1}.'.format(step, errFileName))
                        if tmp[j + k].find('Convergence after') != -1:
                            break
                        if check_maxrot == 0:
                            try:
                                if abs(float(tmp[j + k].split()[6].strip('*'))) > 0.01:
                                    errFileName = os.path.join(QM_DATA_STORAGE, 'molcas_' + str(MDstep) + '.err')
                                    with open(errFileName, "a") as errFile:
                                        errFile.write("Displacement step " + str(step) + ": \n")
                                        with open(logFileName, "r") as logFile:
                                            errFile.write(logFile.read())
                                    logwrt.writewarning('Large orbital rotation during numerical computation for ' +
                                                        'displacement {0}: max ROT greater than 0.01.\n'.format(step) +
                                                        'Check molcas output for step {0} in {1}.'.format(step,
                                                                                                          errFileName))
                                    check_maxrot = 1
                            except:
                                pass

        errFileName = os.path.join(QM_DATA_STORAGE, 'molcas_' + str(MDstep) + '.err')
        itmp = []
        if os.path.exists(errFileName):
            with open(errFileName, "r") as ferr:
                itmp = ferr.readlines()

        # a routine to check for too large gradients
        grad_warn = 0
        if SSPT2 == 'yes' and command[205] == '1' and \
                ((command[1] == 'mdv' and command[14] == '1') or command[1] in ['optxg', 'freqxg', 'irc', 'ts']):
            if abs((sorted_energy[the_state] - float(ref_energy[state])) / (times * disp)) > 0.5: grad_warn = 1
        else:
            if abs((sorted_energy[state] - float(ref_energy[state])) / (times * disp)) > 0.5: grad_warn = 1

        if grad_warn == 1:
            logwrt.writewarning('Large gradient during numerical computation for displacement {0}.\n'.format(step) +
                                'Check molcas output for step {0} in {1}.'.format(step, errFileName))
            # decide whether the error log needs to be printed in the err file
            disp_found = 0
            for j in range(len(itmp)):
                if itmp[j].find('Displacement step ' + str(step) + ':') != -1:
                    disp_found = 1
                    break
            if disp_found == 0:
                with open(errFileName, "a") as errFile:
                    errFile.write("Displacement step " + str(step) + ": \n")
                    with open(logFileName, "r") as logFile:
                        errFile.write(logFile.read())

        # if key 1 'optxg' or 'mdv' with analytical NACs through Molpro we don't need the geomXXX and molpro.log; so,
        # they can be deleted
        if (command[1] == 'mdv' and int(command[85]) > 0 and command[14] == '2') or (
                command[1] == 'mdv' and int(command[85]) < 1) or command[1] in ['optxg', 'irc', 'ts']:
            os.system('rm -rf ' + os.getcwd() + '/geom' + str(step))

        # this is done here, as the state has just been determined
        if command[1] == 'mdv' and command[14] == '0':
            for j in range(nroots):
                if j == state:
                    pass
                else:
                    if abs(DEarray[state][j]) <= float(command[86]):
                        keep = 1
                        break
        if command[1] == 'ci':
            keep = 1
        if keep == 0:
            os.system('rm -rf ' + os.getcwd() + '/geom' + str(step))
        if SSPT2 == 'yes' and command[205] == '1' and \
                ((command[1] == 'mdv' and command[14] == '1') or command[1] in ['optxg', 'freqxg', 'irc', 'ts']):
            if i % 2 == 0:
                if command[10] == '1':
                    forward = sorted_energy[the_state]
                    step += 1
                    continue
                elif command[10] == '0':
                    forward = sorted_energy[the_state]
                    backward = float(ref_energy[state])
            else:
                if command[10] == '1':
                    backward = sorted_energy[the_state]
                elif command[10] == '0':
                    forward = sorted_energy[the_state]
                    backward = float(ref_energy[state])
            if append_x == 1:
                gradient[0].append(-(forward - backward) / (times * disp))
                append_x = 0
                append_y = 1
            elif append_y == 1:
                gradient[1].append(-(forward - backward) / (times * disp))
                append_y = 0
                append_z = 1
            elif append_z == 1:
                gradient[2].append(-(forward - backward) / (times * disp))
                append_z = 0
                append_x = 1
            step += 1
        else:
            if i % 2 == 0:
                if command[10] == '1':
                    forward = []
                    for j in range(nroots):
                        forward.append(sorted_energy[j])
                    step += 1
                    continue
                elif command[10] == '0':
                    forward = []
                    backward = []
                    for j in range(nroots):
                        forward.append(sorted_energy[j])
                        backward.append(float(ref_energy[j]))
            else:
                if command[10] == '1':
                    backward = []
                    for j in range(nroots):
                        backward.append(sorted_energy[j])
                elif command[10] == '0':
                    forward = []
                    backward = []
                    for j in range(nroots):
                        forward.append(sorted_energy[j])
                        backward.append(float(ref_energy[j]))
            # I don't know how to distinguish x,y,z as all folders are identical
            # BUT I know that there is a clear order of reading (x_1, y_1, z_1, x_2, y_2, z_2, etc.)
            # therefore, I have associated the appending with iterating auxiliary labels append_x, append_y, append_z
            for j in range(nroots):
                if append_x == 1:
                    gradient[3 * j].append(-(forward[j] - backward[j]) / (times * disp))
                    if j == nroots - 1:
                        append_x = 0
                        append_y = 1
                elif append_y == 1:
                    gradient[3 * j + 1].append(-(forward[j] - backward[j]) / (times * disp))
                    if j == nroots - 1:
                        append_y = 0
                        append_z = 1
                elif append_z == 1:
                    gradient[3 * j + 2].append(-(forward[j] - backward[j]) / (times * disp))
                    if j == nroots - 1:
                        append_z = 0
                        append_x = 1
            step += 1
    if command[16] == '0':
        for i in range(geometry.NsubH):
            for j in range(len(gradient)):
                gradient[j].append(0.0)
    if command[1] == 'ci':
        sef = shelve.open("cobram-sef")
        sef['gradient2'] = [gradient[3 * (state - 1)], gradient[3 * (state - 1) + 1], gradient[3 * (state - 1) + 2]]
        sef['gradient'] = [gradient[3 * (state)], gradient[3 * (state) + 1], gradient[3 * (state) + 2]]
        sef.close()
    # try is required as in numerical HF/MBPT2/GS-CASSCF/GS-CASPT2 computations there is no gradient[3]...
    try:
        # in SS-CASPT2 with key 206 the grdient is computed only for a reduced number of states, mostly only for the photoactive state...
        if not gradient[3]:
            gradient = [gradient[0], gradient[1], gradient[2]]
        # ...therefore the results are stored in the first three arrays EVEN if gradient[3] exists (it is empty though)
        else:
            gradient = [gradient[3 * (state)], gradient[3 * (state) + 1], gradient[3 * (state) + 2]]
    # ...therefore the first three arrays are read
    except:
        gradient = [gradient[0], gradient[1], gradient[2]]
    # in case of THS (14 == 1) the scaling will be performed along the gradient difference vector (206 == 0)    
    # which must be computed taking into account the redistribution of the link atom gradients on the H and M atoms 
    if int(command[85]) > 0 and command[14] == '1' and command[206] == '0' and os.path.exists(os.getcwd() + '/HOP'):
        sef = shelve.open("cobram-sef")
        grad_state = sef['gradient']
        gradch_state = sef['gradch']
        gradch_newstate = sef['gradch2']
        sef.close()
        computeGD(grad_state, gradient, gradch_state, gradch_newstate, geometry, command)
    sef = shelve.open("cobram-sef")
    if command[1] == 'mdv':
        try:
            sef['gradient2'] = sef['gradient']
        except:
            pass
    sef['gradient'] = gradient
    sef.close()
    if int(command[210])%3 != 0:
        gradient[2].append(0.0)
        if int(command[210])%3 == 1:
            gradient[1].append(0.0)
    
    if logwrt.VERBOSITY_LEVEL > 0:
        logwrt.writelog("\nGRADIENT ALONG COMPUTED COORDINATES:\n\n")
        if int(command[210]) > 0:
            logwrt.writelog('!!! Parallel computation with external coordinates !!!!\n')
            logwrt.writelog('If less then 3*n coordinates are used, the missing ones to fit cartesian format will be set equal to 0.0 for printout\n')
        for at in range(len(gradient[0])):
            logwrt.writelog("{0:16.12f} {1:16.12f} {2:16.12f}\n".format(gradient[0][at], gradient[1][at], gradient[2][at]))    

    return gradient


# this routine computed numerical NACs
@Timer("numer NAC")
def molcasNAC(command, geometry):
    sef = shelve.open("cobram-sef", 'r')
    aux_data = sef['aux_data']
    scal_factors = aux_data[18]
    nroots = sef['nroots']
    NAC = sef['NAC_empty']
    calc_coupl = sef['calc_coupl']
    if os.path.exists(os.getcwd() + '/SS-CASPT2'):
        cas_order2print = []
        cas_order = sef['cas_order']
        for i in range(len(cas_order)):
            cas_order2print.append([])
            cas_order2print[i] = [x + 1 for x in cas_order[i]]
        cas_ref_order = sef['cas_ref_order']
        cas_ref_order2print = [x + 1 for x in cas_ref_order]
    sef.close()

    disp = float(command[12]) / constants.Bohr2Ang
    o = np.zeros((nroots, nroots))
    o_plus = np.zeros((nroots, nroots))
    o_minus = np.zeros((nroots, nroots))
    dc = np.zeros((nroots, nroots))
    o_ref = np.zeros((nroots, nroots))
    append_x, append_y, append_z = 1, 0, 0

    # find out whether +/- or only - displacement was performed
    if command[10] == '0':
        tot_disp = 3
    else:
        tot_disp = 6

    # define the number of atoms to displace
    if geometry.NsubH == 0:
        nAtoms, HAtoms = geometry.NatomQM, geometry.NatomQM
    else:
        if command[16] == '0':
            nAtoms, HAtoms = geometry.NatomQM, geometry.NatomQM
        else:
            nAtoms, HAtoms = geometry.NatomQM + geometry.NsubH, geometry.NatomQM

    if os.path.exists(os.getcwd() + '/SS-CASPT2'):
        if int(command[2]) > 0:
            logwrt.writelog('Order of PT2 states at reference geometry {0}\n'.format(np.array(cas_ref_order2print)))

    step = 1
    if int(command[210]) > 0:
        irange = 2*int(command[210])
    else:
        irange = tot_disp *  nAtoms
    for i in range(irange):

        # first read content of output file for geom+str(step) calculation ...
        with open(os.path.join('geom' + str(step), 'molcas.log')) as enefile:
            tmp = [line.strip() for line in enefile]

        oMats = readOverlap(tmp, nroots, o_ref, o, i, command)
        o_ref = oMats[0]
        o = oMats[1]

        if i % 2 == 0:
            for z in range(len(calc_coupl)):
                j = calc_coupl[z][0]
                k = calc_coupl[z][1]
                o_plus[j][k] = o[j][k]
                o_plus[k][j] = o[k][j]
            if command[10] == '1':
                step += 1
                continue
            elif command[10] == '0':
                for z in range(len(calc_coupl)):
                    j = calc_coupl[z][0]
                    k = calc_coupl[z][1]
                    if o_plus[j][k] < 0:
                        dc[j][k] = -0.5 * (abs((o_plus[j][k] - o_ref[j][k]) / (disp)) + abs(
                            (o_plus[k][j] - o_ref[k][j]) / (disp)))
                    else:
                        dc[j][k] = 0.5 * (abs((o_plus[j][k] - o_ref[j][k]) / (disp)) + abs(
                            (o_plus[k][j] - o_ref[k][j]) / (disp)))
        else:
            if command[10] == '1':
                for z in range(len(calc_coupl)):
                    j = calc_coupl[z][0]
                    k = calc_coupl[z][1]
                    o_minus[j][k] = o[j][k]
                    o_minus[k][j] = o[k][j]
                    if (o_plus[j][k] - o_minus[j][k]) < 0:
                        dc[j][k] = -0.5 * (abs((o_plus[j][k] - o_minus[j][k]) / (2 * disp)) + abs(
                            (o_plus[k][j] - o_minus[k][j]) / (2 * disp)))
                    else:
                        dc[j][k] = 0.5 * (abs((o_plus[j][k] - o_minus[j][k]) / (2 * disp)) + abs(
                            (o_plus[k][j] - o_minus[k][j]) / (2 * disp)))
            elif command[10] == '0':
                for z in range(len(calc_coupl)):
                    j = calc_coupl[z][0]
                    k = calc_coupl[z][1]
                    o_plus[j][k] = o[j][k]
                    o_plus[k][j] = o[k][j]
                    if o_plus[j][k] < 0:
                        dc[j][k] = -0.5 * (abs((o_plus[j][k] - o_ref[j][k]) / (disp)) + abs(
                            (o_plus[k][j] - o_ref[k][j]) / (disp)))
                    else:
                        dc[j][k] = 0.5 * (abs((o_plus[j][k] - o_ref[j][k]) / (disp)) + abs(
                            (o_plus[k][j] - o_ref[k][j]) / (disp)))

        for z in range(len(calc_coupl)):
            j = calc_coupl[z][0]
            k = calc_coupl[z][1]
            dc[j][k] = dc[j][k] * scal_factors[j][k]
            if append_x == 1:
                NAC[j][k][0].append(dc[j][k])
                if z == len(calc_coupl) - 1:
                    append_x = 0
                    append_y = 1
            elif append_y == 1:
                NAC[j][k][1].append(dc[j][k])
                if z == len(calc_coupl) - 1:
                    append_y = 0
                    append_z = 1
            elif append_z == 1:
                NAC[j][k][2].append(dc[j][k])
                if z == len(calc_coupl) - 1:
                    append_z = 0
                    append_x = 1
            if int(command[210])%3 != 0:
                NAC[j][k][2].append(0.0)
                if int(command[210])%3 == 1:
                    NAC[j][k][1].append(0.0) 
            for idc in range(3):
                NAC[k][j][idc] = [-ele for ele in NAC[j][k][idc]]

        step += 1

    # projection of DC orthogonal to the translation and rotational vector seems unnecessary
    # as QM in MM embedding not translationally and rotationally invariant
    if int(command[210]) == 0:
        for i in range(geometry.NsubH):
            for z in range(len(calc_coupl)):
                j = calc_coupl[z][0]
                k = calc_coupl[z][1]
                NAC[j][k][0].append(0.0)
                NAC[k][j][0].append(0.0)
                NAC[j][k][1].append(0.0)
                NAC[k][j][1].append(0.0)
                NAC[j][k][2].append(0.0)
                NAC[k][j][2].append(0.0)

    sef = shelve.open("cobram-sef")
    sef['NAC'] = NAC
    sef.close()

    # and now the directory can be removed (gradient has already been extracted)
    for step in range(1, irange + 1):
        shutil.rmtree('geom' + str(step))
    if logwrt.VERBOSITY_LEVEL > 0:
        logwrt.writelog("\nCOMPUTED NACs:\n")
        for nac in range(len(calc_coupl)):
               j = calc_coupl[nac][0]
               k = calc_coupl[nac][1]
               logwrt.writelog("NAC between  states {0}-{1}\n".format(j+1, k+1))
               for at in range(len(NAC[j][k][0])):
                   logwrt.writelog("{0:16.12f} {1:16.12f} {2:16.12f}\n".format(NAC[j][k][0][at], NAC[j][k][1][at], NAC[j][k][2][at]))
        

    return NAC


# read overlap matrix (for numerical NACs and THS TDNACs)
def readOverlap(log, nroots, o_ref, o, igeom, command):
    sef = shelve.open("cobram-sef")
    nroots = sef['nroots']
    sef.close()
    if os.path.exists(os.getcwd() + '/SS-CASPT2'):
        sef = shelve.open("cobram-sef", 'r')
        # in case of numerical NACs the reference WF is the one at time t
        if command[14] == '0':
            cas_ref_order = sef['cas_ref_order']
            cas_order = sef['cas_order']
        # in case of numerical td-NACs the reference WF is the one at time t-dt
        if command[14] == '1':
            cas_ref_order = sef['cas_ref_order_old']
            cas_order = []
            cas_order.append([])
            cas_order[0] = sef['cas_ref_order']
        sef.close()
    for j in range(len(log)):
        if log[j].find('OVERLAP MATRIX FOR THE ORIGINAL STATES') != -1 and log[j + 2].find(
                'Diagonal, with elements') == -1:
            j += 1
            # routine for reading the lower triangular overlap matrix with arbitrary size
            # <1|1>              # upper triangular matrix corresponds to overlaps at the reference <Psi_i(r)|Psi_j(r)>
            # <2|1> <2|2>        # in the case of SS-PT2 the WFs are not orthogonal, os this should introduce a small correction to SS-PT2 NACs
            # <3|1> <3|2> <3|3>  # values are stored in o_ref
            # <1'|1> <1'|2> <1'|3> <1'|1'>          # lower square matrix (nroots x nroots) corresponds to overlap <Psi_i(r+dr)|Psi_j(r)>
            # <2'|1> <2'|2> <2'|3> <2'|1'> <2'|1'>  # values are stored in o
            # <3'|1> <3'|2> <3'|3> <3'|1'> <3'|2'>  #
            # <3'|3'>
            # currently, SS-PT2 WFs are not stored in the .JobIph or .JobMix, so <Psi_i(r)|Psi_j(r)> overlaps cosists of zeros (no harm)
            # note, that only 5 elements are printed per line
            # m runs over 2 x nroots (.eq raws in the overlap matrix)
            for m in range(2 * nroots):
                # reading next line of the output
                j += 1
                # the first nroot raws correspond to overlaps at the reference <Psi_i(r)|Psi_j(r)>
                if m < nroots:
                    # la is an auxiliary iterator which runs from 0 to 4 (i.e. over the elements in one line)
                    la = 0
                    # l runs over the 0:m+1, i.e. reading only a triangluar matrix (including diagonal elements, therefore m+1)
                    for l in range(m + 1):
                        # if the 5th elements is reached, CR to 0 and continue reading next line from the first element
                        if la >= 5:
                            j += 1
                            la = 0
                        # values are stored in matrix o_ref, which is symmetric
                        ####in case of SS-PT2 the state order can differ from the order of CASSCF states (and hence the order in the overlap matrix)
                        ###if os.path.exists(os.getcwd()+'/SS-CASPT2'):
                        ###    o_ref[cas_ref_order[m]][cas_ref_order[l]] = float(log[j].split()[la])
                        ###    o_ref[cas_ref_order[l]][cas_ref_order[m]] = o_ref[cas_ref_order[m]][cas_ref_order[l]]
                        ###else:
                        o_ref[m][l] = float(log[j].split()[la])
                        o_ref[l][m] = o_ref[m][l]
                        la += 1
                # the second nroot raws correspond to overlaps at the reference <Psi_i(r+dr)|Psi_j(r)>
                elif m >= nroots:
                    # as m runs over 2 x nroots, k maps back to a particular root
                    k = m - nroots
                    la = 0
                    # l runs over the 0:nroots+k+1, i.e. reading only a triangluar matrix (including diagonal
                    # elements, therefore m+1)
                    for l in range(m + 1):
                        if la >= 5:
                            j += 1
                            la = 0
                        if l < nroots:
                            # values are stored in matrix o, which is NOT symmetric
                            # in case of SS-PT2 the state order can differ from the order of CASSCF states
                            # (and hence the order in the overlap matrix)
                            # if os.path.exists(os.getcwd()+'/SS-CASPT2'):
                            # o[cas_order[igeom][k]][cas_ref_order[l]] = float(log[j].split()[la])
                            # else:
                            o[k][l] = float(log[j].split()[la])
                        la += 1
    if command[14] == '1' and os.path.exists(os.getcwd() + '/SS-CASPT2'):
        sef = shelve.open("cobram-sef")
        list_of_signs = sef['list_of_signs']
        #        list_of_signs2=sef['list_of_signs2']
        sef.close()
        for k in range(nroots):
            # first the overlap matrix is modified column-wise to accout for sign changes found in the previous step
            # if list_of_signs2[k] == -1:
            if list_of_signs[k] == -1:
                for m in range(nroots):
                    o[m][k] = -1 * o[m][k]
                # list_of_signs2[k]=1
                list_of_signs[k] = 1
        for k in range(nroots):
            # then we check for negtive diagonal elements and correct overlap matrix row-wise
            if o[k][k] < 0:
                # list_of_signs2[k]=-1
                list_of_signs[k] = -1
                for m in range(nroots):
                    o[k][m] = -1 * o[k][m]
        sef = shelve.open("cobram-sef")
        sef['list_of_signs_tmp'] = list_of_signs
        sef.close()
    elif command[14] == '0':
        for k in range(nroots):
            # we check for negtive diagonal elements and correct overlap matrix row-wise
            if o[k][k] < 0:
                for m in range(nroots):
                    o[k][m] = -1 * o[k][m]
    # in case of SS-PT2 the state order can differ from the order of CASSCF states (and hence the order in the ovlp mat)
    if os.path.exists(os.getcwd() + '/SS-CASPT2'):
        o_tmp = []
        o_tmp2 = []
        for i in range(nroots):
            o_tmp.append([])
            o_tmp2.append([])
        for i in range(nroots):
            for j in range(nroots):
                o_tmp[i].append(0.0)
                o_tmp2[i].append(0.0)
        for i in range(nroots):
            for j in range(i, nroots):
                sign_ij = math.copysign(1.0, o_ref[i][j])
                sign_ji = math.copysign(1.0, o_ref[j][i])
                o_tmp[i][j] = math.copysign(o_ref[cas_ref_order[i]][cas_ref_order[j]], sign_ij)
                o_tmp[j][i] = math.copysign(o_ref[cas_ref_order[j]][cas_ref_order[i]], sign_ji)
                sign_ij = math.copysign(1.0, o[i][j])
                sign_ji = math.copysign(1.0, o[j][i])
                o_tmp2[i][j] = math.copysign(o[cas_order[igeom][i]][cas_ref_order[j]], sign_ij)
                o_tmp2[j][i] = math.copysign(o[cas_order[igeom][j]][cas_ref_order[i]], sign_ji)
        o_ref = copy.deepcopy(o_tmp)
        o = copy.deepcopy(o_tmp2)

    return [o_ref, o]


# in case of 'freqxgp' molcas.molcasEneGradCrg is executed for the equilibrium computation (step = 0, par_num == 1) and once ALL SP have been computated (step > 0, par_num == 2)
# in case of 'optxgp' and 'mdvp' molcas.molcasEneGradCrg is executed in every step at the reference geometry (par_num == 1)
@Timer("EneGradCrg")
def molcasEneGradCrg(command, geometry, charges, step, par_num):
    sef = shelve.open("cobram-sef", 'r')
    DEarray = sef['DEarray']
    NAC = sef['NAC_empty']
    if int(command[85]) > 0 and command[14] == '1':
        TDC = sef['TDC_empty']
    nroots = sef['nroots']
    state = sef['newstate']
    sef.close()

    command = command
    natom = len(geometry.modelH[0])
    SCF = 'no'
    MBPT2 = 'no'
    CAS = 'no'
    SSPT2 = 'no'
    MSPT2 = 'no'
    # relevant for SCF and MBPT2 electric fields and dipole moments
    ef_occur = 0
    dip_occur = 0
    Enumber = 0
    molcasenergyS0 = 0.0
    tmp_root = 0
    EFLD = ''
    # initialize in groups of 5 to easy keeping track
    (nCI, basissetsize, closedorbs, closedels, occorbs) = (0, 0, 0, 0, 0)
    (totalels, charge, numstates, spin, symmetry) = (0, 0, 0, 0, 0)
    tmp, molcasenergy, cas_order, gradient, gradch = [], [], [], [], []
    _charges, dipole, CIhead, CHhead, CHtail = [], [], [], [], []
    casenergy, CInew, allCH, scal_factors = [], [], [], []
    for i in range(nroots):
        gradient.extend([[], [], []])
        gradch.extend([[], [], []])
        dipole.append([])
        CInew.append([])
        scal_factors.append([])
        for j in range(nroots):
            scal_factors[i].append(0.0)

    if command[101] == '1':
        NUMER = 'yes'
        if step == 0:
            logwrt.writelog("Activating internal Molcas numerics for gradients\n")
    elif int(command[9]) > 1:
        NUMER = 'no'
        if step == 0:
            logwrt.writelog("Activating Cobram numerics for gradients\n")
    else:
        logwrt.writelog("Using Molcas analytical gradients\n")
        NUMER = 'no'

    # .log file is fetched from the geomXXX
    # the output with the equilibrium calculation is required twice, therefore "step-1"
    if command[1] == 'freqxg' and par_num == 2 and os.path.exists(os.getcwd() + '/geom' + str(step - 1)):
        logwrt.writelog('Extracting QM results from geom' + str(step - 1) + ' calculation \n')
        logFileName = os.path.join('geom' + str(step - 1), "molcas.log")
    # else use the file in the root directory of the calculation
    else:
        logFileName = "molcas.log"

    # read content of the log file
    with open(logFileName) as enefile:
        for line in enefile:
            tmp.append(line.strip())

    # determine the type of computation and store dummy label on disk
    if step == 0:
        for i in range(len(tmp)):
            if tmp[i].find('MOLCAS executing module SCF') != -1 or (
                    tmp[i].find('&SCF') != -1 and tmp[i - 2].find('()()()()()()') != -1):
                system('touch ' + os.getcwd() + '/SCF')
                SCF = 'yes'
                # state=0 #GS is the default
            if tmp[i].find('MOLCAS executing module RASSCF') != -1 or (
                    tmp[i].find('&RASSCF') != -1 and tmp[i - 2].find('()()()()()()') != -1):
                system('touch ' + os.getcwd() + '/RASSCF')
                CAS = 'yes'
                SCF = 'no'
            if tmp[i].find('MOLCAS executing module MBPT2') != -1 or (
                    tmp[i].find('&MBPT2') != -1 and tmp[i - 2].find('()()()()()()') != -1):
                system('touch ' + os.getcwd() + '/MBPT2')
                SCF = 'no'
                MBPT2 = 'yes'
            if tmp[i].find('MOLCAS executing module CASPT2') != -1 or (
                    tmp[i].find('&CASPT2') != -1 and tmp[i - 2].find('()()()()()()') != -1):
                pt2_prop = 0
                system('rm ' + os.getcwd() + '/RASSCF')
                system('touch ' + os.getcwd() + '/SS-CASPT2')
                SSPT2 = 'yes'
                for ii in range(i + 1, len(tmp)):
                    if tmp[ii].find('MULTI-STATE CASPT2 SECTION') != -1:
                        system('rm ' + os.getcwd() + '/SS-CASPT2')
                        system('touch ' + os.getcwd() + '/MS-CASPT2')
                        SSPT2 = 'no'
                        MSPT2 = 'yes'
                        break
                CAS = 'no'
                SCF = 'no'
                break

    # figure out what type of calculation is done
    else:
        if os.path.exists(os.getcwd() + '/SCF'):
            SCF = 'yes'
        if os.path.exists(os.getcwd() + '/RASSCF'):
            CAS = 'yes'
            SCF = 'no'
        elif os.path.exists(os.getcwd() + '/MBPT2'):
            SCF = 'no'
            MBPT2 = 'yes'
        elif os.path.exists(os.getcwd() + '/SS-CASPT2'):
            SSPT2 = 'yes'
            pt2_prop = 0
        elif os.path.exists(os.getcwd() + '/MS-CASPT2'):
            MSPT2 = 'yes'
            pt2_prop = 0

    tmp1 = []
    # some preliminary things:
    # a) check for large orbital rotations after the last iteration
    # b) in case of 'freqxgp' or internal Molcas numerical computations check for large orbital rotations at displaced geometries
    # c) in the case of SSPT2 figure out which is the CASSCF state whose properties will be used
    # d) check for MCLR failure (SA-CASSCF only) and skip the whole readout procedure (basically perform a velocity Verlet for 2 x time-step)
    non_conv = 0
    for i in range(1, len(tmp)):
        if tmp[len(tmp) - i].find('_NOT_CONVERGED_') != -1:
            break
        if tmp[len(tmp) - i].find('Non-zero return code') != -1:
            for j in range(i + 1, len(tmp)):
                if tmp[len(tmp) - j].find('_NOT_CONVERGED_') != -1:
                    non_conv = 1
                    break
                elif tmp[len(tmp) - j].find('SPACES ARE TOO DISSIMILAR') != -1:
                    logwrt.fatalerror(
                        'Molcas returned an error. Active spaces of current and previous steps too dissimilar. ' +
                        '\nPlease check last entry of molcasALL.log')
            if non_conv == 0:
                logwrt.fatalerror('Molcas returned an error. Please check last entry of molcasALL.log')
            else:
                break

    if CAS == 'yes':
        for i in range(len(tmp)):
            # MCLR is also executed in case of NAC computation and the following line makes the difference when looking for a gradient
            if (tmp[i].find('Lagrangian multipliers are calculated for state no') != -1 or
                    tmp[i].find('Lagrangian multiplier is calculated for root no') != -1):
                k = 0
                for j in range(i + 1, len(tmp)):
                    if tmp[j].find('Iteration       Delta       Res') != -1:
                        while True:
                            k += 1
                            try:
                                if tmp[j + k].strip().split()[0].strip() == 'Warning':
                                    continue
                                if tmp[j + k].strip().split()[0].strip() == 'Perturbation':
                                    break
                                if tmp[j + k].strip().split()[0].strip() == 'No' and tmp[j + k].strip().split()[
                                    1].strip() == 'convergence':
                                    break
                            # sys.exit() raises an Exception (SystemExit), in a try: block it redirects to except!
                            except SystemExit:
                                sys.exit()
                            except:
                                pass
                        break

    if SSPT2 == 'yes' or MSPT2 == 'yes' or CAS == 'yes':
        casenergy = []
        mspt2energy = []
        unsorted_energy = []
        displ_geom = 0
        check_maxrot = 0
        cionly = 0
        for i in range(len(tmp)):
            if tmp[i].find('CI only, no orbital optimization will be done') != -1:
                cionly = 1
            if cionly == 0 and tmp[i].find('max ROT') != -1:
                k = 1
                while True:
                    k += 1

                    if tmp[i + k].find('No convergence after') != -1:
                        errFileName = os.path.join(QM_DATA_STORAGE, 'molcas_' + str(step) + '.err')
                        shutil.copy(logFileName, errFileName)
                        logwrt.fatalerror(
                            'No convergence in CASSCF numerical computation for step {0}.\n'.format(step) +
                            'Check molcas output for step {0} in {1}.'.format(step, errFileName))

                    if tmp[i + k].find('Convergence after') != -1:
                        tmp2 = tmp[i + k + 1].split()
                        try:
                            if abs(float(tmp2[6].strip('*'))) > 0.1:
                                errFileName = os.path.join(QM_DATA_STORAGE, 'molcas_' + str(step) + '.err')
                                shutil.copy(logFileName, errFileName)
                                logwrt.fatalerror(
                                    'Large orbital rotation during numerical computation for step {0}: max ROT greater than 0.1.\n'.format(
                                        step) +
                                    'Check molcas output for step {0} in {1}.'.format(step, errFileName))
                            else:
                                break
                        except SystemExit:
                            sys.exit()
                        except:
                            break
                    if tmp[i + k].find('Wave function printout') != -1:
                        break
                    if check_maxrot == 0:
                        try:
                            tmp2 = tmp[i + k].split()
                            if abs(float(tmp2[6].strip('*'))) > 0.01:
                                errFileName = os.path.join(QM_DATA_STORAGE, 'molcas_' + str(step) + '.err')
                                shutil.copy(logFileName, errFileName)
                                logwrt.writewarning(
                                    'Large orbital rotation during numerical computation for step {0}: max ROT greater than 0.01.\n'.format(
                                        step) +
                                    'Check molcas output for step {0} in {1}.'.format(step, errFileName))
                                check_maxrot = 1
                        except:
                            pass
                    if (command[1] == 'freqxg' and par_num == 2) or (NUMER == 'yes' and displ_geom > 0):
                        tmp2 = tmp[i + k].split()
                        try:
                            if abs(float(tmp2[6].strip('*'))) > 0.01:
                                if command[1] == 'freqxg' and par_num == 2:
                                    errFileName = os.path.join(QM_DATA_STORAGE, 'molcas_' + str(step) + '.err')
                                    shutil.copy(logFileName, errFileName)
                                    logwrt.writewarning(
                                        'Large orbital rotation during numerical computation for step {0}: '.format(
                                            step) +
                                        'max ROT greater than 0.01.\nCheck molcas output for step {0} in {1}.'.format(
                                            step, errFileName))
                                    break
                                elif NUMER == 'yes' and displ_geom > 0:
                                    logwrt.writewarning(
                                        'Large orbital rotation during numerical computation for displaced geometry {0}: '.format(
                                            displ_geom) +
                                        'max ROT greater than 0.01. Check molcasALL.log.')
                                    break
                        except:
                            pass
                displ_geom += 1
            if (SSPT2 == 'yes' or MSPT2 == 'yes') and tmp[i].find('root number ') != -1:
                if len(casenergy) == nroots and NUMER == 'yes':
                    pass
                else:
                    tmp1 = tmp[i].split()
                    casenergy.append(float(tmp1[7]))
                    logwrt.writelog("CASSCF energy for root {0}: {1} Hartree\n".format(tmp1[4], casenergy[-1]))
                    # if CASSCF NACs are computed analytically via Molpro during a parallel numerical run (useful for SSPT2) the GS energy
                    # is stored as reference for determining if the Molcas -> Molpro transformation did go well
                    if tmp1[4] == '1':
                        molcasenergyS0 = float(tmp1[7])
            if SSPT2 == 'yes' and tmp[i].find(' CASPT2 Root') != -1:
                if len(unsorted_energy) == nroots and NUMER == 'yes':
                    pass
                else:
                    tmp1 = tmp[i].split()
                    unsorted_energy.append(float(tmp1[6]))
                    logwrt.writelog(
                        "Unsorted SS-CASPT2 energy for root {0}: {1} Hartree\n".format(tmp1[3], unsorted_energy[-1]))
        # allow to perform computations with a smaller number of CASPT2 states than CASSCF states
        if (SSPT2 == 'yes' or MSPT2 == 'yes') and nroots != len(casenergy):
            if SSPT2 == 'yes':
                logwrt.writelog("Only {0} roots included at SS-CASPT2 level\n".format(nroots))
            elif MSPT2 == 'yes':
                logwrt.writelog("Only {0} roots included at MS-CASPT2 level\n".format(nroots))
            casenergy = casenergy[0:nroots]
        # collect and sort SS-PT2 energies
        # determine which state to read DIPOLES and EFLD for
        if SSPT2 == 'yes':
            molcasenergy = list(unsorted_energy)
            molcasenergy.sort()
            for i in range(len(molcasenergy)):
                for j in range(len(unsorted_energy)):
                    if molcasenergy[i] == unsorted_energy[j]:
                        cas_order.append(j)
            if SSPT2 == 'yes':
                logwrt.writelog(
                    "Will use PT2 corrected electric field for CASSCF state {0}: as it corresponds to SS-CASPT2 state {1}\n".format(
                        cas_order[state] + 1, state + 1))
            elif MSPT2 == 'yes':
                logwrt.writelog(
                    "Will use PT2 corrected electric field for PM-CAS state {0}\n".format(cas_order[state] + 1))
            sef = shelve.open("cobram-sef")
            if step > 0:
                sef['cas_ref_order_old'] = sef['cas_ref_order']
            sef['cas_ref_order'] = cas_order
            sef.close()
        # fill the scal_factor matrix (in SS-PT2 order)
        if SSPT2 == 'yes':
            for i in range(len(casenergy)):
                for j in range(i + 1, len(casenergy)):
                    for k in range(len(molcasenergy)):
                        if molcasenergy[k] == unsorted_energy[i]:
                            for l in range(len(molcasenergy)):
                                if molcasenergy[l] == unsorted_energy[j]:
                                    scal_factors[k][l] = abs(
                                        (casenergy[j] - casenergy[i]) / (unsorted_energy[j] - unsorted_energy[i]))
                                    scal_factors[l][k] = scal_factors[k][l]
                                    if command[1] == 'mdv' or command[1] == 'ci':
                                        logwrt.writelog("SS-PT2 scaling factor for NAC(" + str(k + 1) + "," + str(
                                            l + 1) + ") is " + str(scal_factors[k][l]) + "\n")
    for i in range(len(tmp)):
        tmp1, tmp2, tmp3, tmp4, tmp5 = [], [], [], [], []
        # collect info about the basis set
        # double check against Molcas output regarding the basis set size
        # to enter in this routine one needs to specify each atom's basis set individually
        if tmp[i].find('Basis set information') != -1 and basissetsize == 0:
            for ii in range(i + 1, len(tmp)):
                if tmp[ii].find('Basis set information') != -1:
                    break
                if tmp[ii].find('Shell  nPrim  nBasis') != -1:
                    j = 0
                    while tmp[ii + j].split():
                        j = j + 1
                        try:
                            if tmp[ii + j].split()[0] == 's':
                                basissetsize = basissetsize + int(tmp[ii + j].split()[2]) * 1
                            elif tmp[ii + j].split()[0] == 'p':
                                basissetsize = basissetsize + int(tmp[ii + j].split()[2]) * 3
                            elif tmp[ii + j].split()[0] == 'd' and command[194] == '1':
                                basissetsize = basissetsize + int(tmp[ii + j].split()[2]) * 6
                            elif tmp[ii + j].split()[0] == 'd' and command[194] != '1':
                                basissetsize = basissetsize + int(tmp[ii + j].split()[2]) * 5
                        except:
                            break
                if tmp[i].find('Basis set specifications') != -1:
                    tmp2 = tmp[i + 2].split()
                    if (basissetsize != int(tmp2[2])):
                        logwrt.fatalerror(
                            'The total number of basis functions does not agree ' + str(basissetsize) + ' vs. ' + str(
                                tmp2[2]))
                    else:
                        logwrt.writelog('Total number of basis functions: ' + str(basissetsize) + "\n")
        # collect other info which is necessary to automatically generate the Molpro input in case of a DC computation
        # a detour is necessary as some options are appearing in the CASPT2 section, where core orbitals are frozen and inactive orbitals do not match
        if (tmp[i].find('MOLCAS executing module RASSCF') != -1 or (
                tmp[i].find('&RASSCF') != -1 and tmp[i - 2].find('()()()()()()') != -1)) and closedorbs == 0:
            for j in range(i + 10, len(tmp)):
                if tmp[j].find('()()()()()()') != -1:
                    break
                if tmp[j].find('Number of inactive orbitals') != -1:
                    closedorbs = int(tmp[j].split()[4])
                    closedels = 2 * closedorbs
                    totalels = closedels + activeels
                if tmp[j].find('Number of active orbitals') != -1:
                    activeorbs = int(tmp[j].split()[4])
                    occorbs = closedorbs + activeorbs
                if tmp[j].find('Number of electrons in active shells') != -1:
                    activeels = int(tmp[j].split()[6])
                if tmp[j].find('Total molecular charge') != -1:
                    charge = float(tmp[j].split()[3])
                if tmp[j].find('Number of root(s) required') != -1:
                    numstates = int(tmp[j].split()[4])
                if tmp[j].find('Spin quantum number') != -1:
                    spin = float(tmp[j].split()[3])
                if tmp[j].find('State symmetry') != -1:
                    symmetry = int(tmp[j].split()[2])
        if CAS == 'yes' and tmp[i].find('root number ') != -1:
            if len(molcasenergy) == nroots and (NUMER == 'yes' or command[1] == 'ci'):
                pass
            else:
                tmp2 = tmp[i].split()
                molcasenergy.append(float(tmp2[7]))
                logwrt.writelog("CASSCF energy for root {0}: {1} Hartree\n".format(tmp2[4], tmp2[7]))
                # save energy of the GS (could be any state, but GS is always computed) for comparing Molcas to Molpro
                if tmp2[4] == '1':
                    molcasenergyS0 = float(tmp2[7])
        if MSPT2 == 'yes' and tmp[i].find('MS-CASPT2 Root') != -1:
            if len(molcasenergy) == nroots and NUMER == 'yes':
                pass
            else:
                molcasenergy.append(float(tmp[i].split()[6]))
                logwrt.writelog(
                    "MS-CASPT2 energy for root {0}: {1} Hartree\n".format(tmp[i].split()[3], molcasenergy[-1]))
        if SCF == 'yes' and tmp[i].find('Total SCF energy') != -1:
            if len(molcasenergy) == 1 and NUMER == 'yes':
                pass
            else:
                tmp2 = tmp[i].split()
                molcasenergy.append(float(tmp2[4]))
                logwrt.writelog("SCF energy: {0} Hartree\n".format(molcasenergy[0]))
        if MBPT2 == 'yes' and tmp[i].find('Total MBPT2 energy ') != -1:
            if len(molcasenergy) == 1 and NUMER == 'yes':
                pass
            else:
                tmp2 = tmp[i].split()
                molcasenergy.append(float(tmp2[4]))
                logwrt.writelog("MBPT2 energy: {0} Hartree\n".format(molcasenergy[0]))
        # Gradient format changed in v8.1
        if tmp[i].find('Molecular gradients') != -1 and tmp[i + 6].split()[0] == 'X' and tmp[i + 6].split()[
            1] == 'Y' and tmp[i + 6].split()[2] == 'Z':
            low_state = 0
            if not gradient[3 * state]:
                pass
            else:
                state = state - 1
                low_state = 1
            for j in range(natom):
                gr = tmp[i + j + 8].split()
                gradient[3 * state].append(-float(gr[1]))
                gradient[3 * state + 1].append(-float(gr[2]))
                gradient[3 * state + 2].append(-float(gr[3]))
            if command[120] == '1':
                j = j + 1
                for k in range(len(geometry.list_MEDIUM)):
                    gr = tmp[i + j + k + 8].split()
                    gradch[3 * state].append(float(gr[1]))
                    gradch[3 * state + 1].append(float(gr[2]))
                    gradch[3 * state + 2].append(float(gr[2]))
            if low_state == 1:
                state = state + 1
        elif tmp[i].find('Molecular gradients') != -1:
            low_state = 0
            if not gradient[3 * state]:
                pass
            else:
                state = state - 1
                low_state = 1
            for j in range(natom * 3):
                gr = tmp[i + j + 6].split()
                if gr[1] == 'x':
                    gradient[3 * state].append(-float(gr[2]))
                if gr[1] == 'y':
                    gradient[3 * state + 1].append(-float(gr[2]))
                if gr[1] == 'z':
                    gradient[3 * state + 2].append(-float(gr[2]))
            if command[120] == '1':
                j = j + 1
                for k in range(len(geometry.list_MEDIUM) * 3):
                    gr = tmp[i + j + k + 6].split()
                    if gr[1] == 'x':
                        gradch[3 * state].append(float(gr[2]))
                    if gr[1] == 'y':
                        gradch[3 * state + 1].append(float(gr[2]))
                    if gr[1] == 'z':
                        gradch[3 * state + 2].append(float(gr[2]))
            if low_state == 1:
                state = state + 1
        # from Molcas 8.1
        elif tmp[i].find('Numerical gradient,') != -1:
            tmp_state = int(tmp[i].split()[3]) - 1
            for j in range(natom):
                gr = tmp[i + j + 4].split()
                gradient[3 * tmp_state].append(-float(gr[1]))
                gradient[3 * tmp_state + 1].append(-float(gr[2]))
                gradient[3 * tmp_state + 2].append(-float(gr[3]))
            if command[120] == '1':
                j = j + 1
                for k in range(len(geometry.list_MEDIUM)):
                    gr = tmp[i + j + k + 4].split()
                    gradch[3 * tmp_state].append(float(gr[1]))
                    gradch[3 * tmp_state + 1].append(float(gr[2]))
                    gradch[3 * tmp_state + 2].append(float(gr[3]))
        elif tmp[i].find('Numerical gradient') != -1:
            low_state = 0
            if not gradient[3 * state]:
                pass
            else:
                state = state - 1
                low_state = 1
            for j in range(natom):
                gr = tmp[i + j + 4].split()
                gradient[3 * state].append(-float(gr[1]))
                gradient[3 * state + 1].append(-float(gr[2]))
                gradient[3 * state + 2].append(-float(gr[3]))
            if command[120] == '1':
                j = j + 1
                for k in range(len(geometry.list_MEDIUM)):
                    gr = tmp[i + j + k + 4].split()
                    gradch[3 * state].append(float(gr[1]))
                    gradch[3 * state + 1].append(float(gr[2]))
                    gradch[3 * state + 2].append(float(gr[3]))
            if low_state == 1:
                state = state + 1
        # Time-derivaive NACs according to the Tully-Hammes-Schaefer scheme
        if command[1] == 'mdv' and int(command[85]) > 0 and command[14] == '1' and step > 0 and tmp[i].find(
                'OVERLAP MATRIX FOR THE ORIGINAL STATES') != -1 and tmp[i + 2].find('Diagonal, with elements') == -1:
            sef = shelve.open("cobram-sef")
            # the array calc_coupl is deleted as the NACs are recovered therefore when it comes to computing the TDNACs it is empty
            # therefore a sane calc_coupl is saved as calc_tdc in the shelve
            calc_coupl = sef['calc_tdc']
            list_of_signs = sef['list_of_signs']
            if step == 1:
                o_old = []
            elif step > 1:
                o_old = sef['o_old']
            sef.close()
            if len(calc_coupl) > 0:
                o_ref, o = [], []
                for j in range(nroots):
                    o.append([])
                    o_ref.append([])
                for j in range(nroots):
                    for k in range(nroots):
                        o[j].append(0.0)
                        o_ref[j].append(0.0)
                oMats = readOverlap(tmp, nroots, o_ref, o, 0, command)
                o_ref = oMats[0]
                o = oMats[1]
                # WF propagation is done inside the overlap routine
                if SSPT2 == 'yes':
                    pass
                else:
                    # WF sign propagation through overlap
                    for k in range(nroots):
                        # first the overlap matrix is modified column-wise to accout for sign changes found in the previous step
                        if list_of_signs[k] == -1:
                            for m in range(nroots):
                                o[m][k] = -1 * o[m][k]
                            list_of_signs[k] = 1
                    for k in range(nroots):
                        # then we check for negtive diagonal elements and correct overlap matrix row-wise
                        if o[k][k] < 0:
                            list_of_signs[k] = -1
                            for m in range(nroots):
                                o[k][m] = -1 * o[k][m]
                sef = shelve.open("cobram-sef")
                sef['list_of_signs'] = list_of_signs
                sef['o_old'] = o
                sef.close()
                for k in range(nroots):
                    norm = 0
                    for m in range(nroots):
                        norm += o[k][m] ** 2
                    if norm < 0.95:
                        logwrt.writewarning(
                            "State {0} has an overlap {1} below 0.95 with states at previous step. Possibly a new state!".format(
                                k + 1, norm))
                if step > 0:
                    for z in range(len(calc_coupl)):
                        j = calc_coupl[z][0]
                        k = calc_coupl[z][1]
                        # formulas from Barbatti, Chem. Phys. 2009, 356, 147
                        if not o_old:
                            if z == 0:
                                logwrt.writelog("TDNACs computed using WFs at t and t-dt\n")
                            TDC[j][k] = (o[j][k] - o[k][j]) / (2 * float(command[83]) * float('41.341373337'))
                        else:
                            if z == 0:
                                logwrt.writelog("TDNACs computed using WFs at t, t-dt and t-2dt\n")
                            TDC[j][k] = (3 * o[j][k] - 3 * o[k][j] - o_old[j][k] + o_old[k][j]) / (
                                    4 * float(command[83]) * float('41.341373337'))
                        if SSPT2 == 'yes':
                            TDC[j][k] = TDC[j][k] * scal_factors[j][k]
                        TDC[k][j] = -TDC[j][k]
                    logwrt.writelog("Time-derivative coupling tau(" + str(step) + ")\n")
                    for k in range(nroots):
                        for l in range(nroots):
                            logwrt.writelog('%12.8f\n' % TDC[k][l])
                        logwrt.writelog("\n")
            else:
                sef = shelve.open("cobram-sef")
                o = []
                sef['o_old'] = o
                list_of_signs = []
                # list_of_signs2=[]
                for k in range(nroots):
                    list_of_signs.append(1)
                    # list_of_signs2.append(1)
                sef['list_of_signs'] = list_of_signs
                # sef['list_of_signs2']=list_of_signs2
                sef.close()
        # NACs available since Molcas v8.1.11
        if ((command[1] == 'mdv' and int(command[85]) > 0 and step != 0) or command[1] == 'ci') and (
                tmp[i].find('Non-adiabatic coupling') != -1 or tmp[i].find('Total derivative coupling') != -1):
            # NACs are requested in the input automatically according to their order calc_coupl
            # consequtively, they are read in the same order (as Molcas does not list in the output which NAC is computed)
            # as calc_coupl is re-filled later its elements are evaluated and deleted 1-by-1
            sef = shelve.open("cobram-sef")
            calc_coupl = sef['calc_coupl']
            j = calc_coupl[0][0]
            k = calc_coupl[0][1]
            calc_coupl.pop(0)
            logwrt.writelog("NAC <" + str(j + 1) + "|d/dR|" + str(k + 1) + ">\n")
            for l in range(natom):
                nac = tmp[i + l + 8].split()
                NAC[j][k][0].append(float(nac[1]))
                NAC[k][j][0].append(-float(nac[1]))
                NAC[j][k][1].append(float(nac[2]))
                NAC[k][j][1].append(-float(nac[2]))
                NAC[j][k][2].append(float(nac[3]))
                NAC[k][j][2].append(-float(nac[3]))
            for i in range(geometry.NsubH):
                NAC[j][k][0][-1 * (i + 1)] = 0.0
                NAC[k][j][0][-1 * (i + 1)] = 0.0
                NAC[j][k][1][-1 * (i + 1)] = 0.0
                NAC[k][j][1][-1 * (i + 1)] = 0.0
                NAC[j][k][2][-1 * (i + 1)] = 0.0
                NAC[k][j][2][-1 * (i + 1)] = 0.0
            for l in range(natom):
                logwrt.writelog("%12.6f" % NAC[j][k][0][l] + "%12.6f" % NAC[j][k][1][l] + "%12.6f\n" % NAC[j][k][2][l])
            sef['calc_coupl'] = calc_coupl
            sef.close()
        if CAS == 'yes' and tmp[i].lower().find('mulliken population analysis for root number') != -1:
            if len(CHhead) == nroots and (NUMER == 'yes' or command[1] == 'ci'):
                pass
            else:
                CHhead.append(i)
        # Expectation values for CASSCF
        if CAS == 'yes' and tmp[i].find('Expectation values of various properties for root number') != -1:
            if len(CHtail) == nroots and (NUMER == 'yes' or command[1] == 'ci'):
                pass
            else:
                CHtail.append(i)
                tmp1 = tmp[i].split()
                tmp_state = int(tmp1[8]) - 1
                if int(command[2]) > 0:
                    logwrt.writelog("Found electric field for the CAS state " + str(tmp_state + 1) + "\n")
                for ii in range(i + 1, len(tmp)):
                    if tmp[ii].find('Electric field:') != -1:
                        for j, crg in enumerate(charges.CRG_MEDIUM):
                            gradch[3 * tmp_state].append(float(tmp[ii + j + 2].split()[1]) * (-crg))
                            gradch[3 * tmp_state + 1].append(float(tmp[ii + j + 2].split()[2]) * (-crg))
                            gradch[3 * tmp_state + 2].append(float(tmp[ii + j + 2].split()[3]) * (-crg))
                        break
        # Expectation values for SS-CASPT2
        if (SSPT2 == 'yes') and tmp[i].find('Compute H0 matrices for state') != -1:
            if len(CHhead) == nroots and NUMER == 'yes':
                pass
            else:
                for ii in range(i + 1, len(tmp)):
                    if tmp[ii].lower().find('mulliken population analysis') != -1:
                        CHhead.append(ii)
                        break
                for ii in range(i + 1, len(tmp)):
                    if tmp[ii].find('Expectation values of various properties') != -1:
                        CHtail.append(ii)
                        pt2_prop = 1
                        break
                if pt2_prop == 0:
                    logwrt.fatalerror(
                        'Activate PROP in &CASPT2 (since Molcas 8.0 NOPROP is default) if you perform SS-PT2 dynamics.')
                tmp1 = tmp[i].split()
                # decide if this is the state of interest
                # note that for SSPT2 one has to collect the properties of the CASSCF state "cas_order[state]"
                # which upon PT2 correction becomes the state of interest "state"
                tmp_state = cas_order.index(int(tmp1[5]) - 1)
                for ii in range(i + 1, len(tmp)):
                    if tmp[ii].find('Electric field:') != -1:
                        for j, crg in enumerate(charges.CRG_MEDIUM):
                            gradch[3 * tmp_state].append(float(tmp[ii + j + 2].split()[1]) * (-crg))
                            gradch[3 * tmp_state + 1].append(float(tmp[ii + j + 2].split()[2]) * (-crg))
                            gradch[3 * tmp_state + 2].append(float(tmp[ii + j + 2].split()[3]) * (-crg))
                        break
        #in the case of X/R/MS-CASPT2 the EF for the state of interest is computed in separate SS-PT2 block 
	#only for the state of interst
        if (MSPT2 == 'yes') and tmp[i].find('Type of calculation') != -1 and tmp[i].split()[3] == 'SS-CASPT2':
            pt2_prop = 1
            for j in range(nroots):
                CHhead.append([])
                CHtail.append([])
            for ii in range(i + 1, len(tmp)):
                if tmp[ii].lower().find('mulliken population analysis') != -1:
                    CHhead[state] = ii
                if tmp[ii].find('Expectation values of various properties') != -1:
                    CHtail[state] = ii
                if tmp[ii].find('Electric field:') != -1:
                    logwrt.writelog("Electric field for state {0}\n".format(state+1))
                    for j, crg in enumerate(charges.CRG_MEDIUM):
                        gradch[3 * state].append(float(tmp[ii + j + 2].split()[1]) * (-crg))
                        gradch[3 * state + 1].append(float(tmp[ii + j + 2].split()[2]) * (-crg))
                        gradch[3 * state + 2].append(float(tmp[ii + j + 2].split()[3]) * (-crg))
                    break
        if SCF == 'yes' and tmp[i].find('Electric field:') != -1:
            if ef_occur > 0 and NUMER == 'yes':
                pass
            else:
                for j, crg in enumerate(charges.CRG_MEDIUM):
                    gradch[0].append(float(tmp[i + j + 2].split()[1]) * (-crg))
                    gradch[1].append(float(tmp[i + j + 2].split()[2]) * (-crg))
                    gradch[2].append(float(tmp[i + j + 2].split()[3]) * (-crg))
            ef_occur += 1
        if MBPT2 == 'yes' and tmp[i].find('Electric field:') != -1:
            if ef_occur > 1 and NUMER == 'yes':
                pass
            else:
                # Electric field appears twice in case of MBPT2 computations
                # we need the second output, this is why gradch is re-set
                gradch = []
                gradch.extend([[], [], []])
                for j, crg in enumerate(charges.CRG_MEDIUM):
                    gradch[0].append(float(tmp[i + j + 2].split()[1]) * (-crg))
                    gradch[1].append(float(tmp[i + j + 2].split()[2]) * (-crg))
                    gradch[2].append(float(tmp[i + j + 2].split()[3]) * (-crg))
            ef_occur += 1
        if (SCF == 'yes' or MBPT2 == 'yes') and tmp[i].find('Molecular charges') != -1:
            if len(CHhead) == 1 and NUMER == 'yes':
                pass
            else:
                CHhead.append(i)
        if MBPT2 == 'yes' and tmp[i].find('Molecular properties') != -1:
            if dip_occur > 1 and NUMER == 'yes':
                pass
            else:
                CHtail = []
                CHtail.append(i)
            dip_occur += 1
        if SCF == 'yes' and tmp[i].find('Molecular properties') != -1:
            if dip_occur > 0 and NUMER == 'yes':
                pass
            else:
                CHtail.append(i)
            dip_occur += 1
        if tmp[i].find('It is root nr') != -1:
            if len(CIhead) == 1 and NUMER == 'yes':
                pass
            else:
                CIhead.append(i)
        if nCI == 0 and tmp[i].find('Its length NCI') != -1:
            tmp1 = tmp[i].split('=')
            nCI = int(tmp1[1])
    # finish collecting data from .log file

    # get dipole moment for state of interest
    if SSPT2 == 'yes' or MSPT2 == 'yes':
        pt2_dipole = 0
    for i in range(nroots):
        if CAS == 'yes':
            dipole_line = tmp[CHtail[i] + 11].split()
        if SSPT2 == 'yes' and pt2_prop == 1:
            dipole_line = tmp[CHtail[cas_order[i]] + 10].split()
            pt2_dipole = 1
        if MSPT2 == 'yes' and pt2_prop == 1:
            try:
                dipole_line = tmp[CHtail[i] + 10].split()
                pt2_dipole = 1
            except:
                pass
        if SCF == 'yes':
            dipole_line = tmp[CHtail[i] + 7].split()
        if MBPT2 == 'yes':
            dipole_line = tmp[CHtail[i] + 7].split()
        try:
            dipole[i].append(float(dipole_line[1]) * 0.39342215569939517797)
            dipole[i].append(float(dipole_line[3]) * 0.39342215569939517797)
            dipole[i].append(float(dipole_line[5]) * 0.39342215569939517797)
            dipole[i].append(float(dipole_line[7]) * 0.39342215569939517797)
        except:
            if CAS == 'yes' or SCF == 'yes' or MBPT2 == 'yes':
                logwrt.fatalerror('Dipole moment information could not be recovered')
            elif SSPT2 == 'yes' or MSPT2 == 'yes':
                dipole[i] = [1.0, 1.0, 1.0, 1.0]
    if (SSPT2 == 'yes' or MSPT2 == 'yes') and pt2_dipole == 0 and geometry.NatomMM > 0 and \
                                                            geometry.calculationType not in ['H', 'M', 'ML']:
        logwrt.writewarning('Dipole moment information could not be recovered')
    if MSPT2 == 'yes' or CAS == 'yes' or (SSPT2 == 'yes' and command[1][0:3] != 'mdv' and command[1][0:2] != 'ci'):
        for i in range(nroots):
            for j in range(nroots):
                scal_factors[i][j] = 1.0
                scal_factors[j][i] = 1.0
    for i in range(nroots):
        CH = []
        try:
            for j in range(CHtail[i] - CHhead[i]):
                try:
                    NE = tmp[CHhead[i] + j].split()
                    if NE[0] == 'N-E':
                        for k in range(1, len(NE)):
                            CH.append(float(NE[k]))
                except:
                    pass
        except:
            pass
        _charges.append(np.array(CH))
        allCH.append(CH)

    # fill array with energy gaps
    if command[1] in ['optxg', 'ts', 'ci', 'irc', 'freqxg', 'mdv'] and (par_num == 1 or par_num == 0):
        sef = shelve.open("cobram-sef")
        sef['CInew'] = CInew
        DEarray = sef['DEarray']
        DE_oldarray = copy.deepcopy(DEarray)
        sef['DE_oldarray'] = DE_oldarray
        for i in range(nroots - 1):
            for j in range(i + 1, nroots):
                DEarray[i][j] = 627.51 * (molcasenergy[j] - molcasenergy[i])
                # make array symmetric so doesn't matter if we want ij or ji
                DEarray[j][i] = DEarray[i][j]
        sef['DEarray'] = DEarray
        sef.close()

    # decide how many states to include in the SS-PT2 computation during numerical gradient calculations
    # this is only possible with THS as only gradients are computed at displaced geometries
    # if there is a close lying CASSCF state it may swap with the state of interest along one displacement so both should be considered in the numerics
    ########b) that a state may swap with the state of interest along one displacement upon PT2 correction
    if SSPT2 == 'yes' and int(command[60]) > 1 and command[205] == '1' and \
                ((command[1] == 'mdv' and command[14] == '1') or command[1] in ['optxg', 'freqxg', 'irc', 'ts']):
        sort_states = []
        unsort_states = []
        # the threshold is estimated like this: a change in the energy of 0.002 Hartree for a 0.001A displacement equals to ~1.0 Hartree/Bohr
        # in the worst case scenario the neighbouring state will shift in the opposite direction by the same value
        # thus, if the gap at the reference geometry is > 0.004 we can assume that the states will not swap during the numerical displacement
        # the large the numerical displacement, the large the energy change is, therefore the threshold is adapted to command[12]
        gap_thres = float(command[204]) * float(command[12]) / 0.001
        DEarray2 = []
        for i in range(nroots):
            DEarray2.append([])
        for i in range(nroots):
            for j in range(nroots):
                DEarray2[i].append(1000.0)
        for i in range(nroots - 1):
            for j in range(i + 1, nroots):
                if SSPT2 == 'yes':
                    DEarray2[i][j] = 627.51 * (casenergy[j] - casenergy[i])
                elif MSPT2 == 'yes':
                    DEarray2[i][j] = 627.51 * (mspt2energy[j] - mspt2energy[i])
                DEarray2[j][i] = DEarray2[i][j]
        # check if both lower and both higher states (if available) are < thres
        cas_state = cas_order[state]
        for istate in range(cas_state - 2, cas_state + 3):
            if istate == cas_state:
                unsort_states.append(cas_state)
            elif istate >= 0 and istate < nroots:
                if abs(DEarray2[cas_state][istate]) < gap_thres:
                    unsort_states.append(istate)
        if SSPT2 == 'yes':
            logwrt.writelog("Following CASSCF states will be considered in numerical computations: {0}\n".format(
                [x + 1 for x in unsort_states]))
        elif MSPT2 == 'yes':
            logwrt.writelog("Following MS-CASPT2 states will be considered in numerical computations: {0}\n".format(
                [x + 1 for x in unsort_states]))
        for istate in unsort_states:
            sort_states.append(cas_order.index(istate))
        # sort_states=list(set(sort_states))
        sort_states.sort()
        if SSPT2 == 'yes':
            logwrt.writelog("SS-CASPT2 order of CASSCF states considered for numerical computations: {0}\n".format(
                [x + 1 for x in sort_states]))
        elif MSPT2 == 'yes':
            logwrt.writelog("MS-CASPT2 order states considered for numerical computations: {0}\n".format(
                [x + 1 for x in sort_states]))
        sef = shelve.open("cobram-sef")
        sef['sort_states'] = sort_states
        sef['unsort_states'] = unsort_states
        sef.close()

    if command[1] == 'mdv' and (par_num == 1 or par_num == 0) and int(command[85]) > 0:
        sef = shelve.open("cobram-sef")
        # check energy difference for states below and above the state of interest
        # and estimate which couplings have to be computed
        calc_coupl = []
        calc_coupl2print = []
        logwrt.writelog("\nThe threshold for activating Tully's FSSH is {} kcal/mol\n".format(command[86]))
        if int(command[201]) == 1:
            # offd should not be overwirrten as it holds information which state besides the active state have coefficients != 0
            if step == 0:
                offd = []
            else:
                offd = sef['offd']
        # else use NAC scheme 1
        else:
            offd = []
        for i in range(nroots):
            if i != state and abs(DEarray[i][state]) <= float(command[86]):
                logwrt.writelog("The coupling between states {0} and {1} included (DE = {2:.2f} kcal/mol)\n".format(
                    i + 1, state + 1, DEarray[i][state]))
                calc_coupl.append([i, state])
                calc_coupl2print.append([i + 1, state + 1])
                if i not in offd:
                    offd.append(i)
        offd.sort()
        if int(command[201]) == 1:
            # remove state from offd (happens immediately after hopping)
            try:
                offd.remove(state)
            except:
                pass
            # modify the way couplings between states != state are being computed
            if len(offd) > 0:
                # print "updated offd", offd
                for i in offd:
                    for j in range(i):
                        if j != state and abs(DEarray[i][j]) <= float(command[86]):
                            logwrt.writelog(
                                "The coupling between states {0} and {1} included (DeltaE = {2} kcal/mol)\n".format(
                                    i + 1, j + 1, DEarray[i][j]))
                            calc_coupl.append([i, j])
                            calc_coupl2print.append([i + 1, state + 1])
        else:
            logwrt.writelog("Also computing coupling: ")
            if len(offd) > 1:
                # print "updated offd", offd
                for i in range(len(offd)):
                    for j in range(i + 1, len(offd)):
                        logwrt.writelog("{0} {1}; ".format(offd[i] + 1, offd[j] + 1))
                        calc_coupl.append([offd[i], offd[j]])
                        calc_coupl2print.append([offd[i] + 1, offd[j] + 1])
            logwrt.writelog("\n")
        if int(command[201]) == 1:
            sef['offd'] = offd
        sef['calc_coupl'] = calc_coupl
        sef['calc_tdc'] = calc_coupl
        sef.close()

    if command[1] == 'mdv' and (par_num == 1 or par_num == 0) and int(command[85]) > 0:
        sef = shelve.open("cobram-sef")
        smart_numerics = sef['smart_numerics']
        sef.close()
        if smart_numerics == 1:
            if len(calc_coupl) != 0:
                command[10] = '1'
            else:
                command[10] = '0'
            if command[10] == '0':
                logwrt.writelog("Gradients Will be computed by + displacement only!\n")
            else:
                logwrt.writelog("Gradients and NACs Will be computed by +/- displacements!\n")

    newstate = state
    ref_energy = list(molcasenergy)

    # in case of CASSCF or PT2 computation select the energy of the state of interest
    if CAS == 'yes' or SSPT2 == 'yes' or MSPT2 == 'yes':
        # in case of mdv the values must be stored as the state might change
        # correct way to create a new copy of a list
        energy = list(molcasenergy)  # CBF.Totalenergy requires an array with all energies
    else:
        # CBF.Totalenergy requires an arrays with all energies,
        # so if there is only one energy one should still supply an array
        energy = []
        energy.append(molcasenergy[0])

    if command[1] == 'ci':
        gradient2 = [np.array(gradient[3 * (state - 1)]), np.array(gradient[3 * (state - 1) + 1]),
                     np.array(gradient[3 * (state - 1) + 2])]
        gradch2 = [np.array(gradch[3 * (state - 1)]), np.array(gradch[3 * (state - 1) + 1]),
                   np.array(gradch[3 * (state - 1) + 2])]
        dipole2 = dipole[state - 1]
        charges2 = _charges[state - 1]
    gradient = [np.array(gradient[3 * state]), np.array(gradient[3 * state + 1]), np.array(gradient[3 * state + 2])]
    gradch = [np.array(gradch[3 * state]), np.array(gradch[3 * state + 1]), np.array(gradch[3 * state + 2])]
    dipole = dipole[state]
    _charges = _charges[state]
    with open('STATE', 'w') as eneout:
        eneout.write(str(newstate))

    Selfenergy = 0.0
    # a bunch of data is supplied to molcasMDV: basis set data, data for preparing the Molpro input, etc.
    aux_data = [NAC, CInew, state, newstate, molcasenergy, molcasenergyS0, nCI, CIhead, basissetsize, closedorbs,
                closedels, occorbs, totalels, charge, numstates, spin, symmetry, _charges, scal_factors, ref_energy]
    Results = [energy, gradient, _charges, Selfenergy, dipole, gradch]

    sef = shelve.open("cobram-sef")
    sef['aux_data'] = aux_data
    sef['NAC'] = NAC
    if int(command[85]) > 0 and command[14] == '1':
        sef['TDC'] = TDC
    sef['DEarray'] = DEarray
    sef.close()
    # in case of THS (14 == 1) the scaling will be performed along the gradient difference vector (206 == 0)
    # which must be computed taking into account the redistribution of the link atom gradients on the H and M atoms
    if int(command[85]) > 0 and command[14] == '1' and os.path.exists(os.getcwd() + '/HOP') and command[206] == '0' and par_num == 0:
        sef = shelve.open("cobram-sef")
        grad_state = sef['gradient']
        gradch_state = sef['gradch']
        sef.close()
        computeGD(grad_state, gradient, gradch_state, gradch, geometry, command)
    sef = shelve.open("cobram-sef")
    if step > 0:
        if par_num == 0:
            sef['gradient2'] = sef['gradient']
        sef['gradch2'] = sef['gradch']
        sef['dipole2'] = sef['dipole']
        sef['charges2'] = sef['charges']
    if par_num == 0:
        sef['gradient'] = gradient
    sef['gradch'] = gradch
    sef['dipole'] = dipole
    sef['charges'] = _charges
    if command[1] == 'ci':
        sef['gradient2'] = gradient2
        sef['gradch2'] = gradch2
        sef['dipole2'] = dipole2
        sef['charges2'] = charges2
    sef['molcasenergy'] = molcasenergy
    sef['casenergy'] = casenergy
    sef['newstate'] = newstate
    sef.close()

    logwrt.writelog("\n")
    return Results


def molcasMDV(command, geometry, charges, cobcom, step, par_num):
    sef = shelve.open("cobram-sef", 'r')
    aux_data = sef['aux_data']
    NAC = sef['NAC']
    DEarray = sef['DEarray']
    molcasenergy = sef['molcasenergy']
    state = sef['state']
    newstate = sef['newstate']
    nroots = int(sef['nroots'])
    CIold = sef['CIold']
    CInew = sef['CInew']
    calc_coupl = sef['calc_coupl']
    old_step = sef['old_step']
    sef.close()

    (molcasenergyS0, nCI, CIhead, basissetsize, closedorbs, closedels, occorbs, totalels, charge, numstates, spin,
     symmetry, _charges, scal_factors, ref_energy) = (
        aux_data[5], aux_data[6], aux_data[7], aux_data[8], aux_data[9], aux_data[10], aux_data[11], aux_data[12],
        aux_data[13], aux_data[14], aux_data[15], aux_data[16], aux_data[17], aux_data[18], aux_data[19])
    natom = len(geometry.modelH[0])
    command = command
    CIrot = 'no'
    doCIvec = 'no'
    if step == 0:
        s = []
        ttry = 'long'
        system('echo ' + str(ttry) + '>TSTEP')
        CIold = list(CInew)
        # with methods which do not rely on NACs (e.g. CI vector rotation, dipole moments, charges) there
        # is no need of a short step close to the crossing region
        if int(command[85]) > 0:
            SItime = float(0)
            sef = shelve.open("cobram-sef")
            sef['SItime'] = SItime
            sef.close()
            for i in range(nroots):
                if abs(DEarray[state][i]) < float(command[86]):
                    ttry = 'short'
                    break
            system('echo ' + str(ttry) + '>TSTEP')
            # print 'Initializing Amplitudes'
            AM = []
            for i in range(nroots):
                AM.append(complex(0.0))
            AM[state] = complex(1.0)
            AM = np.array(AM)
            AM1 = copy.deepcopy(AM)
            sef = shelve.open("cobram-sef")
            sef['AM'] = AM
            sef['AM1'] = AM1
            sef['SItime'] = SItime
            sef.close()
            # append amplitudes from this step to list of amplitdues
            # this is required for spectroscopy
            out = open('AmplitudesALL.dat', 'w')
            out.write('%10.4f' % SItime + '   '),
            for i in range(nroots):
                out.write('%16.12f' % AM[i].real + ' %16.12f ' % AM[i].imag + '   '),
            out.write('\n')
            out.close()
        sef = shelve.open("cobram-sef")
        sef['ttry_old'] = ttry
        sef.close()
    else:
        if command[85] == '1':
            if os.path.exists(os.getcwd() + '/SALTO') and int(command[102]) > 0:
                if os.popen('cat SALTO').read().strip().split()[3] == '0':
                    system('rm SALTO')
                    doCIvec = 'yes'
                else:
                    steps_wo_CI = int(os.popen('cat SALTO').read().strip().split()[3])
                    logwrt.writelog("The CI overlap will not be computed for another {0} steps\n".format(steps_wo_CI))
                    system('echo "steps without hop ' + str(steps_wo_CI - 1) + '" >SALTO')
                    CIold = list(CInew)
            elif os.path.exists(os.getcwd() + '/SALTO') and int(command[102]) == 0 and state != 0:
                system('rm SALTO')
                doCIvec = 'yes'
            elif not os.path.exists(os.getcwd() + '/SALTO') and int(command[102]) > 0:
                doCIvec = 'yes'
            elif not os.path.exists(os.getcwd() + '/SALTO') and int(command[102]) == 0 and state != 0:
                doCIvec = 'yes'
            if doCIvec == 'yes':
                CIF1 = 0.0
                CIF2 = 0.0
                CIF3 = 0.0
                CIF4 = 0.0
                # state starts at 0, nroots starts at 1
                if state != 0 and state < nroots - 1:
                    CIF1 = abs(np.add.reduce(np.array(CInew[state]) * np.array(CIold[state + 1])))
                    CIF2 = abs(np.add.reduce(np.array(CInew[state + 1]) * np.array(CIold[state])))
                    CIF3 = abs(np.add.reduce(np.array(CInew[state]) * np.array(CIold[state - 1])))
                    CIF4 = abs(np.add.reduce(np.array(CInew[state - 1]) * np.array(CIold[state])))
                elif state != 0 and state == nroots - 1:
                    CIF3 = abs(np.add.reduce(np.array(CInew[state]) * np.array(CIold[state - 1])))
                    CIF4 = abs(np.add.reduce(np.array(CInew[state - 1]) * np.array(CIold[state])))
                elif state == 0 and DEarray[0][1] <= float(command[86]):
                    CIF1 = abs(np.add.reduce(np.array(CInew[state]) * np.array(CIold[state + 1])))
                    CIF2 = abs(np.add.reduce(np.array(CInew[state + 1]) * np.array(CIold[state])))
                if int(command[102]) > 0:
                    if state < nroots - 1:
                        logwrt.writelog(
                            'CI overlap <' + str(state + 1) + '(t)|' + str(state + 2) + '(t-dt)> = %4.2f\n' % CIF1)
                        logwrt.writelog(
                            'CI overlap <' + str(state + 2) + '(t)|' + str(state + 1) + '(t-dt)> = %4.2f\n' % CIF2)
                        logwrt.writelog('Energy gap: {0}\n'.format(DEarray[state][state + 1]))
                        logwrt.writelog(
                            'CI overlap <' + str(state + 1) + '(t)|' + str(state) + '(t-dt)> = %4.2f\n' % CIF3)
                        logwrt.writelog(
                            'CI overlap <' + str(state) + '(t)|' + str(state + 1) + '(t-dt)> = %4.2f\n' % CIF4)
                        logwrt.writelog('Energy gap: {0}\n'.format(DEarray[state][state - 1]))
                    elif state == 0:
                        logwrt.writelog(
                            'CI overlap <' + str(state + 1) + '(t)|' + str(state + 2) + '(t-dt)> = %4.2f\n' % CIF1)
                        logwrt.writelog(
                            'CI overlap <' + str(state + 2) + '(t)|' + str(state + 1) + '(t-dt)> = %4.2f\n' % CIF2)
                        logwrt.writelog('Energy gap: {0}\n'.format(DEarray[state][state + 1]))
                    elif state == nroots - 1:
                        logwrt.writelog(
                            'CI overlap <' + str(state + 1) + '(t)|' + str(state) + '(t-dt)> = %4.2f\n' % CIF3)
                        logwrt.writelog(
                            'CI overlap <' + str(state) + '(t)|' + str(state + 1) + '(t-dt)> = %4.2f\n' % CIF4)
                        logwrt.writelog('Energy gap: {0}\n'.format(DEarray[state][state - 1]))
                else:
                    logwrt.writelog('CI overlap <' + str(state + 1) + '(t)|' + str(state) + '(t-dt)> = %4.2f\n' % CIF3)
                    logwrt.writelog('CI overlap <' + str(state) + '(t)|' + str(state + 1) + '(t-dt)> = %4.2f\n' % CIF4)
                    logwrt.writelog('Energy gap: {0}\n'.format(DEarray[state][state - 1]))
                CIold = list(CInew)
                # if int(command[102]) > 0 and abs(CIF1) >= 0.25 and abs(CIF2) >= 0.25 :
                # Thiel's conditions (deltaE < 30 kcal/mol, averaged overlap > 0.5)
                if int(command[102]) > 0 and (CIF1 + CIF2) / 2 > 0.5 and DEarray[state][state + 1] < float(
                        command[86]):
                    logwrt.writelog('Up-hop from state ' + str(state + 1) + ' to state ' + str(
                        state + 2) + ' using CI coeff. rotation\n')
                    newstate = state + 1
                    logwrt.writelog('-----------------------------\n')
                    logwrt.writelog('  !!!!Gimme Hop Joanna!!!\n')
                    logwrt.writelog('-----------------------------\n')
                    system('echo "steps without hop ' + command[102] + '" >SALTO')
                    system('touch HOP')
                # Thiel's conditions (deltaE < 30 kcal/mol, averaged overlap > 0.5)
                # elif abs(CIF3) >= 0.25 and abs(CIF4) >= 0.25 :
                elif (CIF3 + CIF4) / 2 > 0.5 and DEarray[state][state - 1] < float(command[86]):
                    logwrt.writelog('Down-hop from state ' + str(state + 1) + ' to state ' + str(
                        state) + ' using CI coeff. rotation\n')
                    newstate = state - 1
                    logwrt.writelog('-----------------------------\n')
                    logwrt.writelog('  !!!!Gimme Hop Joanna!!!\n')
                    logwrt.writelog('-----------------------------\n')
                    if int(command[102]) > 0:
                        system('echo "steps without hop ' + command[102] + '" >SALTO')
                    elif int(command[102]) == 0 and state != 0:
                        system('echo "hopping to a lower state" >SALTO')
                    elif int(command[102]) == 0 and state == 0:
                        pass
                    system('touch HOP')
                else:
                    logwrt.writelog("No CI rotation, WF is the same!\n")
                    newstate = state
                system('echo ' + str(newstate) + ' > states')
        # after a HOP molcasMDV is executed once more in case of Molpro couplings and state and newstate differ
        elif step == old_step + 1:
            newstate = state

        if command[14] == '2' and len(
                calc_coupl) > 0:  # when NACs are computed numerically key 14 == '0' and this procedure is skipped
            logwrt.writelog(
                'The energy gap threshold of ' + str(float(command[86])) + ' was crossed\nTully is activated\n' +
                'NACs will be computed with Molpro,Molcas MOs will be provided as initial guess\n')
            molpro_ele_order = ['O', 'N', 'C', 'H', 'F', 'P', 'S']  # atom order defined in molpro.py
            elements = geometry.modelH[0]
            # the 2D list tmp_orbitallist contains the MO coefficients in the Molcas order
            # the 2D list orbitallist contains the MO coefficients in the Molpro order
            tmp_orbitallist = [[0 for cols in range(basissetsize)] for rows in range(basissetsize)]
            orbitallist = [[0 for cols in range(basissetsize)] for rows in range(basissetsize)]

            # read MOs from molcas.RasOrb to a temporary array tmp_orbitallist
            enefile = open('molcas.RasOrb', 'r')
            tmp = []
            for line in enefile:
                tmp.append(line.strip('\n'))
            enefile.close()
            icount = -1
            # hard-coded format of molcas.RasOrb differs in versions INPORB 1.1 (18 char) and INPORB 2.0 (22 char)
            orblength = 22
            for iline in range(len(tmp)):
                if tmp[iline].find('#INPORB') != -1 and tmp[iline].split()[1].strip() == '1.1':
                    orblength = 18
                if tmp[iline].find('#ORB') != -1:
                    for jline in range(iline + 1, len(tmp)):
                        if tmp[jline].find('ORBITAL') != -1:
                            icount = icount + 1
                            start = 0
                        elif tmp[jline].find('#OCC') != -1:
                            break
                        else:
                            tmp2 = tmp[jline]
                            for i in range(0, len(tmp2), orblength):
                                tmp_orbitallist[icount][start + i / orblength] = float(
                                    tmp2[i:i + orblength])  # hard coded format of the molcas.RasOrb
                            start = start + (len(tmp2) / orblength)

            # obtain basis sets in use
            basisset_info = readbassisset(geometry, command, cobcom)
            for i in range(natom):
                if basisset_info[i][0] == 'ANO-L' or basisset_info[i][0] == 'ANO-S':
                    logwrt.fatalerror("Basis set " + str(basisset_info[i][0]) + " not supported for NAC computation")

            # reorder the atoms according to molpro_ele_order and store order in the list molpro_at_order
            molcas_at_order = []
            molpro_at_order = []
            for iel in range(len(elements)):
                molcas_at_order.append(iel)
            for iel in range(len(molpro_ele_order)):
                for i, j in enumerate(elements):
                    if j == molpro_ele_order[iel]:
                        molpro_at_order.append(i)

            # reorder atoms and AO coefficients using a) molpro_at_order and b) the predefined basis
            # set map for each row; save to orbitallist
            for kk in range(basissetsize):
                for ii in range(len(elements)):
                    ll = molpro_at_order[ii]
                    basis_label = basisset_info[ll][0]
                    # obtain the postition the first AO coefficient of element 'ii' HAS in the molcasorb
                    start1 = determine_position(molcas_at_order, elements, basisset_info, int(command[194]), ll, 0)
                    # obtain the postition the first AO coefficient of element 'ii' WILL HAVE in the molproorb
                    start2 = determine_position(molpro_at_order, elements, basisset_info, int(command[194]), ii, 0)
                    map_array = basisset_mapping(basis_label, int(command[194]))
                    map_1st = map_array[0]
                    map_2nd = map_array[1]
                    map_3rd = map_array[2]
                    if elements[ll] == 'H':
                        for jj in range(len(map_1st)):
                            # reorder (for each orbital kk and for each coefficient jj)
                            orbitallist[kk][start2 + jj] = tmp_orbitallist[kk][start1 + map_1st[jj] - 1]
                    elif elements[ll] == 'C' or elements[ll] == 'N' or elements[ll] == 'O' or elements[ll] == 'F':
                        for jj in range(len(map_2nd)):
                            orbitallist[kk][start2 + jj] = tmp_orbitallist[kk][start1 + map_2nd[jj] - 1]
                            # With Cartesian basis sets, the elements xx,yy and zz have different definitions
                            # in Molpro and Molcas, multiplication with sqrt(3) is required
                            if ((basis_label == '6-31Gp' or basis_label == '6-31Gpp') and (
                                    jj == 9 or jj == 10 or jj == 11) and (int(command[194]) == 1)):
                                orbitallist[kk][start2 + jj] = orbitallist[kk][start2 + jj] * math.sqrt(3)
                    elif (elements[ll] == 'Si' or elements[ll] == 'P' or elements[ll] == 'S' or elements[ll] == 'Cl'):
                        for jj in range(len(map_3rd)):
                            orbitallist[kk][start2 + jj] = tmp_orbitallist[kk][start1 + map_3rd[jj] - 1]
                            if ((basis_label == '6-31Gp' or basis_label == '6-31Gpp') and (
                                    jj == 14 or jj == 17 or jj == 19) and (int(command[194]) == 1)):
                                orbitallist[kk][start2 + jj] = orbitallist[kk][start2 + jj] * math.sqrt(3)

            # create molproorb.dat
            # the MO array has to be transposed!!!
            fileout = open('molcasorb.dat', 'w')
            fileout.write('BEGIN_DATA,\n')
            for kk in range(basissetsize):
                for ll in range(basissetsize):
                    if (ll + 1) % 5 == 0:
                        tmp = '%15.8f,' % orbitallist[ll][kk]
                        fileout.write(tmp + '\n')
                    else:
                        tmp = '%15.8f,' % orbitallist[ll][kk]
                        fileout.write(tmp)
                if basissetsize % 5 != 0:
                    fileout.write('\n')
            fileout.write('END_DATA,\n')
            fileout.close()

            threshold = float(command[86])
            newcoupl = molpro.chk_coupl(newstate, threshold, DEarray, calc_coupl)
            if newcoupl != []:
                calc_coupl = copy.deepcopy(newcoupl)
            # prepare the molpro input for CIonly run
            # use all the additional info collected previously to construct the input automatically
            gen = []
            keys = [closedorbs, closedels, occorbs, totalels, charge, nroots, spin, symmetry]
            molpro.makeinp(command, geometry, charges, step)
            molpro.launch(command, step)

            # read NACs and re-scale them in case of SS-PT2
            molproenergyS0 = molpro.readDC(command, geometry, step)
            # compare Molcas and Mopro energies for S0
            if abs(molproenergyS0 - molcasenergyS0) > 1e-5:
                logwrt.fatalerror('Molpro and Molcas energies of GS differ by ' +
                                  str(abs(
                                      molproenergyS0 - molcasenergyS0)) + ' which is above the internal threshold of 1e-5')
            else:
                logwrt.writelog(
                    "Transfer of WF from Molcas to Molpo was successful, Molpro and Molcas energies of GS differ by " +
                    str(abs(molproenergyS0 - molcasenergyS0)) + ' (threshold 1e-5)')
        # do not enter Tully directly after HOP when molcasMDV is repeated
        if int(command[85]) > 0 and step == old_step + 1:
            newstate = tully(command, geometry, state, cobcom, step)
            logwrt.writelog("State after propagating the coefficients c_i in Tully FSSH {0}\n".format(newstate + 1))

    # save CIold to sef
    sef = shelve.open("cobram-sef")
    sef['CIold'] = CIold
    sef['newstate'] = newstate
    sef.close()

    # try NOT to write new seward file, instead changing only cobram.command
    if str(newstate) != str(state) and step == old_step + 1:
        enefile = input('cobram.command')
        tmp = []
        for line in enefile:
            tmp.append(line.strip())
        enefile.close()
        for i in range(len(tmp)):
            # add some more flexibility to recognize user input
            if tmp[i].strip().lower()[0:4] == 'rlxr':
                try:
                    int(tmp[i + 1])
                    tmp[i + 1] = str(newstate + 1)
                except:
                    tmp[i] = 'rlxr=' + str(newstate + 1)
            if tmp[i].strip().lower()[0:4] == 'sala':
                try:
                    int(tmp[i + 1])
                    tmp[i + 1] = str(newstate + 1)
                except:
                    tmp[i] = 'sala=' + str(newstate + 1)
        filein = open('cobram.command', 'w')
        for i in range(len(tmp)):
            filein.write(str(tmp[i]) + '\n')
        filein.close()

    if command[1] == 'mdv' and command[203] != '0' and float(command[203]) < abs(
            DEarray[0][1]) and state == 0 and not os.path.exists(os.getcwd() + '/MBPT2'):
        logwrt.writelog("\n#############\nMD in the GS! Switching to MP2!\n#############\n")
        if os.path.exists(os.getcwd() + '/RASSCF'):
            system('rm ' + os.getcwd() + '/RASSCF')
        elif os.path.exists(os.getcwd() + '/SS-CASPT2'):
            system('rm ' + os.getcwd() + '/SS-CASPT2')
        elif os.path.exists(os.getcwd() + '/MS-CASPT2'):
            system('rm ' + os.getcwd() + '/MS-CASPT2')
        system('touch ' + os.getcwd() + '/MBPT2')
        sef = shelve.open("cobram-sef")
        sef['nroots'] = 1
        sef['MP2active'] = 1
        sef.close()
        enefile = input('cobram.command')
        tmp = []
        for line in enefile:
            tmp.append(line.strip())
        enefile.close()
        for i in range(len(tmp)):
            try:
                if tmp[i].strip().lower()[0:7] == '!molcas':
                    startline = i
                if tmp[i].strip().lower()[0:7] == '?molcas':
                    endline = i
                    break
            except:
                pass
        filein = open('cobram.command', 'w')
        for i in range(startline):
            filein.write(str(tmp[i]) + '\n')
        for i in range(endline + 1, len(tmp)):
            filein.write(str(tmp[i]) + '\n')
        filein.write('\n!molcas\n&SCF\nCHARGE=' + str(int(charge)) + '\n\n&MBPT2\n?molcas\n')
        filein.close()
    cobcom = CBF.getCobrammCommand('cobram.command')

    # Brute SA switch off
    if (command[93] != '0') and (state == 0) and step == old_step + 1:
        try:
            wchk = open('wchk', 'r')
            wstep = wchk.read()
        except:
            wchk = open('wchk', 'w')
            wstep = 1
            wchk.write(str(wstep) + '\n')
        wchk.close()
        if int(wstep) > 0:
            if abs(float(DEarray[0][1])) > float(command[93]):
                enefile = open('cobram.command', 'r')
                tmp = []
                tmp1 = []
                for line in enefile:
                    tmp.append(line.strip())
                enefile.close()
                wstep = '0'
                wchk = open('wchk', 'w')
                wchk.write(str(wstep) + '\n')
                wchk.close()
                i = 0
                while (tmp[i].find('?molcas') != 0):
                    i = i + 1
                    # print i,tmp[i]
                    try:
                        tmp1 = tmp[i].split()
                        tmp2 = tmp1[0].lower()
                        if tmp2.find('ciro') == 0:
                            # print "FOUND CIRO"
                            tmp[i + 1] = '1 1'
                            tmp[i + 2] = '1'
                            del tmp[i + 3]
                        if tmp2.find('&rassi') == 0:
                            # print "FOUND RASSI SECTION, deleting!"
                            while (tmp[i].find('End') != 0):
                                del tmp[i]
                            tmp[i] = ''
                            i = i - 1
                        if tmp2.find('&mclr') == 0:
                            # print "FOUND MCLR SECTION, deleting!"
                            while (tmp[i].find('End') != 0):
                                del tmp[i]
                            tmp[i] = ''
                            i = i - 1
                    except:
                        pass
                # change command file
                filein = open('cobram.command', 'w')
                for i in range(len(tmp)):
                    filein.write(str(tmp[i]) + '\n')
                filein.close()
                # we need to clean up here before we go on
                os.system('cp ' + os.getcwd() + '/MOLCAS/molcas.RasOrb' + os.getcwd() + '/INPORB ')
                os.system('rm -r ' + os.getcwd() + '/MOLCAS/*')
                os.system('cp ' + os.getcwd() + '/INPORB ' + os.getcwd() + '/MOLCAS/molcas.RasOrb')

    # if hopping occured and state changed we need to recompute GRADs
    if str(newstate) != str(state) and step == old_step + 1:
        logwrt.writelog("A hop occured, gradients and NACs need to be recomputed for the new state {0}\n".format(
            newstate + 1))
        # use DEarray with newstate to define new coupl which will be computed either numerically or ananlytically
        if int(command[85]) > 0:
            sef = shelve.open("cobram-sef")
            # check energy difference for states below and above the state of interest
            # and estimate which couplings have to be computed
            calc_coupl = []
            calc_coupl2print = []
            logwrt.writelog("\nThe threshold for activating Tully's FSSH is {0} kcal/mol\n".format(command[86]))
            if int(command[201]) == 1:
                # offd should not be overwirrten as it holds information which state besides the active state have coefficients != 0
                if step == 0:
                    offd = []
                else:
                    offd = sef['offd']
            # else use NAC scheme 1
            else:
                offd = []
            for i in range(nroots):
                if i != newstate and abs(DEarray[i][newstate]) <= float(command[86]):
                    logwrt.writelog("The coupling between states {0} and {1} included (DE = {2:.2f} kcal/mol)\n".format(
                        i + 1, newstate + 1, DEarray[i][newstate]))
                    calc_coupl.append([i, newstate])
                    calc_coupl2print.append([i + 1, newstate + 1])
                    if i not in offd:
                        offd.append(i)
            offd.sort()
            if int(command[201]) == 1:
                # remove newstate from offd (happens immediately after hopping)
                try:
                    offd.remove(newstate)
                except:
                    pass
                # modify the way couplings between states !=new state are being computed
                if len(offd) > 0:
                    for i in offd:
                        for j in range(i):
                            if j != newstate and abs(DEarray[i][j]) <= float(command[86]):
                                logwrt.writelog(
                                    "The coupling between states {0} and {1} included (DeltaE = {2} kcal/mol)\n".format(
                                        i + 1, j + 1, DEarray[i][j]))
                                calc_coupl.append([i, j])
                                calc_coupl2print.append([i + 1, newstate + 1])
            else:
                logwrt.writelog("Also computing coupling: ")
                if len(offd) > 1:
                    # print "updated offd", offd
                    for i in range(len(offd)):
                        for j in range(i + 1, len(offd)):
                            logwrt.writelog("{0} {1}; ".format(offd[i] + 1, offd[j] + 1))
                            calc_coupl.append([offd[i], offd[j]])
                            calc_coupl2print.append([offd[i] + 1, offd[j] + 1])
            logwrt.writelog("\n")
            if int(command[201]) == 1:
                sef['offd'] = offd
            sef['calc_coupl'] = calc_coupl
            sef.close()

        sef = shelve.open("cobram-sef")
        sef['old_step'] = step
        sef.close()
        # numerically (if par_num == 3) ...
        if par_num == 3:
            # practically parallel_numerics.run is repeated up to the MDV part
            # right now redundancies: reference computation is repeated (unnecessary), as well as displacements
            # In the case of MS-PT2 and CASSCF the data is already available BUT in the case of PT2 one might need to compute the PT2 energies of the new states
            par_num = 1
            # print "Getting QM results inside MDV"
            QM_Results = CBF.QM(command, cobcom, charges, geometry, step, par_num)
            command[1] = 'mdvp'
            step_old = step
            step = 0
            parallel_numerics.finite_displ_calcs(geometry, command, cobcom, charges, par_num)
            command[1] = 'mdv'
            step = 1
            par_num = 2
            Results = CBF.QM(command, cobcom, charges, geometry, step, par_num)
            QM_Results[1] = Results[0]
            par_num = 3
            step = step_old
        # ... or analytically
        else:
            QM_Results = CBF.QM(command, cobcom, charges, geometry, step, par_num)
        sef = shelve.open("cobram-sef")
        sef['old_step'] = step - 1
        sef.close()
        # in the THS scheme NACs are not available, but the GD is computed and stored in the previous block after hop
        # therefore the check for furstrated hops is performed here and not in tully.py
        if int(command[85]) > 0 and newstate > state and command[14] == '1' and os.path.exists(
                os.getcwd() + '/HOP') and command[206] == '0':
            logwrt.writelog("#############\nHOP to a higher state, check if enough kinetic energy available!\n")
            xvel, yvel, zvel = [], [], []
            velinp = open('velocity.dat')
            vel = velinp.read().split('\n')
            velinp.close()
            for i in range(len(vel) - 1):
                el = vel[i].split()
                xvel.append(float(el[0]))
                yvel.append(float(el[1]))
                zvel.append(float(el[2]))
            xvel = np.array(xvel)
            yvel = np.array(yvel)
            zvel = np.array(zvel)
            sef = shelve.open("cobram-sef", 'r')
            GD = sef['GD']
            xdc = np.array(GD[0])
            ydc = np.array(GD[1])
            zdc = np.array(GD[2])
            sef.close()
            # print 'Gradient difference for scaling velocites after hopping:'
            # for i in range(len(xdc)):
            # print "%12.6f" %xdc[i]+"%12.6f" %ydc[i]+"%12.6f" %zdc[i]
            tstep = float(command[84]) * float('41.341373337')
            condition = product.velscale(geometry, command, xdc, ydc, zdc, xvel, yvel, zvel, DEarray[state][newstate],
                                         tstep, cobcom)
            if condition == 'true':
                logwrt.writelog("HOP accepted!\n#############\n")
            else:
                logwrt.writelog("HOP rejected, not enough kinetic energy to redistribute!\n#############\n")
                newstate = state
                system('rm HOP')
                sef = shelve.open("cobram-sef")
                sef['newstate'] = newstate
                sef['gradient'] = sef['gradient2']
                sef['gradch'] = sef['gradch2']
                sef['dipole'] = sef['dipole2']
                sef['charges'] = sef['charges2']
                QM_Results[1] = sef['gradient2']
                QM_Results[2] = sef['charges2']
                QM_Results[4] = sef['dipole2']
                QM_Results[5] = sef['gradch2']
                sef.close()
                enefile = input('cobram.command')
                tmp = []
                for line in enefile:
                    tmp.append(line.strip())
                for i in range(len(tmp)):
                    if tmp[i].strip().lower()[0:4] == 'rlxr':
                        try:
                            int(tmp[i + 1])
                            tmp[i + 1] = str(newstate + 1)
                        except:
                            tmp[i] = 'rlxr=' + str(newstate + 1)
                    if tmp[i].strip().lower()[0:4] == 'sala':
                        try:
                            int(tmp[i + 1])
                            tmp[i + 1] = str(newstate + 1)
                        except:
                            tmp[i] = 'sala=' + str(newstate + 1)
                filein = open('cobram.command', 'w')
                for i in range(len(tmp)):
                    filein.write(str(tmp[i]) + '\n')
                filein.close()
                cobcom = CBF.getCobrammCommand('cobram.command')
    # if there was not hopping nothing is done
    else:
        QM_Results = 0

    if command[1] == 'mdv' and step == old_step + 1:
        sef = shelve.open("cobram-sef")
        sef['old_step'] = step
        sef.close()

    return QM_Results


def basisset_mapping(basis_label, cartesian):
    ### Molcas to Molpro basis set mapping: one needs to reorder the MOs ###
    ### defined up to 3rd row ###
    ### different maps for each basis set and for each row ###
    ### generally contracted basis sets cannot be used as CADPAC cannot compute DC for them ###
    ### 3rd row not tested yet (extrapolated from 2nd row) ###
    ###STO-3G (the same basis function order in both programs)
    map_STO3G_1st = [1]
    map_STO3G_2nd = [1, 2, 3, 4, 5]
    map_STO3G_3rd = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    ###6-31G (Molcas stores 2x,3x,2y,3y,2z,3z, Molpro stores 2x,2y,2z,3x,3y,3z)
    map_631G_1st = [1, 2]
    map_631G_2nd = [1, 2, 3, 4, 6, 8, 5, 7, 9]
    map_631G_3rd = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 9, 11, 13]
    ###6-31G* (6d)
    map_631Gp6d_1st = [1, 2]
    map_631Gp6d_2nd = [1, 2, 3, 4, 6, 8, 5, 7, 9, 10, 13, 15, 11, 12, 14]
    map_631Gp6d_3rd = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 9, 11, 13, 14, 17, 19, 15, 16, 18]
    ###6-31G* (5d)
    map_631Gp5d_1st = [1, 2]
    map_631Gp5d_2nd = [1, 2, 3, 4, 6, 8, 5, 7, 9, 12, 10, 13, 14, 11]
    map_631Gp5d_3rd = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 9, 11, 13, 16, 14, 17, 18, 15]
    ###6-31G** (6d)
    map_631Gpp6d_1st = [1, 2, 3, 4, 5]
    map_631Gpp6d_2nd = [1, 2, 3, 4, 6, 8, 5, 7, 9, 10, 13, 15, 11, 12, 14]
    map_631Gpp6d_3rd = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 9, 11, 13, 14, 17, 19, 15, 16, 18]
    ###6-31G** (5d)
    map_631Gpp5d_1st = [1, 2, 3, 4, 5]
    map_631Gpp5d_2nd = [1, 2, 3, 4, 6, 8, 5, 7, 9, 12, 10, 13, 14, 11]
    map_631Gpp5d_3rd = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 9, 11, 13, 16, 14, 17, 18, 15]

    # choose basis set and return map
    if basis_label == 'STO-3G':
        map_1st = map_STO3G_1st
        map_2nd = map_STO3G_2nd
        map_3rd = map_STO3G_3rd
    elif basis_label == '6-31G':
        map_1st = map_631G_1st
        map_2nd = map_631G_2nd
        map_3rd = map_631G_3rd
    elif (basis_label == '6-31Gp' or basis_label == '6-31G*') and cartesian == 1:
        map_1st = map_631Gp6d_1st
        map_2nd = map_631Gp6d_2nd
        map_3rd = map_631Gp6d_3rd
    elif (basis_label == '6-31Gp' or basis_label == '6-31G*') and cartesian != 1:
        map_1st = map_631Gp5d_1st
        map_2nd = map_631Gp5d_2nd
        map_3rd = map_631Gp5d_3rd
    elif (basis_label == '6-31Gpp' or basis_label == '6-31G**') and cartesian == 1:
        map_1st = map_631Gpp6d_1st
        map_2nd = map_631Gpp6d_2nd
        map_3rd = map_631Gpp6d_3rd
    elif (basis_label == '6-31Gpp' or basis_label == '6-31G**') and cartesian != 1:
        map_1st = map_631Gpp5d_1st
        map_2nd = map_631Gpp5d_2nd
        map_3rd = map_631Gpp5d_3rd
    else:
        logwrt.fatalerror('Basis set basis set is not supported')

    return [map_1st, map_2nd, map_3rd]


def determine_position(order, elements, basisset_info, cartesian, i, pos):
    ### recursively determine the position 'pos' of the first AO coefficient (i.e. 1s) of atom 'i' in an orbital
    i = i - 1
    # if 'i' is 0, i.e. the first atom in the Molcas/Molpro atom list, return 0 as the starting position, otherwise enter the recursive procedure
    if i == -1:
        return pos
    else:
        j = order[i]
        # find out what is the basis set of the specific atom 'i'
        map_array = basisset_mapping(basisset_info[j][0], cartesian)
        if elements[j] == 'H':
            pos = pos + len(map_array[0])
        elif elements[j] == 'C' or elements[j] == 'N' or elements[j] == 'O' or elements[j] == 'F':
            pos = pos + len(map_array[1])
        elif elements[j] == 'Si' or elements[j] == 'P' or elements[j] == 'S' or elements[j] == 'Cl':
            pos = pos + len(map_array[2])
        pos = determine_position(order, elements, basisset_info, cartesian, i, pos)
        return pos


# a routine for the flexible reading of basis set definitions: !basisset 1,3-5 6-31G* 3s2p1d\n 2,6-7 6-31G 3s2p\n ?basisset
# previous basis set definition via keyword 197 still holds but is overwritten by !basisset ?basisset
def readbassisset(geometry, command, cobcom):
    natom = geometry.NatomH + geometry.NsubH
    basissetkeys = CBF.ReadCobramCommand(cobcom, 'basisset', 'basisset.input')
    if not basissetkeys:
        basisset_info = [[0 for cols in range(2)] for rows in range(natom)]  # a 2D list of length natom with
        # first column basis set label and second column the basis set contraction; basis set contraction is
        # useful for generally contracted basis sets used by Molcas, it is not used with Molpro
        for i in range(natom):
            basisset_info[i][0] = str(command[197])
    else:
        basisset_info = [[0 for cols in range(2)] for rows in range(natom)]
        for i in range(len(basissetkeys)):
            tmp_elements = []
            try:
                tmp_elements = basissetkeys[i].split()[0].strip().split(',')
                for k in range(len(tmp_elements)):
                    tmp2_elements = []
                    tmp2_elements = tmp_elements[k].split('-')
                    if len(tmp2_elements) == 1:
                        basisset_info[int(tmp2_elements[0]) - 1][0] = basissetkeys[i].split()[1].strip()
                        if basissetkeys[i].split()[1].strip() not in ['6-31G', '6-31G*', '6-31Gp', '6-31G**', '6-31Gpp',
                                                                      'STO-3G']:
                            basisset_info[int(tmp2_elements[0]) - 1][1] = basissetkeys[i].split()[2].strip()
                        else:
                            basisset_info[int(tmp2_elements[0]) - 1][1] = 'dummy'
                    else:
                        for l in range(int(tmp2_elements[0]) - 1, int(tmp2_elements[1])):
                            basisset_info[l][0] = basissetkeys[i].split()[1].strip()
                            if basissetkeys[i].split()[1].strip() not in ['6-31G', '6-31G*', '6-31Gp', '6-31G**',
                                                                          '6-31Gpp', 'STO-3G']:
                                basisset_info[l][1] = basissetkeys[i].split()[2].strip()
                            else:
                                basisset_info[l][1] = 'dummy'
            except:
                pass
        for i in range(natom):
            if not basisset_info[i][0] or not basisset_info[i][1]:
                logwrt.fatalerror("Basis set definition for atom " + str(i) + " is missing")
    logwrt.writelog("basisset_info: " + str(basisset_info)+"\n")
    return basisset_info


def save_QM_molcas_step(step, command, copyLog, copyOrb):
    """ save mocals output and orbital files for last single point QM calculation
        in a common directory where QM data is stored """

    # name of the molcas calculation output
    logName = 'molcas.log'
    # file where to store QM results
    allName = os.path.join(QM_DATA_STORAGE, 'molcasALL.log')
    savname = allName = os.path.join(MOLCAS_WRKDIR, 'molcas.log')

    if copyLog:
        # check if the file exists, otherwise print a warning to screen
        if not os.path.exists(logName):
            logwrt.writewarning(
                "molcas file {0} cannot be found: the {1} file will not be updated".format(logName, allName))

        else:
            # write the content of the QM single point calc in the ALL file, decorated with
            # comments that highligh the step number of the QM single point
            with open(allName, 'a') as molcastot:

                molcastot.write('=' * 80 + '\n')
                molcastot.write("Start SP Molcas calculation of STEP : " + str(step) + "\n")
                molcastot.write(' || ' * 20 + '\n')
                molcastot.write(' \/ ' * 20 + '\n')
                with open(logName, "r") as molcasOut:
                    molcastot.write(molcasOut.read())
                molcastot.write(' /\ ' * 20 + '\n')
                molcastot.write(' || ' * 20 + '\n')
                molcastot.write("End SP Molcas calculation of STEP : " + str(step) + "\n")
                molcastot.write('=' * 80 + '\n')
                molcastot.write('\n')

    # save orbital file and molden file once you know what type of calculation you have
    SCF = 'no'
    MBPT2 = 'no'
    CAS = 'no'
    SSPT2 = 'no'
    MSPT2 = 'no'
    if os.path.exists(os.getcwd() + '/SCF'):
        SCF = 'yes'
    if os.path.exists(os.getcwd() + '/RASSCF'):
        CAS = 'yes'
        SCF = 'no'
    elif os.path.exists(os.getcwd() + '/MBPT2'):
        SCF = 'no'
        MBPT2 = 'yes'
    elif os.path.exists(os.getcwd() + '/SS-CASPT2'):
        SSPT2 = 'yes'
    elif os.path.exists(os.getcwd() + '/MS-CASPT2'):
        MSPT2 = 'yes'

    # copy orbitals and molden file when required
    if copyOrb:
        if CAS == 'yes' or SSPT2 == 'yes' or MSPT2 == 'yes':
            srcOrbFile = "molcas.RasOrb"
            targetOrbFile = os.path.join(QM_DATA_STORAGE, 'molcas_' + str(step) + '.RasOrb')
            srcMoldenFile = "molcas.rasscf.molden"
            targetMoldenFile = os.path.join(QM_DATA_STORAGE, 'step_' + str(step) + '_rasscf.molden')
        elif SCF == 'yes':
            srcOrbFile = "molcas.ScfOrb"
            targetOrbFile = os.path.join(QM_DATA_STORAGE, 'molcas_' + str(step) + '.ScfOrb')
            srcMoldenFile = "molcas.scf.molden"
            targetMoldenFile = os.path.join(QM_DATA_STORAGE, 'step_' + str(step) + '_scf.molden')
        elif MBPT2 == 'yes':
            srcOrbFile = "molcas.ScfOrb"
            targetOrbFile = os.path.join(QM_DATA_STORAGE, 'molcas_' + str(step) + '.ScfOrb')
            srcMoldenFile = "molcas.scf.molden"
            targetMoldenFile = os.path.join(QM_DATA_STORAGE, 'step_' + str(step) + '_scf.molden')

        # check if the orbital file exists and save it and zip it, otherwise print a warning to screen
        if not os.path.exists(os.path.join("MOLCAS", srcOrbFile)):
            logwrt.writewarning("molcas file {0} cannot be found: it will not be stored".format(srcOrbFile))
        else:
            shutil.copyfile(os.path.join("MOLCAS", srcOrbFile), targetOrbFile)
            CBF.GZIP(targetOrbFile)
        # check if the molden file exists and save it and zip it, otherwise print a warning to screen
        if not os.path.exists(os.path.join("MOLCAS", srcMoldenFile)):
            logwrt.writewarning("molcas file {0} cannot be found: it will not be stored".format(srcMoldenFile))
        else:
            shutil.copyfile(os.path.join("MOLCAS", srcMoldenFile), targetMoldenFile)
            CBF.GZIP(targetMoldenFile)


def clean_QM_molcas():
    """ clean up the run directory from all the files that have been used to run MOLCAS and that
        are no longer needed at the end of the calculation """

    # list of files/directory to remove
    toRemove = ["SCF", "RASSCF", "MBPT2", "SS-CASPT2", "MS-CASPT2", "STATE", "MOLCAS",
                "molcas.log", "molcas.status"]

    for f in toRemove:
        # if f is an existing file, remove it
        if os.path.isfile(f):
            os.remove(f)
        # if f is an existing directory, clean tree
        elif os.path.isdir(f):
            shutil.rmtree(f)

    # done
    return


def computeGD(grad_state, grad_newstate, gradch_state, gradch_newstate, geometry, command):
    sef = shelve.open("cobram-sef")
    Fxyz_modelH = sef['Fxyz_modelH']
    sef.close()
    GD = []
    fx_state, fy_state, fz_state = [], [], []
    fx_newstate, fy_newstate, fz_newstate = [], [], []
    GD.extend([[], [], []])
    # add MM part to QM gradient
    fx, fy, fz = [], [], []
    if geometry.NatomMM > 0:
        for j in range(len(Fxyz_modelH[1])):
            fx.append(float(-grad_state[0][j] + Fxyz_modelH[0][j]))
            fy.append(float(-grad_state[1][j] + Fxyz_modelH[1][j]))
            fz.append(float(-grad_state[2][j] + Fxyz_modelH[2][j]))
    else:
        fx = copy.deepcopy(-grad_state[0])
        fy = copy.deepcopy(-grad_state[1])
        fz = copy.deepcopy(-grad_state[2])
    correctgrad_state = [np.array(fx), np.array(fy), np.array(fz)]
    fx, fy, fz = [], [], []
    if geometry.NatomMM > 0:
        for j in range(len(Fxyz_modelH[1])):
            fx.append(float(-grad_newstate[0][j] + Fxyz_modelH[0][j]))
            fy.append(float(-grad_newstate[1][j] + Fxyz_modelH[1][j]))
            fz.append(float(-grad_newstate[2][j] + Fxyz_modelH[2][j]))
    else:
        fx = copy.deepcopy(-grad_newstate[0])
        fy = copy.deepcopy(-grad_newstate[0])
        fz = copy.deepcopy(-grad_newstate[0])
    correctgrad_newstate = [np.array(fx), np.array(fy), np.array(fz)]
    if geometry.NsubH > 0:
        # redistribute H-link atom gradient over QM and MM atoms
        coulomb_state = [0.0, np.array(gradch_state[0]), np.array(gradch_state[1]), np.array(gradch_state[2])]
        coulomb_newstate = [0.0, np.array(gradch_newstate[0]), np.array(gradch_newstate[1]), np.array(gradch_newstate[2])]
        high = list(geometry.list_HIGH)
        medium = list(geometry.list_MEDIUM)
        B = geometry.atomLink_BA[0]
        A = geometry.atomLink_BA[1]
        NatomQM = geometry.NatomQM
        hl_count = 0
        ml_count = 0
        for i in range(geometry.atomNum):
            if (high.count(i + 1) == 1 and A.count(i + 1) == 0) or (
                    high.count(i + 1) == 1 and command[61] == '0' and A.count(i + 1) > 0):
                fx_state.append(correctgrad_state[0][hl_count])
                fy_state.append(correctgrad_state[1][hl_count])
                fz_state.append(correctgrad_state[2][hl_count])
                fx_newstate.append(correctgrad_newstate[0][hl_count])
                fy_newstate.append(correctgrad_newstate[1][hl_count])
                fz_newstate.append(correctgrad_newstate[2][hl_count])
                hl_count = hl_count + 1
            if high.count(i + 1) == 1 and command[61] == '1' and A.count(i + 1) > 0:
                # redistribute gradient of H-link over HL atom for first gradient
                aa = A.index(i + 1)
                atype = geometry.atomLabel[A[aa] - 1]
                if atype == 'C':
                    g_jac = float(command[21]) / float(command[24])
                if atype == 'O':
                    g_jac = float(command[22]) / float(command[25])
                if atype == 'N':
                    g_jac = float(command[23]) / float(command[26])
                la = aa + NatomQM
                QMFx, QMFy, QMFz = [], [], []
                QMFx = correctgrad_state[0][hl_count] + correctgrad_state[0][la] * (1 - g_jac)
                QMFy = correctgrad_state[1][hl_count] + correctgrad_state[1][la] * (1 - g_jac)
                QMFz = correctgrad_state[2][hl_count] + correctgrad_state[2][la] * (1 - g_jac)
                if (A.count(i + 1) == 2):
                    aa = aa + 1
                    la = aa + NatomQM
                    QMFx = QMFx + correctgrad_state[0][la] * (1 - g_jac)
                    QMFy = QMFy + correctgrad_state[1][la] * (1 - g_jac)
                    QMFz = QMFz + correctgrad_state[2][la] * (1 - g_jac)
                if (A.count(i + 1) == 3):
                    aa = aa + 1
                    la = aa + NatomQM
                    QMFx = QMFx + correctgrad_state[0][la] * (1 - g_jac)
                    QMFy = QMFy + correctgrad_state[1][la] * (1 - g_jac)
                    QMFz = QMFz + correctgrad_state[2][la] * (1 - g_jac)
                fx_state.append(float(QMFx))
                fy_state.append(float(QMFy))
                fz_state.append(float(QMFz))
                # redistribute gradient of H-link over HL atom for second gradient
                aa = A.index(i + 1)
                la = aa + NatomQM
                QMFx, QMFy, QMFz = [], [], []
                QMFx = correctgrad_newstate[0][hl_count] + correctgrad_newstate[0][la] * (1 - g_jac)
                QMFy = correctgrad_newstate[1][hl_count] + correctgrad_newstate[1][la] * (1 - g_jac)
                QMFz = correctgrad_newstate[2][hl_count] + correctgrad_newstate[2][la] * (1 - g_jac)
                if (A.count(i + 1) == 2):
                    a = a + 1
                    la = a + NatomQM
                    QMFx = QMFx + correctgrad_newstate[0][la] * (1 - g_jac)
                    QMFy = QMFy + correctgrad_newstate[1][la] * (1 - g_jac)
                    QMFz = QMFz + correctgrad_newstate[2][la] * (1 - g_jac)
                if (A.count(i + 1) == 3):
                    a = a + 1
                    la = a + NatomQM
                    QMFx = QMFx + correctgrad_newstate[0][la] * (1 - g_jac)
                    QMFy = QMFy + correctgrad_newstate[1][la] * (1 - g_jac)
                    QMFz = QMFz + correctgrad_newstate[2][la] * (1 - g_jac)
                fx_newstate.append(float(QMFx))
                fy_newstate.append(float(QMFy))
                fz_newstate.append(float(QMFz))
                hl_count = hl_count + 1
            if (medium.count(i + 1) == 1 and B.count(i + 1) == 0) or (
                    medium.count(i + 1) == 1 and command[61] == '0' and B.count(i + 1) > 0):
                fx_state.append(coulomb_state[1][ml_count])
                fy_state.append(coulomb_state[2][ml_count])
                fz_state.append(coulomb_state[3][ml_count])
                fx_newstate.append(coulomb_newstate[1][ml_count])
                fy_newstate.append(coulomb_newstate[2][ml_count])
                fz_newstate.append(coulomb_newstate[3][ml_count])
                ml_count = ml_count + 1
            if medium.count(i + 1) == 1 and command[61] == '1' and B.count(i + 1) == 1:
                # redistribute gradient of H-link over ML atom for first gradient
                bb = B.index(i + 1)
                atype = geometry.atomLabel[A[bb] - 1]
                if atype == 'C':
                    g_jac = float(command[21]) / float(command[24])
                if atype == 'O':
                    g_jac = float(command[22]) / float(command[25])
                if atype == 'N':
                    g_jac = float(command[23]) / float(command[26])
                la = bb + NatomQM
                MMFx, MMFy, MMFz = [], [], []
                MMFx = coulomb_state[1][ml_count] + correctgrad_state[0][la] * g_jac
                MMFy = coulomb_state[2][ml_count] + correctgrad_state[1][la] * g_jac
                MMFz = coulomb_state[3][ml_count] + correctgrad_state[2][la] * g_jac
                fx_state.append(float(MMFx))
                fy_state.append(float(MMFy))
                fz_state.append(float(MMFz))
                # redistribute gradient of H-link over ML atom for second gradient
                MMFx, MMFy, MMFz = [], [], []
                MMFx = coulomb_newstate[1][ml_count] + correctgrad_newstate[0][la] * g_jac
                MMFy = coulomb_newstate[2][ml_count] + correctgrad_newstate[1][la] * g_jac
                MMFz = coulomb_newstate[3][ml_count] + correctgrad_newstate[2][la] * g_jac
                fx_newstate.append(float(MMFx))
                fy_newstate.append(float(MMFy))
                fz_newstate.append(float(MMFz))
                ml_count = ml_count + 1
    else:
        fx_state = correctgrad_state[0]
        fy_state = correctgrad_state[1]
        fz_state = correctgrad_state[2]
        fx_newstate = correctgrad_newstate[0]
        fy_newstate = correctgrad_newstate[1]
        fz_newstate = correctgrad_newstate[2]

    for i in range(len(fx_newstate)):
        GD[0].append(fx_state[i] - fx_newstate[i])
        GD[1].append(fy_state[i] - fy_newstate[i])
        GD[2].append(fz_state[i] - fz_newstate[i])

    sef = shelve.open("cobram-sef")
    sef['GD'] = GD
    sef.close()

