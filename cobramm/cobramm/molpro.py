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

from math import *  # mathematical functions
import os  # operating system utilities
import sys  # system-specific parameters and functions
import shutil  # high-level file operations
import shelve  # python object persistence
import copy  # shallow and deep copy operations

# imports of local modules

import CBF
import logwrt  # manages log file output + start/end procedures
import cobrammenv  # environmental variable for COBRAMM and 3rd-party software
import molcas
from tullyNEW import Tully

# math libraries

import numpy as np  # numpy library for scientific computation

#####################################################################################################

QM_DATA_STORAGE = 'QM_data'
MOLDEN_FILE = "cobramm.molden"
SINGLEPOINT_DIR = "SINGLEPOINT"

resortat = []

#####################################################################################################


def prepare(command, step):
    """ prepare the directory for molpro calculation """

    # check that the environment to run molcas is properly defined
    envDefined, errorMsg = cobrammenv.checkMolproEnv()
    if not envDefined: logwrt.fatalerror(errorMsg)

    # stuff to do at the beginning of a molpro calculation
    if step == 0:

        # print message to log about Molpro executable that is used in the current environment
        molproexe = os.getenv('MOLPRO_EXE')  # path of the Molpro executable
        logwrt.writelog('using Molpro executable {}\n'.format(molproexe))

        # create temporary directories to store Molpro calculation files
        if os.path.exists(QM_DATA_STORAGE): shutil.rmtree(QM_DATA_STORAGE)
        os.mkdir(QM_DATA_STORAGE)

        # if a user supply the INP.wfu input file, it will be used as initial wavefunction for the calculation
        if os.path.exists('INP.wfu'):
            shutil.copy("INP.wfu", "work.wfu")

        # create directory to store SP calculations done on top of trajectories/paths
        if command[200] != '0':
            if os.path.exists(SINGLEPOINT_DIR): shutil.rmtree(SINGLEPOINT_DIR)
            os.mkdir(SINGLEPOINT_DIR)

    gen, keys = [], []
    return [keys, gen]


def makeinp(command, geometry, CRG_emb, step):
    """create the input file of molpro which can be wariable for each step"""

    covcost = float('0.52917720830008323378')
    guess = command[6]
    Nproc = command[7]
    zero = float('0.0')
    ciro_key = 0
    n = 1

    sef = shelve.open("cobram-sef", 'r')
    NAC = sef['NAC']
    if command[51] == '6':
        aux_data = sef['aux_data']
    NAC_old = sef['NAC_old']
    DEarray = sef['DEarray']
    nroots = sef['nroots']
    sef.close()

    if step == 0:
        # convention for state labeling in cobram:  0=GS, 1=1st ES, 2=2nd ES"
        # Molcas convention is different: 1=GS, 2=1st ES, 3=2nd ES"
        cobcom = CBF.getCobrammCommand('cobram.command')
        keys = CBF.ReadCobramCommand(cobcom, 'molpro', 'molpro.input')
        for i in range(len(keys)):
            ele = keys[i].strip().split(';')
            for j in range(len(ele)):
                try:
                    if ele[j].lower()[0:5] == 'state':
                        nroots = int(ele[j].split(",")[1])
                        ciro_key = 1
                        break
                except:
                    pass
            if ciro_key == 1:
                break

        # generate empty array NAC to take the couplings
        if ciro_key == 1:
            if int(command[85]) > 3 or command[1] == 'ci':
                logwrt.writelog("Initializing NAC matrix for {0} states\n".format(nroots))
            NAC = []
            for i in range(nroots):
                NAC.append([])
            for i in range(nroots):
                for j in range(nroots):
                    NAC[i].append([])
            for i in range(nroots):
                for j in range(nroots):
                    NAC[i][j].extend([[], [], []])
            NAC_old = copy.deepcopy(NAC)
            # We need an empty array to re-initialise the NAC array
            # every time if a value is given in command 86!
            NAC_empty = copy.deepcopy(NAC)

            sef = shelve.open("cobram-sef")
            sef['NAC_old'] = NAC_old
            sef['NAC_empty'] = NAC_empty
            sef.close()

            # logwrt.writelog( "Initializing transition energy matrix for {0} states\n".format(nroots) )
            for i in range(nroots):
                DEarray.append([])
            for i in range(nroots):
                for j in range(nroots):
                    DEarray[i].append(1000.0)

            # Generate an empty nroots x nroots matrix for usage in tully.py
            # maybe we should shift all this to cobram.py
            zeromat = np.zeros((nroots, nroots))
            sef = shelve.open("cobram-sef")
            sef['zeromat'] = zeromat
            sef.close()

        else:
            DEarray.append([])
            DEarray[0].append(1000.0)
            DEarray[0].append(1000.0)  # element [0][1]

        # determine the state of interest
        if os.path.exists('STATE'):
            with open("STATE") as f:
                state = int(f.read())
        else:
            state = 0

    elif step != 0 and command[51] == '6':
        sef = shelve.open("cobram-sef")
        state = sef['state']
        sef.close()
    else:
        state = 0
    if state != 0:
        DE = DEarray[0][1]
    else:
        DE = 1000.0
    if command[199] != '0' and DE < -float(command[199]) and step != 0 and state == 0:
        cobcom = CBF.getCobrammCommand('cobram.command')
        keys = CBF.ReadCobramCommand(cobcom, 'molproSS', 'molproinSS.input')
    else:
        cobcom = CBF.getCobrammCommand('cobram.command')
        keys = CBF.ReadCobramCommand(cobcom, 'molpro', 'molproin.input')

    molproin = []
    lattice = []
    logwrt.writelog('\nConstructing the Molpro input file...\n')

    # get elements,coordinate in Angstrom
    elements = [geometry.atomLabel[i - 1] for i in geometry.list_QM] + ["H"] * geometry.NsubH
    X = geometry.modelH[0]
    Y = geometry.modelH[1]
    Z = geometry.modelH[2]

    global resortat
    resortat = []
    # resort atoms for molpro, atoms of the same elements in one block:
    # modify if you need other atom type
    for i in range(len(elements)):
        if elements[i] == 'O':
            resortat.append(i)
    for i in range(len(elements)):
        if elements[i] == 'N':
            resortat.append(i)
    for i in range(len(elements)):
        if elements[i] == 'C':
            resortat.append(i)
    for i in range(len(elements)):
        if elements[i] == 'H':
            resortat.append(i)
    for i in range(len(elements)):
        if elements[i] == 'F':
            resortat.append(i)
    for i in range(len(elements)):
        if elements[i] == 'P':
            resortat.append(i)
    for i in range(len(elements)):
        if elements[i] == 'S':
            resortat.append(i)
    for i in range(len(elements)):
        if elements[i] == 'S' or elements[i] == 'P' or elements[i] == 'F' or elements[i] == 'H' or \
                elements[i] == 'C' or elements[i] == 'N' or elements[i] == 'O':
            pass
        else:
            logwrt.fatalerror("NO ATOM TYPE DEFINED IN MOLPRO RESORT ROUTINE: " + elements[i] + " not found")
    for i in range(len(elements)):
        b = resortat[i]

    # construct molpro input
    # check if basis is not already contained in the molpro section of cobram.command
    do_basis = 'yes'
    for i in range(len(keys)):
        if keys[i].find('basis') != -1 or keys[i].find('BASIS') != -1:
            do_basis = 'no'
    for i in range(0, len(elements)):
        n = n + 1
        j = resortat[i]  # use new sorting for input
        # add an integer to the element label for atom specific basis set definition
        # but only when basissetkeys is defines, otherwise use without numbering for compatibility
        if do_basis == 'yes':
            molproin.append(' ' + elements[j] + str(i + 1) + ',   ')
        else:
            molproin.append(' ' + elements[j] + ',   ')

        molproin.append('%14.8f' % X[j] + ',')
        molproin.append('%14.8f' % Y[j] + ',')
        molproin.append('%14.8f' % Z[j] + '\n')
    molproin.append('}\n')

    # write external charges
    if len(CRG_emb) != 0 and geometry.calculationType not in ['H', 'M', 'ML']:
        # logwrt.writelog('Appending the atomic point charges to molpro-lattice file\n')
        molproin.append('LATTICE,infile=lattice,NUCONLY;')
        lattice.append('Point charges\n' + str(len(CRG_emb)) + '\n')
        for i in range(len(CRG_emb)):
            if geometry.list_MM[i] in geometry.list_MEDIUM:
                lattice.append('%14.8f,' % (geometry.pod[0][i]) + '%14.8f,' % (geometry.pod[1][i]) + '%14.8f,' % (
                    geometry.pod[2][i]) + '%12.8f,' % CRG_emb[i] + '%5i' % 1 + '\n')
            else:
                lattice.append('%14.8f,' % (geometry.pod[0][i]) + '%14.8f,' % (geometry.pod[1][i]) + '%14.8f,' % (
                    geometry.pod[2][i]) + '%12.8f,' % CRG_emb[i] + '%5i' % 0 + '\n')

    # write input molpro file
    filein = open('molpro.input', 'w')
    filein.write('#step number')  # Olli debug
    filein.write(str(step))
    filein.write('\n')
    # Molcas reads MB, Molpro MW
    if command[53][-2:].lower() == 'mw':
        filein.write("Memory," + command[53][0:-1] + ";\n")
    elif command[53][-2:].lower() == 'gw':
        filein.write("Memory," + str(int(command[53][0:-2]) * 1000) + ",M;\n")
    elif command[53][-2:].lower() == 'mb':
        filein.write("Memory," + str(int(command[53][0:-2]) / 8) + ",M;\n")
    elif command[53][-2:].lower() == 'gb':
        filein.write("Memory," + str(int(command[53][0:-2]) * 1000 / 8) + ",M;\n")
    # assuming memory is in MW
    elif command[53][-1].lower() == 'm':
        filein.write("Memory," + command[53] + ";\n")
    # assuming memory is in GW
    elif command[53][-1].lower() == 'g':
        filein.write("Memory," + str(int(command[53][0:-1]) * 1000) + ",M;\n")
    else:
        # if non of the above try to convert to an integer
        try:
            filein.write("Memory," + str(int(command[53])) + ",M;\n")
        except:
            logwrt.writelog("Memory card 53 could not be recognized. Using a default memory of 2048MB\n")
            filein.write("Memory,256,M;\n")
    filein.write('Gprint,basis;\n')
    filein.write('Gprint,civector;\n')
    # this affects RS2 computations and is required in order to prevent Molpro crashing
    filein.write('GTHRESH,ENERGY=1.d-10,GRADIENT=1.d-10;\n')
    filein.write('Gdirect;\n')
    # CI coeffs are not required during Molcas MD
    if command[51] == '7' and command[1] in ['mdt', 'mdo', 'mdv'] and command[85] == '1':
        filein.write('Gthresh,printci=0.0;\n')
    # filein.write('Gprint,ref;\n')
    filein.write('Gprint,orbital;\n')
    filein.write('file,2,WORK.wfu;\n')
    # filein.write('geomtyp=xyz\n')               # OW Adapted input to molpro2009.2
    # explicit specification of the d-functions
    if command[194] == '1':
        filein.write('CARTESIAN\n')
    else:
        filein.write('SPHERICAL\n')
    filein.write('ANGSTROM\n')
    filein.write('NOSYM\n')
    filein.write('geometry={\n')
    filein.write(str(len(elements)) + '\n')
    filein.write('GeomXYZ\n')
    for i in range(len(molproin)):
        filein.write(str(molproin[i]))

    # flexible basis set definition
    if do_basis == 'yes':
        # the varible basis has a max length 1024: we cannot specify the atom basis sets 1-by-1 like in Molcas
        basisset_info = molcas.readbassisset(geometry, command, cobcom)
        filein.write('\nbasis\n')
        # count basis set definitions, the basis set with the most definitions becomes default
        basisset_list = []
        basisset_count = [0 for rows in
                          range(len(elements))]  # maximal possible length (if every atom has a different basis set)
        for i in range(len(elements)):
            if basisset_info[i][0] == '6-31Gpp':
                basisset_info[i][0] = '6-31G**'
            elif basisset_info[i][0] == '6-31Gp':
                basisset_info[i][0] = '6-31G*'
            if basisset_info[i][0] in basisset_list:
                j = basisset_list.index(str(basisset_info[i][0]))
                basisset_count[j] = basisset_count[j] + 1
            else:
                basisset_list.append(str(basisset_info[i][0]))
                j = basisset_list.index(str(basisset_info[i][0]))
                basisset_count[j] = basisset_count[j] + 1
        default_bs = basisset_list[basisset_count.index(max(basisset_count))]
        logwrt.writelog('Default basis set: {0}\n'.format(default_bs))
        filein.write('default=' + default_bs)
        # count basis set definitions per atom type, the basis set with the most definitions
        # (if different from default) becomes atom type default
        # for all other atom from the same atom type define basis sets 1-by-1
        # does not work for ANO-L with contractions as they require special format
        unique_elements = []
        appearances = elements.count(elements[0])
        # unique_elements is a 1D list with even elements the atom type and odd elements a list of
        # the positions of the atom type: ['N',[4],'O',[5,6],'C',[1,2,3]]
        unique_elements.append(elements[0])
        unique_elements.append([0 for rows in range(appearances)])
        for i in range(len(elements)):
            if elements[i] not in unique_elements:
                appearances = elements.count(elements[i])
                unique_elements.append(elements[i])
                unique_elements.append([0 for rows in range(appearances)])
        for i in range(len(unique_elements) // 2):
            index = 0
            for j in range(len(elements)):
                if unique_elements[2 * i] == elements[j]:
                    unique_elements[2 * i + 1][index] = j
                    index = index + 1

        for i in range(len(unique_elements) // 2):
            basisset_list = []
            basisset_count = [0 for rows in range(len(unique_elements[2 * i + 1]))]
            for j in range(len(unique_elements[2 * i + 1])):
                if basisset_info[unique_elements[2 * i + 1][j]][0] in basisset_list:
                    k = basisset_list.index(str(basisset_info[unique_elements[2 * i + 1][j]][0]))
                    basisset_count[k] = basisset_count[k] + 1
                else:
                    basisset_list.append(str(basisset_info[unique_elements[2 * i + 1][j]][0]))
                    k = basisset_list.index(str(basisset_info[unique_elements[2 * i + 1][j]][0]))
                    basisset_count[k] = basisset_count[k] + 1
            default_ibs = basisset_list[basisset_count.index(max(basisset_count))]
            if default_ibs != default_bs:
                filein.write(',' + unique_elements[2 * i] + '=' + default_ibs)
            for j in range(len(unique_elements[2 * i + 1])):
                if basisset_info[unique_elements[2 * i + 1][j]][0] != default_ibs:
                    # for atom specific basis set definitions 'resortat' must be used as the atom order
                    # in the Molpro file is different
                    filein.write(',' + str(unique_elements[2 * i]) + str(
                        resortat.index(unique_elements[2 * i + 1][j]) + 1) + '=' + str(
                        basisset_info[unique_elements[2 * i + 1][j]][0]))
        filein.write('\nend\n')

    # add the matrop section for reading the Molcas MOs
    if command[51] == '6':
        filein.write('{matrop;\n')
        filein.write('read,molcasorb,type=orb,subtype=natural,file=molcasorb.dat;\n')
        filein.write('print,molcasorb;\n')
        filein.write('save,molcasorb,2140.2,orbitals;\n')
        filein.write('}\n')

    # automatically create the molpro input using Molcas data (only when no user input in cobram.command)
    if command[51] == '6':
        closedorbs = aux_data[9]
        closedels = aux_data[10]
        occorbs = aux_data[11]
        totalels = aux_data[12]
        charge = int(aux_data[13])
        numstates = aux_data[14]
        spin = int(aux_data[15])
        symmetry = aux_data[16]
        filein.write('{multi;\nstart,2140.2,natural;\ndont,orbital;\nocc,' + str(occorbs) + ';\nclosed,' + str(
            closedorbs) + ';\nwf,' + str(totalels) + ',' + str(symmetry) + ',' + str(spin) + ';\nstate,' + str(
            numstates) + ';\n')
        rec = 1

        sef = shelve.open("cobram-sef")
        calc_coupl = sef['calc_coupl']
        newstate = sef['newstate']
        threshold = float(command[86])
        sef.close()

        tmp_coupl = chk_coupl(newstate, threshold, DEarray, calc_coupl)
        calc_coupl = list(tmp_coupl)
        for z in range(len(calc_coupl)):
            j = calc_coupl[z][0]
            k = calc_coupl[z][1]
            filein.write('CPMCSCF,NACM,' + str(j + 1) + '.1,' + str(k + 1) + '.1,record=510' + str(rec) + '.1;\n')
            rec = rec + 1
        filein.write('}\n')
    filein.close()

    # state=str('1')
    # print "In the current implementation a gradient is required even in case of a SP."
    # print "If you have not specified it manually Cobram will try to include the missing commands in the molpro.input."
    grad_check = 1
    dm_check = 0
    pop_check = 0
    molden_check = 0
    multi_check = 0
    fileout = open('molpro.input', 'a')
    # check in which state we start from CPMCSCF,GRAD (OW)
    # set 5100.1 as the default record for gradients
    # decide whether to add a {force} section
    for i in range(len(keys)):
        ele = keys[i].strip().split(',')
        try:
            if ele[0].lower()[0:6] == '{multi' or ele[0].lower()[0:5] == 'multi':
                multi_check = 1
                for ii in range(i + 1, len(keys)):
                    ele1 = keys[ii].strip().split(',')
                    if ele1[0].lower() in ['dm;}', 'dm}', 'dm;', 'dm']:
                        dm_check = 1
                    try:
                        if ele1[-1][-1] == '}' and dm_check == 0:
                            # logwrt.writelog("Adding keyword DM to input\n")
                            dm_check = 1
                            ele1[-1] = ele1[-1][:-1]
                            try:
                                if ele1[-1][-1] == ';':
                                    ele1[-1] = ele1[-1] + '\nDM}'
                                else:
                                    ele1[-1] = ele1[-1] + ';\nDM}'
                            except:
                                ele1[-1] = 'DM}'
                            keys[ii] = ','.join(ele1)
                    except:
                        pass
                    try:
                        if ele1[1].lower() == 'grad':
                            state = int(ele1[2].split('.')[0]) - 1
                            logwrt.writelog("State taken from the MULTI block in the input {0}\n".format(state + 1))
                            # check will be set to 1 if the user has included a CPMCSCF,GRAD section in multi,
                            # this would tell Cobram whether a force section is necessary
                            # grad_check=1
                            # set the record for GRAD to 5100.1 independently from user's input
                            for j in range(len(ele1)):
                                if ele1[j][0:6] == 'record':
                                    if ele1[j][7:11] != '5100' and command[1] != 'ci':
                                        tmp = ele1[j][7:11]
                                        logwrt.writewarning(
                                            'modifying user input for GRAD record from {0}.1 to 5100.1'.format(
                                                ele1[j][7:11]))
                                        ele1[j] = ele1[j].replace(tmp, '5100')
                                        keys[ii] = ','.join(ele1)
                        if command[1] == 'ci' and ele1[1].lower() == 'nacm':
                            state = int(ele1[3].split('.')[0]) - 1
                    except:
                        pass
            if ele[0].lower()[0:4] == '{rs2' or ele[0].lower()[0:3] == 'rs2':
                for j in range(len(ele)):
                    if ele[j].lower()[0:4] == 'root':
                        state = int(ele[j].strip(';').split('=')[1]) - 1
                        logwrt.writelog("state taken from the RS2 block in the input {0}\n".format(state + 1))
                        pt2_check = 1
                        multi_check = 0
            if (ele[0].lower()[0:6] == '{force' or ele[0].lower()[0:5] == 'force') and grad_check == 1:
                # if the user has defined the {force} section by himself Cobram sets check to 0;
                # additionally the record is synchronized to the record set in multi
                grad_check = 0
                if multi_check == 1:
                    ele1 = keys[i + 1].strip().split(',')
                    if ele1[1].lower()[0:4] == 'samc' and ele1[1][0:4] != '5100' and command[1] != 'ci':
                        logwrt.writewarning(
                            'modifying user input for FORCE record from {0}.1 to 5100.1'.format(ele1[1][0:4]))
                        tmp = ele1[1][0:4]
                        # set the record for GRAD to 5100.1 independently from user's input
                        ele1[1] = ele1[1].replace(tmp, '5100')
                        keys[i + 1] = ','.join(ele1)
            if ele[0].lower()[0:4] == '{pop' or ele[0].lower()[0:3] == 'pop':
                pop_check = 1
            if ele[0].lower()[0:4] == '{put' or ele[0].lower()[0:3] == 'put':
                molden_check = 1
        except:
            pass
        fileout.write(keys[i] + '\n')

    # if check was set to 1 above a {force} section is added (hard-coded record 5100)
    if command[51] == '6':
        grad_check = 0
        pop_check = 1
        molden_check = 1
    if grad_check == 1 and multi_check == 1 and int(command[60]) > 1:
        # logwrt.writelog( "Adding a FORCE section for state {0} to input\n".format(state+1) )
        fileout.write('\n{force;\nSAMC,5100.1;\n}\n')
    elif grad_check == 1 and int(command[60]) > 1:
        # logwrt.writelog( "Adding a FORCE section\n" )Adding
        fileout.write('\n{force;\n}\n')
    if command[51] == '6':
        for i in range(1, rec):
            fileout.write('\n{force;\nSAMC,510' + str(i) + '.1;\n}')
        fileout.write('\n')
    # add a POP section
    if pop_check == 0 and int(command[60]) > 1:
        # logwrt.writelog( "Adding POP section to input for each state\n" )
        for i in range(nroots):
            fileout.write('{POP;\nDENSITY,STATE=' + str(i + 1) + '.1;\n}\n')
    # adding a MOLDEN section
    if molden_check == 0:
        # logwrt.writelog( "Adding MOLDEN section to input\n" )
        fileout.write('put,molden,{0}\n'.format(MOLDEN_FILE))

    if int(command[60]) > 1 and len(CRG_emb) != 0 and geometry.calculationType not in ['H', 'M', 'ML'] \
            and command[51] == '7':
        a = 0
        chMED = []
        fileout.write('{Property;\nDENSITY,STATE=' + str(state + 1) + '.1;\n')
        for i in range(len(CRG_emb)):
            if geometry.list_MM[i] in geometry.list_MEDIUM:
                fileout.write('EF,0,%14.8f,' % (geometry.pod[0][i] / covcost) + '%14.8f,' % (
                        geometry.pod[1][i] / covcost) + '%14.8f;\n' % (geometry.pod[2][i] / covcost))
                chMED.append(CRG_emb[i])
                if a == 25:
                    fileout.write('}\n{Property;\nDENSITY,STATE=' + str(state + 1) + '.1;\n')
                a = a + 1
        fileout.write('}')
        if command[1] == 'ci':
            a = 0
            chMED = []
            fileout.write('{Property;\nDENSITY,STATE=' + str(state) + '.1;\n')
            for i in range(len(CRG_emb)):
                if geometry.list_MM[i] in geometry.list_MEDIUM:
                    fileout.write('EF,0,%14.8f,' % (geometry.pod[0][i] / covcost) + '%14.8f,' % (
                            geometry.pod[1][i] / covcost) + '%14.8f;\n' % (geometry.pod[2][i] / covcost))
                    chMED.append(CRG_emb[i])
                    if a == 25:
                        fileout.write('}\n{Property;\nDENSITY,STATE=' + str(state) + '.1;\n')
                    a = a + 1
            fileout.write('}')
    else:
        chMED = []
    fileout.close()

    # write the file with the charges
    filelat = open('lattice', 'w')
    for i in range(len(lattice)):
        filelat.write(str(lattice[i]))
    filelat.close()

    # shelve the runtime data in cobram system-exchange file!
    sef = shelve.open("cobram-sef")
    sef['DEarray'] = DEarray
    sef['NAC'] = NAC
    sef['nroots'] = nroots
    sef['state'] = state
    if step == 0:
        sef['newstate'] = state
    sef.close()

    return chMED


def launch(command, step):
    command = command
    guess = command[6]
    Nproc = command[7]
    os.putenv('step', str(step))

    # define molpro executable from the environmental variable MOLPRO_EXE
    molpro_exe = os.getenv('MOLPRO_EXE')
    logwrt.writelog('Launching Molpro job ... ')
    os.system(molpro_exe + ' -d$PWD -I$PWD -W$PWD < molpro.input > molpro.log 2>&1')
    os.system('rm -f molpro.xml*;mv molpro.log_1 molpro.log')
    os.system('if [ -f molpro.log_2 ] ;then mv molpro.log_2 molpro.log ;fi')
    logwrt.writelog(' Done!\n')


# *******************************************************************
#
# now all tools for switching on an off the couplings computations
# (OW 09/2015)
#
# *******************************************************************
def chk_coupl(newstate, threshold, DEarray, calc_coupl):
    """check which couplings to compute in molpro and delete entries if necessary"""

    tmp_coupl = []
    if (len(calc_coupl)) <= 4:
        tmp_coupl = calc_coupl
    elif (len(calc_coupl)) > 4:
        logwrt.writewarning("""There are more than 4 couplings to compute!
This apparently does not work with the current version of molpro.
Trying to eliminate OFF-diagonal couplings of far separated states""")
        for i in range(len(calc_coupl)):
            if abs(DEarray[calc_coupl[i][0]][calc_coupl[i][1]]) >= threshold:
                logwrt.writelog("the states {0} and {1} are separated by {2} kcal/mol, deleting this entry\n".format(
                    calc_coupl[i][0] + 1, calc_coupl[i][1] + 1,
                    DEarray[calc_coupl[i][0]][calc_coupl[i][1]]))
            else:
                tmp_coupl.append(calc_coupl[i])
        #       calc_coupl=copy.deepcopy(tmp_coupl)
        logwrt.writelog("\nThe new couplings to compute are: {0}\n".format(tmp_coupl))
        if (len(tmp_coupl)) > 4:
            logwrt.fatalerror("""
********************************************************
**  FATAL - can not eliminate enough couplings, ABORT **
**      You can still try to lower the value of       **
**            command 86 and restart                  **
********************************************************""")
    return tmp_coupl


def del_coupl(tmp):
    """delete all couplings from tmp list"""
    tmp2 = []
    # first run, delete all couplings and force entries so we do not get a mess
    for i in range(len(tmp)):
        try:
            if tmp[i].split(',')[0].lower() == 'cpmcscf' and (tmp[i].split(',')[1].lower() == 'nacm'):
                logwrt.writelog("eliminating coupling line {0}\n".format(tmp[i]))
                tmp[i] = ''
        except:
            pass
        try:
            if tmp[i].split()[0].lower()[0:6] == "{force" and tmp[i + 1].split(',')[1][0:6] != "5100.1":
                # tmp[i-1]=''
                tmp[i] = ''
                tmp[i + 1] = ''
                tmp[i + 2] = ''
        except:
            pass
        try:
            if tmp[i].strip().lower()[0:9] == "{property":
                first_line = i
                for j in range(i + 1, len(tmp)):
                    if tmp[j].strip().lower()[0:1] == "}":
                        last_line = j
                        break
                for j in range(first_line, last_line + 1):
                    tmp[j] = ''
        except:
            pass
    # now remove blank lines so we have a clean input
    for i in range(len(tmp)):
        if tmp[i] != '':
            tmp2.append(tmp[i])
    tmp = copy.deepcopy(tmp2)
    return tmp


def add_coupl(tmp, calc_coupl, newstate):
    """add CPMCSCF entries into molpro.input"""

    logwrt.writelog('changing molpro input to calc derivative coupling\n')
    modi = 'F'
    rec = 0
    for i in range(len(tmp)):
        tmp1 = tmp[i].split(',')
        if tmp[i].split(',')[0].lower() == 'cpmcscf' and tmp[i].split(',')[1].lower() == 'grad':
            tmp[i] = 'CPMCSCF,GRAD,' + str(newstate + 1) + '.1,record=5100.1\n'
            rec = 1
            # print "FOUND entry, now changing, calc_coupl is:",calc_coupl
            for j in range(len(calc_coupl)):
                tmp[i] = tmp[i] + 'CPMCSCF,NACM,' + str(calc_coupl[j][0] + 1) + '.1,' + str(
                    calc_coupl[j][1] + 1) + '.1,record=510' + str(rec) + '.1;\n'
                rec = rec + 1
                # print "inserting into molpro.input: compute couplings for",calc_coupl[j][0],calc_coupl[j][1]
            # print "**********"
            logwrt.writelog("adding coupling line {0}".format(tmp[i]))
            modi = 'T'
            # print "modi has been set to T now",modi
    #    if modi == 'T' :
    for j in range(1, rec):
        tmp.append('{Force\nSAMC,510' + str(j) + '.1\n}')
        logwrt.writelog("Adding line: " + '{Force\nSAMC,510' + str(j) + '.1\n}\n')
    molpr = copy.deepcopy(tmp)
    return molpr


def write_molpro(molpr):
    filein = open('molpro.input', 'w')
    for i in range(len(molpr)):
        filein.write(str(molpr[i]) + '\n')
    filein.write('\n')
    filein.close()


def write_cobram(molpr):
    """adapt the cobram.command to have the new molpro entries also there"""
    tmp = []
    tmp1 = []
    modi = 'F'
    posi = 'F'
    # now do the same in cobram.command
    with open('cobram.command') as enefile:
        for line in enefile:
            tmp.append(line.strip())
    for i in range(len(tmp)):
        if tmp[i] == '!molpro':
            begin_molpro = i
            posi = 'T'
        if posi == 'T' and tmp[i].split(',')[0].lower() == 'cpmcscf' and tmp[i].split(',')[1].lower() == 'grad':
            begin_modi = i
        if tmp[i] == '?molpro':
            end_molpro = i

    for i in range(len(molpr)):
        if molpr[i].split(',')[0].lower() == 'cpmcscf' and molpr[i].split(',')[1].lower() == 'grad':
            molpr_start = i
    for i in range(begin_modi):
        tmp1.append(tmp[i])
    for i in range(molpr_start, len(molpr)):
        tmp1.append(molpr[i])
    for i in range(end_molpro, len(tmp)):
        tmp1.append(tmp[i])

    filein = open('cobram.command', 'w')
    for i in range(len(tmp1)):
        filein.write(str(tmp1[i]) + '\n')
    filein.close()
    # cobcom=CBF.getCobrammCommand('cobram.command')


def read_molpro():
    tmp = []
    with open('molpro.input') as enefile:
        for line in enefile:
            tmp.append(line.strip())
    return tmp


def couplings_on(DEarray, threshold, newstate, calc_coupl):
    logwrt.writelog("Coupling computations switched on\n")
    newcoupl = chk_coupl(newstate, threshold, DEarray, calc_coupl)
    if newcoupl:
        calc_coupl = copy.deepcopy(newcoupl)
    worklist = read_molpro()
    del_list = del_coupl(worklist)
    molpr = add_coupl(del_list, calc_coupl, newstate)
    write_molpro(molpr)
    # print "****************************"
    # print "writing new cobram.command"
    # print "****************************"
    write_cobram(molpr)


def couplings_off():
    logwrt.writelog("Switching off coupling computations in molpro\n")
    worklist = read_molpro()
    molpr = del_coupl(worklist)
    write_molpro(molpr)
    # print "****************************"
    # print "writing new cobram.command"
    # print "****************************"
    write_cobram(molpr)


# ******************************************************************
# End of tools for couplings
# ******************************************************************


def molproEneGradCrg(command, geometry, CRG_emb, step, chMED, cobcom):
    logwrt.writelog("\n")
    natom = len(geometry.modelH[0])
    CIrot = 'no'

    sef = shelve.open("cobram-sef", 'r')
    # take care we do not carry old entries if command 86 is active
    if float(command[86]) < float(1000):
        NAC = sef['NAC_empty']
    else:
        NAC = sef['NAC']
    NAC_old = sef['NAC_old']
    nroots = sef['nroots']
    newstate = sef['newstate']
    sef.close()

    gradient = []
    dipole = []
    pt2dipole = []
    for i in range(nroots):
        pt2dipole.append([])
    for i in range(nroots):
        for j in range(4):
            pt2dipole[i].append(0.0)
    gradch = []
    charges = []
    for i in range(nroots):
        gradient.extend([[], [], []])
        gradch.extend([[], [], []])
    tmp = []
    molproenergy = []
    pt2energy = []
    for i in range(nroots):
        pt2energy.append(0.0)
    mp2energy = []
    xdc, ydc, zdc = [], [], []
    CHhead, CHtail = [], []
    allCH = []
    state = 0
    a = 0
    CIhead, CHhead, = [], []
    kkk = []
    for i in range(len(resortat)):
        for j in range(len(resortat)):
            if i == resortat[j]:
                kkk.append(j)
    with open('molpro.log') as enefile:
        for line in enefile:
            tmp.append(line.strip())

    mp2check = 0
    pt2check = 0
    for i in range(len(tmp)):
        if tmp[i] == '** WVFN ****  MAXIMUM NUMBER OF ITERATIONS REACHED':
            if int(command[198]) == 1:
                logwrt.fatalerror('NO CONVERGENCY IN THE MOLPRO WF')
            else:
                logwrt.writelog('***********************************************************************')
                logwrt.writelog('**   WARNING!!!! No Convergency in molpro WF, but I will continue!!  **')
                logwrt.writelog('**                 Carefully check the result!!!                     **')
                logwrt.writelog('***********************************************************************')

        tmp1, tmp2, tmp3, tmp4 = [], [], [], []

        try:
            tmp1 = tmp[i].split('=')
            tmp2 = tmp[i].split()
            tmp3 = tmp[i].split(':')
            # if tmp2[0]== 'mp2;' or tmp2[0]== 'MP2' or tmp2[0]== 'Mp2' :
            #      print "This is an MP2 calculation, taking the MP2 energy as QM energy"
            #      mp2check=1   # skip RHF section
            if tmp2[0] == 'CI' and tmp2[1] == 'vector':
                # print 'found CI-Vector'
                CIhead.append(i)
                # print 'CIhead(0)=',CIhead
            if (tmp2[0] == '!MCSCF' or tmp2[0] == '!RKS') and tmp2[3].lower() == 'energy':
                molproenergy.append(float(tmp2[4]))
                logwrt.writelog("Energy of CASSCF state {0}: {1}\n".format(len(molproenergy), molproenergy[-1]))
            if tmp2[0] == '!RSPT2' and tmp2[3].lower() == 'energy':
                pt2check = 1
                istate = int(tmp2[2].split('.')[0]) - 1
                pt2energy[istate] = float(tmp2[4])
                logwrt.writelog("Energy of CASPT2 state {0}: {1}\n".format(istate + 1, pt2energy[istate]))
            if (tmp2[0] == '!MP2') and tmp2[2].lower() == 'energy':
                mp2check = 1
                mp2energy.append(float(tmp2[3]))
                logwrt.writelog("MP2 energy: {0}\n".format(mp2energy[0]))
            if (tmp2[0] == '!RHF') and tmp2[3].lower() == 'energy':
                molproenergy.append(float(tmp2[4]))
                logwrt.writelog("SCF energy: {0}\n".format(molproenergy[0]))
            # get dipole moments in atomic units (OW Jan. 2012)
            if tmp2[0] == '!RSPT2' and tmp2[3].strip().lower()[0:6] == 'dipole':
                istate = int(tmp2[2].split('.')[0]) - 1
                dix = float(tmp2[5])
                diy = float(tmp2[6])
                diz = float(tmp2[7])
                pt2dipole[istate][0] = dix
                pt2dipole[istate][1] = diy
                pt2dipole[istate][2] = diz
                pt2dipole[istate][3] = float(sqrt(dix * dix + diy * diy + diz * diz))
            elif tmp2[3].strip().lower()[0:6] == 'dipole':
                dipole.append([])
                dix = float(tmp2[5])
                diy = float(tmp2[6])
                diz = float(tmp2[7])
                dipole[-1].append(dix)
                dipole[-1].append(diy)
                dipole[-1].append(diz)
                dipole[-1].append(float(sqrt(dix * dix + diy * diy + diz * diz)))
            if (tmp2[0] == 'RSPT2' or tmp2[0] == 'SA-MC' or tmp2[0] == 'MC' or tmp2[0] == 'SCF' or tmp2[0] == 'RKS' or
                tmp2[0] == 'MP2' or tmp2[0] == 'RHF') and (tmp2[1] == 'GRADIENT'):
                state = int(tmp2[4].split('.')[0]) - 1
                if command[1] != 'ci' and state != newstate:
                    logwrt.fatalerror(
                        'The state ' + str(newstate + 1) + ' Cobram expects to find gradient for and the state ' + str(
                            state + 1) + ' for which Molpro has computed the gradient do not coincide in step' + str(
                            step) + '!')
                jb = 0
                blanks = 'no'
                while jb < natom:
                    if tmp[i + jb + 4] == '':
                        blanks = 'yes'
                        break
                    jb = jb + 1
                jb = 0
                while jb < natom:
                    k = kkk[jb]
                    if blanks == 'yes':
                        k1 = k / 50
                    elif (blanks == 'no'):
                        k1 = 0
                    gr = tmp[i + k + k1 + 4].split()
                    gradient[3 * state].append(-float(gr[1]))
                    gradient[3 * state + 1].append(-float(gr[2]))
                    gradient[3 * state + 2].append(-float(gr[3]))
                    jb = jb + 1
                # logwrt.writelog( "Force for state {0}\n".format(state+1 ) )
                for iatom in range(len(geometry.modelH[0])):
                    if abs(gradient[3 * state][iatom]) > 1.0 or abs(gradient[3 * state + 1][iatom]) > 1.0 or abs(
                            gradient[3 * state + 2][iatom]) > 1.0:
                        logwrt.writewarning('Large gradient during numerical computation detected. Possible problems!')
                        break
                # for iatom in range(len(geometry.modelH[0])):
                # logwrt.writelog( "{0:12.6f} {1:12.6f} {2:12.6f}\n".format( gradient[3*state][iatom],
                # gradient[3*state+1][iatom], gradient[3*state+2][iatom]) )
                if command[1] == 'ci' and state != newstate:
                    sef = shelve.open("cobram-sef")
                    sef['gradient2'] = [gradient[3 * (state)], gradient[3 * (state) + 1], gradient[3 * (state) + 2]]
                    sef.close()
                elif command[1] != 'ci' or (command[1] == 'ci' and state == newstate):
                    gradient = [gradient[3 * (state)], gradient[3 * (state) + 1], gradient[3 * (state) + 2]]
                # logwrt.writelog( "****************************\n" )

            if ((command[1] == 'mdv' and int(command[85]) > 3) or command[1] == 'ci') and (
                    tmp2[0] == 'SA-MC' or tmp2[0] == 'MC') and tmp2[1] == 'NACME':
                state1 = int(tmp2[4].split('.')[0])
                state2 = int(tmp2[6].split('.')[0])
                xdc = []
                ydc = []
                zdc = []
                jb = 0
                blanks = 'no'
                while jb < natom:
                    if tmp[i + jb + 4] == '':
                        blanks = 'yes'
                        break
                    jb = jb + 1
                jb = 0
                while jb < natom:
                    k = kkk[jb]
                    if blanks == 'yes':
                        k1 = k / 50
                    elif blanks == 'no':
                        k1 = 0
                    # logwrt.writelog('jb is '+str(jb)+' k is '+str(k)+' k1 is '+str(k1));
                    gr = tmp[i + k + k1 + 4].split()
                    # logwrt.writelog('gr is '+str(gr))
                    xdc.append(float(gr[1]))
                    ydc.append(float(gr[2]))
                    zdc.append(float(gr[3]))
                    jb = jb + 1
                NAC[state1 - 1][state2 - 1] = [xdc, ydc, zdc]
                for dc in range(3):
                    NAC[state2 - 1][state1 - 1][dc] = [-ele for ele in NAC[state1 - 1][state2 - 1][dc]]
                logwrt.writelog("NAC < {0} |d/dR| {1} >\n".format(state1, state2))
                for iatom in range(len(geometry.modelH[0])):
                    logwrt.writelog("{0:12.6f} {1:12.6f} {2:12.6f}\n".format(NAC[state1 - 1][state2 - 1][0][iatom],
                                                                             NAC[state1 - 1][state2 - 1][1][iatom],
                                                                             NAC[state1 - 1][state2 - 1][2][iatom]))
                logwrt.writelog("****************************\n")

            if tmp2[0] == 'OPERATOR' and tmp2[1] == 'EF':
                try:
                    tmp_state = int(tmp[i + 2].strip().split()[8].split('.')[0]) - 1
                    if a == 0:
                        logwrt.writelog("Found electric fields for state {0}\n".format(tmp_state + 1))
                except:
                    logwrt.fatalerror('Could not identify state in the PROPERTY routine')
                tmp2 = tmp[i + 12].split()
                gradch[3 * tmp_state].append(float(tmp2[2]) * (-chMED[a]))
                gradch[3 * tmp_state + 1].append(float(tmp2[3]) * (-chMED[a]))
                gradch[3 * tmp_state + 2].append(float(tmp2[4]) * (-chMED[a]))
                if a == len(chMED) - 1:
                    if int(command[2]) > 0:
                        logwrt.writelog("Maximal number of movable point charges reached\n")
                    a = 0
                else:
                    a = a + 1
            if tmp3[0] == 'Number of CSFs':
                nCI = int(tmp2[6])
                # print 'found CSFs, nCI=',nCI
            if tmp2[0] == '1PROGRAM' and tmp2[2] == 'POP' and tmp2[4] == 'population':
                for ii in range(i, len(tmp)):
                    try:
                        tmp2 = tmp[ii].split()
                        if tmp2[6] == 'Type=MCSCF/CHARGE':
                            if int(command[2]) > 0:
                                logwrt.writelog('Found MCSCF charges\n')
                            CHhead.append(ii + 5)
                            break
                        elif tmp2[6] == 'Type=RHF/CHARGE':
                            if int(command[2]) > 0:
                                logwrt.writelog('Found SCF charges\n')
                            CHhead.append(ii + 5)
                            break
                        elif tmp2[6] == 'Type=RKS/CHARGE':
                            if int(command[2]) > 0:
                                logwrt.writelog('Found RKS charges\n')
                            CHhead.append(ii + 5)
                            break
                        else:
                            logwrt.fatalerror('POP section found but charges are neither MCSCF, nor RHF nor RKS')
                    except SystemExit:
                        sys.exit()
                    except:
                        pass
        except:
            pass

    # selecting the right point charge gradient
    if command[1] == 'ci':
        sef = shelve.open("cobram-sef")
        sef['gradch2'] = [np.array(gradch[3 * (state - 1)]), np.array(gradch[3 * (state - 1) + 1]),
                          np.array(gradch[3 * (state - 1) + 2])]
        sef.close()
        # logwrt.writelog( 'Force of point charges for lower state:\n' )
        # for i in range(len(gradch[3*(state-1)])):
        # logwrt.writelog( "{0:12.6f} {1:12.6f} {2:12.6f}\n".format(gradch[3*(state-1)][i],
        # gradch[3*(state-1)+1][i],gradch[3*(state-1)+2][i]) )
    gradch = [np.array(gradch[3 * (state)]), np.array(gradch[3 * (state) + 1]), np.array(gradch[3 * (state) + 2])]

    if mp2check == 1:
        molproenergy = mp2energy
    if pt2check == 1:
        molproenergy = pt2energy
        dipole = pt2dipole

    # selecting the right set of dipole moments
    if command[1] == 'ci':
        sef = shelve.open("cobram-sef")
        sef['dipole2'] = dipole[state - 1]
        sef.close()
        # logwrt.writelog( 'Dipole moments for lower state:\n' )
        # for i in range(len(dipole[state-1])):
        # logwrt.writelog( "{0:6.2f}\n".format(dipole[state-1][i]) )
    dipole = dipole[state]

    for i in range(len(CHhead)):
        CH = []
        jb = 0
        while jb < natom:
            k = kkk[jb]
            NE = tmp[CHhead[i] + k].split()
            OLLI = NE[8] + NE[9]
            CH.append(float(OLLI))
            jb = jb + 1
        if i == state:
            charges = np.array(CH)
        if command[1] == 'ci' and i == state - 1:
            # logwrt.writelog( 'Charges for lower state:\n' )
            # for j in range(len(CH)):
            # logwrt.writelog( "{0:6.2f}".format(CH[j]) )
            sef = shelve.open("cobram-sef")
            sef['charges2'] = np.array(CH)
            sef.close()
        allCH.append(CH)

    newstate = state

    string1 = [3, 4, 5, 6, 7, 8, 9, 14, 18, 19, 20, 21, 25, 26, 27, 33, 34, 35]
    string2 = [1, 2, 10, 11, 12, 13, 15, 16, 17, 22, 23, 24, 28, 29, 30, 31, 32]

    if command[1] == 'ci':
        # actually only one value is required (DE[state][state-1])
        sef = shelve.open("cobram-sef")
        DEarray = sef['DEarray']
        for i in range(nroots - 1):
            for j in range(i + 1, nroots):
                DEarray[i][j] = 627.51 * (molproenergy[j] - molproenergy[i])
                # OW make array symmetric so doesn't matter if we want ij or ji
                DEarray[j][i] = DEarray[i][j]
        sef['DEarray'] = DEarray
        sef.close()

    if command[1] == 'mdv':

        # fill DE array
        sef = shelve.open("cobram-sef")
        DEarray = sef['DEarray']
        DE_oldarray = copy.deepcopy(DEarray)
        sef['DE_oldarray'] = DE_oldarray
        for i in range(nroots - 1):
            for j in range(i + 1, nroots):
                DEarray[i][j] = 627.51 * (molproenergy[j] - molproenergy[i])
                # OW make array symmetric so doesn't matter if we want ij or ji
                DEarray[j][i] = DEarray[i][j]
        sef['DEarray'] = DEarray
        sef.close()

        if int(command[85]) > 3:
            sef = shelve.open("cobram-sef")
            # check energy difference for states below and above the state of interest
            # and estimate which couplings have to be computed
            calc_coupl = []
            calc_coupl2print = []

            logwrt.writelog("Threshold for activating Tully FSSH {0}\n".format(command[86]))
            # activate NAC scheme 2
            if int(command[201]) == 1:
                # offd should not be overwirrten as it holds information which state besides
                # the active state have coefficients != 0
                if step == 0:
                    offd = []
                else:
                    offd = sef['offd']
            # else use NAC scheme 1
            else:
                offd = []
            for i in range(nroots):
                if abs(DEarray[i][state]) <= float(command[86]):
                    logwrt.writelog("{0} is smaller than {1}: compute couplings for {2} {3}\n".format(DEarray[i][state],
                                                                                                      command[86],
                                                                                                      i + 1, state + 1))
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
                                    "{0} is smaller than {1}: compute couplings for {2} {3}\n".format(DEarray[i][j],
                                                                                                      command[86],
                                                                                                      i + 1, j + 1))
                                calc_coupl.append([i, j])
                                calc_coupl2print.append([i + 1, state + 1])
                    sef['offd'] = offd
            else:
                if len(offd) > 1:
                    # print "updated offd", offd
                    for i in range(len(offd)):
                        for j in range(i + 1, len(offd)):
                            logwrt.writelog("also computing coupling {0} {1}\n".format(offd[i] + 1, offd[j] + 1))
                            calc_coupl.append([offd[i], offd[j]])
                            calc_coupl2print.append([offd[i] + 1, offd[j] + 1])
            logwrt.writelog("Final set of NACs to compute {0}\n".format(calc_coupl2print))
            sef['calc_coupl'] = calc_coupl
            sef.close()

        if int(step) == int('0'):
            CI22 = []
            CI12 = []
            s = []

            DE = DEarray[state - 1][state]
            if abs(DE) < float(command[86]):
                ttry = 'short'
            else:
                ttry = 'long'
            sef = shelve.open("cobram-sef")
            sef['ttry_old'] = ttry
            sef.close()
            os.system('echo ' + str(ttry) + '>TSTEP')
            try:
                for j in range(nCI):
                    s.append(tmp[j + CIhead[0] + 3].split())
                    s[j][1] = float(s[j][1])
                    s[j][2] = float(s[j][2])
                s.sort(key=lambda x: x[0])
                CIout = open('ci-coeff-1.dat', 'w')
                for i in range(len(s)):
                    CIout.write(' %14.10f' % s[i][1] + ' %14.10f\n' % s[i][2])
                CIout.close()
            except:
                pass

            if command[88] == str('1'):
                half1 = []
                half2 = []
                for i in range(1, len(allCH)):
                    sum1 = 0.0
                    sum2 = 0.0
                    for j in range(len(allCH[0])):
                        if j + 1 in string1:
                            sum1 = sum1 + allCH[i][j]
                        if j + 1 in string2:
                            sum2 = sum2 + allCH[i][j]
                    ##print sum1,sum2,sum1+sum2
                    half1.append(sum1)
                    half2.append(sum2)
                newstate = half1.index(max(half1)) + 1
                os.system('echo ' + str(newstate) + ' > states')
            else:  # Olli debug
                newstate = state  #

        else:
            if command[85] == '1' and ((state != int('0') and (
                    abs(DEarray[state][state - 1]) <= float(command[86]) or abs(DEarray[state][state + 1]) <= float(
                command[86]))) or (state == int('0') and abs(DEarray[state][state + 1]) <= float(command[86]))):

                CI11 = []
                CI22 = []
                CI12 = []
                CI21 = []
                s = []
                for j in range(nCI):
                    s.append(tmp[j + CIhead[0] + 3].split())
                    s[j][1] = float(s[j][1])
                    s[j][2] = float(s[j][2])

                s.sort(key=lambda x: x[0])

                CIin = open('ci-coeff-1.dat')
                vec = CIin.read().split('\n')
                CIin.close()
                for i in range(0, len(vec) - 1):
                    line = vec[i].split()
                    CI21.append(float(line[0]))
                    CI11.append(float(line[1]))
                    CI22.append(s[i][1])
                    CI12.append(s[i][2])
                CI21 = np.array(CI21)
                CI11 = np.array(CI11)
                CI22 = np.array(CI22)
                CI12 = np.array(CI12)

                CIF1 = (abs(np.add.reduce(CI11 * CI22)))
                CIF2 = (abs(np.add.reduce(CI21 * CI22)))
                CIF3 = (abs(np.add.reduce(CI21 * CI12)))
                CIF4 = (abs(np.add.reduce(CI12 * CI11)))

                CIout = open('ci-coeff-1.dat', 'w')
                for i in range(len(CI22)):
                    CIout.write(' %14.10f' % CI22[i] + ' %14.10f\n' % CI12[i])
                CIout.close()
                # set CIrot only when "orbital rotations" is requested, otherwise interferences
                # with Tully when DC direction is computed
                if command[85] == '1':
                    if abs(CIF1) >= 0.25 and abs(CIF3) >= 0.25:
                        if command[85] == str('1'):
                            newstate = int('0')
                        else:
                            CIrot = 'yes'
                    else:
                        newstate = state
                else:
                    newstate = state

                if command[88] == str('1'):
                    DE = DEarray[state - 1][state]
                    if float(command[89]) <= DE:
                        if str(newstate) == str(state):
                            half1 = []
                            half2 = []
                            for i in range(1, len(allCH)):
                                sum1 = 0.0
                                sum2 = 0.0
                                for j in range(len(allCH[0])):
                                    if j + 1 in string1:
                                        sum1 = sum1 + allCH[i][j]
                                    if j + 1 in string2:
                                        sum2 = sum2 + allCH[i][j]
                                half1.append(sum1)
                                half2.append(sum2)
                            newstate = half1.index(max(half1)) + 1
                os.system('echo ' + str(newstate) + ' > states')

            else:
                newstate = state

            # MASSEY PARAMETER and TULLY FEWEST SWICTHING METHOD
            if int(command[85]) > 3:
                newstate = tully(command, geometry, state, cobcom, step)
                logwrt.writelog("State after propagating the coefficients c_i in Tully FSSH {0}\n".format(newstate + 1))
        #                print "NAC vectors after tully:"
        #                #save NAC vector from previous step
        #                NAC_old=copy.deepcopy(NAC)
        #                print "NAC_old vectotrs after coyping",NAC_old
        #                sef=shelve.open("cobram-sef")
        #                sef['NAC_old']=NAC_old
        #                sef.close()

        # AN start 01.16
        if int(command[85]) > 3:
            threshold = float(command[86])
            if len(calc_coupl) > 0:
                if int(command[2]) > 0:
                    logwrt.writelog('Changing derivative coupling section in molpro.input '
                                    '(adding / modifing coupling lines if necessary)\n')
                couplings_on(DEarray, threshold, newstate, calc_coupl)
                cobcom = CBF.getCobrammCommand('cobram.command')
            else:
                if int(command[2]) > 0:
                    logwrt.writelog(
                        'Changing derivative coupling section in molpro.input (removing coupling lines if necessary)\n')
                couplings_off()
                cobcom = CBF.getCobrammCommand('cobram.command')
        ########################################
        # Is the code below equal the code above?
        ########################################
        #            print "state=", int(state)
        #            print "newstate=", int(newstate)
        #            print "DEarray[",int(state)-1,"][",int(state),"]=", DEarray[state-1][state]
        #            for k in range(nroots):
        #                for l in range(nroots):
        #                    print "DE[",k,"][",l,"]=", DEarray[k][l]
        #            DE=DEarray[state-1][state]
        #            threshold=float(command[86])
        #
        #            ## activate derivative coupling calculation on S1
        #            if abs(DE) < float(command[86]) and newstate !=0 :
        #                if int(command[85]) > 3:
        #                  print 'changing molpro input, calc derivative coupling'
        #                  couplings_on(DEarray,threshold,newstate,calc_coupl)
        #                  cobcom=CBF.getCobrammCommand('cobram.command')
        #
        #	    # activate derivative coupling calculation on S0
        #            if abs(DE) < float(command[86]) and newstate ==0 :
        #                if int(command[85]) > 3:
        #                  print "in S0 routine!"
        #                  print "newstate here is:", newstate
        #                  print 'changing molpro input, calc derivative coupling'
        #                  couplings_on(DEarray,threshold,newstate,calc_coupl)
        #                  cobcom=CBF.getCobrammCommand('cobram.command')
        #
        #	    # deactivate derivative coupling calculation
        #            if abs(DE) > float(command[86]):
        #                if int(command[85]) > 3 :
        #                  print "Deactivate DC"
        #                  print 'changing molpro input, deactivate derivative coupling calculation '
        #                  couplings_off()
        #                  cobcom=CBF.getCobrammCommand('cobram.command')

        # newstate must be saved to sef before the re-launch routine is called
        sef = shelve.open("cobram-sef")
        sef['newstate'] = newstate
        sef.close()
        # AN end 01.16

        if str(newstate) != str(state):
            if int(command[2]) > 1:
                logwrt.writelog('Changing gradient section in molpro.input after hopping\n')
            tmp = []
            tmp1 = []
            with open('molpro.input') as enefile:
                for line in enefile:
                    tmp.append(line.strip())
            for i in range(len(tmp)):
                try:
                    tmp1 = tmp[i].split(',')
                    if tmp[i].split(',')[0].lower() == 'cpmcscf' and tmp[i].split(',')[1].lower() == 'grad':
                        tmp[i] = 'CPMCSCF,GRAD,' + str(newstate + 1) + '.1,record=5100.1'
                    if tmp[i].strip() == '{Property;':
                        tmp[i + 1] = 'DENSITY,STATE=' + str(newstate + 1) + '.1;'
                except:
                    pass
            filein = open('molpro.input', 'w')
            for i in range(len(tmp)):
                filein.write(str(tmp[i]) + '\n')
            filein.close()

            if int(command[2]) > 0:
                logwrt.writelog('Changing gradient section in cobram.command after hopping\n')
            tmp = []
            tmp1 = []
            with open('cobram.command') as enefile:
                for line in enefile:
                    tmp.append(line.strip())
            for i in range(len(tmp)):
                try:
                    tmp1 = tmp[i].split(',')
                    if tmp[i].split(',')[0].lower() == 'cpmcscf' and tmp[i].split(',')[1].lower() == 'grad':
                        tmp[i] = 'CPMCSCF,GRAD,' + str(newstate + 1) + '.1,record=5100.1'
                except:
                    pass
            filein = open('cobram.command', 'w')
            for i in range(len(tmp)):
                filein.write(str(tmp[i]) + '\n')
            filein.close()
            cobcom = CBF.getCobrammCommand('cobram.command')

    # dynamical weighting (OW 05/2009)
    if (command[92] != '0') and (state == 0):
        logwrt.writelog("#####################################\n")
        logwrt.writelog("starting dynamical weighting routine \n")
        logwrt.writelog("reducing S1 weight to zero after hop \n")
        logwrt.writelog("works only with 2 roots at the moment\n")
        logwrt.writelog("#####################################\n")
        logwrt.writelog("Energy difference for start={0}\n".format(command[92]))
        try:
            wchk = open('wchk', 'r')
            wstep = wchk.read()
        except:
            wchk = open('wchk', 'w')
            wstep = 1
            wchk.write(str(wstep) + '\n')
        wchk.close()
        if int(wstep) < 21:
            try:
                if DE > float(command[92]):
                    enefile = open('cobram.command', 'r')
                    tmp = []
                    tmp1 = []
                    newstr = ''
                    for line in enefile:
                        tmp.append(line.strip())
                    enefile.close()
                    nweigh = abs((exp(-(int(wstep) * 0.12 - 2.52) ** 2) - 1))
                    wstep = int(wstep) + 1
                    wchk = open('wchk', 'w')
                    wchk.write(str(wstep) + '\n')
                    if wstep == 21:
                        logwrt.writelog("we have arrived in S0, switching off\n")
                        nweigh = 0
                    logwrt.writelog('new weight for S1: {0}'.format(nweigh))
                    if int(command[2]) > 0:
                        logwrt.writelog('Changing weigths in the CASSCF routine in cobram.command')
                    for i in range(len(tmp)):
                        try:
                            tmp1 = tmp[i].split(',')
                            if tmp1[0] == 'Weight':
                                tmp[i] = tmp1[0] + ',' + tmp1[1] + ',' + str(nweigh) + ';'
                            if (wstep == 21) and (tmp1[0] == 'CPMCSCF'):
                                tmp[i] = '}'
                                ck = i + 1
                            if (wstep == 21) and (tmp1[0] == 'FORCE;'):
                                tmp[i] = ''
                                tmp[ck] = 'FORCE;'
                        except:
                            pass
                    # change command file
                    filein = open('cobram.command', 'w')
                    for i in range(len(tmp)):
                        filein.write(str(tmp[i]) + '\n')
                    filein.close()
                    wchk.close()
            except:
                pass
    # end dynamically switch off SA

    # Brute SA switch off
    if (command[93] != '0') and (state == 0):
        try:
            wchk = open('wchk', 'r')
            wstep = wchk.read()
        except:
            wchk = open('wchk', 'w')
            wstep = 1
            wchk.write(str(wstep) + '\n')
        wchk.close()
        if int(wstep) > 0:
            try:
                # DE=(molproenergy[1]-molproenergy[0])*float('627.51')
                logwrt.writelog('DE in SA check={0} and wstep={1}\n'.format(DE, wstep))
                if float(DE) > float(command[93]):
                    logwrt.writelog('I am in SA-if clause\n')
                    enefile = open('cobram.command', 'r')
                    tmp = []
                    tmp1 = []
                    for line in enefile:
                        tmp.append(line.strip())
                    enefile.close()
                    logwrt.writelog('switching off stateaveraging at dE= {0} kcal/mol\n'.format(command[93]))
                    logwrt.writelog('changing cobram input\n')
                    wstep = '0'
                    wchk = open('wchk', 'w')
                    wchk.write(str(wstep) + '\n')
                    wchk.close()
                    for i in range(len(tmp)):
                        try:
                            tmp1 = tmp[i].split(',')
                            if tmp1[0] == 'Weight':
                                tmp[i] = tmp1[0] + ',' + tmp1[1] + ',' + str('0') + ';'
                            if tmp1[0] == 'CPMCSCF':
                                logwrt.writelog("Adding a bracket in line 1353\n")
                                tmp[i] = '}'
                                ck = i + 1
                            if tmp1[0] == 'FORCE;':
                                tmp[i] = ''
                                tmp[ck] = 'FORCE;'
                        except:
                            pass
                    # change command file
                    filein = open('cobram.command', 'w')
                    for i in range(len(tmp)):
                        filein.write(str(tmp[i]) + '\n')
                    filein.close()
            except:
                pass
    # Intruder state detection (OW 05/2009)
    if (command[94] != '0') and (state != 0):
        logwrt.writelog("###################################\n")
        logwrt.writelog("starting intruder detection routine\n")
        logwrt.writelog("###################################\n")
        enefile = open('cobram.command', 'r')
        tmp = []
        tmp1 = []
        nstates = 0
        newstr = ''
        deintr = 0
        for line in enefile:
            tmp.append(line.strip())
        enefile.close()
        # read how many states we have
        for i in range(len(tmp)):
            try:
                tmp1 = tmp[i].split(',')
                if tmp1[0] == 'State':
                    nstates = int(tmp1[1].strip(';'))
            except:
                pass
        deintr = abs(molproenergy[nstates - 1] - molproenergy[nstates - 2]) * float('627.51') * (-1)
        logwrt.writelog("deintr={0}\n".format(deintr))
        if deintr > float('-30'):
            #             nweigh=float(exp(-0.005*deintr**2))*float(command[92])
            #             nweigh=float((1.1/tanh(1.1*(deintr+20))-0.1/tanh(0.1*(deintr+20))+1)/4)*float(command[92])
            #           use SIGMOIDAL
            nweigh = float((1 / (1 + exp(-0.4 * (deintr + 18))) / 2) * float(command[94]))
            if (command[95] and command[96] != ''):
                p1 = float(command[95])
                p2 = float(command[96])
                nweigh = float((1 / (1 + exp(p1 * (deintr + p2))) / 2) * float(command[94]))
            if nweigh < float(0.01):
                nweigh = 0
            logwrt.writelog('new weight for state {0}: {1}\n'.format(nstates, nweigh))
            logwrt.writelog('changing cobram.command file\n')
            for i in range(len(tmp)):
                try:
                    tmp1 = tmp[i].split(',')
                    if tmp1[0] == 'Weight':
                        for j in range(nstates):
                            newstr = newstr + tmp1[j] + ','
                        tmp[i] = newstr + str(nweigh) + ';'
                except:
                    pass
            # change command file
            filein = open('cobram.command', 'w')
            for i in range(len(tmp)):
                filein.write(str(tmp[i]) + '\n')
            filein.close()

    # end intruder state detection

    # newstate, 'jfejffpowefref'
    # CBF.Totalenergy requires an arrays with all energies
    energy = list(molproenergy)
    # energy=molproenergy[int(newstate)]
    Selfenergy = 0.0
    #    dipole=[0.0,0.0,0.0,0.0]

    #    grdQ=input('lattice.grad')
    #    a,b=0,0
    #    gradch=[[],[],[]]
    #    for line in grdQ :
    #        if b == 0 or b == 1 :
    #           pass
    #        else:
    #          if geometry.list_MM[a] in geometry.list_MEDIUM:
    #            ele=line.split()
    #            gradch[0].append(float(ele[0]))
    #            gradch[1].append(float(ele[1]))
    #            gradch[2].append(float(ele[2]))
    #            print a, 'gradch = ',  gradch[0]
    #          a=a+1
    #        b=b+1
    #    grdQ.close()

    # in case of a hop Cobram re-launches the QM computation to obtain the new gradient
    # printing to the dynamics.out should be controlled, otherwise some data is listed twice
    # which affects cobram-mdv-grep.py
    # controll is achieved via storing the step and comparing to the previous run
    sef = shelve.open("cobram-sef")
    old_step = sef['old_step']
    sef.close()
    if command[1] == 'mdv' and step == old_step + 1:
        sef = shelve.open("cobram-sef")
        sef['old_step'] = step
        sef.close()
        space = []
        # print 'newstate=',str(newstate),'current state=',str(state) #Olli debug
        if str(newstate) == str(state):
            R = [energy, gradient, charges, Selfenergy, dipole, gradch]
            Results = R
        else:
            launch(command, step)
            Results = molproEneGradCrg(command, geometry, CRG_emb, step, chMED, cobcom)
            # this one is for fans of "Inception": when a HOP occurs, Cobram re-runs the last QM computation with a
            # new gradient, and collects all the data for the new state (gradient, EFs, charges, dipole moments);
            # so Cobram goes once more through a nested molcasEneGrdCrg routine, this time avoiding this section and
            # stores them in sef after that Cobram returns here to exit the parent molcasEneGrdCrg routine, where
            # gradients, charges, etc. are assigned to the previous state as at the end of the routine
            # everything has to be stored in sef we need to load all up-to-date values from sef
            # as they have been modified in the nested molcasEneGrdCrg routine
            sef = shelve.open("cobram-sef", 'r')
            NAC = sef['NAC']
            gradient = sef['gradient']
            gradch = sef['gradch']
            dipole = sef['dipole']
            charges = sef['charges']
            molproenergy = sef['molproenergy']
            sef.close()
    else:
        R = [energy, gradient, charges, Selfenergy, dipole, gradch]
        Results = R

    # define a logical variable to compute single point along the path/trajectory
    # default is non to compute the point
    computeSP = False

    # decide whther to perform a single point computation in case of IRC
    if command[200] == '-1':
        if step == 1:
            computeSP = True
        elif step > 1:
            for line in reversed(open('geometry.log').readlines()):
                tmp.append(line.strip())
            for i in range(len(tmp)):
                try:
                    if tmp[i].find('Converged?') != -1:
                        for j in (range(1, 5)):
                            if tmp[i - j].split()[4] == 'YES':
                                computeSP = True
                            elif tmp[i - j].split()[4] == 'NO':
                                computeSP = False
                                break
                        break
                except:
                    pass
        if computeSP: logwrt.writelog("\nComputing a single point at previous step {0}\n".format(step - 1))

    # perform a single point computation every N steps
    if int(command[200]) > 0 and step % int(command[200]) == 0:
        computeSP = True
        logwrt.writelog("Computing a single point at step {0}\n".format(step))
        for fileName in ["lattice", "molpro.input", "work.wfu", "WORK.wfu"]:
            if os.path.exists(fileName): shutil.copy(fileName, os.path.join(SINGLEPOINT_DIR, fileName))

    # prepare and execute single point calculation
    if computeSP:
        basis_key = 0
        section = CBF.ReadCobramCommand(cobcom, 'singlepoint', 'singlepoint.input')
        for i in range(len(section)):
            if section[i].strip().lower()[0:5] == 'basis':
                logwrt.writelog("Basis set provided by the user\n")
                basis_key = 1
        filein = open(os.path.join(SINGLEPOINT_DIR, "molpro.input"))
        tmp = []
        for line in filein:
            tmp.append(line)
        filein.close()
        conto = 'yes'
        newinp = []
        for i in range(len(tmp)):
            if basis_key == 1 and tmp[i].strip().lower()[0:5] == 'basis':
                conto = 'no'
            elif basis_key == 0 and tmp[i].strip().lower()[0:7] == 'basis={':
                for j in range(i, len(tmp)):
                    newinp.append(tmp[j])
                    if '}' in tmp[j]:
                        break
                conto = 'no'
            elif basis_key == 0 and tmp[i].strip().lower()[0:6] == 'basis=':
                newinp.append(tmp[i])
                conto = 'no'
            elif basis_key == 0 and tmp[i].strip().lower()[0:6] == 'basis,':
                newinp.append(tmp[i])
                conto = 'no'
            elif basis_key == 0 and tmp[i].strip().lower()[0:5] == 'basis':
                for j in range(i, len(tmp)):
                    newinp.append(tmp[j])
                    if 'end' in tmp[j]:
                        break
                conto = 'no'
            if conto == 'yes':
                newinp.append(tmp[i])
        filein.close()
        filein = open(os.path.join(SINGLEPOINT_DIR, "molpro.input"), 'w')
        for i in range(len(newinp)):
            filein.write(newinp[i])
        for i in range(len(section)):
            filein.write(section[i] + '\n')
        filein.close()
        guess = command[6]
        Nproc = command[7]

        # define molpro executable from the environmental variable MOLPRO_EXE
        molpro_exe = os.getenv('MOLPRO_EXE')
        os.system(molpro_exe + ' -d$PWD/{0} -I$PWD/{0} -W$PWD/{0} < {0}/molpro.input > {0}/molpro.log 2>&1'.format(
            SINGLEPOINT_DIR))

        if os.path.exists(os.path.join(SINGLEPOINT_DIR, "molpro.log_1")):
            shutil.move(os.path.join(SINGLEPOINT_DIR, "molpro.log_1"), os.path.join(SINGLEPOINT_DIR, "molpro.log"))
        if os.path.exists(os.path.join(SINGLEPOINT_DIR, "molpro.log_2")):
            shutil.move(os.path.join(SINGLEPOINT_DIR, "molpro.log_2"), os.path.join(SINGLEPOINT_DIR, "molpro.log"))

        # write the content of the QM single point calc in the ALL file, decorated with
        # comments that highligh the step number of the QM single point
        # name of the output file
        allName = os.path.join(SINGLEPOINT_DIR, "molproALL.log")
        # decide the step number depending on the type of calculation (IRC or not)
        if command[200] == '-1':
            spStepNr = step - 1
        else:
            spStepNr = step
        # name of the Molpro log file to cat into the allName file
        logName = os.path.join(SINGLEPOINT_DIR, "molpro.log")

        # open outfile and append Molpro log
        with open(allName, 'a') as molprotot:
            molprotot.write('=' * 80 + '\n')
            molprotot.write("Start SP Molpro calculation of STEP : " + str(spStepNr) + "\n")
            molprotot.write(' || ' * 20 + '\n')
            molprotot.write(' \/ ' * 20 + '\n')
            with open(logName, "r") as molproOut:
                molprotot.write(molproOut.read())
            molprotot.write(' /\ ' * 20 + '\n')
            molprotot.write(' || ' * 20 + '\n')
            molprotot.write("End SP Molpro calculation of STEP : " + str(spStepNr) + "\n")
            molprotot.write('=' * 80 + '\n')
            molprotot.write('\n')

        logwrt.writelog("Output stored in {0}/molproALL.log\n".format(SINGLEPOINT_DIR))

    if command[200] == '-1':
        for fileName in ["lattice", "molpro.input", "work.wfu", "WORK.wfu"]:
            if os.path.exists(fileName): shutil.copy(fileName, os.path.join(SINGLEPOINT_DIR, fileName))

    os.system('echo ' + str(newstate) + ' > STATE')

    # save all data to sef
    # any newly added sef-store command must be added above to the sef-restore for the case of a HOP
    sef = shelve.open("cobram-sef")
    sef['NAC'] = NAC
    sef['gradient'] = gradient
    sef['gradch'] = gradch
    # if len(gradch[0]) > 0 and int(command[60]) > 1:
    # logwrt.writelog( 'Force of point charges for upper state:\n' )
    # for i in range(len(gradch[0])):
    # logwrt.writelog( "{0:12.6f} {1:12.6f} {2:12.6f}\n".format(gradch[0][i], gradch[1][i], gradch[2][i]) )
    sef['dipole'] = dipole
    # if int(command[60]) > 1:
    # logwrt.writelog( 'Dipole moments for upper state:\n' )
    # for i in range(len(dipole)):
    # logwrt.writelog( "{0:6.2f}\n".format( dipole[i] ) )
    sef['charges'] = charges
    # if len(charges) > 0 and int(command[60]) > 1:
    # logwrt.writelog( 'Charges for upper state:\n' )
    # for i in range(len(charges)):
    # logwrt.writelog( "{0:6.2f}\n".format( charges[i] ) )
    sef['molproenergy'] = molproenergy
    sef.close()

    return Results


# read DC and GS energy after DC computation during Molcas MD
def readDC(command, geometry, step):
    sef = shelve.open("cobram-sef", 'r')
    aux_data = sef['aux_data']
    scal_factors = aux_data[18]
    NAC = sef['NAC_empty']
    sef.close()

    natom = len(geometry.modelH[0])
    tmp = []
    err_state = 1
    enefile = open('molpro.log', 'r')
    for line in enefile:
        tmp.append(line.strip())
    enefile.close()
    os.system('cat molpro.log >> molproALL.log')
    os.system('rm molpro.log')
    kkk = []
    for i in range(len(resortat)):
        for j in range(len(resortat)):
            if i == resortat[j]:
                kkk.append(j)
    for i in range(len(tmp)):
        if tmp[i].find('Variable memory released') != -1:
            err_state = 0
        try:
            tmp2 = tmp[i].split()
            if tmp[i].find('MCSCF STATE 1.1 Energy') != -1:
                molproenergyS0 = float(tmp[i].split()[4])
            if (tmp2[0] == 'SA-MC' and tmp2[1] == 'NACME') or (tmp2[0] == 'MC' and tmp2[1] == 'NACME'):
                xdc, ydc, zdc = [], [], []
                state1 = int(tmp2[4].split('.')[0])
                state2 = int(tmp2[6].split('.')[0])
                jb = 0
                blanks = 'no'
                while jb < natom:
                    if tmp[i + jb + 4] == '':
                        blanks = 'yes'
                        break
                    jb = jb + 1
                jb = 0
                while jb < natom:
                    k = kkk[jb]
                    if blanks == 'yes':
                        k1 = k / 50
                    elif blanks == 'no':
                        k1 = 0
                    logwrt.writelog('jb is ' + str(jb) + ' k is ' + str(k) + ' k1 is ' + str(k1))
                    gr = tmp[i + k + k1 + 4].split()
                    logwrt.writelog('gr is ' + str(gr))
                    xdc.append(float(gr[1]))
                    ydc.append(float(gr[2]))
                    zdc.append(float(gr[3]))
                    jb = jb + 1
                if os.path.exists(os.getcwd() + '/SS-CASPT2'):
                    xdc[:] = [j * scal_factors[state1 - 1][state2 - 1] for j in xdc]
                    ydc[:] = [j * scal_factors[state1 - 1][state2 - 1] for j in ydc]
                    zdc[:] = [j * scal_factors[state1 - 1][state2 - 1] for j in zdc]
                    logwrt.writelog(
                        "Rescaling CASSCF NAC {0} {1} with a scaling factor (relevant only for SSPT2): {2} ".format(
                            state1, state2, scal_factors[state1 - 1][state2 - 1]))
                NAC[state1 - 1][state2 - 1] = [xdc, ydc, zdc]
                for dc in range(3):
                    NAC[state2 - 1][state1 - 1][dc] = [-ele for ele in NAC[state1 - 1][state2 - 1][dc]]

                logwrt.writelog("NAC < {0} |d/dR| {1} >\n".format(state1, state2))
                for iatom in range(len(geometry.modelH[0])):
                    logwrt.writelog("{0:12.6f} {1:12.6f} {2:12.6f}\n".format(NAC[state1 - 1][state2 - 1][0][iatom],
                                                                             NAC[state1 - 1][state2 - 1][1][iatom],
                                                                             NAC[state1 - 1][state2 - 1][2][iatom]))
                logwrt.writelog("****************************\n")
        except:
            pass
    if err_state == 1:
        logwrt.fatalerror('Something went wrong with the Molpro computation!')
    sef = shelve.open("cobram-sef")
    sef['NAC'] = NAC
    return molproenergyS0


def save_QM_molpro_step(step, command, copyLog, copyOrb):
    """ save Molpro output and orbital files for last single point QM calculation
        in a common directory where QM data is stored """

    # name of the gaussian calculation (common to chk and log file)
    logName = 'molpro.log'
    # file where to store QM results
    allName = os.path.join(QM_DATA_STORAGE, 'molproALL.log')

    if copyLog:
        # check if the file exists, otherwise print a warning to screen
        if not os.path.exists(logName):
            logwrt.writewarning(
                "Molpro file {0} cannot be found: the {1} file will not be updated".format(logName, allName))

        else:
            # write the content of the QM single point calc in the ALL file, decorated with
            # comments that highligh the step number of the QM single point
            with open(allName, 'a') as molprotot:

                molprotot.write('=' * 80 + '\n')
                molprotot.write("Start SP Molpro calculation of STEP : " + str(step) + "\n")
                molprotot.write(' || ' * 20 + '\n')
                molprotot.write(' \/ ' * 20 + '\n')
                with open(logName, "r") as molproOut:
                    molprotot.write(molproOut.read())
                molprotot.write(' /\ ' * 20 + '\n')
                molprotot.write(' || ' * 20 + '\n')
                molprotot.write("End SP Molpro calculation of STEP : " + str(step) + "\n")
                molprotot.write('=' * 80 + '\n')
                molprotot.write('\n')

    # names of the wavefunction files
    orbFile = "work.wfu"
    targetOrbFile = os.path.join(QM_DATA_STORAGE, 'work' + str(step) + '.wfu')
    # names of the molden files
    srcMoldenFile = MOLDEN_FILE
    targetMoldenFile = os.path.join(QM_DATA_STORAGE, 'step_' + str(step) + '_scf.molden')

    # copy orbitals and molden file when required
    if copyOrb:

        # check if the orbital file exists and save it and zip it, otherwise print a warning to screen
        if not os.path.exists(orbFile):
            logwrt.writewarning("Molpro file {0} cannot be found: it will not be stored".format(orbFile))
        else:
            shutil.copyfile(orbFile, targetOrbFile)
            CBF.GZIP(targetOrbFile)

        # check if the molden file exists and save it and zip it, otherwise print a warning to screen
        if not os.path.exists(srcMoldenFile):
            logwrt.writewarning("Molpro file {0} cannot be found: it will not be stored".format(srcMoldenFile))
        else:
            shutil.copyfile(srcMoldenFile, targetMoldenFile)
            CBF.GZIP(targetMoldenFile)


def clean_QM_molpro():
    """ clean up the run directory from all the files that have been used to run MOLPRO and that
        are no longer needed at the end of the calculation """

    # list of files/directory to remove
    toRemove = ["ci-coeff-1.dat", "lattice", "STATE", "work.wfu", MOLDEN_FILE,
                "molpro.input", "molpro.log", "WFU", "AmplitudesALL.dat", ]

    for f in toRemove:
        # if f is an existing file, remove it
        if os.path.isfile(f):
            os.remove(f)
        # if f is an existing directory, clean tree
        elif os.path.isdir(f):
            shutil.rmtree(f)

    # done
    return
