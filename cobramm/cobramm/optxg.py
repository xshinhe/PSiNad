#! /usr/bin/env python3
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

# TODO: the infinite while loops should have a safety mechanism that kill the
#  running job after some time, otherwise there is the possibility that
#  the COBRAMM jobs hangs indefinitely in case of error

# TODO: the conversion factor bohr -> ang should be defined in a separate module
#  with all physical constants

#####################################################################################################

# import statments of module from python standard library

import glob  # unix style pathname pattern expansion
import subprocess  # run external program as child process
import re  # regular expression operations
import os  # operating system utilities
import time  # time functions

# imports of local modules

import CBF
import logwrt  # manages log file output + start/end procedures
import cobrammenv  # environmental variable for COBRAMM and 3rd-party software
import constants  # physical constants and conversion factors

# imports of user-defined classes

from timer import Timer  # keep timings of the different sections of the code
from extSignal import extSignal  # communication protocol between processes

# math libraries

import numpy as np  # numpy library for scientific computation


#####################################################################################################

@Timer("gau optimizer")
def GAUSSIAN_optimizator_X(step, geometry, command, cobcom, E_tot, Tot_Grad, Tot_dipole):

    # waiting time between each iteration in infinite loops (in seconds)
    wait = 0.1
    # maximum waiting time in infinite loops (in minute)
    maxwait = 100.

    # preparation is done only at step 0 (write input file, launch gaussian optimizer)
    if step == 0:

        # check that the environment for gaussian is actually defined
        envDefined, errorMsg = cobrammenv.checkGaussianOptEnv()
        if not envDefined: logwrt.fatalerror(errorMsg)

        if command[1] == 'freqxg':
            logwrt.writelog('COBRAMM is now setting up the frequency calculation\n')
        elif command[1] in ['optxg', 'irc', 'ci', 'ts']:
            logwrt.writelog('COBRAMM is now setting up the geometry optimization\n')
        logwrt.writelog('The "external" module of Gaussian will be used\n\n')

        # get from environment the variables for the GAUSSIAN program
        gaussian_exe = os.getenv('GAUSSIAN_EXE')
        gpath = os.getenv('GAUSSIAN_DIR')
        # define profile for running the version of gaussian defined by gaussian_exe and gpath
        cobrammenv.setGaussianProfile(gaussian_exe, gpath)

        # prepare input file for Gaussian optimization
        geometry_optxg_inp(geometry, command, cobcom, E_tot, Tot_Grad)

        # check if the cobramext command is available
        if not cobrammenv.which("cobramext"):
            logwrt.fatalerror("cobramext command is not available!\n Please add the COBRAM_PATH/cobramm directory to " +
                              "the PATH env variable")

            # write messages to log
        logwrt.writelog("Launching Gaussian optimizer in the background!\n")
        logwrt.writelog("Optimization is being performed by " + gaussian_exe + "\n\n")

        # launch Gaussian optimization with external option
        subprocess.Popen("nice -20 {0} < geometry.com 1> geometry.log 2> /dev/null &".format(gaussian_exe), shell=True)


        # wait until the file geometry.dat is present in the working dir (only in the first step,
        # in subsequent step the geometry.dat is always present)
        while True:
            if os.path.exists('geometry.dat'): break
            time.sleep(wait)

            # it may happen that the geometry.dat has not been written, because the geometry calculation terminated
            # with error check if the geometry.log contains error termination, in this case stop execution
            if os.path.exists('geometry.log'):
                gauTerm = checkGauTermination('geometry.log')
                if gauTerm == 2: logwrt.fatalerror('Gaussian geometry optimization terminated with error, please ' +
                                                   'check the file geometry.log for details')

    # now read the names of external communication files of Gaussian from the geometry.dat file
    gaussian_exe = os.getenv('GAUSSIAN_EXE')
    time.sleep(wait)
    with open('geometry.dat', "r") as f:
        q = f.readline().split()
    # and store the names in the gnames list
    gnames = []
    if gaussian_exe == 'g09' or gaussian_exe == 'g16':
        gnames = [q[0], q[1], q[2], q[3]]
    elif gaussian_exe == 'g03':
        gnames = ['R', q[0], q[1], 'N/A']

    # when this is not the step 0, we need to give the QMMM forces to Gau, and let it read them and use
    # them to do an optimization step
    if step != 0:
        # open the external output / gaussian input file and write QMMM energy and gradient
        with open(gnames[2], "w") as fileout:
            fileout.write(
                '%20.12e' % E_tot + '%20.12e' % Tot_dipole[0] + '%20.12e' % Tot_dipole[1] + '%16.8e\n' % Tot_dipole[
                    2])  # write energy and dipole
            for i in range(geometry.NatomHM):  # lists[7] is list of atoms belonging to the H and M layers
                fileout.write('%20.12e' % Tot_Grad[0][i] + '%20.12e' % Tot_Grad[1][i] + '%20.12e' % Tot_Grad[2][
                    i] + '\n')  # write gradient

        # define object for communicating with COBRAM
        comm = extSignal()
        # send signal to cobramext on the channel 0
        comm.sendSignal(0)
        # now wait till get signal from cobramext on the channel 1
        comm.waitTillSignal(1)
        # now gaussian can proceed: it will read the forces from gnames[2] and perform one optimization step

    # now wait until one of two conditions are met: a) gaussian has stopped or b) it has written a gnames[1] file
    wait_tot = 0.0
    while True:

        # attempt to read the new geometry from the gnames[1] file
        try:
            # read the content of the file and store the list of lines
            with open(gnames[1], "r") as inp:
                fileLines = inp.readlines()
            # check the number of moving atoms in the first lines
            natoms = float(fileLines[0].split()[0])
            # if the file is long enough and the number of atoms corresponds, break the loop
            if len(fileLines) >= natoms + 1 and geometry.NatomHM == natoms: break
        except:
            pass

        # check if the gaussian geometry.log file has the termination line (because of convergence
        # reached or no more iteration available)
        gauTerm = checkGauTermination('geometry.log')
        if gauTerm != 0: 
            break

        # if the waiting time has been too much, abort COBRAMM calculation
        if wait_tot > 1. * maxwait:
            logwrt.fatalerror("Something went wrong with the Gaussian optimizer. " +
                              "Waited for {} minutes and still can't find a new geometry".format(maxwait) +
                              "\nUsually this means that there is a ghost external optimizer process (g16 -> l402.exe -> cobramext) from a previous COBRAMM run in the same working directory." + 
                              "\nPlease check the presence of ghost 'cobramext' processes and kill them before trying again.\n")

        # wait until the gnames[1] file is written by gaussian
        time.sleep(wait)
        wait_tot += wait

    gauTerm = checkGauTermination('geometry.log')

    # if the optimization is converged, return None
    if gauTerm != 0:
        newGeom = None

    # else read the geometry from the gnames[1] file, the input of the Gau_External
    # (it is the output of gaussian and contains the new geometry)
    else:
        with open(gnames[1], "r") as inp:
            fileLines = inp.readlines()

        # the first line of the gnames[1] file contains some general information
        _splitted = fileLines[0].split()
        natoms = int(_splitted[0])  # number of atoms
        # devreq = _splitted[1]  # requested derivative
        # Chm = int(_splitted[2])  # charge of high-medium
        # Mhm = int(_splitted[3])  # spin of high-medium

        # the rest of the file has the geometry
        an, x, y, z, c = [], [], [], [], []
        for i in range(1, natoms + 1):
            _splitted = fileLines[i].split()
            an.append(int(_splitted[0]))
            x.append(float(_splitted[1]) * constants.Bohr2Ang)
            y.append(float(_splitted[2]) * constants.Bohr2Ang)
            z.append(float(_splitted[3]) * constants.Bohr2Ang)
            c.append(float(_splitted[4]) * constants.Bohr2Ang)
        newGeom = [np.array(x), np.array(y), np.array(z)]

        # remove gnames[1] file, so at the next iteration this one won't be there
        os.remove(gnames[1])

    # now extract from the geometry.log file the information about convergence criteria
    convValues, convTable = getConvCriteria("geometry.log")

    return convValues, convTable, newGeom


#####################################################################################################

def geometry_optxg_inp(geometry, command, cobcom, energy, gradient):
    # extract from the input directives the relevant sections
    redun = CBF.ReadCobramCommand(cobcom, 'redundant', 'redundant.dat')
    if command[1].strip() == 'irc':
        inputkey = CBF.ReadCobramCommand(cobcom, 'irc', 'irc.dat')
    elif command[1].strip() == 'ts':
        inputkey = CBF.ReadCobramCommand(cobcom, 'ts', 'ts.dat')
    elif command[1].strip() == 'ci':
        inputkey = CBF.ReadCobramCommand(cobcom, 'ci', 'ci.dat')
    elif command[1].strip() == 'freqxg':
        inputkey = CBF.ReadCobramCommand(cobcom, 'freqxg', 'freqxg.dat')
    else:
        inputkey = CBF.ReadCobramCommand(cobcom, 'optxg', 'optxg.dat')

    high = list(geometry.list_HIGH)
    medium = list(geometry.list_MEDIUM)
    at = []
    x = []
    y = []
    z = []
    # grep medium-high geometry
    for i in range(geometry.atomNum):
        if high.count(i + 1) == 1:
            x.append(geometry.cartesian[0][i])
            y.append(geometry.cartesian[1][i])
            z.append(geometry.cartesian[2][i])
            at.append(geometry.atomLabel[i])
        if medium.count(i + 1) == 1:
            x.append(geometry.cartesian[0][i])
            y.append(geometry.cartesian[1][i])
            z.append(geometry.cartesian[2][i])
            at.append(geometry.atomLabel[i])

    mem = command[5]  # allocated memory for geometry optimization
    Nproc = '1'

    if command[60].strip() == 'sp' or command[60].strip() == '1':
        maxoptstep = '1'
    else:
        maxoptstep = str(int(command[60].strip()) - 1)  # maximum number of opt step

    foundredundant = redundantfinder(geometry, command, redun)
    redundant = foundredundant[0]
    F = foundredundant[1]

    # define charge and multeplicity for the geometry calculation
    Chm, Mhm = getCM(at, x, y, z, mem)

    new_inputkey = []
    coboptkey = []
    for i in range(len(inputkey)):
        if inputkey[i].startswith('coboptkey') == 1:
            coboptkey.append(inputkey[i].split('=')[1].lower())
        else:
            new_inputkey.append(inputkey[i])
    inputkey = new_inputkey

    usesymm = 0
    if "symm" in coboptkey:
        usesymm = 1
        logwrt.writelog('The optimization will use symmetry: WARNING possible reorientation')
    # if the user specifies nosymm in the !optxg ?optxg (or similar) inputs Cobram must register this
    elif "nosymm" in inputkey:
        usesymm = 1

    geom = open('geometry.com', 'w')
    geom.write('%chk=geometry\n')
    geom.write('%int=geometry\n')
    geom.write('%d2e=geometry\n')
    geom.write('%rwf=geometry\n')
    geom.write('%Nproc=' + Nproc + '\n')
    geom.write('%mem=' + mem + '\n')
    if command[1] == 'freqxg' and len(inputkey) == 0:
        inputkey = ['#P freq=(numer,step=' + str(int(float(command[12]) / 0.0001)) + ',HPModes)']
    elif command[1] in ['optxg', 'ci'] and len(inputkey) == 0:
        inputkey = ['#p opt(noeigentest,nolinear,nomicro,maxstep=' + command[68] + ')']
    elif command[1] in ['ts'] and len(inputkey) == 0:
        inputkey = ['#p opt(ts,NewEstmFC,noeigentest,nolinear,nomicro,maxstep=' + command[68] + ') iop(1/9=1,5/6=6)']
    elif command[1] in ['irc'] and len(inputkey) == 0:
        inputkey = ['#p IRC=(FCCards,Readvector,MaxPoints=' + command[65] + ',StepSize=' + command[66] + ') Use=L115']

    for i in range(len(inputkey)):
        geom.write(inputkey[i] + '\n')

    if command[1] in ['bomdxg', 'fullbomdxg']:
        # steptraj = int(maxoptstep) - (geometry.NatomHM * 6 + 2)
        ntraj = str(1)
        geom.write('iop(1/6=' + ntraj + ')\n')
    else:
        geom.write('iop(1/6=' + maxoptstep + ')\n')
        if command[67] != '':
            geom.write('iop(1/7=' + command[67] + ')\n')

    # get name of gaussian executable (no need to check... see at the beginning of GAUSSIAN_optimizator_X!)
    gaussian_executable = os.getenv('GAUSSIAN_EXE')

    if gaussian_executable == 'g09' or gaussian_executable == 'g16':
        if usesymm == 0:
            geom.write('external=cobramext nosymm iop(3/5=30)\n')
        else:
            geom.write('external=cobramext iop(3/5=30)\n')
    elif gaussian_executable == 'g03':
        if usesymm == 0:
            geom.write('external nosymm\n')
        else:
            geom.write('external \n')

    geom.write('\n')
    geom.write('Optimization requested by COBRAM\n')
    geom.write('\n')
    geom.write(Chm + ' ' + Mhm + '\n')

    at2 = []
    for i in range(geometry.NatomHM):
        if at[i].lower() == 'ag':
            at2.append('cu')
        else:
            at2.append(at[i])

    for i in range(geometry.NatomHM):
        j = i + 1
        if j in F:
            label_opt = '-1'
        else:
            label_opt = '0'
        if command[1] == 'freqxg':
            label_opt = '0'
        geom.write(at2[i] + '\t' + '%2s' % label_opt + '%16.8f' % x[i] + '%16.8f' % y[i] + '%16.8f\n' % z[i])
    geom.write(' \n')

    # so far this has been adapted for following the gradient on the ES with a unit force matrix
    # modifications are required if the user wants to run IRC from a TS (can be done also within optxg)
    if command[1] == 'irc':
        # write energy
        geom.write("%24.16f\n" % energy)
        # write gradient
        for i in range(geometry.NatomHM):
            geom.write("{0:12.8f}{1:12.8f}{2:12.8f}".format(-gradient[0][i], -gradient[1][i], -gradient[2][i]))
            if i % 2 == 1:
                geom.write("\n")
        if geometry.NatomHM % 2 == 1:
            geom.write("\n")
        # create a 3*N x 3*N unit matrix of forces
        m = np.eye(3 * geometry.NatomHM)
        # write the matrix in the correct format to input (six values per line)
        k = 0
        for i in range(3 * geometry.NatomHM):
            for j in range(i + 1):
                geom.write("%12.8f" % (m[i][j]))
                k = k + 1
                if k == 6:
                    geom.write("\n")
                    k = 0
        # write an additional blank line
        if (((geometry.NatomHM * 3) * (geometry.NatomHM * 3 + 1)) / 2) % 6 != 0:
            geom.write("\n")
        geom.write("\n")
        # now we read the gradient one more time (this time it gives the direction to follow in the IRC)
        k = 0
        for i in range(geometry.NatomHM):
            for f in (-gradient[0][i], -gradient[1][i], -gradient[2][i]):
                geom.write("{0:10.6f}".format(f))
                k += 1
                if k == 8:
                    geom.write("\n")
                    k = 0

    # if command[3] == '3' or command[4] == '3':
    #     for i in range(len(redundant_A)):
    #         geom.write(redundant_A[i])

    if command[1] != 'freqxg':
        for i in range(len(redundant)):
            geom.write(str(redundant[i]) + '\n')
    geom.write(' \n')

    if len(redun) == 0:
        pass
        # logwrt.writelog('No redundant section was given in the cobram.command')
    else:
        logwrt.writelog('Following Coordinates will be frozen according to the ModRedundant section: \n')
        for i in range(len(redun)):
            logwrt.writelog(redun[i])
        logwrt.writelog('\n')
    if len(F) != 0:
        logwrt.writelog('Also the following atoms will be frozen \n')
        for i in range(len(F)):
            logwrt.writelog(str(F[i]) + '\n')
    else:
        pass
        # logwrt.writelog('No other atoms to be frozen ')
    geom.write(' \n')
    geom.close()
    logwrt.writelog('\n')


#####################################################################################################

def redundantfinder(geometry, command, redun):
    list_HIGH = geometry.list_HIGH
    list_MEDIUM = geometry.list_MEDIUM
    list_BA = geometry.atomLink_BA
    list_MEDIUM_HIGH = geometry.list_MEDIUM_HIGH
    calculation_type = geometry.calculationType
    list_DC = geometry.atomLink_DC
    A = list_BA[1]
    C = list_DC[1]

    mmsteps = command[74]
    k = mmsteps.split('/')
    mmstepstot = k[0]
    # if you optimize the mm part it makes no sense to freeze any border!
    if (command[3] != '0' or command[4] != '0') and int(mmstepstot) >= 0:
        command[3] = '0'
        command[4] = '0'
        logwrt.writelog(' part of the LOW layer is optimized command(74)=' + command[
            74] + ' (total step to perform = ' + mmstepstot + '):\n')
        logwrt.writelog(' thus it makes no sense to freeze the QM-MM border(s): all constraints are released!\n')
        # logwrt.writelog( 'setting command(3)=0' )
        # logwrt.writelog( 'setting command(4)=0' )
        # logwrt.writelog( '\n\n' )
    if int(mmstepstot) >= 0 and calculation_type in ['HL', 'HM', 'ML', 'HML', ]:
        command[3] = '0'
        command[4] = '0'
        # logwrt.writelog(' there is no LOW layer or the LOW layer is NOT optimized command(74)='+command[74]+' :')
        # logwrt.writelog(' thus it is a sensible choice to set the QM-MM border(s) free')
        # logwrt.writelog( 'setting command(3)=0' )
        # logwrt.writelog( 'setting command(4)=1' )
        # logwrt.writelog( '\n\n' )

    if command[4] == '1' or command[4] == '3':
        if command[4] == '1':
            # bord = 'F'
            # logwrt.writelog('command(4)=1 [default is 0]')
            logwrt.writelog('Atoms of outer border will be kept frozen\n')
        elif command[4] == '3':
            # bord = 'A'
            # logwrt.writelog('command(4)=3 [default is 0]')
            logwrt.writelog('Atoms of outer border will be active!\n')
    elif command[4] == '0':
        # logwrt.writelog('command(4)= '+str(command[4])+' [default is 0]')
        logwrt.writelog('Atoms of outer border will be free!\n')
    else:
        logwrt.fatalerror('command(4)= ' + str(command[4]) + ' is not valid')

    # find outer border atoms and assign the correct index
    list_border_in_HL = []
    list_border_in_HM = []
    list_border_in_ML = []
    list_border_in_HML = []
    list_border = []
    if calculation_type == 'HL':
        for i in range(len(list_HIGH)):
            if list_HIGH[i] in A:
                list_border_in_HL.append(i + 1)
        list_border = list_border_in_HL[:]
        # logwrt.writelog(' so the outer border is between the HIGH and the LOW layers')
    elif calculation_type == 'HM':
        for i in range(len(list_HIGH)):
            if list_HIGH[i] in A:
                list_border_in_HM.append(i + 1)
        list_border = list_border_in_HM[:]
        # logwrt.writelog(' so the outer border is between the HIGH and the MEDIUM layers')
    elif calculation_type == 'HML':
        for i in range(len(list_MEDIUM_HIGH)):
            if list_MEDIUM_HIGH[i] in C:
                list_border_in_HML.append(i + 1)
        list_border = list_border_in_HML[:]
        # logwrt.writelog(' so the outer border is between the MEDIUM and the LOW layers')
    elif calculation_type == 'ML':
        for i in range(len(list_MEDIUM)):
            if list_MEDIUM[i] in C:
                list_border_in_ML.append(i + 1)
        list_border = list_border_in_ML[:]
        # logwrt.writelog(' so the outer border is between the MEDIUM and the LOW layers')
    elif calculation_type == 'L':
        list_border = []
        # logwrt.writelog(' so there is no outer border')
    elif calculation_type == 'M':
        list_border = []
        # logwrt.writelog(' so there is no outer border')
    elif calculation_type == 'H':
        list_border = []
        # logwrt.writelog(' so there is no outer border')

    # create a list of redundant atom to keep frozen by label
    list_M_in_HM = []
    list_H_in_HM = []
    for i in range(len(list_MEDIUM_HIGH)):
        if list_MEDIUM_HIGH[i] in list_MEDIUM:
            list_M_in_HM.append(i + 1)
        elif list_MEDIUM_HIGH[i] in list_HIGH:
            list_H_in_HM.append(i + 1)

    bF, bA, bN = [], [], []
    mF, mA, mN = [], [], []
    if calculation_type == 'HML' or calculation_type == 'HM' or calculation_type == 'HL' or calculation_type == 'ML':
        if command[3] == '0':
            mF = []
            mA = []
            mN = list_H_in_HM + list_M_in_HM
        elif command[3] == '1':
            mF = list_M_in_HM[:]
            mA = []
            mN = list_H_in_HM[:]
        elif command[3] == '3':
            mF = []
            mA = list_M_in_HM[:]
            mN = list_H_in_HM[:]
        if command[4] == '0':
            bF = []
            bA = []
            bN = list_border[:]
        elif command[4] == '1':
            bF = list_border[:]
            bA = []
            bN = []
        elif command[4] == '3':
            bF = []
            bA = list_border[:]
            bN = []
    else:
        mF, mA, mN = [], [], []
        bF, bA, bN = [], [], []

    f_QM1, f_QM2, f_QM3, f_QM4 = [], [], [], []
    addredu = []
    for i in range(len(redun)):
        if 'F' or 'f' in redun[i].strip().split():
            if len(redun[i].strip().split()) == 2:
                f_QM1.append(int(redun[i].split()[-2]))  # 1 F
            else:
                f_QM2.append(redun[i])  # 1 2 F; 1 2 3 F; 1 2 3 F
        else:
            if len(redun[i].strip().split()) == 2:
                f_QM3.append(int(redun[i].split()[-2]))  # 1 F
            else:
                f_QM4.append(redun[i])  # 1 2 F; 1 2 3 F; 1 2 3 F
        addredu.append(redun[i])

    f, a, n = [], [], []
    for i in range(len(list_MEDIUM_HIGH)):
        j = i + 1
        if (j in mF or j in bF) and j not in f_QM1 and j not in f_QM3:
            f.append(j)
        if (j in mA or j in bA) and j not in f_QM1 and j not in f_QM3:
            a.append(j)
        if (j in mN or j in bN) and j not in f_QM1 and j not in f_QM3:
            n.append(j)
    for i in range(len(bF)):
        j = bF[i]
        if j in a:
            a.remove(j)
    F = []  # frozen atoms
    # A = []  # active atoms
    # N = []
    for i in range(len(f)):
        if f[i] not in a:
            F.append(f[i])
    A = a[:]
    # N = n[:]
    redundant_F = []
    redundant_A = []
    redundant_bF = []
    for i in range(len(F)):
        redundant_F.append(str(F[i]) + ' F\n')
    for i in range(len(A)):
        redundant_A.append(str(A[i]) + ' A\n')
    for i in range(len(bF)):
        redundant_bF.append(str(bF[i]) + ' F\n')
    redundant = addredu
    foundredundant = [redundant, F]
    return foundredundant


#####################################################################################################


def getCM(at, x, y, z, mem):
    # prepare an input file for a fake calculation to check charge and moltepl of the optimization run
    with open('pseudogeometry.com', 'w') as geom:
        geom.write('%chk=pseudogeometry\n%int=pseudogeometry\n%d2e=pseudogeometry\n%rwf=pseudogeometry\n')
        geom.write('%Nproc=1\n%mem=' + mem + '\n')
        geom.write('#p UFF iop(1/6=1)\n\n')
        geom.write('Optimization requested by COBRAM  to find the correct charge-multiplicity set\n\n')
        geom.write('0 1\n')
        for i in range(len(at)):
            label_opt = '0'
            if at[i].lower() == 'ag':
                geom.write('cu' + '\t' + '%2s' % label_opt + '%16.8f' % x[i] + '%16.8f' % y[i] + '%16.8f\n' % z[i])
            else:
                geom.write(at[i] + '\t' + '%2s' % label_opt + '%16.8f' % x[i] + '%16.8f' % y[i] + '%16.8f\n' % z[i])
        geom.write(' \n')

    # get from environment the variables for the GAUSSIAN program
    gaussian_exe = os.getenv('GAUSSIAN_EXE')
    gpath = os.getenv('GAUSSIAN_DIR')
    # define profile for running the version of gaussian defined by gaussian_exe and gpath
    cobrammenv.setGaussianProfile(gaussian_exe, gpath)

    # now run the pseudogeometry calculation with gaussian
    fin = open("pseudogeometry.com", "r")
    fout = open("pseudogeometry.log", "w")
    ferr = open("pseudogeometry.err", "w")
    subprocess.call(["nice", "-1", gaussian_exe], stdin=fin, stderr=ferr, stdout=fout)
    fin.close()
    fout.close()
    ferr.close()

    # check the termination of the pseudogeometry run
    gauTerm = checkGauTermination("pseudogeometry.log")

    # now check if pseudogeometry terminated normally or with error
    C, M = None, None
    if gauTerm == 1:  # in the case of normal termination, charge 0 and multeplicity 1 works fine
        C, M = '0', '1'
    elif gauTerm == 2:  # with error termination, we will pretend that the system is a duplet
        C, M = '0', '2'
    else:  # if no termination line is found, abort with error
        logwrt.fatalerror(" Error in charge and multeplicity assignment of geometry optimization... "
            "please inspect the pseudogeometry.log file ")

    # clean up and return charge and multeplicity values
    for fl in glob.glob("pseudogeometry.*"):
        os.remove(fl)
    return [C, M]


#####################################################################################################


def getConvCriteria(fileName):
    # get rms and max forces and displacement from gaussian optimization output

    # extract the text of the gaussian output file
    with open(fileName, "r") as f:
        text = f.read()

    try:
        # from the output text, find the last table with the convergence values
        convTable = re.findall("( *Item *Value *Threshold.*?Converged\?.*?)Predicted change", text, re.DOTALL)[-1]

        # now get Maximum/RMS Force Displacement from each line
        tableLines = convTable.split("\n")
        try:
            FMax = float(tableLines[1].split()[2])
        except ValueError:
            FMax = None
        try:
            FRMS = float(tableLines[2].split()[2])
        except ValueError:
            FRMS = None
        try:
            DisplMax = float(tableLines[3].split()[2])
        except ValueError:
            DisplMax = None
        try:
            DisplRMS = float(tableLines[4].split()[2])
        except ValueError:
            DisplRMS = None

        convTable = "\n" + convTable

    except:  # the file does not contain any convergence information
        convTable = """
         Item               Value     Threshold  Converged?
 Maximum Force               N.A.        N.A.       NO
 RMS     Force               N.A.        N.A.       NO
 Maximum Displacement        N.A.        N.A.       NO
 RMS     Displacement        N.A.        N.A.       NO
 """
        FMax = None
        FRMS = None
        DisplMax = None
        DisplRMS = None

    # return both the list of values and the convergence table string
    return [FMax, FRMS, DisplMax, DisplRMS], convTable


#####################################################################################################


def checkGauTermination(fileName):
    # extract the last ten lines of the gaussian output file
    with open(fileName, "r") as f:
        fileTail = " ".join(f.readlines()[-10:])

    # now check if gaussian terminated normally or with error
    if "Normal termination" in fileTail:  # in the case of normal termination, return status = 1
        return 1
    elif "Error termination" in fileTail:  # with error termination
        # return 3 if optimization has not converged
        with open(fileName, 'r') as f:
            logFile = f.read()
            if re.findall('Number of steps exceeded', logFile):
                return 3
        # for other errors return 2
        return 2
    else:  # if no termination line is found, return 0
        return 0


#####################################################################################################


def getNormalModes(fileName, high_precision=False):
    # read gaussian optimization file content
    # if high_precision is True, return dict (key: freq, value: list of normal modes)
    # if high precision is False, return string
    with open(fileName, "r") as f:
        fileText = f.read()

    # extract section on Harmonic frequencies with regular expression
    base_regex = '( Harmonic frequencies.*?)\n'
    if high_precision:
        regex = base_regex + ' Harmonic frequencies'
    else:
        regex = '\.\d{5}\n' + base_regex + ' -------------------\n - Thermochemistry'

    normModesString = re.findall(regex, fileText, re.DOTALL)

    try:
        if high_precision:
            freq_blobs = normModesString[0].split('Frequencies ---')[1:]
            freq = []
            NM = []
            for blob in freq_blobs:
                _freq = re.findall(r'-?\d+\.\d{4}', blob.split('\n')[0])
                nm = [[] for i in range(len(_freq))]
                freq += _freq
                lines = blob.split('Coord Atom Element:')[1].split('\n')[1:]
                for line in lines:
                    linenm = re.findall(r'-?\d+\.\d{5}', line)
                    for i in range(len(linenm)):
                        nm[i].append(float(linenm[i]))
                NM += nm

            normModes = {}
            for f, n in zip(freq, NM):
                normModes[float(f)] = n

        else:
            normModes = normModesString[0]

        return normModes
    except:
        return None

#####################################################################################################

