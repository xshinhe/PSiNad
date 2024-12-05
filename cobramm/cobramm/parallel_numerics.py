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

import time
import sys
import os
import stat
import shutil
import subprocess
import shelve

# imports of local modules

import CBF
import logwrt
import inpdata
import constants

# import of local classes

from QM import QM  # QM class controls the QM calculation, define its input and stores the output
from timer import Timer  # keep timings of the different sections of the code

# math libraries

import numpy as np  # numpy library for scientific computation


#####################################################################################################

@Timer("parall. freq")
def newFreqParallelRun(cobcom, command, geometry, charges, equilCalc, setState=None, forceState=None):
    """ Driver function for a parallel calculations of the finite difference calculations that
    are necessary for a frequency calculation. This driver works with the new QM interface
    defined by the QM class.
    This functions executes all the required SP calculations, and then returns the results
    as QM class instances that are subsequently read during the main COBRAMM loop. """

    # initial message
    logwrt.writelog("Entering parallel execution of the finite differences for the frequency computation.\n")

    # define dictionary that stores the option for atom displacements
    # define "length" and "linkfollow" options according to command[12] and commmand[16]
    Displ = {"iAt": None,
             "iCoord": None,
             "iDir": None,
             "length": float(command[12]),
             "linkfollow": True if command[16] == "0" else False}

    QMcalculations = []  # list to store the QM objects for the calculations with displaced geometries
    QMtorun = []  # list to store the QM objects of the calculations that needs to be run

    # the first element of the calculation list should be the equilibrium calculation, which is repeated in step 1
    QMcalculations.append(equilCalc)

    # loop over atoms, atom coordinates and +/- directions
    for iAtHM in range(geometry.NatomHM):
        for iCoord in range(3):
            for iDir in [+1, -1]:

                # get identifier of the atom
                idAt = geometry.list_MEDIUM_HIGH[iAtHM]

                # the current atom is one of the High layer ones
                if idAt in geometry.list_HIGH:

                    # get the corresponding index of the atom in the HIGH LAYER list
                    iAtH, = np.where(geometry.list_HIGH == idAt)
                    # store the values in the displacement dictionary
                    Displ["iAt"], Displ["iCoord"], Displ["iDir"] = iAtH, iCoord, iDir
                    # construct the instance of the QM object
                    QMcalc = QM(cobcom, command, geometry, charges, 1, equilCalc.saveRestartFile(), displacement=Displ, setState=setState, forceState=forceState)
                    # store the instance in the full list of QM calc, and in the list of calc.s to be run
                    QMcalculations.append(QMcalc), QMtorun.append(QMcalc)

                # the current atom is one of the High layer ones
                elif idAt in geometry.list_MEDIUM:

                    # just append the equilibrium QM calculation in the full list of calculations
                    QMcalculations.append(equilCalc)

    # the number of available processors for parallel execution is defined by command[9]
    nproc, memory, nprocQM = int(command[9]), command[53], int(command[7])
    # message about parallel execution parameters: number of calculations to execute and number of processors
    logwrt.writelog("Running {0} calculations distibuted on {1} available cores\n".format(len(QMtorun), nproc))
    logwrt.writelog('Avg number of SP computed sequentially on each core is between {0} and {1}\n'.format(
        len(QMtorun) // nproc, len(QMtorun) // nproc + 1))

    # print message with time at start of this function
    logwrt.writelog("Start numerical computations at: "
                    "{0}\n".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())))

    # run the calculation with the specific QM function
    QM.runQM(QMtorun, memory=memory, nprocQM=nprocQM, nprocParall=nproc)

    logwrt.writelog("Finish numerical computations at: "
                    "{0}\n\n".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())))

    # and now return the full list of calculations
    return QMcalculations

#####################################################################################################


@Timer("parall. deriv")
def newGradientParallerRun(cobcom, command, geometry, charges, equilCalc):
    """ Driver function for a parallel calculations of the gradient with finite difference calculations.
    This driver works with the new QM interface defined by the QM class.
    This functions executes all the required SP calculations, and then returns the results
    by updating equilCalc QM class instance that is passes as input. """

    # initial message
    logwrt.writelog("Entering parallel execution of the finite differences for the gradient computation.\n")

    # define dictionary that stores the option for atom displacements
    # define "length" and "linkfollow" options according to command[12] and commmand[16]
    Displ = {"iAt": None,
             "iCoord": None,
             "iDir": None,
             "length": float(command[12]),
             "linkfollow": True if command[16] == "0" else False}

    QMcalculations = []  # list to store the QM objects for the calculations with displaced geometries
    QMtorun = []  # list to store the QM objects of the calculations that needs to be run

    # print a message to inform about the behaviour with respect to H atoms of atom links
    if geometry.NsubH != 0:
        if command[16] == '0':
            print("Numerical differentiation over link atoms is NOT performed. " +
                  "They will be displaced together with the neighboring QM atom.")
        else:
            print("Including link atoms in numerical differentiation.")
            print("Number of link atoms "+str(geometry.NsubH)+"\n")

    # the final dimension of the gradient will be in any case equal to the size of modelH
    nAtoms = geometry.NatomQM + geometry.NsubH

    # define the list of directions to consider, and for the future define the denominator of the fin diff formula
    if command[10] == '0':
        directList = [+1]
        denominator = 1.0 * float(command[12]) / constants.Bohr2Ang
    else:
        directList = [+1, -1]
        denominator = 2.0 * float(command[12]) / constants.Bohr2Ang

    # store with a list the atom, coordinate and direction definition
    displDefs = []

    # loop over atoms, atom coordinates and +/- directions
    if int(command[210]) > 0:
        nnAtoms=int(command[210])
        irange=1
    else:
        nnAtoms=nAtoms
        irange=3
    for iAt in range(nnAtoms):
        for iCoord in range(irange):
            for iDir in directList:

                # this is a High layer atom atom, or the finite difference is computed for all the atoms
                if iAt < geometry.NatomQM or command[16] != "0":

                    # store the values in the displacement dictionary
                    Displ["iAt"], Displ["iCoord"], Displ["iDir"] = iAt, iCoord, iDir
                    # construct the instance of the QM object
                    QMcalc = QM(cobcom, command, geometry, charges, equilCalc.saveRestartFile(), displacement=Displ)
                    # store the instance in the full list of QM calc, and in the list of calc.s to be run
                    QMcalculations.append(QMcalc), QMtorun.append(QMcalc)

                # this is the H of an atom link, and according to command[16] we should skip this derivative
                else:

                    # store the equilibium geometry in the full calculation list
                    QMcalculations.append(equilCalc)

                # store the atom, coordinate and direction in a list to be used later to reference the QM calc
                displDefs.append([iAt, iCoord, iDir])

    # the number of available processors for parallel execution is defined by command[9]
    nproc, memory, nprocQM = int(command[9]), command[53], int(command[7])
    # message about parallel execution parameters: number of calculations to execute and number of processors
    logwrt.writelog("Running {0} calculations distibuted on {1} available cores\n".format(len(QMtorun), nproc))
    logwrt.writelog('Avg number of SP computed sequentially on each core is between {0} and {1}\n'.format(
        len(QMtorun) // nproc, len(QMtorun) // nproc + 1))

    # print message with time at start of this function
    logwrt.writelog("Start numerical computations at: "
                    "{0}\n".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())))

    # run the calculation with the specific QM function
    QM.runQM(QMtorun, memory=memory, nprocQM=nprocQM, nprocParall=nproc)

    # first prepare a dictionary that stores the forward and backward energy for each electronic state
    forward, backward = {}, {}
    for state, energy in equilCalc.outputData.get("energy").items():
        forward[state] = np.empty((3, nAtoms))
        backward[state] = np.empty((3, nAtoms))
        # when command[10] is "0", we have only a forward difference, thus the backward value is the equilibrium E
        if command[10] == '0':
            backward[state].fill(energy)

    # loop over the computed QM displaced calculations and fill backward
    # also forward is filled when command[10] != '0', since QM calcs are present for iDir == -1
    for qm, disp in zip(QMcalculations, displDefs):
        iAt, iCoord, iDir = disp
        for state, energy in qm.outputData.get("energy").items():
            if iDir == +1:
                forward[state][iCoord, iAt] = energy
            else:
                backward[state][iCoord, iAt] = energy

    # now we can construct a dictionary with the gradients for all the states for which the energy is available
    gradients = {}
    for state in equilCalc.outputData.get("energy"):
        gradients[state] = (forward[state] - backward[state]) / denominator

    # and finally save the gradient in the equilCalc object
    equilCalc.outputData.set("gradient", gradients)

    logwrt.writelog("Finish numerical computations at: "
                    "{0}\n\n".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())))


#####################################################################################################

@Timer("old parall.")
def run(QM_Results, step, geometry, command, cobcom, charges, par_num):
    """ Driver function of the parallel calculations in which the finite differentiation
        is parallelized by distributing the single finite difference calculations
        on the different available nodes """

    ##########################################################
    # parallel numerical FREQs
    ##########################################################

    # all SP calculations are performed during step 0 and stored in directories geomXXX
    # in the subsequent steps if a QM computation should be done, Cobram just fetches the
    # QM .log file (molcas.log, gaussian-QM.log, etc.) from the geomXXX directory to the
    # main directory and continues normally

    if command[1] == 'freqxg':

        # perform the parallel numerical computation, par_num == 1 at this point

        logwrt.writelog("Entering parallel execution of the finite differences for the frequency computation.\n")

        # the number of available processors is defined by command[9]
        nProcs = int(command[9])

        tot_disp = 6
        nRepeats = tot_disp * geometry.NatomQM // nProcs
        logwrt.writelog('Number of SP computed sequentially on each core: ' +
                        'between {0} and {1}\n'.format(nRepeats, nRepeats + 1))

        # initialize a multi-dimensional array, where the dimension is equal to the number of requested processors
        list_steps = []
        for iProc in range(nProcs):
            list_steps.append([])

        # this is the id label of the processor that is used to distribute the calculations
        iProc = 0
        # while freqstep runs over the steps of the frequency calculation that involves displacements
        # of all the H and M atoms, we need another counter for the H atoms alone that labels the QM calculations
        QMstep = 1

        # loop over all the displacements to consider
        for freqstep in range(1, tot_disp * geometry.NatomHM + 1):

            # an iterator over the requested processors
            if iProc == nProcs: iProc = 0

            # skip this value of QMstep if the atom we are dealing with at this step is in the M layer
            if freq_QM_skipper(freqstep, geometry, command): continue

            # if the directory "geom"+str(QMstep) already exists AND there is a file SUCCESS therein,
            # this means that the SP is already computed and will be skipped (this is for restarts)
            stepPath = "geom" + str(freqstep)
            if os.path.isdir(stepPath):
                if os.path.exists(os.path.join(stepPath, "SUCCESS")):  # allow for restart
                    logwrt.writewarning("SP geom{0} is already present: skipping this QM calculation.".format(freqstep))
                    continue
                else:  # clean up previous run
                    shutil.rmtree(stepPath)
            # create a directory geomX with X = QMstep
            os.mkdir(stepPath)

            # distribute the SP over the requested processors in a round-robin way
            list_steps[iProc].append(freqstep)

            # create the input files BUT DO NOT launch yet
            # it should be determined which atom is to be moved according to 'QMstep'
            CBF.QM(command, cobcom, charges, geometry, QMstep, par_num)

            # copy input to geomX
            if command[51] == '6':
                shutil.move("seward.input", stepPath)
            else:
                logwrt.fatalerror('Parallel finite-difference frequencies are not implemented yet for {0}'.format(
                    inpdata.getQMCode(command)))

            # assign next calculation to another processor
            iProc += 1
            # move to the next QM calculation
            QMstep += 1

        logwrt.writelog('\n')

        # the actual parallel computation
        QM2(command, list_steps, nProcs)

        # par_num is set to 2, so in the next execution of CBF.QM, the energies and gradient are read from the
        # calculations that have been already done
        par_num = 2

    ##########################################################
    # parallel numerical GRADs & NACs
    ##########################################################

    elif command[1] in ['optxg', 'mdv', 'irc', 'ci', 'ts']:

        # because of issues with the MOLCAS functions for input writing,
        # temporarily restore the command[1] to the version with "p"
        command[1] += "p"
        # freqxgp performs the parallel numerical computation, par_num == 1 at this point
        finite_displ_calcs(geometry, command, cobcom, charges, par_num)
        # return to the standard command without "p" at the end
        command[1] = command[1][:-1]

        # par_num is set to 2 so that in the execution of CBF.QM in the following steps only molcas.molcasEne
        # (and possible molcas.molcasNAC) is executed but not molcas.launch
        par_num = 2
        # GRADs (and possibly NACs) are computed by processing the indivual SP computations
        # the returned tmp_Results is used to update QM_Results
        Gradient, NAC = CBF.QM(command, cobcom, charges, geometry, step, par_num)

        # print warning when computed gradient is too large
        if int(command[210]) == 0:
            for iatom in range(geometry.NatomQM + geometry.NsubH):
                if abs(Gradient[0][iatom]) > 1.0 or abs(Gradient[1][iatom]) > 1.0 or abs(Gradient[2][iatom]) > 1.0:
                    logwrt.writewarning('Large gradient during numerical computation detected. Possible problems!')
                    break

        # when the calculation is a CoIn optimization, then check also the gradient on the lower state
        if command[1] == 'ci':
            sef = shelve.open("cobram-sef")
            gradient2 = sef['gradient2']
            sef.close()
            for iatom in range(geometry.NatomQM + geometry.NsubH):
                if abs(gradient2[0][iatom]) > 1.0 or abs(gradient2[1][iatom]) > 1.0 or abs(gradient2[2][iatom]) > 1.0:
                    logwrt.writewarning('Large gradient during numerical computation detected. Possible problems!')
                    break

        # store the gradient in QM_Results
        QM_Results[1] = Gradient

        if command[1] == 'mdv':
            # par_num is set to 3 so that only the MDV routine can be executed
            par_num = 3
            # the molcas.molcasMDV routine is called once the GRADs and NACs have been computed
            tmp_Results = CBF.QM(command, cobcom, charges, geometry, step, par_num)
            if tmp_Results != 0:
                # QM_Results points to tmp_Results
                QM_Results = tmp_Results

        par_num = 1

    return QM_Results, par_num


def finite_displ_calcs(geometry, command, cobcom, charges, par_num):
    """ This function is entered in a parallel environment to compute the first derivative by finite differences
    in the case of a optxg, mdv, irc, ci or ts calculation"""

    # the number of runs is computed according to key 9 and the number of atoms
    nProcs = int(command[9])

    # define the number of atoms to displace
    if geometry.NsubH == 0:
        nAtoms, HAtoms = geometry.NatomQM, geometry.NatomQM
    else:
        if command[16] == '0':
            print("Numerical differentiation over link atoms (if present) is NOT performed. " +
                  "They will be displaced together with the neighboring QM atom.")
            nAtoms, HAtoms = geometry.NatomQM, geometry.NatomQM
        else:
            print("Including link atoms in numerical differentiation.")
            nAtoms, HAtoms = geometry.NatomQM + geometry.NsubH, geometry.NatomQM
            print("Number of link atoms "+str(geometry.NsubH)+"\n")

    # a possibility to displace in +/- (key 10 .ne 0) or + direction (key 10 .eq 0)
    if command[10] == '0':
        tot_disp = 3
    else:
        tot_disp = 6
    nRepeats = tot_disp * nAtoms // nProcs
    print('Number of SP computed sequentially on each core: between ' + str(nRepeats) + ' and ' + str(nRepeats + 1))

    # initialize a multi-dimensional array, where the dimension is equal to the number of requested processors
    list_steps = []
    for iProc in range(nProcs):
        list_steps.append([])

    # initialize counters
    iProc = 0

    # loop over all the displacements
    if int(command[210]) > 0 :
        irange = 2*int(command[210])
    else:
        irange = tot_disp * nAtoms
    for step in range(1, irange+1):

        if iProc == nProcs: iProc = 0

        # if the directory geomX already exists AND there is a file SUCCESS therein, this means that the
        # SP is already computed and will be skipped (this is for restarts)
        stepPath = "geom" + str(step)
        if os.path.isdir(stepPath):
            if os.path.exists(os.path.join(stepPath, "SUCCESS")):  # allow for restart
                logwrt.writewarning("SP geom{0} is already present: skipping this QM calculation.".format(step))
                continue
            else:  # clean up previous run
                shutil.rmtree(stepPath)
        # create a directory geomX with X = QMstep
        os.mkdir(stepPath)

        # distribute the SP over the requested processors in a round-robin way
        list_steps[iProc].append(step)

        # create the input files BUT DO NOT launch yet
        # it should be determined which atom is to be moved according to 'step'
        CBF.QM(command, cobcom, charges, geometry, step, par_num)
        # copy input to geomX
        if command[51] == '6':
            shutil.move("seward.input", stepPath)
        else:
            logwrt.fatalerror('Parallel finite-difference derivatives are not implemented for {0}'.format(
                inpdata.getQMCode(command)))

        iProc += 1

    # the actual parallel computation
    QM2(command, list_steps, nProcs)


def QM2(command, list_steps, nProcs):
    """ in this routine an external script 'freqxt' (system independent) or 'freqxtBO' (infrastructure dependent)
    creates files submitter files 'run' in each directory geomX, the submitter script is then executed """

    # print message with time at start of this function
    logwrt.writelog("Start numerical computations at: "
                    "{0}\n".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())))

    logwrt.writelog("Running single-point QM calculations in parallel with {0}\n".format(inpdata.getQMCode(command)))

    for i in range(nProcs):
        if list_steps[i]:  # this if-statement accounts for the fact that there might be less jobs than processors

            # the submitter file for each set of SP jobs will be stored in the first folder of each sub-list
            sub_dir = list_steps[i][0]

            # out of each sub-list create a string with the geometry numbers
            ilist = " ".join([str(igeom) for igeom in list_steps[i]])
            logwrt.writelog("List of geometries @ node {0} : {1}\n".format(i, ilist))

            # create submitter "run" script, depending on the type of QM calculation
            commandline = []
            if command[51] == '1':  # QM calulation by Gaussian
                commandline = ["freqext", command[130], os.getcwd(), "gaussian", command[53], str(sub_dir), ilist]
            elif command[51] == '6':  # QM calculation by MOLCAS
                commandline = ["freqext", command[130], os.getcwd(), "molcas", command[53], str(sub_dir), ilist]
            elif command[51] == '7':  # QM calculation by MOLPRO
                commandline = ["freqext", command[130], os.getcwd(), "molpro", command[53], str(sub_dir), ilist]
            else:
                logwrt.fatalerror('Parallel submitter is not yet implemented ' +
                                  'for {0}'.format(inpdata.getQMCode(command)))
            subprocess.call(commandline)

            # name and path of the run script
            runScript = os.path.join("geom" + str(sub_dir), "run")

            # give exec permission to the run script
            st = os.stat(runScript)
            os.chmod(runScript, st.st_mode | stat.S_IEXEC)

            # python executes each submitter 'run' BUT DOES NOT wait for it to finish
            subprocess.Popen(["nice", "-2", runScript], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # flush standard output
    sys.stdout.flush()

    # now wait that all submitter script have finished and SUCCESS or FAIL files are present in the SP directories
    for i in range(nProcs):
        for j in list_steps[i]:

            # name of the SUCCESS and FAIL files
            successFile = os.path.join("geom" + str(j), "SUCCESS")
            failFile = os.path.join("geom" + str(j), "FAIL")

            # wait until calculation j of node i has finisced
            while True:
                time.sleep(1.0)

                # if the SP was a success the computation continues
                if os.path.exists(successFile):
                    break
                # if the SP has failed, stop COBRAMM execution
                elif os.path.exists(failFile):
                    logwrt.fatalerror(
                        'Finish step {0} @ node {1} unsuccessfully. Cobram is going to die now.'.format(j, i))

    logwrt.writelog(
        "Finish numerical computations at: {0}\n\n".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())))


def freq_QM_skipper(step, geometry, command):
    """ this function returns true when the QM calculation for a given atom displacement can be skipped """

    # the step 0 is the equilibrium configuration, it should NOT be skipped!
    if step == 0: return False

    # define atom and coordinate to displace
    if command[1] in ['optxg', 'mdv', 'irc', 'ci', 'ts'] and command[10] == '0':
        iAt = (step - 1) // 3
    else:
        iAt = (step - 1) // 6

    # decide whether to skip the QM calculation at the given step:
    if geometry.list_MEDIUM_HIGH[iAt] in geometry.list_HIGH:
        return False
    else:  # return True (=skip step) only when the atom is in the M layer
        return True
