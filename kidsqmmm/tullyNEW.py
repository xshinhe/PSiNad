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


from math import *
import numpy as np
import math
import os
from os import system
import CBF
import random
import product
import shelve
import copy
import logwrt
import constants
from QMCalc import QM
from typing import Union


class Tully:

    def __init__(self, QMPrev: Union[QM,None], QMNow: Union[QM,None], actualstate: int, actualtime: float, command, cobcom, amplitudes, geometry, step):
        self.actualstate = actualstate
        self.newstate = actualstate
        self.NACPrev = None
        self.NACNow = None
        self.TDCPrev = None
        self.TDCNow = None
        self.DEarrayNow = None
        self.DEarrayPrev = None
        self.tstepNow = None
        self.tstepPrev = None
        self.HT = None
        self.H0 = None
        self.DV0 = None
        self.DVT = None
        self.AM = amplitudes
        self.SItime = actualtime
        self.geometry = geometry
        self.P = None
        self.rnum = None
        # velocity attribute will be a dictonary of arrays: each item is the velocity at t, t-dt/2, t-dt and t-dt*3/2
        self.velocity = {}



        # execute only if we have both instances for current and previous timestep (i.e., if step > 0)
        if step == 0:
            #if step == 0 we only need to load energy gaps and determine the tstep

            #initialize required variables and set attributes
            nroots = QMNow.outputData.get('nroots')
            self.DEarrayNow = self.buildDEarray(QMNow, nroots)
    
            # in case of deltaE hop, we always use long timestep
            if command[85] == '2':
                logwrt.writelog("Time step is set to {0} in this step\n".format(command[83]))
                ttry = 'long'
                QMNow.setTstep(ttry)
                system('echo ' + str(ttry) + '>TSTEP')
            # in case of TDC, DV0 and DVT are directly initialized with TDC at previous and current step
            elif command[14] == '1':
                # tstep for next step is set and saved in ttry variable
                # for compatibility with rest of the code (old VVerlet modules) the ttry string is saved into TSTEP file
                logwrt.writelog("Time step is set to {0} in this step\n".format(command[84]))
                ttry = 'short'
                QMNow.setTstep(ttry)
                system('echo ' + str(ttry) + '>TSTEP')
            # mdv using NACs
            else:
                # the energy difference between state and state+1 and state-1, respectively, is fetched
                # by comparison to keyword/key ediff (86) the time step (short or long) is dtermined
                # NOTE: also with TDNACs it is possbile to have short and long steps!
                if self.actualstate != 0:
                    #if actual state id GS, we cannot fetch deltaE with lower state
                    DE1 = self.DEarrayNow[self.actualstate][self.actualstate - 1]
                else:
                    DE1 = float(command[86])
                if int(command[81]) == 0:
                    highest_root = nroots - 1
                else:
                    highest_root = int(command[81]) - 1
                if self.actualstate < highest_root:
                    #if actual state id the highest root, we cannot fetch deltaE with higher state
                    DE2 = self.DEarrayNow[self.actualstate][self.actualstate + 1]
                else:
                    DE2 = float(command[86])
    
                # type of tstep for next step is determined based on deltaE and saved in ttry variable
                # for compatibility with rest of the code (old VVerlet modules) the ttry string is saved into TSTEP file
                if abs(DE1) < float(command[86]) or abs(DE2) < float(command[86]):
                    logwrt.writelog("Time step is set to {0} in this step (SHORT)\n".format(command[84]))
                    ttry = 'short'
                else:
                    logwrt.writelog("Time step is set to {0} in this step (LONG)\n".format(command[83]))
                    ttry = 'long'
                QMNow.setTstep(ttry)
                system('echo ' + str(ttry) + '>TSTEP')
        #step > 0
        else:
            logwrt.writelog("\n")
            logwrt.startSubSection("TULLY'S FSSH")
    
            #initialize required variables and set attributes
            nroots = QMNow.outputData.get('nroots')
            self.DEarrayNow = self.buildDEarray(QMNow, nroots)
            #variable comp_state is initialized to None for compatibility with later use of determineCopulingScheme module
            comp_state = None
            # following additional attributes are needed only for hop based on DC (either TDC or NAC)
            if command[85] != '2':
                zeromat = np.zeros((nroots, nroots))
                self.P = copy.deepcopy(zeromat)
                self.HT = np.zeros(zeromat.shape, dtype=np.complex128)
                self.H0 = np.zeros(zeromat.shape, dtype=np.complex128)
                self.DV0 = np.zeros(zeromat.shape, dtype=np.complex128)
                self.DVT = np.zeros(zeromat.shape, dtype=np.complex128)
                
                self.NACNow = QMNow.outputData.get("nac")
                self.TDCNow = QMNow.outputData.get("tdc")
    
                self.DEarrayPrev = self.buildDEarray(QMPrev, nroots)
                self.NACPrev = QMPrev.outputData.get("nac")
                self.TDCPrev = QMPrev.outputData.get("tdc")
                self.tstepPrev = QMPrev.tstep
    
                # following section contains the initialization of derivative couplings (DV0 and DVT) and the determination of
                # time step for next step (by comparison of deltaE with command[86])...these are also to be performed only if hop is NOT based on deltaE
                if command[14] == '1':
                    # in case of TDC: DV0 and DVT (i.e. velocity x NAC projection at previous and current step) are directly initialized with TDC at previous and current step
                    self.DVT = np.array(self.TDCNow)
                    self.DV0 = np.array(self.TDCPrev)
                    # tstep for next step is set and saved in ttry variable
                    # for compatibility with rest of the code (old VVerlet modules) the ttry string is saved into TSTEP file
                    logwrt.writelog("Time step is set to {0} in this step\n".format(command[84]))
                    tstep = float(command[83]) * 41.341373337
                    ttry = 'short'
                    QMNow.setTstep(ttry)
                    self.tstepNow = ttry
                    system('echo ' + str(ttry) + '>TSTEP')
                    comp_state = None
                # mdv using NACs
                else:
                    # the energy difference between state and state+1 and state-1, respectively, is fetched
                    # by comparison to keyword/key ediff (86) the time step (short or long) is dtermined
                    # NOTE: also with TDNACs it is possbile to have short and long steps!
                    if self.actualstate != 0:
                        #if actual state is GS, we cannot fetch deltaE with lower state
                        DE1 = self.DEarrayNow[self.actualstate][self.actualstate - 1]
                    else:
                        DE1 = float(command[86])
                    if int(command[81]) == 0:
                        highest_root = nroots - 1
                    else:
                        highest_root = int(command[81]) - 1
                    if self.actualstate < highest_root:
                        #if actual state is the highest root, we cannot fetch deltaE with higher state
                        DE2 = self.DEarrayNow[self.actualstate][self.actualstate + 1]
                    else:
                        DE2 = float(command[86])
        
                    # type of tstep for next step is determined based on deltaE and saved in ttry variable
                    # for compatibility with rest of the code (old VVerlet modules) the ttry string is saved into TSTEP file
                    if abs(DE1) < float(command[86]) or abs(DE2) < float(command[86]):
                        logwrt.writelog("Time step is set to {0} in this step (SHORT)\n".format(command[84]))
                        tstep = float(command[84]) * 41.341373337
                        ttry = 'short'
                    else:
                        logwrt.writelog("Time step is set to {0} in this step (LONG)\n".format(command[83]))
                        tstep = float(command[83]) * 41.341373337
                        ttry = 'long'
                    QMNow.setTstep(ttry)
                    self.tstepNow = ttry
                    system('echo ' + str(ttry) + '>TSTEP')
                    #set variavle comp_state to closest state (needed?)
                    if abs(DE1) <= abs(DE2):
                        comp_state = self.actualstate - 1
                    else:
                        comp_state = self.actualstate + 1

            TULLY, THS = self.determineCopulingScheme(command, step, comp_state)
    
            if TULLY or THS:
                logwrt.writelog("\nEntering Tully's FSSH algorithm\n")

                if command[85] == '2':
                    #i.e., hopping based on deltaE
                    hop = False
                    for i in range(nroots):
                        if i != self.actualstate:
                            if not hop:
                                if self.DEarrayNow[self.actualstate][i] < float(command[207]):   
                                    logwrt.writelog('-----------------------------\n')
                                    logwrt.writelog('  !!!!Gimme Hop Joanna!!!\n')
                                    logwrt.writelog('-----------------------------\n')
                                    logwrt.writelog('hopping from state ' + str(self.actualstate + 1) + ' to --> ' + str(i + 1) + "\n")
                                    logwrt.writelog('hop ONLY based on energy difference\n')
                                    logwrt.writelog('-----------------------------\n')
                                    self.newstate = i
                                    hop = True
                                    #creation of HOP file is left for compatibility with old style VVerlet
                                    system('touch HOP')
                                    #CHECK FOR BACK-HOP
                                    if i > self.actualstate and command[87] == '0':
                                        logwrt.writelog("HOP rejected because backhop is not allowed (key 87)!\n")
                                        self.newstate = self.actualstate
                                        hop = False
                                        #again, HOP file is still needed by VVerlet at present
                                        system('rm HOP')
                else:
                    #i.e. hop using TDC or NAC
                    #
                    # LOADING VELOCITIES
                    # initialize empty velocity array at t-dt (xvel0, yvel0, zvel0) and t (xvel, yvel, zvel) ...
                    # for now I leave it like this (read/write velocities in files) for compatibility with old-style VVerlet
                    # in the future, it is better to add velocities as attributes of calculation (QMNow and QMPrev), like other quantities
                    xvel, yvel, zvel = [], [], []
                    xvel0, yvel0, zvel0 = [], [], []
                    # ... and load velocities from velocity.dat and velocityOLD.dat
                    velinp = open('velocity.dat')
                    vel = velinp.read().split('\n')
                    velinp.close()
                    velinp = open('velocityOLD.dat')
                    vel0 = velinp.read().split('\n')
                    velinp.close()
    
                    # assign velocities of QM atoms
                    for i in range(len(self.geometry.list_MEDIUM_HIGH)):
                        if self.geometry.list_MEDIUM_HIGH[i] in self.geometry.list_HIGH:
                            el = vel[i].split()
                            xvel.append(float(el[0]))
                            yvel.append(float(el[1]))
                            zvel.append(float(el[2]))
                            el = vel0[i].split()
                            xvel0.append(float(el[0]))
                            yvel0.append(float(el[1]))
                            zvel0.append(float(el[2]))
                    # assign velocities of atom-links (zeros!)
                    for _ in range(self.geometry.NsubH):
                            xvel.append(0.0)
                            yvel.append(0.0)
                            zvel.append(0.0)
                            xvel0.append(0.0)
                            yvel0.append(0.0)
                            zvel0.append(0.0)
                    xvel = np.array(xvel)
                    yvel = np.array(yvel)
                    zvel = np.array(zvel)
                    xvel0 = np.array(xvel0)
                    yvel0 = np.array(yvel0)
                    zvel0 = np.array(zvel0)
                    self.velocity['t-dt/2'] = np.array([xvel, yvel, zvel])
                    self.velocity['t-dt*3/2'] = np.array([xvel0, yvel0, zvel0])
    
                    if command[14] != '1':
                        if int(command[81]) == 0:
                            dynroots = nroots
                        else:
                            dynroots = int(command[81])
                        self.DVT, self.DV0 = self.computeDV(command, xvel, yvel, zvel, xvel0, yvel0, zvel0, dynroots, QMNow, QMPrev)
            
                        # print the time-derivative couplings matrices
                        logwrt.writelog("\nTime-derivative couplings at this step\n", 1)
                        logwrt.writelog(logwrt.matrix_prettystring(self.DVT, ".8f"), 1)
                        logwrt.writelog("\n", 1)
                        logwrt.writelog("\nTime-derivative couplings at previous step\n", 1)
                        logwrt.writelog(logwrt.matrix_prettystring(self.DV0, ".8f"), 1)
                        logwrt.writelog("\n", 1)
    
                    # following lines are executed by both Tully and THS
                    # LOADING AMPLITUDES
                    # initialize amplitudes every time we enter in Tully after a LONG time time step (?)
                    if not list(self.AM):
                        AM = list(self.AM)
                        logwrt.writelog('Initializing Amplitudes\n')
                        for i in range(nroots):
                            self.AM.append(complex(0.0))
                        self.AM[self.actualstate] = complex(1.0)
                        for i in range(nroots):
                            logwrt.writelog("state {0} = {1:8.6f} + {2:8.6f}*i\n".format(i+1, self.AM[i].real, self.AM[i].imag))
                        logwrt.writelog("\n")
                        self.AM = np.array(self.AM)
                        AM1 = copy.deepcopy(self.AM)
                    else:
                        # otherwise Amplitudes were loaded from the shelve
                        AM1 = copy.deepcopy(self.AM)
    
                    # CREATE HAMILTONIAN AT TIME t-dt (previous step, H0) and t (present step, HT)
                    for i in range(nroots):
                        for j in range(i, nroots):
                            if i == j:
                                self.HT[i][j] = complex(self.DEarrayNow[0][i] / 627.51)
                                self.H0[i][j] = complex(self.DEarrayPrev[0][i] / 627.51)
                            else:
                                self.HT[i][j] = complex(self.DVT[i][j] * 1j)
                                self.H0[i][j] = complex(self.DV0[i][j] * 1j)
                                if command[14] != '1':
                                    self.HT[j][i] = complex(-self.DVT[i][j] * 1j)
                                    self.H0[j][i] = complex(-self.DV0[i][j] * 1j)
                                elif command[14] == '1':
                                    self.HT[j][i] = complex(self.DVT[j][i] * 1j)
                                    self.H0[j][i] = complex(self.DV0[j][i] * 1j)
    
                    # the GS energy is set as reference 0.0
                    self.HT[0][0] = complex(0.0)
                    self.H0[0][0] = complex(0.0)
            
                    logwrt.writelog("The Hamilton matrix at the present step\n", 1)
                    logwrt.writelog(logwrt.matrix_prettystring(self.HT), 1)
                    logwrt.writelog("\n", 1)
                    logwrt.writelog("The Hamilton matrix at the previous step\n", 1)
                    logwrt.writelog(logwrt.matrix_prettystring(self.H0), 1)
                    logwrt.writelog("\n", 1)
    
                    self.calculateProbability(tstep, nroots, AM1, command)
                
                    if command[80] != '0':
                        # for testing purposes
                        # a user defined value for the random number can be specified
                        self.rnum = float(command[80])
                        logwrt.writelog('User-defined "random" number (for testing purposes): {0:10.6f}\n'.format(self.rnum))
                    else:
                        # get a random number between 0 and 1
                        self.rnum = float(random.randrange(0, 1000000, 1)) / 1000000
                        logwrt.writelog('Random number selected for the hopping algorithm: {0:10.6f}\n'.format(self.rnum))
    
                    # APPLY DECOHERENCE CORRECTION
                    if command[85] == '1':
                        logwrt.writelog('\nUse decoherence correction from Granucci & Persico\n')
                        self.decoherence(tstep)
                        logwrt.writelog("Electronic states amplitudes after correction\n")
                        for i in range(nroots):
                            logwrt.writelog("state {0} = {1:8.6f} + {2:8.6f}*i\n".format(i + 1, self.AM[i].real, self.AM[i].imag))
                        logwrt.writelog("\n")
                    hop = False
    
                    # MAKE DECISION FOR HOPPING
                    # for HOPs from state 0:
                    # to 1: rnum < g_01
                    # to 2: g_01 < rnum < g_01 + g_02
                    # to 3: g_01 + g_02 < rnum < g_01 + g_02 + g_03
                    # for HOPs from state 1:
                    # to 0: rnum < g_10
                    # to 2: g_10 < rnum < g_10 + g_12
                    # to 3: g_10 + g_12 < rnum < g_10 + g_12 + g_13
                    # we need Psum and Psum_old
                    Psum = 0
                    Psum_old = 0
                    for i in range(nroots):
                        if i != self.actualstate:
                            Psum = Psum + self.P[self.actualstate][i]
                            if not hop:
                                # hopping condition is evaluated here
                                if Psum_old < self.rnum < Psum:
                                    logwrt.writelog('-----------------------------\n')
                                    logwrt.writelog('  !!!!Gimme Hop Joanna!!!\n')
                                    logwrt.writelog('-----------------------------\n')
                                    logwrt.writelog('hopping from state ' + str(self.actualstate + 1) + ' to --> ' + str(i + 1) + "\n")
                                    logwrt.writelog('-----------------------------\n')
                                    self.newstate = i
                                    hop = True
                                    #again, HOP file is still needed by VVerlet at present
                                    system('touch HOP')
                                    #CHECK FOR ALLOWED/FORBIDDEN BACK-HOP
                                    if i > self.actualstate and command[87] == '0':
                                        logwrt.writelog("HOP rejected because backhop is not allowed (key 87)!\n")
                                        self.newstate = self.actualstate
                                        hop = False
                                        #again, HOP file is still needed by VVerlet at present
                                        system('rm HOP')
                                else:
                                    self.newstate = self.actualstate
                            # print probabilities for states != active state
                            logwrt.writelog(' The Tully Hopping Probability to state ' + str(
                                i + 1) + ' is :  %10.6f ' % Psum_old + ' (total - current) %10.6f ' % self.P[self.actualstate][i] + ' (current) %10.6f ' % Psum + ' (total)' + "\n")
                        else:
                            # if i == active state
                            logwrt.writelog(' You are in this state ' + "\n")
                        # update cumulative probability
                        Psum_old = Psum
        
                logwrt.writelog("EXIT from Tully's FSSH algorithm\n\n")
            else:
                # if not TULLY or THS
                logwrt.writelog("Nothing to be done in Tully FSSH, resetting Amplitudes\n")
                self.AM = []
        
                # OW 7/18 print populations and occupations even when tully is inactive
        
                occup = []
                popul = []
        
                for i in range(nroots):
                    self.AM.append(complex(0.0))
                    occup.append(float(0.0))
                    popul.append(float(0.0))
                self.AM[self.actualstate] = complex(1.0)
                occup[self.actualstate] = float(1.0)
                popul[self.actualstate] = float(1.0)
                sef = shelve.open("cobram-sef")
        
                # OW introduce stopping criterion
        
                if occup[0] == 1 and float(command[122]) > 0:
                    stop = float(sef['stop'])
                    stop = stop + tstep / float('41.341373337')
                    if float(command[122]) > 0 and (stop >= float(command[122])):
                        sef.close()
                        logwrt.writelog('The trajectory has been in the ground state for '
                                        + str(command[122]) + ' fs, calculation stopped')
                        if not logwrt.DEBUG_COBRAMM_RUN:  CBF.garbager(self.geometry, command)
                        logwrt.cobramend()
                    else:
                        logwrt.writelog('the trajectory is in S0 for ' + str(stop) +
                                        ' fs, continuing until we reach ' + str(command[122]) + ' fs')
                        sef['stop'] = stop
                else:
                    sef['stop'] = '0.0'
        
                logwrt.writelog("Electronic states amplitudes\n")
                for i in range(nroots):
                    logwrt.writelog("state {0} = {1:8.6f} + {2:8.6f}*i\n".format(i + 1, self.AM[i].real, self.AM[i].imag))
                logwrt.writelog("\n")
        
                self.AM = np.array(self.AM)
                AM1 = copy.deepcopy(self.AM)
        
#                sef['AM'] = self.AM
                sef.close()
        
                with open('AmplitudesALL.dat', 'a') as out:
                    out.write('%10.4f' % self.SItime + '   '),
                    for i in range(nroots):
                        out.write('%16.12f' % self.AM[i].real + ' %16.12f ' % self.AM[i].imag + '   '),
                    out.write('\n')

        # NACs and TDNACs are stored in the shelve for compatibility with old VVerlet and product modules
        # !!!DEBUG !!! TO BE REVISED WHEN NEW VVERLET IS WRITTEN
        sef = shelve.open("cobram-sef")
        #update also old DEarray for new interface
        sef['DE_oldarray'] = copy.deepcopy(self.DEarrayNow)
        sef['newstate'] = self.newstate
        sef['SItime'] = self.SItime
        if command[85] == '1' and command[14] == '1':
            sef['TDC_old'] = copy.deepcopy(self.DVT)
        elif command[85] == '1' and command[14] == '0':
            NAC_old_format = []
            for root1 in range(nroots):
                NAC_old_format.append([])
                for root2 in range(nroots):
                    try:
                        NAC_old_format[root1].append([self.NACNow[root1][root2][0], self.NACNow[root1][root2][1], self.NACNow[root1][root2][2]])
                    except:
                        NAC_old_format[root1].append([[],[],[]])
            sef['NAC'] = NAC_old_format
        sef.close()
     
    # ===============================================================================================================  

    def buildDEarray(self, QMCalc, nroots):
        #builds a matrix of energy differences
        DEarray = [[1000.0 for j in range(nroots)] for i in range(nroots)]
        energies = QMCalc.outputData.get("energy")
        for i in range(nroots):
            for j in range(nroots):
                if i != j:
                    DEarray[i][j] = abs(energies[i] - energies[j]) * 627.51
        return np.array(DEarray)

    # =============================================================================================================== 

    def determineCopulingScheme(self, command, step, comp_state=None):
        # returns TULLY and THS variables (Booleans) that will determine the coupling used

        # in order to evaluate Tully one needs NACs at two consequtive time steps (for interpolation)
        # at step t-2dt: the energy gap between two states goes below threshold set in ediff, NACs are assigned for comput
        # -> in tully ttry is set to 'short' and stored in TSTEP, tstepNow and tstepPrev are 'long'
        # at step t-dt : an empty NAC(t-2dt) is stored to NACPrev, the NACs are computed and stored in NAC(t-dt) (NACNow)
        # -> in tully ttry is set to 'short' and stored in TSTEP, tstepPrev (from step t-2dt) is 'short'
        # at step t    : NAC(t-dt) are stored to NACPrev, the NACs are computed and stored in NAC(t) (NACNow)
        # -> in tully ttry is set to 'short' and stored in TSTEP, tstepPrev (from step t-dt) is 'short'
        # ENTER Tully and compute hopping probability between t-dt and t
        TULLY = False
        # with THS one can use longer time steps (i.e. both 'long' and 'short' should be set to the value of 'long')
        # in order to evaluate Tully one needs TDNACs at two consequtive time steps (for interpolation)
        # at step t-2dt: the energy gap between two states goes below threshold set in ediff, TDNACs are assigned for comput
        # at step t-dt : an empty TDNAC(t-2dt) is stored in TDCPrev, TDNAC(t-dt) is computed from WF(t-2dt) and WF(t-dt) (at step 1 or after hopping) or from WF(t-3dt), WF(t-2dt) and WF(t-dt)
        # at step t : TDNAC(t-dt) is stored as TDCPrev,TDNAC(t) is computed from WF(t-2dt), WF(t-dt) and WF(t)
        # ENTER Tully and compute hopping probability between t-dt and t
        THS = False

        if command[85] == '2':
            TULLY = True
            logwrt.writelog("Hopping scheme based only on energy difference\n")
        elif command[14] == '1' and step > 1:
            #if we use TDC and we are not in the first two steps
            THS = True
            logwrt.writelog("Solving the TDSE numerically\n")
        elif command[14] != '1' and self.NACNow[comp_state][self.actualstate] != [[], [], []] and \
            self.NACPrev[comp_state][self.actualstate] != [[], [], []] and self.tstepPrev == 'short' and self.tstepNow == 'short':
            #if we use NACs and we have NACs available from current and previous step and we have equal time steps (short) for current and previous step
            TULLY = True
            logwrt.writelog("Solving the TDSE numerically\n")

        return TULLY, THS
    
    # =============================================================================================================== 

    def computeDV(self, command, xvel, yvel, zvel, xvel0, yvel0, zvel0, nroots, QMNow, QMPrev):
        DVT, DV0 = copy.deepcopy(self.DVT), copy.deepcopy(self.DV0)
        # DETERMINE IF NACs CHANGED SIGN
        for i in self.NACNow:
            for j in self.NACNow[i]:
                self.NACNow[i][j] = np.array(self.NACNow[i][j])
                self.NACPrev[i][j] = np.array(self.NACPrev[i][j])
        for state_i in range(nroots):
            for state_j in range(state_i+1, nroots):
                # set the sign of the NAC(t) so that the angle with the NAC(t-dt) is < 90 deg
                # compute angle between the NAC(t-dt) and the NAC(t)
                angleDCp = product.angle(self.NACNow[state_i][state_j][0], self.NACNow[state_i][state_j][1], self.NACNow[state_i][state_j][2],
                                         self.NACPrev[state_i][state_j][0], self.NACPrev[state_i][state_j][1], self.NACPrev[state_i][state_j][2])[0]
                logwrt.writelog("Angle between  NACs at this time step and at previous time step < " + str(state_i + 1) + " |dR| " + str(state_j + 1) + " > is " + str(angleDCp) + "\n")
                # if the angle between the NAC(t) and NAC(t-dt) is > 90, change the sign
                # note that CIrot, which was taking care of this by comparing the WFs is no longer needed
                if abs(angleDCp) > 90:
                    if int(command[2]) > 0:
                        logwrt.writelog('The DC vector has changed sign! Sign correction applied!\n')
                    self.NACNow[state_i][state_j] = np.array(self.NACNow[state_i][state_j]) * float(-1.0)
                    self.NACNow[state_j][state_i] = np.array(self.NACNow[state_j][state_i]) * float(-1.0)
        
        # COMPUTE TIME-DERIVATIVE NACs (i.e. PROJECTIONS NAC x velocity) at
        # time t-dt (past step, H0) and t (present step, HT)

        #STEP 1) PROPAGATE VELOCITIES FROM t-0.5*dt (current) and t-1.5*dt (prev. step) to t and t-dt, respectively
        # this is identical to the first part of VVerlet propagation, except for the gradient used (in this case, we always use grdient of
        # current activa state as we do not know yet if a hop will occurr)
        # the velocities produced here will never be saved, so that the VVerlet is completrely not affected by this preliminary propagation
        # we always use SHORT timestep (command[84]) for propagation, because the condition to enter the computeDV module is that TULLY variable is true (i.e., we have used short timestep and
        # computed NACs sinbce at least two time steps ago)
        at = [self.geometry.atomLabel[i - 1] for i in self.geometry.list_QM]
        masses = []
        for i in range(len(at)):
            atommass = constants.atommass(at[i])
            if atommass is None:
                logwrt.fatalerror('no mass for ' + at[i] + ', please modify the code')
            else:
                masses.append(atommass)
        tstep = float(command[84]) * 41.341373337
        self.velocity['t'] = self.velocity['t-dt/2'] - 0.5 * tstep * np.array(QMNow.gradient(self.actualstate)) / masses
        self.velocity['t-dt'] = self.velocity['t-dt*3/2'] - 0.5 * tstep * np.array(QMPrev.gradient()[0]) / masses
        #for prev timestep, the state for which the gradient was computed might be different from self.actualstate (if a hop occurred at t-dt)
        #so we try to extract whatever gradient is saved in QMPrev (only gradient of previous active state should be in the library)
        if len(QMPrev.gradient()) > 1:
            logwrt.fatalerror("More than one gradient available when trying to propagate velocity of previous time step: I do not know what to choose!\n")

        #STEP 2) COMPUTE DC (DV)
        for i in range(nroots):
            for j in range(nroots):
                if i != j:
                    DVT[i][j] = np.dot(np.reshape(np.array(self.NACNow[i][j]), 3*len(at)), np.reshape(self.velocity['t'], 3*len(at)))
                    DV0[i][j] = np.dot(np.reshape(np.array(self.NACPrev[i][j]), 3*len(at)), np.reshape(self.velocity['t-dt'], 3*len(at)))

        return np.array(DVT), np.array(DV0)

    # ===============================================================================================================  

    def decoherence(self, tstep):
        #F: EKIN file must be removed after rewiting VVerlet 
        ekin=float(os.popen('cat EKIN').read())
        ASum=0
        #F: following assignment is probably needed to avoid division by zero
        if ekin == 0.0 :
           ekin=0.8
        '''following lines apply decoherence to all states.
        elev is the exponent for decoherence correction
        DEarray contans the delta E in kcal/mol (by construction), so iuìt must be translated to Hartree (1/627.51)
        Formula is taken from Granucci, Persico, JCP 126 (2007), 134114
        the 0.5 factor was added later (see Sharc code)'''
        for i in range(len(self.AM)):
            if self.actualstate != i:
               elev = -0.5 * tstep/((1.0/(abs(self.DEarrayNow[self.actualstate][i]/627.51)))*(1.0+0.1/ekin))
               self.AM[i]=self.AM[i]*math.exp(elev)
               ASum += abs(self.AM[i])**2
        #update amplitude of active state
        self.AM[self.actualstate] = self.AM[self.actualstate]*math.sqrt((1.0-ASum)/(abs(self.AM[self.actualstate])**2))
        #update sef
        sef=shelve.open("cobram-sef")
        sef['AM']=self.AM
        SItime=sef['SItime']
        sef.close()
        #update Amplitudes and AmpliutudesALL
        #TO REMEMBER !!!DEBUG:!!!! THERE IS A TIME MISMATCH BETWEEN AMPLITUDES (REAL) AND AmplitudesALL.dat FILE
        #I THINK THE ORIGIN IS HERE...CHECK!
        out=open('Amplitudes.dat','w')
        for i in range(len(self.AM)):
           out.write('%16.12f' %self.AM[i].real+' %16.12f ' %+self.AM[i].imag)
        out.close()
        out=open('AmplitudesALL.dat','a')
        out.write('%10.4f' %SItime+'   '),
        for i in range(len(self.AM)):
            out.write('%16.12f' %self.AM[i].real+' %16.12f ' %self.AM[i].imag+'   '),
        out.write('\n')
        out.close()

    # ===============================================================================================================

    def calculateProbability(self, tstep, nroots, AM1, command, cycle=100):
        # this function propagates the self.AM attributes (amplitudes) interpolating between t-dt and t with n=cycle (100) intermediate steps
        # and then calculates hopping probabilities (self.P)

        # set micro time step for electronic TDSE dt'
        dt = tstep / cycle
        # set Hamiltonian increment
        dH = (self.HT - self.H0) / cycle
        
        # CREATE AUXILIARY MATRICES
        b = np.zeros((nroots, nroots))
        g = np.zeros((nroots, nroots))
        
        He2 = np.zeros(b.shape, dtype=np.complex128)
        one = copy.deepcopy(He2)
        for i in range(nroots):
            one[i][i] = complex(1.0)
        # zero matrix
        He2 = np.array(He2)
        HeOLD = He2
        # unit matrix
        one = np.array(one)

        # START PROPAGATION
        for i in range(cycle):
            Ht = (self.H0 + (.5 + i) * dH) * dt
            # BUILT UP EXPONENTIAL e^Ht THROUGH SERIES EXPANSION
            He1 = one
            He = one
            fact = complex(1.0)
            for lop in range(1, 100):
                for j in range(nroots):
                    for k in range(nroots):
                        He2[j][k] = 0.0
                        # compute the lop-1 power of Ht in every term of the expansion (Ht^(lop-1))
                        for l in range(nroots):
                            He2[j][k] = He2[j][k] + Ht[l][k] * He1[j][l]
                # compute the factorial (lop-1)!
                fact = fact * (1j / lop)
                # He1 contains the lop-1 power of Ht for use in the next term
                He1 = He2
                # obtain the term Ht^(lop-1)/(lop-1)!
                He2 = He2 * fact
                # add term to expansion
                He = He + He2
                # check for convergence
                summ = sqrt(abs(np.add.reduce(np.add.reduce((HeOLD - He) ** 2)).real))
                if summ < 1.0e-60:
                    break
                if lop == 99:
                    logwrt.fatalerror('no convergency in the unitary propagator')
                HeOLD = He
            # propagate the Amplitudes
            AAA = product.product(self.AM, He)
            # Amplitudes at the end of the micro step t'+dt'
            self.AM = np.array(AAA)
            # Amplitdues at the center of the micro step as the average between
            # the Amplitudes at micro time step t' and t'+dt'
            actA = (AM1 + self.AM) * .5

            # COMPUTE HOPPING PROBABILITY
            # interpolate between initial TDNAC(t-dt) and final TDNAC(t) to get TDNAC(t-dt+i*dt')
            DVtmp = self.DV0 + (self.DVT - self.DV0) * (i + .5) / cycle
            DVtmp = np.array(DVtmp)
            for k in range(nroots):
                for j in range(k, nroots):
                    b[k][j] = float(-2.0 * ((actA[k] * np.conjugate(actA[j])) * DVtmp[k][j]).real)
                    b[j][k] = float(-2.0 * ((actA[k] * np.conjugate(actA[j])) * -DVtmp[k][j]).real)
                    g[k][j] = g[k][j] + b[k][j] * dt
                    g[j][k] = g[j][k] + b[j][k] * dt
            AM1 = copy.deepcopy(self.AM)

        # STORE AND WRITE AMPLITUDES
        #sef = shelve.open("cobram-sef")
        #sef['AM'] = self.AM
        #sef.close()

        # I have removed the writing of amplitudes for current step in Amplitudes.dat file (not needed I think)
        # reset to previous behavior if I am wrong :)
        logwrt.writelog("\nElectronic states amplitudes\n")
        for i in range(nroots):
            logwrt.writelog("state {0} = {1:8.6f} + {2:8.6f}*i\n".format(i+1, self.AM[i].real, self.AM[i].imag))
        logwrt.writelog("\n")
        
        #self.SItime = self.SItime + tstep / 41.341373337

        # CALCULATE PROBABILITY
        # indices are reversed here to indicate hopping probability from state k to l
        for k in range(nroots):
            for j in range(nroots):
                if abs(self.AM[k])**2 != 0:
                    self.P[k][j] = g[j][k] / abs(self.AM[k])**2
                else:
                    self.P[k][j] = 0.0
