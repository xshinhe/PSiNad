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

from typing import Union
import shelve  # python object persistence
import copy  # shallow and deep copy operations
import time # timing
import sys  # system-specific parameters and functions
import multiprocessing as mp # run child processes in parallel

# imports of local objects

from QMCalc import QM  # QM class controls the QM calculation, define its input and stores the output
from orbitals import AtomicOrbital  # object to parse and store orbital information from QM output files

# imports of local modules

import logwrt  # manages log file output + start/end procedures

# math libraries

import numpy as np  # numpy library for scientific computation


#####################################################################################################

class TDCoupling:

    def __init__(self, QMPrev: Union[QM, None], QMNow: Union[QM, None], threshold: float, 
            deltaT: float, nprocParall: int = 1, ediff: float = 1000.0, 
            couplingwithGS: bool = False, standardHST: bool = False, transformS: bool = False):
        """ Compute the time derivative coupling with finite difference formula according to the
         method described in Ryabinkin, Nagesh, & Izmaylov (2015). "Fast numerical evaluation of time-derivative
         nonadiabatic couplings for mixed quantum–classical methods". J.Phys.Chem.Lett, 6(21), 4200-4203.
         The method requires the CIS / TD(Tamm-Damkov) expansion of the wavefunction at the present time
         and at the previous time-step, then can be applied only to TD-DFT or CIS calculations.

        :param QMPrev: QM instance with the results from the QM calculation at previous time-step
        :param QMNow: QM instance with the results of the QM calculation at present time-step
        """

        # first initialize everything to None
        self.nroots = None
        self.DEarray = None
        self.AOovlp = None
        self.MOoverlap = None
        self.MOoccup = None
        self.statesNow = None
        self.statesPrev = None
        self.psioverlap = None
        self.psioverlapPrev = None
        self.tdcouplings = None
        self.signs = None

        # set threshold for discarding overlap of Slater determinant; 
        # if a diagonal of a SD has not element with magnitude > SDthres the overlap of this SD with other SDs is set to 0.0
        SDthres = 1.e-2

        # when the given QM is defined, we can compute some information on the states at present time
        if QMNow:

            # store the number of states that is considered here, that is equal to the number of
            # excitations + 1 (the ground state)
            self.nroots = len(QMNow.outputData.get("cis_coeffs")) + 1

            # initialize an array with a square matrix with the nr of states as dimension
            self.DEarray = np.zeros((self.nroots, self.nroots))
            # fill array with deltaE of the relevant states (all numbers should be positive)
            for i in range(self.nroots):
                for j in range(self.nroots):
                    if j <= i: continue
                    self.DEarray[i][j] = 627.51 * (QMNow.energydict[j] - QMNow.energydict[i])
                    self.DEarray[j][i] = self.DEarray[i][j]

        # to construct overlap, both QMNow and QMPrevious need to be defined
        if QMNow and QMPrev:

            if standardHST:
                logwrt.writelog("\nTime-derivative couplings computed with Tully Hammes-Schiffer formula\n", 1)
            elif transformS:
                logwrt.writelog("\nTime-derivative couplings computed by transformation of the MO overlap matrix\n", 1)
            else:
                logwrt.writelog("\nTime-derivative couplings computed with Ryabinkin, Nagesh, & Izmaylov formula\n", 1)
            if couplingwithGS:
                logwrt.writelog("Couplings with GS will be computed (be careful!)\n", 1)
            else:
                logwrt.writelog("Couplings with GS will be set to zero\n", 1)

            # time derivative couplings can be computed only if the basis func and the cis coeffs are available
            self.statesNow = copy.deepcopy(QMNow.outputData.get("cis_coeffs"))
            self.statesPrev = copy.deepcopy(QMPrev.outputData.get("cis_coeffs"))
            self.signs = copy.deepcopy(QMPrev.outputData.get("signs"))
            # at time step 1 set all signs to +1
            if self.signs == []:
                for iroot in range(self.nroots):
                    self.signs.append(1)
            self.psioverlapPrev = copy.deepcopy(QMPrev.outputData.get("psioverlap"))
            if not self.statesNow or not self.statesPrev or not QMNow.orbitals or not QMPrev.orbitals:
                logwrt.writewarning("missing wavefunction data, cannot compute td couplings")

            # check that the occupations of the two QM calculations are identical, and store the value
            if QMNow.orbitals.MOoccup != QMPrev.orbitals.MOoccup:
                logwrt.writewarning("occupations of the two QM differ, something might be wrong!")
                self.MOoccup = None
            else:
                self.MOoccup = QMNow.orbitals.MOoccup

            # as a first control check, check the orthonormalization of the CIS solutions
            for K, leftState in self.statesNow.items():
                for J, rightState in self.statesNow.items():
                    if J < K: continue  # skip triangle and do only unique elements
                    ov = self._CIScoeffSum(leftState, rightState)
                    if K == J and not np.isclose(ov, 1.0, atol=1.e-5):
                        logwrt.writewarning("normalization of state {0} is {1}".format(K, ov))
                    elif K != J and not np.isclose(ov, 0.0, atol=2.e-5):
                        logwrt.writewarning("overlap between states {0} and {1} is {2}".format(K, J, ov))

            # compute the overlap between the AOs of the two functions
            AO1 = QMNow.orbitals.atomicOrbitals
            AO2 = QMPrev.orbitals.atomicOrbitals
            self.AOovlp = np.zeros([len(AO1), len(AO2)])
            for ileft, ao_left in enumerate(AO1):
                for iright, ao_right in enumerate(AO2):
                    self.AOovlp[ileft, iright] = AtomicOrbital.AOoverlap(ao_left, ao_right)

            # now compute the overlap between the MOs
            MO1 = np.array(QMNow.orbitals.MOcoeffs)
            MO2 = np.array(QMPrev.orbitals.MOcoeffs)

            if standardHST or transformS:
                self.MOoverlap = np.matmul(MO1, np.matmul(self.AOovlp, MO2.T))
                # print the new MO overlap matrix, after the transformation
                logwrt.writelog("\nOverlap matrix\n", 3)
                logwrt.writelog(logwrt.matrix_prettystring(self.MOoverlap, fmt=".6f"), 3)
                logwrt.writelog("\n", 3)
            else:

                self.orbitalTransform = np.matmul(MO1, np.matmul(self.AOovlp, MO2.T)).round().T
                Overlap = np.matmul(MO1, np.matmul(self.AOovlp, MO2.T)).T
             
                # print the MO overlap matrix, rounded up to define the orbital transformation
                logwrt.writelog("\nRounded overlap matrix - ORBITAL TRANSFORMATION\n", 3)
                logwrt.writelog(logwrt.matrix_prettystring(self.orbitalTransform, fmt=".0f"), 3)
                logwrt.writelog("\n", 3)
             
                #check if there is more than one element in orbitalTransform = 1 or -1
                self.orbitalTransform = self._removeAmbiguousValues(self.orbitalTransform, Overlap)
             
                # print the MO overlap matrix, rounded up to define the orbital transformation
                logwrt.writelog("\nRounded overlap matrix - ORBITAL TRANSFORMATION after removing ambiguities\n", 3)
                logwrt.writelog(logwrt.matrix_prettystring(self.orbitalTransform, fmt=".0f"), 3)
                logwrt.writelog("\n", 3)
             
                # use the overlap to transform the orbitals in phase and order, to get the a continous basis set in dt
                MO1 = np.matmul(self.orbitalTransform, MO1)
                # use the same transformation for the CI coefficients
                self.statesNow = self._transformCICoefficients(self.statesNow, self.orbitalTransform)
             
                logwrt.writelog("\nAFTER TRANSFORMATION OF THE CURRENT WF\n", 3)
                logwrt.writelog("*** Molecular orbitals ***\n", 3)
                for i, mo1, mo2 in zip(range(len(MO1)), MO1, MO2):
                    logwrt.writelog("MO number {0}\n".format(i+1), 3)
                    logwrt.writelog("prev: " + str(mo2) + "\n", 3)
                    logwrt.writelog("now: " + str(mo1) + "\n", 3)
                logwrt.writelog("*** CI coefficients ***\n", 3)
                for state in self.statesNow:
                    logwrt.writelog("CI root number {0}\n".format(state), 3)
                    logwrt.writelog("prev:" + str(self.statesPrev[state]) + "\n", 3)
                    logwrt.writelog("now:" + str(self.statesNow[state]) + "\n", 3)
             
                # now compute the new overlap
                self.MOoverlap = np.matmul(MO1, np.matmul(self.AOovlp, MO2.T))
             
                # print the new MO overlap matrix, after the transformation
                logwrt.writelog("\nOverlap matrix after the transformation\n", 3)
                logwrt.writelog(logwrt.matrix_prettystring(self.MOoverlap, fmt=".6f"), 3)
                logwrt.writelog("\n", 3)

                # now compute the overlap between CI vectors, assuming no MO change for simplicity
                civecoverlap = np.zeros((self.nroots, self.nroots))
                for i, stateleft in self.statesPrev.items():
                    for j, stateright in self.statesNow.items():
                        # compute overlap (if the MO are assumed equal, then the Slater determinants are orthonormal)
                        # loop over all the excitations contained in the left and right expansion
                        for exL in stateleft:
                            if exL in stateright:
                                civecoverlap[i, j] += stateleft[exL] * stateright[exL]
                # the ground state is by construction always well identified, so we can set the overlap = 1
                civecoverlap[0, 0] = 1.0
             
                # print the new CI eigenvectors overlap
                logwrt.writelog("\nCI roots overlap after the MO transformation (assuming no MO change)\n", 0)
                logwrt.writelog(logwrt.matrix_prettystring(civecoverlap, fmt=".6f"), 0)
                logwrt.writelog("\n", 0)

            # extract the list of excitations to include in the coupling matrix
            exstatelist = list(set(self.statesPrev.keys()) & set(self.statesNow.keys()))
            exstatelist.sort()

            # initialize the overlap as a matrix of zeros
            self.tdcouplings = np.zeros((self.nroots, self.nroots))

            if transformS:
                # initialize the overlap as a matrix of zeros
                self.psioverlap = np.zeros((self.nroots, self.nroots))

                # find out number of occ. and virt. orbitals
                occMO=self.MOoccup.count("O")
                virtMO=self.MOoccup.count("V")

                start = time.time()
                # extract the sub-matrices of occ and virt orbtial overlap
                OccS = self.MOoverlap[0:occMO,0:occMO]
                VirtS = self.MOoverlap[occMO:,occMO:]

                # find the inverse which is also the transformation to a unit matrix
                transOccS = np.linalg.inv(OccS)
                transVirtS = np.linalg.inv(VirtS)

                # put the dict statesPrev into a list of occMO x virtMO matrices
                CISarrayPrev = []
                for i, state in self.statesPrev.items():
                    CISarray = np.zeros((occMO,virtMO))
                    for ex, coeff in state.items():
                        CISarray[ex[0]-1][ex[1]-occMO-1] = coeff
                    CISarrayPrev.append(CISarray)

                # put the dict statesNow into a list of occMO x virtMO matrices
                CISarrayNowT = []
                for i, state in self.statesNow.items():
                    CISarray = np.zeros((occMO,virtMO))
                    for ex, coeff in state.items():
                        CISarray[ex[0]-1][ex[1]-occMO-1] = coeff
                    # apply transformatinal matrices from left and right
                    CISarrayNowT.append(np.matmul(transOccS.T, np.matmul(CISarray, transVirtS)))

                # compute WF overlap S_ij = sum_a C_ia * C_ja (all cross terms cancel out as MO overlap is zero)
                for i in range(1,self.nroots):
                    for j in range(1,self.nroots):
                        self.psioverlap[i, j] = np.sum(np.multiply(CISarrayPrev[i-1],CISarrayNowT[j-1]))

                end = time.time()
                logwrt.writelog("Time spent to compute WF overlap: {0}s\n".format(end-start), 2)

                # print the new CI eigenvectors overlap
                logwrt.writelog("\nCI roots overlap after the MO transformation (with exact formula)\n", 2)
                logwrt.writelog(logwrt.matrix_prettystring(self.psioverlap, fmt=".6f"), 2)
                logwrt.writelog("\n", 2)

                # correct sign of the WF to assure continuity
                self.psioverlap, self.signs = self._correctSigns(self.psioverlap, self.signs, self.nroots)

                # save overlap and signs for evaluation in next step
                QMNow.outputData.set('signs', self.signs)
                QMNow.outputData.set('psioverlap', self.psioverlap)

                # now construct the td couplings
                if self.psioverlapPrev is None:
                    # if there is only one overlap matrix (that is at step 1) compute couplings from
                    # the antisymmetrized projection using hammes-schiffer-tully formula
                    logwrt.writelog("Compute couplings using overlap matrix of WFs at the previous 2 steps\n", 2)
                    self.tdcouplings = 0.5 * (self.psioverlap - self.psioverlap.T) / deltaT
                else:
                    # else use finite differences to compute the couplings using overlap matrices at two consecutive points
                    logwrt.writelog("Compute couplings using overlap matrices of WFs at the previous 3 steps\n", 2)
                    self.tdcouplings = 0.25 * (3*(self.psioverlap - self.psioverlap.T) - (self.psioverlapPrev - self.psioverlapPrev.T)) / deltaT

            elif standardHST:
                # initialize the overlap as a matrix of zeros
                self.psioverlap = np.zeros((self.nroots, self.nroots))

                # build MO overlap matrix for closed shell SD 
                occMO=self.MOoccup.count("O")
                closedShellOverlapMat = self.MOoverlap[0:occMO,0:occMO]

                start = time.time()
                # store overlaps of determinants as they can be re-used
                # keys are 4-index identifiers occ->vit (left), occ->virt (right) 
                detoverlap = {}

                # analyze which TD couplings will be computed according to the energy gap threshold
                compStates = []
                skipStates = [] 
                for i in self.statesPrev.keys():
                    for j in self.statesNow.keys():
                        if j > i and self.DEarray[i][j] > ediff:
                            skipStates.append((i, j))
                        elif j > i and self.DEarray[i][j] <= ediff:
                            compStates.append((i, j))
                if compStates == []:
                    logwrt.writelog("No ES-ES couplings will be computed in this step\n", 1)
                elif compStates != []:
                    logwrt.writelog("Compute ES-ES couplings for: ", 1)
                    for el in compStates:
                        logwrt.writelog("{0}-{1} ".format(el[0]+1, el[1]+1), 1)
                    logwrt.writelog("\n", 1)
                if skipStates != []:
                    logwrt.writelog("Skip ES-ES couplings for: ", 1)
                    for el in skipStates:
                        logwrt.writelog("{0}-{1} ".format(el[0]+1, el[1]+1), 1)
                    logwrt.writelog("\n", 1)
                # now compute the overlap between CI vectors
                expRtruncInfo = {}
                for i, stateleft in self.statesPrev.items():
                    # truncate the length of CI expansion according to user-specified threshold
                    fullLength = len(stateleft)
                    if threshold < 1.:
                        stateleft = self._truncateExpansion(stateleft, threshold)
                        logwrt.writelog("Length of truncated expansion for state {0:d} from previous step: {1:d} out of {2:d}\n".format(i+1, len(stateleft), fullLength), 2)
                        # compute norm of the tuncated expansion
                        stateleft_norm = np.sqrt(sum([val**2 for val in stateleft.values()]))
                    else: 
                        stateleft_norm = 1.0
                    for j, stateright in self.statesNow.items():
                        if self.DEarray[i][j] <= ediff:
                            # truncate the length of CI expansion according to user-specified threshold
                            fullLength = len(stateright)
                            if threshold < 1.:
                                stateright = self._truncateExpansion(stateright, threshold)
                                if j not in expRtruncInfo.keys():
                                    expRtruncInfo[j] = (len(stateright), fullLength)
                                # compute norm of the tuncated expansion
                                stateright_norm = np.sqrt(sum([val**2 for val in stateright.values()]))
                            else:
                                stateright_norm = 1.0
                            # parallelize the calculation of the SD overlaps 
                            #if nprocParall > 1:
                            #    with mp.Pool(processes=nprocParall) as pool:
                            #        for exL in stateleft.keys():
                            #            result = pool.apply_async(self._slateroverlapP, (exL, list(stateright.keys()), detoverlap, self.MOoverlap, copy.deepcopy(closedShellOverlapMat), occMO))
                            #            detoverlap.update(result.get())
                         
                            # loop over all the excitations contained in the left and right expansion
                            for exL, coeffL in stateleft.items():
                                # check if the column of the SD containing the virtual orbital exL[1] has an element with magnitude > SDthres
                                if abs(sorted(self.MOoverlap[exL[1]-1,0:occMO], key=lambda coeff: abs(coeff),reverse=True)[0]) < SDthres:
                                    for exR, coeffR in stateright.items():
                                        # check that the element with the virtual orbitals exL[1] and exR[1] is > SDthres
                                        if abs(self.MOoverlap[exL[1]-1,exR[1]-1]) < SDthres:
                                            # otherwise set to 0.0
                                            detoverlap[exL+exR] = 0.0
                                        else:
                                            # if the element of the SD overlap was alrady computed take it from detoverlap
                                            try:
                                                # add contribution to the i,j element of the overlap for this couple of slater determ.
                                                self.psioverlap[i, j] += coeffL * detoverlap[exL+exR] * coeffR
                                            # othewise compute and store in detoverlap
                                            except:
                                                #print("Compute determinant for ", exL, "and ", exR)
                                                detoverlap[exL+exR] = self._slateroverlap(self.MOoverlap, copy.deepcopy(closedShellOverlapMat), occMO, exL, exR)
                                                self.psioverlap[i, j] += coeffL * detoverlap[exL+exR] * coeffR
                                else:
                                    for exR, coeffR in stateright.items():
                                        # check if the column of the SD containing the virtual orbital exR[1] has an element with magnitude > SDthres
                                        if abs(sorted(self.MOoverlap[0:occMO,exR[1]-1], key=lambda coeff: abs(coeff),reverse=True)[0]) < SDthres and abs(self.MOoverlap[exL[1]-1,exR[1]-1]) < SDthres:
                                            # otherwise set to 0.0
                                            detoverlap[exL+exR] = 0.0
                                        else:
                                            # if the element of the SD overlap was alrady computed take it from detoverlap
                                            try:
                                                # add contribution to the i,j element of the overlap for this couple of slater determ.
                                                self.psioverlap[i, j] += coeffL * detoverlap[exL+exR] * coeffR
                                            # othewise compute and store in detoverlap
                                            except:
                                                #print("Compute overlap for ", exL, "and ", exR)
                                                detoverlap[exL+exR] = self._slateroverlap(self.MOoverlap, copy.deepcopy(closedShellOverlapMat), occMO, exL, exR)
                                                self.psioverlap[i, j] += coeffL * detoverlap[exL+exR] * coeffR
                            self.psioverlap[i, j] /= (stateleft_norm*stateright_norm) 
                        else:
                            self.psioverlap[i, j] = 0.0
                logwrt.writelog("\n")
                for i, lengths in expRtruncInfo.items():
                    logwrt.writelog("Length of truncated expansion for state {0:d} from current step: {1:d} out of {2:d}\n".format(i+1, lengths[0], lengths[1]), 2)
                end = time.time()
                logwrt.writelog("Time spent to compute WF overlap: {0}s\n".format(end-start), 2)

                # when one wants to compute explicitely the coupling with the GS, compute the relevant overlap elements
                # otherwise they will be left = 0
                if couplingwithGS:
                    for i, state in self.statesPrev.items():
                        if self.DEarray[0][i] <= ediff:
                            for ex, coeff in state.items():
                                detoverlap = self._slateroverlap(self.MOoverlap, closedShellOverlapMat, occMO, ex, (0, 0))
                                self.psioverlap[i, 0] += coeff * detoverlap * 1.0
                                detoverlap = self._slateroverlap(self.MOoverlap, closedShellOverlapMat, occMO, (0, 0), ex)
                                self.psioverlap[0, i] += 1.0 * detoverlap * coeff
                        else:
                            self.psioverlap[i, 0] = 0.0
                            self.psioverlap[0, i] = 0.0
                    # compute overlap for GS
                    self.psioverlap[0, 0] = self._slateroverlap(self.MOoverlap, closedShellOverlapMat, occMO, (0, 0), (0, 0))

                # print the new CI eigenvectors overlap
                logwrt.writelog("\nCI roots overlap after the MO transformation (with exact formula)\n", 2)
                logwrt.writelog(logwrt.matrix_prettystring(self.psioverlap, fmt=".6f"), 2)
                logwrt.writelog("\n", 2)

                # correct sign of the WF to assure continuity
                self.psioverlap, self.signs = self._correctSigns(self.psioverlap, self.signs, self.nroots)

                # save overlap and signs for evaluation in next step
                QMNow.outputData.set('signs', self.signs)
                QMNow.outputData.set('psioverlap', self.psioverlap)

                # now construct the td couplings 
                if self.psioverlapPrev is None:
                    # if there is only one overlap matrix (that is at step 1) compute couplings from 
                    # the antisymmetrized projection using hammes-schiffer-tully formula
                    logwrt.writelog("Compute couplings using overlap matrix of WFs at the previous 2 steps\n", 2)
                    self.tdcouplings = 0.5 * (self.psioverlap - self.psioverlap.T) / deltaT
                else:
                    # else use finite differences to compute the couplings using overlap matrices at two consecutive points
                    logwrt.writelog("Compute couplings using overlap matrices of WFs at the previous 3 steps\n", 2)
                    self.tdcouplings = 0.25 * (3*(self.psioverlap - self.psioverlap.T) - (self.psioverlapPrev - self.psioverlapPrev.T)) / deltaT
                
            else:
                # compute the couplings between ES among themselves
                for K in exstatelist:
                    for J in exstatelist:
                        # skip diagonal elements
                        if K == J: continue
                        # for the others, compute the relevant terms
                        term1 = self._CIScoeffSum(self.statesPrev[K], self.statesNow[J])
                        term23 = self._CISsumWithOvlp(self.statesPrev[K], self.statesPrev[J], self.MOoverlap)
                        # add the terms in an antisymmetrized formula
                        self.tdcouplings[K, J] += 0.5 * (term1 + term23) / deltaT
                        self.tdcouplings[J, K] -= 0.5 * (term1 + term23) / deltaT

                # when requested, the coupling with GS is computed using the GS/K-state td derivative expression
                # this option should be used with caution: standard TD-DFT is known to fail for GS-ES crossings
                # since the WF for the GS is intrinsically single-reference
                if couplingwithGS:
                    for K in exstatelist:
                        kzero = self._CISgroundexcitedOverlap(self.statesPrev[K], self.statesNow[K],
                                                              self.MOoverlap) / deltaT
                        self.tdcouplings[K, 0] = kzero
                        self.tdcouplings[0, K] = - kzero

            # print the time-derivative couplings matrix
            logwrt.writelog("\nTime-derivative couplings at this step\n", 1)
            logwrt.writelog(logwrt.matrix_prettystring(self.getTDMatrix(), ".8f"), 1)
            logwrt.writelog("\n", 1)

        # now, for compatibility with what is already implemented in COBRAMM,
        # call the function that updates some relevant quantities saved in the shelve
        self.updateShelve()

    # ===============================================================================================================

    def updateShelve(self):

        sef = shelve.open("cobram-sef")

        # store the number of roots of the electronic structure
        sef['nroots'] = self.nroots

        # store a matrix of zeros, with nroots as size
        if self.nroots:
            sef['zeromat'] = np.zeros((self.nroots, self.nroots))
        else:
            sef['zeromat'] = None

        # store the new array with DeltaE in sef['DEarray'], and backup the old one in sef['DE_oldarray']
        sef['DE_oldarray'] = copy.deepcopy(sef['DEarray'])
        sef['DEarray'] = self.DEarray

        # store the new array with TDC in sef['TDC'], and backup the old one in sef['TDC_old']
        sef['TDC_old'] = copy.deepcopy(sef['TDC'])
        sef['TDC'] = self.getTDMatrix()

        sef.close()

    # ===============================================================================================================

    def getTDMatrix(self):
        """ Returns the TD couplings in form of a numpy 2D array """

        if self.tdcouplings is not None:  # when the TD couplings are defined, give the matrix with the values
            return self.tdcouplings
        else:  # otherwise return an array of zeros
            return np.zeros((self.nroots, self.nroots))

    # ===============================================================================================================

    @staticmethod
    def _CISgroundexcitedOverlap(expansionLeft: dict, expansionRight: dict, overlap: np.ndarray) -> float:
        """
        Compute an antisymmetrized summation of the form:
        Sum_{i,a} Cleft_{ia} overlap{i,a} - Sum_{i,a} overlap{a, i} Cright_{ia}
        with i being an index for the occupied orbital and a for the virtual orbitals.
        The expansions of Cleft and Cright are defined as dictionary, with the coefficients labelled
        by the integer tuples (i, a). The dictionaries do not need to list all the possible excitations,
        it is assumed that the excitations that are not listed in the dictionary are equal to 0.
        The overlap is defined as a numpy array with two dimensions, labelled by the ordinal number of the orbital
        (attention! orbital indices start from 1, while array indices start from 0!)

        :param expansionLeft: dictionary listing the coefficients of a CIS expansion Cleft_{ia}
        :param expansionRight: dictionary listing the coefficients of a CIS expansion Cright_{ia}
        :param overlap: numpy 2D array with the overlap between the orbitals overlap{i/a,j/b}
        :return: floating point number with sum computed by this function
        """

        # initialize value of the summation
        s = 0.0

        # loop over all the excitations contained in the left and right expansion
        for exL, cL in expansionLeft.items():
            i, a = exL
            s += cL * overlap[a - 1, i - 1]
        # and now subtract the symmetric element, to get antisymmetric sum
        for exR, cR in expansionRight.items():
            i, a = exR
            s -= overlap[i - 1, a - 1] * cR
        # divide by two (antisymmetric sum: 0.5 * [A_K0 - A_OK] )
        s = s / 2.0

        # return the summation
        return s

    # ===============================================================================================================

    @staticmethod
    def _CIScoeffSum(expansionLeft: dict, expansionRight: dict) -> float:
        """
        Compute a summation of the form Sum_{i,a} Cleft_{ia} Cright_{ia} where i-a are single
        excitations labelling the coefficients. The expansions of Cleft and Cright are defined
        as dictionaries, with the coeffiecients labelled by the integer tuples (i, a).
        The dictionaries do not need to list all the possible excitations, it is assumed that the
        excitations that are not listed in the dictionary are equal to 0.

        :param expansionLeft: dictionary listing the coefficients of a CIS expansion Cleft_{ia}
        :param expansionRight: dictionary listing the coefficients of a CIS expansion Cright_{ia}
        :return: floating point number with Sum_{i,a} Cleft_{ia} Cright_{ia}
        """

        # initialize value of the summation
        s = 0.0

        # loop over all the excitations contained in the left expansion
        for ex in expansionLeft.keys():
            # only if an excitation is in both expansions we need to add something, otherwise contribution is null
            if ex in expansionRight:
                s += expansionLeft[ex] * expansionRight[ex]

        # return the summation
        return s

    # ===============================================================================================================

    @staticmethod
    def _CISsumWithOvlp(expansionLeft: dict, expansionRight: dict, overlap: np.ndarray) -> float:
        """
        Compute a summation of the form:
        Sum_{i,a,b} Cleft_{ia} overlap{a,b} Cright_{ib} + Sum_{i,j,a} Cleft_{ia} overlap{j,i} Cright_{ja}
        with i, j being indices for the occupied orbital and a,b indices of virtual orbitals.
        The expansions of Cleft and Cright are defined as dictionaries, with the coefficients labelled
        by the integer tuples (i, a). The dictionaries do not need to list all the possible excitations,
        it is assumed that the excitations that are not listed in the dictionary are equal to 0.
        The overlap is defined as a numpy array with two dimensions, labelled by the ordinal number of the orbital
        (attention! orbital indices start from 1, while array indices start from 0!)

        :param expansionLeft: dictionary listing the coefficients of a CIS expansion Cleft_{ia}
        :param expansionRight: dictionary listing the coefficients of a CIS expansion Cright_{ia}
        :param overlap: numpy 2D array with the overlap between the orbitals overlap{i/a,j/b}
        :return: floating point number with sum computed by this function
        """

        # initialize value of the summation
        s = 0.0

        # loop over all the excitations contained in the left and right expansion
        for exL in expansionLeft.keys():
            i, a = exL
            for exR in expansionRight.keys():
                j, b = exR
                # when the excitations start from the same orbital, there is a positive contribution to add
                if i == j:
                    s += expansionLeft[exL] * overlap[a - 1, b - 1] * expansionRight[exR]
                # when the excitations end in the same orbital, there is a negative contribution to add
                if a == b:
                    s -= expansionLeft[exL] * overlap[j - 1, i - 1] * expansionRight[exR]

        # return the summation
        return s

    # ===============================================================================================================

    @staticmethod
    def _removeAmbiguousValues(orbTransf: np.ndarray, overlap: np.ndarray) -> np.ndarray:
        """
        Remove ambiguity in the transformation by assuring that there is only one element per row which is 1 or -1
        """

        for istart, overlapstate in enumerate(orbTransf):
            if len(np.where(np.isclose(abs(overlapstate), 1.0))[0]) > 1:
                logwrt.writelog("Warning! Strong orbital mixing for reference orbital {:d} with orbitals: ".format(istart+1), 2)
                iend, = np.where(np.isclose(abs(overlapstate), 1.0))
                for iorb in iend:
                    logwrt.writelog("{0:d} ({1:.2f}) ".format(iorb+1, overlap[istart][iorb]), 2)
                logwrt.writelog("\n", 2)
                if istart in iend:
                    orbTransf[istart] = [(0 if i != istart else val) for i,val in enumerate(overlapstate)]
                    logwrt.writelog("Diagonal value retained {0:d},{1:d} = {2:d}\n".format(istart+1,istart+1,int(orbTransf[istart][istart])), 2)
                else:
                    ind, = np.where(np.equal(abs(overlap[istart]),np.max(abs(overlap[istart]))))
                    orbTransf[istart] = [(0 if i != ind[0] else val) for i,val in enumerate(overlapstate)]
                    logwrt.writelog("Off-diagonal value retained {0:d},{1:d} = {2:d}\n".format(istart+1,ind[0]+1,int(orbTransf[istart][ind[0]])), 2)
            elif len(np.where(np.isclose(abs(overlapstate), 1.0))[0]) == 0:
                logwrt.writelog("Warning! Strong orbital mixing for reference orbital {:d} with orbitals: ".format(istart+1), 2)
                for i,val in enumerate(overlap[istart]):
                    thres = np.max(abs(overlap[istart]))/2
                    if abs(val) > thres:
                        logwrt.writelog("{0:d} ({1:.2f}) ".format(i+1, val), 1)
                logwrt.writelog("\n", 1)
                ind, = np.where(np.equal(abs(overlap[istart]),np.max(abs(overlap[istart]))))
                orbTransf[istart] = [(0 if i != ind[0] else (1 if orbTransf[istart][ind[0]] > 0 else -1)) for i in range(len(overlapstate))]
                if ind[0] == istart:
                    logwrt.writelog("Diagonal value retained {0:d},{1:d} = {2:d}\n".format(istart+1,ind[0]+1,int(orbTransf[istart][ind[0]])), 2)
                else:
                    logwrt.writelog("Off-diagonal value retained {0:d},{1:d} = {2:d}\n".format(istart+1,ind[0]+1,int(orbTransf[istart][ind[0]])), 2)
        for istart,overlapstate in enumerate(orbTransf):
            if len(np.where(np.isclose(abs(overlapstate), 1.0))[0]) == 0:
                 logwrt.fatalerror("Warning! The transformation matrix has an empty row {:d}\n".format(istart+1))
        for istart,overlapstate in enumerate(orbTransf.T):
            if len(np.where(np.isclose(abs(overlapstate), 1.0))[0]) == 0:
                 logwrt.writelog("Warning! The transformation matrix has an empty column {:d}\n".format(istart+1), 2)
        for istart,overlapstate in enumerate(orbTransf):
            if len(np.where(np.isclose(abs(overlapstate), 1.0))[0]) > 1:
                 logwrt.fatalerror("Warning! Ambiguity could not be resolved for row {:d}\n".format(istart+1))
        for istart,overlapstate in enumerate(orbTransf.T):
            if len(np.where(np.isclose(abs(overlapstate), 1.0))[0]) > 1:
                 logwrt.writelog("Warning! Ambiguity could not be resolved for column {:d}\n".format(istart+1), 2)

        return orbTransf

    # ===============================================================================================================

    @staticmethod
    def _transformCICoefficients(CIcoeffs: dict, orbTransf: np.ndarray) -> dict:
        """
        Transform the coefficients of a CI expansion contained in CIcoeffs, but applying the
        transformation of orbTransf. orbTransf is a matrix of 1., 0. and -1., which describes an
        MO transformation of phase and labelling that relates the basis set of the old QM and the
        new one. For conveniency, this matrix is recasted to a mapping dictionary and then used to
        relabel and change sign of the CI coefficients.

        :param CIcoeffs: dictionary with the CI coefficients to transform
        :param orbTransf: matrix that contains the MO transformation that we need to apply to the CI coeff.s
        :return: a new dictionary with the modified CI coeffiecients
        """

        mo2mo = TDCoupling._overlap2statemapping(orbTransf)

        CIcoeffsNew = {}
        for istate, ciexpansion in CIcoeffs.items():
            newexpansion = {}
            for exc, coeff in ciexpansion.items():
                # make a copy of the excitation data
                newexc = [exc[0], exc[1]]
                newcoeff = coeff
                # now check if occupied orbital is transformed
                if exc[0] in mo2mo:
                    newexc[0] = abs(mo2mo[exc[0]])
                    newcoeff *= np.sign(mo2mo[exc[0]])
                # and then do the same for virtual orbital
                if exc[1] in mo2mo:
                    newexc[1] = abs(mo2mo[exc[1]])
                    newcoeff *= np.sign(mo2mo[exc[1]])
                # store the transformed element in the data dictionary for this excited state
                newexpansion[tuple(newexc)] = newcoeff
            # store the excited state in the final dictionary
            CIcoeffsNew[istate] = newexpansion

        return CIcoeffsNew

    # ===============================================================================================================

    @staticmethod
    def _overlap2statemapping(overlapmatrix: np.ndarray) -> dict:
        """
        Take as input an overlap matrix, and based on whether the values are close to 1, 0 or -1
        decide if between the left and right state there has been a swap or a phase change.
        Returns a dictionary, which contains only those state for which there has been an actual
        change. The key of the dictionary is the initial state, the value is the final state
        with + sign when the phase is the same and - sign when the phase has changed.

        :param overlapmatrix: numpy array with a square overlap matrix
        :return: dictionary that maps left states to right states, with information on phase change
        """

        # initialize a dictionary to store the mapping {istate: +/- jstate} (+ is saved only when istate != jstate)
        statemap = {}

        for istart, overlapstate in enumerate(overlapmatrix.round()):
            # find the element in a row that is close to 1
            iend, = np.where(np.isclose(overlapstate, 1.0))
            try:
                # save the state mapping only when there is an actual change of state labelling
                if iend[0] != istart:
                    statemap[istart + 1] = iend[0] + 1
            # when there is no element close to one, it means that there is an element close to -1
            except IndexError:
                iend, = np.where(np.isclose(overlapstate, -1.0))
                # in this case always save the results, because in any case there is a phase change
                statemap[istart + 1] = -(iend[0] + 1)

        # retunr the dictionary with the mapping
        return statemap

    # ===============================================================================================================

    @staticmethod
    def _slateroverlap(fullMOoverlap: np.ndarray, MOoverlap: np.ndarray, occ: int,
                       exleft: tuple, exright: tuple) -> float:
        """
        Compute the a overlap between two slater determinants as determinant of the overlap sub-matrix,
        constructed by selecting the elements of the MO overlap according to slater-condon generalized rules
        works only for singly excited SD

        @param fullMOoverlap: overlap matrix of all the MO, occupied and virtual
        @param MOoverlap: overlap matrix of occupied orbitals in the current SD
        @param exleft: (i, a) tuple labelling the excitation on the left of the matrix element
        @param exright: (i, a) tuple labelling the excitation on the right of the matrix element
        @return: the matrix
        """

        # replace row and column in closed shell MO overlap according to matrix elements 
        MOoverlap[exleft[0]-1,:] = fullMOoverlap[exleft[1]-1,0:occ] 
        MOoverlap[:,exright[0]-1] = fullMOoverlap[0:occ,exright[1]-1]
        MOoverlap[exleft[0]-1,exright[0]-1] = fullMOoverlap[exleft[1]-1,exright[1]-1]

        # compute the determinant of the matrix
        SDoverlap = np.linalg.det(MOoverlap)

        return SDoverlap

    # ===============================================================================================================

    def _slateroverlapP(self, exL: tuple, exRList: list, detoverlap: dict, 
                        fullMOoverlap: np.ndarray, MOoverlap: np.ndarray, occ: int) -> dict:
        """ Compute a sub-set of Slater determinants"""

        exLdetoverlap = {}
        for exR in exRList:
            if exL+exR not in detoverlap:
                #exLdetoverlap[exL+exR] = self._slateroverlap(fullMOoverlap, copy.deepcopy(MOoverlap), occ, exL, exR)
                MOoverlap[exL[0]-1,:] = fullMOoverlap[exL[1]-1,0:occ]
                MOoverlap[:,exR[0]-1] = fullMOoverlap[0:occ,exR[1]-1]
                MOoverlap[exL[0]-1,exR[0]-1] = fullMOoverlap[exL[1]-1,exR[1]-1]
                exLdetoverlap[exL+exR] = np.linalg.det(MOoverlap)
        
        return exLdetoverlap

    # ===============================================================================================================

    @staticmethod
    def _truncateExpansion(expansion: dict, threshold: float) -> dict:
        """
        Truncate the CI expansion so that the total sum_i |CI_i|^2 > threshold
        """

        # sort CI coeff. in decreasing order
        expansion_sort = sorted(expansion.items(), key=lambda coeff: abs(coeff[1]),reverse=True)

        maxel = len(expansion_sort)
        # determine the length of the truncated expansion
        cumulative_weight = 0
        for el in range(len(expansion_sort)):
            cumulative_weight += expansion_sort[el][1]**2
            if cumulative_weight > threshold:
                maxel = el+1
                break

        # truncate CI expansion
        expansion_trunc = expansion_sort[0:maxel]

        # turn list into a dictionary
        expansion={}
        for el in expansion_trunc:
            expansion[el[0]] = el[1]

        return expansion

    # ===============================================================================================================

    @staticmethod
    def _correctSigns(psioverlap: np.ndarray, signs: list, nroots: int) -> np.ndarray:

        # first the overlap matrix is modified row-wise to accout for sign changes found in the previous step
        for i in range(nroots):
            if signs[i] == -1:
                for j in range(nroots):
                    psioverlap[i][j] *= -1.
                    signs[i] = 1
        
        # then we check for negtive diagonal elements and correct overlap matrix column-wise
        for i in range(nroots):
            if psioverlap[i][i] < 0.:
                for j in range(nroots):
                    psioverlap[j][i] *= -1.
                    signs[i] = -1
            else:
                signs[i] = 1

        return psioverlap, signs

# ===============================================================================================================

class MolcasTDCoupling:

    def __init__(self, QMPrev: Union[QM, None], QMNow: Union[QM, None], deltaT: float):
        """  Time-derivaive NACs according to the Tully-Hammes-Schaefer scheme.
        formulas from Barbatti, Chem. Phys. 2009, 356, 147

        :param QMPrev: QM instance with the results from the QM calculation at previous time-step
        :param QMNow: QM instance with the results of the QM calculation at present time-step
        """

        # first initialize everything to None
        self.nroots = None
        self.DEarray = None
        self.psioverlap = None
        self.psioverlapPrev = None
        self.tdcouplings = None
        self.signs = None

        # when the given QM is defined, we can compute some information on the states at present time
        if QMNow:

            # store the number of states that is considered here, that is equal to the number of
            # excitations + 1 (the ground state)
            self.nroots = QMNow.outputData.get("nroots")

            # initialize an array with a square matrix with the nr of states as dimension
            self.DEarray = np.zeros((self.nroots, self.nroots))
            energydict = QMNow.outputData.get("energy")
            # fill array with deltaE of the relevant states (all numbers should be positive)
            for i in range(self.nroots):
                for j in range(self.nroots):
                    if j <= i: continue
                    self.DEarray[i][j] = 627.51 * (energydict[j] - energydict[i])
                    self.DEarray[j][i] = self.DEarray[i][j]

        # to construct overlap, both QMNow and QMPrevious need to be defined
        if QMNow and QMPrev:
            self.signs = copy.deepcopy(QMPrev.outputData.get("signs"))
            # at time step 1 set all signs to +1
            if self.signs == []:
                for iroot in range(self.nroots):
                    self.signs.append(1)

            self.psioverlap = np.array(QMNow.outputData.get("psioverlap"))
            self.psioverlapPrev = QMPrev.outputData.get("psioverlap")
            try:
                self.psioverlapPrev == None
            except ValueError:
                self.psioverlapPrev = np.array(self.psioverlapPrev)

            # correct sign of the WF to assure continuity
            self.psioverlap, self.signs = self._correctSigns(self.psioverlap, self.signs, self.nroots)

            # save overlap and signs for evaluation in next step
            QMNow.outputData.set('signs', self.signs)
            QMNow.outputData.set('psioverlap', self.psioverlap)

            # now construct the td couplings 
            if self.psioverlapPrev is None:
                # if there is only one overlap matrix (that is at step 1) compute couplings from 
                # the antisymmetrized projection using hammes-schiffer-tully formula
                logwrt.writelog("Compute couplings using overlap matrix of WFs at the previous 2 steps\n", 2)
                self.tdcouplings = 0.5 * (self.psioverlap - self.psioverlap.T) / deltaT
            else:
                # else use finite differences to compute the couplings using overlap matrices at two consecutive points
                logwrt.writelog("Compute couplings using overlap matrices of WFs at the previous 3 steps\n", 2)
                self.tdcouplings = 0.25 * (3*(self.psioverlap - self.psioverlap.T) - (self.psioverlapPrev - self.psioverlapPrev.T)) / deltaT

            # print the time-derivative couplings matrix
            logwrt.writelog("\nTime-derivative couplings at this step\n", 1)
            logwrt.writelog(logwrt.matrix_prettystring(self.getTDMatrix(), ".8f"), 1)
            logwrt.writelog("\n", 1)

#        # now, for compatibility with what is already implemented in COBRAMM,
#        # call the function that updates some relevant quantities saved in the shelve
#        self.updateShelve()

#    # ===============================================================================================================
#
#    def updateShelve(self):
#
#        sef = shelve.open("cobram-sef")
#
#        # store the number of roots of the electronic structure
#        sef['nroots'] = self.nroots
#
#        # store a matrix of zeros, with nroots as size
#        if self.nroots:
#            sef['zeromat'] = np.zeros((self.nroots, self.nroots))
#        else:
#            sef['zeromat'] = None
#
#        # store the new array with DeltaE in sef['DEarray'], and backup the old one in sef['DE_oldarray']
#        sef['DE_oldarray'] = copy.deepcopy(sef['DEarray'])
#        sef['DEarray'] = self.DEarray
#
#        # store the new array with TDC in sef['TDC'], and backup the old one in sef['TDC_old']
#        sef['TDC_old'] = copy.deepcopy(sef['TDC'])
#        sef['TDC'] = self.getTDMatrix()
#
#        sef.close()

    # ===============================================================================================================

    def getTDMatrix(self):
        """ Returns the TD couplings in form of a numpy 2D array """

        if self.tdcouplings is not None:  # when the TD couplings are defined, give the matrix with the values
            return self.tdcouplings
        else:  # otherwise return an array of zeros
            return np.zeros((self.nroots, self.nroots))

    # ===============================================================================================================

    @staticmethod
    def _correctSigns(psioverlap: np.ndarray, signs: list, nroots: int) -> np.ndarray:

        # first the overlap matrix is modified column-wise to accout for sign changes found in the previous step
        for k in range(nroots):
            if signs[k] == -1:
                for m in range(nroots):
                    psioverlap[m][k] *= -1.
                    signs[k] = 1
        
        # then we check for negtive diagonal elements and correct overlap matrix row-wise
        for k in range(nroots):
            if psioverlap[k][k] < 0.:
                for m in range(nroots):
                    psioverlap[k][m] *= -1.
                    signs[k] = -1
            else:
                signs[k] = 1

        return psioverlap, signs  

