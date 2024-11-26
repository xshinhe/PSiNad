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
import shelve
import os           # filesystem utilities
import numpy as np  # numpy library for scientific computation
import kids_log       # manages log file output + start/end procedures
from Layers import Layers  # object to store information on the geometry and the layers definition
# from harmonicSampling import HarmonicSampling  # class to handle Wigner sampling
# from cobrammDriver import CobrammOutput  
import constants
from pprint import pprint



class QMMM:
    """Object that generate, store and process data for the QMMM energy and gradient"""

    def __init__(self, ks_config, geometry, QM_Results, MM_Results, step, prevx1=None, prevx2=None):
        """Constructor of the QMMM class. The inputs are the data for the QM and MM calculations,
         that are then combined to construct the QM gradients and energies. """

        # initialize attributes of the QMMM class instance
        # QM/MM energies: it is a dictionary with {nstate: energy of the state nstate}
        self.energies = None
        # QM/MM gradient: it is a dictionary with {nstate: gradient of the state nstate}
        self.gradient = None
        # projected gradient for CI optimization
        self.cigradient = None
        # x1 and x2 vectors defining branching plane for CI optimization (only for new branching plane)
        self.x1 = None
        self.x2 = None

        # extract the different parts of the MM_Results
        # Fxyz_second, E_MM_real, E_MM_real_nocharge, Fxyz_modelnoc, Fxyz_modelH, E_MM_modelH = MM_Results
        
        # decide whether to use the correction scheme for the H atom gradient from the QM calculation
        if ks_config.get_nested('QMMM.do_correct', True):
            doCorrect = True
        else:
            doCorrect = False
        # we treat H as virtual atom, so always set doCorrect true
        doCorrect = True

        # =======================================
        # compute the QMMM energies and gradients
        # =======================================

        if "H" not in geometry.calculationType:  # MM only calculation: extract MM results
            gradMM_modelH = MM_Results.grad_modelH
            gradMM_real_nocharge = MM_Results.grad_real_modelnoc
            gradMM_real = MM_Results.grad_real 

            self.energies = {0: MM_Results.energy_real}
            G = [], [], []
            for i, gx, gy, gz in zip(range(len(gradMM_real)), *gradMM_real):
                if i + 1 in geometry.list_MEDIUM_HIGH:
                    G[0].append(gx), G[1].append(gy), G[2].append(gz)
            self.gradient = {0: G}

        elif geometry.calculationType == "H":  # QM only calculation: copy QM results
            self.energies = QM_Results.energydict
            self.gradient = QM_Results.gradientdict
            self.nac = QM_Results.nacdict

        else:  # real QMMM calculation, build energy and gradient for each electronic state with subtractive scheme
            # here E_MM_real_nocharge contains coulumb interaction between M&L, so should be subtracted with QM's selfenergy
            gradMM_modelH = MM_Results.grad_modelH
            gradMM_real_nocharge = MM_Results.grad_real_modelnoc
            gradMM_real = MM_Results.grad_real 

            self.energies = {nstate: eqm + MM_Results.energy_real_modelnoc - MM_Results.energy_modelH - QM_Results.selfenergy
                             for nstate, eqm in QM_Results.energydict.items()}
            self.gradient = {}
            for nstate in QM_Results.gradientdict:
                # when some piece is missing, do not compute the QMMM, instead just set it to None
                if QM_Results.gradient(nstate) is None or \
                        (QM_Results.gradcharges(nstate) is None and "M" in geometry.calculationType):
                    gradQMMM = None
                else:
                    gradQMMM = subtractiveScheme(ks_config, geometry, QM_Results.gradient(nstate), gradMM_modelH,
                                                 gradMM_real_nocharge, QM_Results.gradcharges(nstate),
                                                 correctAtomLink=doCorrect)
                # store gradient
                self.gradient[nstate] = gradQMMM

            # add QMMM nac = QM nac + fullnaccharge
            self.nac = {}
            for istate in QM_Results.nacdict:
                if not istate in self.nac:
                    self.nac[istate] = {}
                for kstate in QM_Results.nacdict[istate]:
                    nacQM = QM_Results.nacdict[istate][kstate]
                    if QM_Results.nacCRGdict is not None:
                        nacCRG = QM_Results.nacCRGdict[istate][kstate]
                        nacQMMM = ReduceNacScheme(ks_config, geometry, nacQM, nacCRG, correctAtomLink=doCorrect)
                    else:
                        nacQMMM = nacQM
                    
                    self.nac[istate][kstate] = nacQMMM

        # now we can proceed with additional operations for CI gradient computation

        if ks_config.args.type == 'ci':

            # extract data from shelve
            sef = shelve.open("cobram-sef", 'r')
            state = sef['state']
            DEarray = sef['DEarray']
            #NAC = sef['NAC']
            sef.close()
            #if we use traditional branching plane
            if ks_config.get_nested('QM.bp_type', 'gmean') == 'nacs':
                NAC = QM_Results.nac(state-1, state)
    
                # reshape NAC and add zeros for the MEDIUM atoms
                nac = []
                for i in range(geometry.NatomHM):
                    if geometry.list_MEDIUM_HIGH[i] in geometry.list_HIGH:
                        ind = list(geometry.list_HIGH).index(geometry.list_MEDIUM_HIGH[i])
                        nac.append([NAC[0][ind], NAC[1][ind], NAC[2][ind]])
                    elif geometry.list_MEDIUM_HIGH[i] in geometry.list_MEDIUM:
                        nac.append([0.0, 0.0, 0.0])
                    else:
                        Log.fatalError('Atom ' + str(geometry.list_MEDIUM_HIGH[i]) + ' not found among HL and ML atoms!')
                nac = np.array(nac)

            # compute energy difference between states in atomic units
            deltaE = self.energies[state] - self.energies[state - 1]

            # compute gradient with Bearpark's Gradient Projection method
            gradientlow = self.gradient[state-1]
            gradienthigh = self.gradient[state]
            #if we use traditional branching plane
            if ks_config.get_nested('QM.bp_type', 'gmean') == 'nacs':

                self.cigradient = BearparkGradientProjection(geometry.NatomHM, gradientlow, gradienthigh, nac, deltaE)
            else:
                self.cigradient, self.x1, self.x2 = BearparkGradientProjection_newBP(geometry.NatomHM, gradientlow, gradienthigh, prevx1, prevx2, deltaE, step)

    # =============================================================================================================

    def getenergy(self, nstate=None):
        """Return a list containing the QMMM electronic state energy (or energies)"""

        try:  # exception handler, return None if that state requested from input is not available
            # when nstate is None, no state is requested... return the full list in order of increasing energy
            if nstate is None:
                return [self.energies[n] for n in sorted(self.energies.keys())]
            # otherwise return the energy of the state requested as input
            else:
                return self.energies[nstate]
        except (KeyError, TypeError):
            return None

    # =============================================================================================================

    def getgradient(self, nstate=None):
        """Return the gradient on the QM part (or gradients) computed with the QM calculation """

        try:  # exception handler, return None if that state requested from input is not available
            # when nstate is None, no state is requested... return the full list in order of increasing energy
            if nstate is None:
                return [self.gradient[n] for n in sorted(self.gradient.keys())]
            # otherwise return the gradient of the state requested as input
            else:
                return self.gradient[nstate]
        except (KeyError, TypeError):
            return None

    # =============================================================================================================

    def sharcQMMMoutfile(self, QMoutText, geometry):
        """After computing the QMMM energies, this function uses the information from a QM.out file (sharc)
        given as input and updates it with the results obtained from QMMM. """

        def FORTRANformat(f, prec, exp_digits):
            """This function prints numbers in FORTRAN scientific notation"""
            s = ("{0:."+str(prec)+"E}").format(f)
            mantissa, exp = s.split('E')
            if mantissa[0] not in ["+", "-"]:
                mantissa = " " + mantissa
            return "{0}E".format(mantissa)+('{0:+0'+str(exp_digits+1)+'d}').format(int(exp))

        # initialize the text of the QMMM.out file
        QMMMoutText, i = "", 0

        # loop over the lines of the QM.out file
        linebyline = QMoutText.splitlines()
        while i < len(linebyline):

            # the line is the starting line of a Hamiltonian section
            if linebyline[i].split() and linebyline[i].split()[0] == "!" and linebyline[i].split()[1] == "1":
                QMMMoutText += linebyline[i] + "\n"

                nstates = int(linebyline[i+1].split()[0])
                QMMMoutText += linebyline[i+1] + "\n"

                for istate in range(nstates):
                    energies = [float(e) for e in linebyline[i+2+istate].split()]
                    energies[2*istate] = self.getenergy(istate)
                    for e in energies:
                        QMMMoutText += "{0} ".format(FORTRANformat(e, 9, 3))
                    QMMMoutText += "\n"

                i += nstates + 2

            # the line is the starting line of a gradient section
            if linebyline[i].split() and linebyline[i].split()[0] == "!" and linebyline[i].split()[1] == "3":
                # write header line with different number of gradient size
                QMMMoutText += linebyline[i].split("x")[0] + "x{0}x".format(geometry.atomNum) + \
                               linebyline[i].split("x")[2] + "\n"
                i += 1
                # loop over the states
                for istate in self.gradient:
                    noldgradient = int(linebyline[i].split()[0])
                    QMMMoutText += "{0} ".format(geometry.atomNum) + linebyline[i].split(" ", 1)[1] + "\n"
                    if self.gradient[istate] is not None:
                        for iatom, atomgrad in enumerate(zip(*self.gradient[istate])):
                            if iatom >= geometry.NatomHM: break
                            for g in atomgrad:
                                QMMMoutText += "{0} ".format(FORTRANformat(g, 9, 3))
                            QMMMoutText += "\n"
                        for iatom in range(geometry.atomNum-geometry.NatomHM):
                            QMMMoutText += ("{0} ".format(FORTRANformat(0.0, 9, 3)))*3 + "\n"
                    else:
                        for iatom in range(geometry.atomNum):
                            QMMMoutText += ("{0} ".format(FORTRANformat(0.0, 9, 3))) * 3 + "\n"
                    i += noldgradient + 1

            # in all the other cases, just copy the line
            else:
                QMMMoutText += linebyline[i] + "\n"
                i += 1

        return QMMMoutText


###################################################################################################################

def ReduceNacScheme(ks_config, geometry, nacQM, nacCRG, correctAtomLink=True):
    x, y, z = [], [], []
    jMedium, jHigh = 0, 0  # counters for the list of MEDIUM and HIGH layer atoms

    # loop over the full real geometry, HIGH, MEDIUM and LOW atoms
    for i in range(geometry.atomNum):

        # the atom (index i+1) is in the HIGH LAYER
        if i+1 in geometry.list_HIGH:
            # gradient of the QM part
            if ks_config.get_nested('QMMM.free_QM', True): 
                grad = [nacQM[ix][jHigh] for ix in range(3)]

                # when requested, correct the gradient of the link atoms by distributing the H atom gradient
                if correctAtomLink:
                    # select the atom links in which the a atom is the i+1 HIGH atom
                    for ilink, b, a in zip(range(len(geometry.atomLink_BA[0])), *geometry.atomLink_BA):
                        if i+1 == a:
                            # compute the bond length ration that is necessary to evaluate dR(L)/d(R(Q))
                            g_jac = geometry.hbondratio(geometry.atomLabel[a-1])
                            # correct Q1 by dR(L)/d(R(Q))
                            jHAtom = geometry.NatomQM + ilink  # position of the H atom in the modelH gradient
                            grad = [grad[ix] + nacQM[ix][jHAtom] * (1 - g_jac) for ix in range(3)]
            else:
                #freeze HighLayer (experimental)
                grad = [0, 0, 0]
            # store the gradient computed for this HIGH atom
            x.append(grad[0]), y.append(grad[1]), z.append(grad[2])
            # increment the counter of the HIGH atoms
            jHigh += 1

        # the atom (index i+1) is in the MEDIUM LAYER
        if i+1 in geometry.list_MEDIUM:

            # gradient of the MM part
            grad = [nacCRG[ix][jMedium] for ix in range(3)]

            # when requested, correct the gradient of the link atoms by distributing the H atom gradient
            if ks_config.get_nested('QMMM.free_QM', True) and correctAtomLink:
                # select the atom links in which the a atom is the i+1 HIGH atom
                for ilink, b, a in zip(range(len(geometry.atomLink_BA[0])), *geometry.atomLink_BA):
                    if i+1 == b:
                        # compute the bond length ration that is necessary to evaluate dR(L)/d(R(Q))
                        g_jac = geometry.hbondratio(geometry.atomLabel[a-1])
                        # correct Q1 by dR(L)/d(R(Q))
                        jHAtom = geometry.NatomQM + ilink  # position of the H atom in the modelH gradient
                        grad = [grad[ix] + nacQM[ix][jHAtom] * g_jac for ix in range(3)]

            # store the gradient computed for this HIGH atom
            x.append(grad[0]), y.append(grad[1]), z.append(grad[2])
            # increment the counter of the HIGH atoms
            jMedium += 1

    # now x,y and z contain the components of the HIGH-MEDIUM QM-MM gradient
    return x, y, z


def subtractiveScheme(ks_config, geometry, gradQM_modelH, gradMM_modelH, gradMM_real_nocharge, gradQM_coulomb,
                      correctAtomLink=True):

    # compute the QM components for the modelH part, by taking the difference between the
    # QM gradient and the MM gradient

    grad_correct_modelH = np.array(gradQM_modelH) - np.array(gradMM_modelH)

    # now process atom links to redistribute the component of the derivative according to the chain rule dR(L)/d(R(Q))
    # the resulting gradient will be defined only for HIGH and MEDIUM atoms, the gradient of LOW atoms is discarded

    x, y, z = [], [], []
    jMedium, jHigh = 0, 0  # counters for the list of MEDIUM and HIGH layer atoms

    # loop over the full real geometry, HIGH, MEDIUM and LOW atoms
    for i in range(geometry.atomNum):

        # the atom (index i+1) is in the HIGH LAYER
        if i+1 in geometry.list_HIGH:

            # gradient of the QM part
            if ks_config.get_nested('QMMM.free_QM', True): 
                grad = [gradMM_real_nocharge[ix][i] + grad_correct_modelH[ix][jHigh] for ix in range(3)]

                # when requested, correct the gradient of the link atoms by distributing the H atom gradient
                if correctAtomLink:
                    # select the atom links in which the a atom is the i+1 HIGH atom
                    for ilink, b, a in zip(range(len(geometry.atomLink_BA[0])), *geometry.atomLink_BA):
                        if i+1 == a:
                            # compute the bond length ration that is necessary to evaluate dR(L)/d(R(Q))
                            g_jac = geometry.hbondratio(geometry.atomLabel[a-1])
                            # correct Q1 by dR(L)/d(R(Q))
                            jHAtom = geometry.NatomQM + ilink  # position of the H atom in the modelH gradient
                            grad = [grad[ix] + grad_correct_modelH[ix][jHAtom] * (1 - g_jac) for ix in range(3)]
            else:
                #freeze HighLayer (experimental)
                grad = [0, 0, 0]
            # store the gradient computed for this HIGH atom
            x.append(grad[0]), y.append(grad[1]), z.append(grad[2])
            # increment the counter of the HIGH atoms
            jHigh += 1

        # the atom (index i+1) is in the MEDIUM LAYER
        if i+1 in geometry.list_MEDIUM:

            # gradient of the MM part
            grad = [gradMM_real_nocharge[ix][i] + gradQM_coulomb[ix][jMedium] for ix in range(3)]

            # when requested, correct the gradient of the link atoms by distributing the H atom gradient
            if ks_config.get_nested('QMMM.free_QM', True) and correctAtomLink:
                # select the atom links in which the a atom is the i+1 HIGH atom
                for ilink, b, a in zip(range(len(geometry.atomLink_BA[0])), *geometry.atomLink_BA):
                    if i+1 == b:
                        # compute the bond length ration that is necessary to evaluate dR(L)/d(R(Q))
                        g_jac = geometry.hbondratio(geometry.atomLabel[a-1])
                        # correct Q1 by dR(L)/d(R(Q))
                        jHAtom = geometry.NatomQM + ilink  # position of the H atom in the modelH gradient
                        grad = [grad[ix] + grad_correct_modelH[ix][jHAtom] * g_jac for ix in range(3)]

            # store the gradient computed for this HIGH atom
            x.append(grad[0]), y.append(grad[1]), z.append(grad[2])
            # increment the counter of the HIGH atoms
            jMedium += 1

    # now x,y and z contain the components of the HIGH-MEDIUM QM-MM gradient
    return x, y, z


###################################################################################################################

def BearparkGradientProjection(ncoords, gradientlow, gradienthigh, nac, deltaE):
    """Computing the effective gradient for optimizing CIs"""

    # initialize gradient components
    x, y, z = [], [], []

    # reshape gradientlow and gradienthigh
    s0 = np.array(list(zip(*gradientlow))[0:ncoords])
    s1 = np.array(list(zip(*gradienthigh))[0:ncoords])

    # computation gradient difference (GD, x1)
    x1 = s1 - s0
    # compute the norm of the x1 vector
    x1v = np.reshape(x1, 3 * ncoords)
    x1v = x1v / np.linalg.norm(x1v)

    # computation of the scaled derivative coupling (DC, x2)
    # scaling of the DC with the energy difference (E_high - E_low)
    x2 = deltaE * nac

    # compute g1 = 2 * (E_high - E_low) * GD(normalized)
    g1 = 2 * deltaE * x1v

    # orthogonalize x1 and x2
    identity = np.eye(3 * ncoords)
    P = identity - np.outer(x1v, x1v)
    x2v = np.reshape(x2, 3 * ncoords)
    x2v = x2v / np.linalg.norm(x2v)
    Px2v = np.dot(P, x2v)
    Px2v = Px2v / np.linalg.norm(Px2v)

    # compute g2 = P(grad_high) with P the projector orthogonal to the branching space
    # P = I - A * (A_T * A)^(-1) * A_T
    A = np.column_stack((x1v, Px2v))
    A_T = np.transpose(A)
    ATA = np.dot(A_T, A)
    ATAinv = np.linalg.inv(ATA)
    AATAinv = np.dot(A, ATAinv)
    P = identity - np.dot(AATAinv, A_T)
    s1v = np.reshape(s1, 3 * ncoords)
    g2 = np.dot(P, s1v)

    # final gradient G = g1 + g2
    G = g1 + g2

    # store final gradient G to x,y,z
    G_reshaped = np.reshape(G, (-1, 3))
    for i in range(ncoords):
        x.append(G_reshaped[i][0])
        y.append(G_reshaped[i][1])
        z.append(G_reshaped[i][2])

    # return the gradient to output
    return x, y, z

###################################################################################################################

def BearparkGradientProjection_newBP(ncoords, gradientlow, gradienthigh, prevx1, prevx2, deltaE, step):
    """Computing the effective gradient for optimizing CIs with NAC-free branching plane (J. Chem. Theory Comput., Vol. 6, No. 5, 2010)"""

    # initialize gradient components
    x, y, z = [], [], []

    # reshape gradientlow and gradienthigh
    s0 = np.array(list(zip(*gradientlow))[0:ncoords])
    s1 = np.array(list(zip(*gradienthigh))[0:ncoords])

    # computation gradient difference (GD, x1) 
    x1 = s1 - s0
    # compute the norm of the x1 vector
    x1v = np.reshape(x1, 3 * ncoords)
    x1v = x1v / np.linalg.norm(x1v)

    # computation of the x2 vector:
    if step == 0:
        #x2 = mean of gradients for step 0
        x2 = (s1 + s0)/2
    else:
        x2 = ((np.dot(prevx2,x1v)*prevx1) - (np.dot(prevx1,x1v)*prevx2)) / (np.sqrt((np.dot(prevx2,x1v))**2 + (np.dot(prevx1,x1v))**2))

    # compute g1 = 2 * (E_high - E_low) * GD(normalized)
    g1 = 2 * deltaE * x1v

    identity = np.eye(3 * ncoords)

    if step == 0:
        # orthogonalize x1 and x2
        P = identity - np.outer(x1v, x1v)
        x2v = np.reshape(x2, 3 * ncoords)
        x2v = x2v / np.linalg.norm(x2v)
        Px2v = np.dot(P, x2v)
        Px2v = Px2v / np.linalg.norm(Px2v)
    else:
        Px2v = x2

    # compute g2 = P(grad_high) with P the projector orthogonal to the branching space
    # P = I - A * (A_T * A)^(-1) * A_T
    A = np.column_stack((x1v, Px2v))
    A_T = np.transpose(A)
    ATA = np.dot(A_T, A)
    ATAinv = np.linalg.inv(ATA)
    AATAinv = np.dot(A, ATAinv)
    P = identity - np.dot(AATAinv, A_T)
    s1v = np.reshape(s1, 3 * ncoords)
    g2 = np.dot(P, s1v)

    # final gradient G = g1 + g2
    G = g1 + g2

    # store final gradient G to x,y,z
    G_reshaped = np.reshape(G, (-1, 3))
    for i in range(ncoords):
        x.append(G_reshaped[i][0])
        y.append(G_reshaped[i][1])
        z.append(G_reshaped[i][2])

    # return the gradient to output as well as x1 and x2 vectors to save for next step
    return [x, y, z], x1v, Px2v

###################################################################################################################

def GradientProjectionOnNormalModes(geometry, ks_config, gradient):
    """Compute the projection of the gradient on the set of normal modes below a frequency threshold defined by ks_config[41] """
    
    freqdir = os.getcwd() 
    # parse the results of the frequency calculation run using CobrammOutput
    freqresults = None #CobrammOutput(freqdir, store_files=True)

    # extract equilibrium geometry and format as a simple 1D vector of 3N elements
    geomvector_1D = []
    #for x, y, z in zip(*geometry.model):
    for x, y, z in zip(*geometry.MEDIUM_HIGH):
        geomvector_1D.append(x), geomvector_1D.append(y), geomvector_1D.append(z)

    # calculate center of mass
    massvec = np.array(freqresults.coord_masses[::3])
    com = np.array([0, 0, 0])
    mass = sum(freqresults.coord_masses)/3
    for i in range(geometry.NatomHM):
        com[0] = sum(np.array([x * m for x, m in zip(geomvector_1D, freqresults.coord_masses)])[0::3])
        com[1] = sum(np.array([x * m for x, m in zip(geomvector_1D, freqresults.coord_masses)])[1::3])
        com[2] = sum(np.array([x * m for x, m in zip(geomvector_1D, freqresults.coord_masses)])[2::3])
    com = com / mass

    # shift coordinates to center of mass (3xN array)
    geomvector_com = []
    geomvector_com.append(np.array(geomvector_1D)[0::3] - com[0])
    geomvector_com.append(np.array(geomvector_1D)[1::3] - com[1])
    geomvector_com.append(np.array(geomvector_1D)[2::3] - com[2])
    geomvector_com_1D = []
    for x, y, z in zip(*geomvector_com):
        geomvector_com_1D.append(x), geomvector_com_1D.append(y), geomvector_com_1D.append(z)

    # calculate moment of inertia
    moimat = np.zeros((3,3))
    moimat[0,0] = np.dot(massvec,geomvector_com[1]**2 + geomvector_com[2]**2)
    moimat[1,1] = np.dot(massvec,geomvector_com[0]**2 + geomvector_com[2]**2)
    moimat[2,2] = np.dot(massvec,geomvector_com[0]**2 + geomvector_com[1]**2)
    for i in range(3):
        for j in range(i+1,3):
            moimat[i,j] = -np.dot(massvec,geomvector_com[i]*geomvector_com[j])
            moimat[j,i] = moimat[i,j]

    # diagonalize moment of inertia matrix
    eigenVal, eigenVec = np.linalg.eig(moimat)
    eigenVec = np.transpose(eigenVec)

    # create translation and rotation vectors in mass-weighted coordinates
    # the vectors should be orthonormalized
    transX,transY,transZ,rot1,rot2,rot3 = [np.zeros(3*geometry.NatomHM) for i in range(6)]
    transX[0::3] = np.sqrt(freqresults.coord_masses[0::3])
    transY[1::3] = np.sqrt(freqresults.coord_masses[1::3])
    transZ[2::3] = np.sqrt(freqresults.coord_masses[2::3])
    norm = np.linalg.norm(transX)
    transX = transX/norm; transY = transY/norm; transZ = transZ/norm

    for i in range(geometry.NatomHM):
        for j in range(3):
            rot1[3*i+j] = ( np.dot(geomvector_com_1D[3*i:3*i+3],eigenVec[1]) * eigenVec[j][2] - np.dot(geomvector_com_1D[3*i:3*i+3],eigenVec[2]) * eigenVec[j][1] ) * np.sqrt(massvec[i])
            rot2[3*i+j] = ( np.dot(geomvector_com_1D[3*i:3*i+3],eigenVec[2]) * eigenVec[j][0] - np.dot(geomvector_com_1D[3*i:3*i+3],eigenVec[0]) * eigenVec[j][2] ) * np.sqrt(massvec[i])
            rot3[3*i+j] = ( np.dot(geomvector_com_1D[3*i:3*i+3],eigenVec[0]) * eigenVec[j][1] - np.dot(geomvector_com_1D[3*i:3*i+3],eigenVec[1]) * eigenVec[j][0] ) * np.sqrt(massvec[i])
    norm = np.linalg.norm(rot1); rot1 = rot1/norm
    norm = np.linalg.norm(rot2); rot2 = rot2/norm
    norm = np.linalg.norm(rot3); rot3 = rot3/norm

    # generate 3N-6 vectors orthogonal to transX-Z and rot1-3 through Gram-Schmidt orthogonalization
    cart2intmat = []
    cart2intmat.append(transX); cart2intmat.append(transY); cart2intmat.append(transZ);
    cart2intmat.append(rot1); cart2intmat.append(rot2); cart2intmat.append(rot3);
    while len(cart2intmat) < 3*geometry.NatomHM:
        # generate random vector
        randvec = np.random.rand(3*geometry.NatomHM)
        norm = np.linalg.norm(randvec); randvec = randvec/norm
        # subtract iteratively the projection of the previous vectors from randvec
        for j in range(len(cart2intmat)):
            proj = np.dot(cart2intmat[j],randvec)
            randvec = randvec - proj*cart2intmat[j]
        # the residual is the new orthogonal vector
        norm = np.linalg.norm(randvec); randvec = randvec/norm
        if norm < 0.005:
            continue
        # add new vector to list of arrays
        cart2intmat.append(randvec)

    # define harmonic sampling of the molecular oscillations
    mol_oscillator = None # HarmonicSampling(geometry, geomvector_1D, freqresults.coord_masses, freqresults.force_matrix, cart2intmat, 0, -1)

    projGradient = np.zeros(3*geometry.NatomHM) 
    # QMMM gradient is 3xN matrix where x = gradient[0], y = gradient[1], z = gradient[2]
    # eigenVecCart is a 3xN array where x = eigenVecCart[0::3], y = eigenVecCart[1::3], z = eigenVecCart[2::3]
    # reshape QMMM gradient to eigenVecCart shape
    gradient = np.reshape(gradient, 3*geometry.NatomHM, order='F')
    for i in range(len(mol_oscillator.eigenVal)):
        if mol_oscillator.eigenVal[i] > 0. and mol_oscillator.eigenVal[i] < ( float(ks_config[40]) * constants.wavnr2au ) ** 2 and i == 4:
            #print("Project mode ", i)
            a = np.dot(gradient,mol_oscillator.eigenVecCart[i])/np.linalg.norm(mol_oscillator.eigenVecCart[i])
            #print("projection = ", a)
            projGradient = projGradient + a*mol_oscillator.eigenVecCart[i]/np.linalg.norm(mol_oscillator.eigenVecCart[i])
    a = np.dot(projGradient,gradient)
    norm = np.linalg.norm(gradient)
    #print("Projection on sub-space below {0:6.1f} = {1:5.1f}%\n".format(float(ks_config[40]),100*a/norm))
    # reshape projGradient to shape of gradient
    projGradient = np.reshape(projGradient, (3,geometry.NatomHM), order='F' )
    return projGradient 
