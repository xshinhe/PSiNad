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

import math  # mathematical functions
import sys  # system-specific parameters and functions
import os  # miscellaneous operating system interfaces
import re  # regular expressions

# math libraries

import numpy as np  # numpy library for scientific computation

# import local modules and classes

import constants  # values of physical constants and conversion factors
import logwrt  # manages log file output + start/end procedures
import subprocess  # run external program as child process

#####################################################################################################

# dictionaries to map the S, P, D, F notation to the values of the exponents
# in the cartesian functions

CARTESIAN_SHELLS = {"S": [(0, 0, 0)],
                    "P": [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
                    "D": [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)],
                    "F": [(3, 0, 0), (0, 3, 0), (0, 0, 3), (1, 2, 0), (2, 1, 0), (2, 0, 1),
                          (1, 0, 2), (0, 1, 2), (0, 2, 1), (1, 1, 1)]}

CARTESIAN_LABELS = {"S": (0, 0, 0), "PX": (1, 0, 0), "PY": (0, 1, 0), "PZ": (0, 0, 1),
                    "DXX": (2, 0, 0), "DYY": (0, 2, 0), "DZZ": (0, 0, 2),
                    "DXY": (1, 1, 0), "DXZ": (1, 0, 1), "DYZ": (0, 1, 1),
                    "FXXX": (3, 0, 0), "FYYY": (0, 3, 0), "FZZZ": (0, 0, 3),
                    "FXYY": (1, 2, 0), "FXXY": (2, 1, 0), "FXXZ": (2, 0, 1),
                    "FXZZ": (1, 0, 2), "FYZZ": (0, 1, 2), "FYYZ": (0, 2, 1), "FXYZ": (1, 1, 1)}

# dictionaries to map the S, P, D, F notation to the values of the angular momentum
# of the spherical harmonics

SPHARM_SHELLS = {"S": [(0, 0)],
                 "P": [(1, 0), (1, 1), (1, -1)],
                 "D": [(2, 0), (2, 1), (2, -1), (2, 2), (2, -2)],
                 "F": [(3, 0), (3, 1), (3, -1), (3, 2), (3, -2), (3, 3), (3, -3)]}

SPHARM_LABELS = {"S": (0, 0), "P0": (1, 0), "P+1": (1, 1), "P-1": (1, -1),
                 "D0": (2, 0), "D+1": (2, 1), "D-1": (2, -1), "D+2": (2, 2), "D-2": (2, -2),
                 "F0": (3, 0), "F+1": (3, 1), "F-1": (3, -1), "F+2": (3, 2), "F-2": (3, -2),
                 "F+3": (3, 3), "F-3": (3, -3)}

# dictionaries with the mapping from cartesian to spherical GTOs

CARTESIAN_SPHARM_TRANS = {
    # D FUNCTIONS - sherical label (l, m): {cartesian label (lx, ly, lz): coeff of the expansions)}
    (2, 0):  {(2, 0, 0): -0.5, (0, 2, 0): -0.5, (0, 0, 2): 1.0},
    (2, +1): {(1, 0, 1): 1.0},
    (2, -1): {(0, 1, 1): 1.0},
    (2, +2): {(2, 0, 0): math.sqrt(3./4.), (0, 2, 0): -math.sqrt(3./4.)},
    (2, -2): {(1, 1, 0): 1.0},
    # F FUNCTIONS - sherical label (l, m): {cartesian label (lx, ly, lz): coeff of the expansions)}
    (3, 0):  {(0, 0, 3): 1., (2, 0, 1): -3./(2.*math.sqrt(5.)), (0, 2, 1): -3./(2.*math.sqrt(5.))},
    (3, +1): {(3, 0, 0): -math.sqrt(3./2.)/2., (1, 2, 0): -math.sqrt(3./10.)/2., (1, 0, 2): math.sqrt(6./5.)},
    (3, -1): {(0, 3, 0): -math.sqrt(3./2.)/2., (2, 1, 0): -math.sqrt(3./10.)/2., (0, 1, 2): math.sqrt(6./5.)},
    (3, +2): {(2, 0, 1): math.sqrt(3./4.), (0, 2, 1): -math.sqrt(3./4.)},
    (3, -2): {(1, 1, 1): 1.},
    (3, +3): {(3, 0, 0): math.sqrt(5./2.)/2., (1, 2, 0): -3./(2.*math.sqrt(2.))},
    (3, -3): {(0, 3, 0): -math.sqrt(5./2.)/2., (2, 1, 0): 3./(2.*math.sqrt(2.))},
}


#####################################################################################################

class AtomicOrbital:
    """ This class defines the structure that holds the data of a single contracted atomic orbital,
    with given: atomic center, radial part (either cartesian function or spherical harmonic),
    lists of exponential and contraction coefficients
    The class also defines a method to compute the overlap integral between two given generic orbitals """

    def __init__(self, center, alphaList, coeffList, cartesian, radialNumbers):

        # store a copy of the input values
        self.center = center
        self.alphaList = alphaList
        self.coeffList = coeffList
        self.cartesian = cartesian
        self.radialNumbers = radialNumbers

        # compute normalization factors
        if self.cartesian:  # when the angular part is cartesian, it is required to compute normaliz coeffs
            self.normalization = [normalization(alpha, self.radialNumbers) for alpha in self.alphaList]

    ###########################################################################################

    def sphericalExpansion(self):
        """If the AO defined in self is spherical, return the expansion of the orbital
        in terms of the corresponding cartesian orbitals, with the same type of contraction scheme"""

        if self.cartesian:
            # in this case the expansion is trivial, it is just the function itself
            expcoeffs, cartesorbs = [1.0], [self]

        else:
            # in this case construct the lists according to the appropriate linear transformation
            expcoeffs, cartesorbs = [], []
            # define the letter for this type of orbital
            orbtype = AtomicOrbital._AOLabel(self.radialNumbers, self.cartesian)[0]
            # now create the lists of coefficients and atomic orbitals, with a loop on a single row of the transf.
            for cartNr in CARTESIAN_SHELLS[orbtype]:
                # since the matrix is sparse, only nonzero elements are defined, so check when the element is present
                if cartNr in CARTESIAN_SPHARM_TRANS[self.radialNumbers]:
                    expcoeffs.append(CARTESIAN_SPHARM_TRANS[self.radialNumbers][cartNr])
                    cartesorbs.append(AtomicOrbital(self.center, self.alphaList, self.coeffList,
                                                    True, radialNumbers=cartNr))

        # return the expansion coefficients and the cartesian orbitals
        return zip(expcoeffs, cartesorbs)

    ###########################################################################################

    def prettyLog(self):
        """Return a string with a short formatted description of the atomic orbital data"""

        string = "Atomic Orbital {0}\n".format(AtomicOrbital._AOLabel(self.radialNumbers, self.cartesian))
        string += "center: {0:6.3f} {1:6.3f} {2:6.3f}".format(*self.center)
        string += "   contract. degree: {0}\n".format(len(self.coeffList))
        string += "contraction coeff.s: {0}\n".format(self.coeffList)
        string += "contraction alpha's: {0}\n".format(self.alphaList)
        if self.cartesian:
            string += "cartesian function with exponents: {0} {1} {2}\n".format(*self.radialNumbers)
            string += "normalization coeff.s: {0}\n".format(self.normalization)

        return string

    ###########################################################################################

    @staticmethod
    def AOoverlap(AO1, AO2):
        """Compute the overlap between two atomic orbitals defined as AO1 and AO2, instantes
        of the AtomicOrbital object. The function combines the overlap integrals between the
        primitive gaussian functions, computed via the primitiveOverlap() function, using
        the coefficients given by the contraction scheme. """

        # initialize the overlap between the orbitals
        orbOvlp = 0.0

        # for spherical harmonics, expand in terms of cartesian and call this function recursively
        if not AO1.cartesian or not AO2.cartesian:
            for cartCoeff1, cartAO1 in AO1.sphericalExpansion():
                for cartCoeff2, cartAO2 in AO2.sphericalExpansion():
                    orbOvlp += cartCoeff1 * AtomicOrbital.AOoverlap(cartAO1, cartAO2) * cartCoeff2

        else:  # when the AOs are cartesian, simply use the primitive overlap formula taking care of the contraction
            # loop over the contractions of AO1 and AO2
            for c1, alpha1, N1 in zip(AO1.coeffList, AO1.alphaList, AO1.normalization):
                for c2, alpha2, N2 in zip(AO2.coeffList, AO2.alphaList, AO2.normalization):
                    # compute the overlap between the two primitive functions
                    primOvlp = _primitiveOverlap(alpha1, alpha2, AO1.center, AO2.center,
                                                 AO1.radialNumbers, AO2.radialNumbers)
                    # now add to the overlap summation
                    orbOvlp += N1 * c1 * primOvlp * c2 * N2

        return orbOvlp

    ###########################################################################################

    @staticmethod
    def _AOLabel(inpNumbers, cartesian):
        """Given a set of numbers labelling a single basis functions, returns the letter that labels
        the shell and the type of that basis function (i.e. "S","PX","D-1", ... )
        """

        if cartesian:  # for a cartesian function, search in the CARTESIAN_LABELS dictionary
            labDict = CARTESIAN_LABELS
        else:  # for spherical harmonics, use the SPHARM_LABELS dictionary
            labDict = SPHARM_LABELS

        # search inpNumbers inside the dictionary
        for name, labelNumbers in labDict.items():
            if labelNumbers == inpNumbers:
                return name  # return the label

        # if still inside the function, it means that the numbers were not found, return None
        return None


#####################################################################################################

class Shell:
    """ This class defines the structure that holds the data for defining a shell of atomic orbitals (AO).

    The attributes of this objects are:
    - idAtom: the indentifier of the atom of the guassian function shell
    - center: the cartesian coordinates (in a.u.) of the center of the gaussian functions, usually
              equal to the coordinates of the atom labelled by idAtom
    - nContract: contraction degree of the gaussian basis function shell
    - alphaList: list of the exponential parameters (in a.u.) for the gaussian contraction
    - coeffList: list of the expansion parameters for the gaussian contraction
    - cartesian: logical variable, True for cartesian functions, False for spherical harmonics
    - nFunctions: number of gaussian basis functions contained in the shell
    - radNumbers: list of the numbers labelling the radial part of the basis function, i.e.
                  (l, ml) for spherical harmonics, (nx,ny,nz) for cartesian functions
    - atomicOrbitals: list of the AOs of the shell, defined each as an instance of the AtomicOrbital class

    The class is initialized by giving:
    - idAtom: the indentifier of the atom of the guassian function shell
    - center:  the cartesian coordinates (in a.u.) of the center of the gaussian functions
    - angMom: a label that defines the angular momentum of the shell function ("S","P","D" ... )
    - cartesian: logical variable, True for cartesian functions, False for spherical harmonics
    - alphaList: list of the exponential parameters (in a.u.) for the gaussian contraction
    - coeffList: list of the expansion parameters for the gaussian contraction

    """

    def __init__(self, idAtom, center, angMom, cartesian, alphaList, coeffList):
        """construct the shell"""

        # these two numbers identify the atom of the shell
        self.idAtom = idAtom
        # store the center of the AO shell, that should be the center of the Atom idAtom
        self.center = center

        # store the degree of contraction of the AO shell
        self.nContract = len(alphaList)

        # store the list of exponential and expansion coefficients of the shell contraction
        self.alphaList = alphaList
        if len(coeffList) == self.nContract:
            self.coeffList = coeffList
        else:
            logwrt.fatalerror("contraction degree mismatch in atomic orbital shell definition\n"
                              "Please check the lists of exponential and expansion coeff.s of the shell contraction")

        # store the type of radial function (cartesian/spherica harm)
        # for S and P functions, cartesian functions are always used, regardless of the input cartesian value
        if angMom.upper() in ["S", "P"]:
            self.cartesian = True
        else:
            self.cartesian = cartesian

        # define according to angMom the number of functions and the list of "quantum numbers" for the radial part
        # Note that the order of the functions matters! We have assumed the same ordering that is
        # used in Gaussian output:
        # - for spherical harmonics: 0, +1, -1, +2, -2, +3, -3 .. etc etc
        # - for cartesian functions: X,Y,Z, | XX,YY,ZZ,XY,XZ,YZ | XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ
        try:
            if self.cartesian:
                self.radNumbers = CARTESIAN_SHELLS[angMom.upper()]
            else:
                self.radNumbers = SPHARM_SHELLS[angMom.upper()]
        except:
            logwrt.fatalerror("Cannot construct a shell of type {}.\n".format(angMom.upper()) +
                              "Only S,P,D,F functions are currently implemented.")
        # define the number of atomic orbitals contained in the shell
        self.nFunctions = len(self.radNumbers)

    ###########################################################################################

    @staticmethod
    def __angMomLabel__(inpNumbers, cartesian):
        """Given a set of numbers labelling a single basis functions, returns the letter that labels
        the shell of that basis function (i.e. "S","P","D" ... )
        """

        if cartesian:  # for a cartesian function, search in the CARTESIAN_LABELS dictionary
            labDict = CARTESIAN_LABELS
        else:  # for spherical harmonics, use the SPHARM_LABELS dictionary
            labDict = SPHARM_LABELS

        # search inpNumbers inside the dictionary
        for name, labelNumbers in labDict.items():
            if labelNumbers == inpNumbers:
                return name[0]  # return the first letter of the function name, that labels the shell

        # if still inside the function, it means that the numbers were not found, return None
        return None

    ###########################################################################################

    @property
    def atomicOrbitals(self):
        """This function returns the atomic orbitals, as a list of instances of the class AtomicOrbital,
        extracted from the set of shell functions defined in the instance of Shell. """

        return [AtomicOrbital(center=self.center, alphaList=self.alphaList, coeffList=self.coeffList,
                              cartesian=self.cartesian, radialNumbers=self.radNumbers[idAO])
                for idAO in range(self.nFunctions)]

    ###########################################################################################

    def prettyLog(self):
        """Return a string with a short formatted description of the shell data"""

        string = "Shell of atomic orbitals {0} on atom {1}\n".format(
            Shell.__angMomLabel__(self.radNumbers[0], self.cartesian), self.idAtom)
        string += "center: {0:6.3f} {1:6.3f} {2:6.3f}".format(*self.center)
        string += "   contract. degree: {0}\n".format(self.nContract)
        string += "contraction coeff.s: {0}\n".format(self.coeffList)
        string += "contraction alpha's: {0}\n".format(self.alphaList)
        if self.cartesian:
            string += "the shell contains {0} cartesian orbitals\n".format(self.nFunctions)

        return string

#####################################################################################################


class Orbitals:
    """" This class defines a set of molecular orbitals, by getting from input the definition
    of the atomic orbitals (given as a list of instances of the class Shell) and
    the coefficients of the molecular orbitals expansion.
    The order of the expansion coefficients is important! the coefficients should follow:
    1) the order of the AO shells that are read from input
    2) within each shell, the functions are ordered as follow:
       for spherical harmonics: 0, +1, -1, +2, -2, +3, -3 .. etc etc
       for cartesian functions: X,Y,Z, | XX,YY,ZZ,XY,XZ,YZ | XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ """

    def __init__(self, AOshells, MOcoeffs, MOoccup = None):

        # store the definition of the shells
        self.AOshells = AOshells

        # create a list that maps a progressive integer index to the specific atomic orbital in the full basis
        self.ID2AtomicOrb = []
        for iShell, shell in enumerate(self.AOshells):
            for iOrb in range(shell.nFunctions):
                self.ID2AtomicOrb.append([iShell, iOrb])

        # store the total number of atomic orbitals
        self.AOnumber = len(self.ID2AtomicOrb)

        # store the coefficients of the atomic orbitals
        self.MOcoeffs = MOcoeffs
        # store the labels for the occupation of the MOs (when available)
        self.MOoccup = MOoccup

    ###########################################################################################

    @property
    def atomicOrbitals(self):
        """This function returns the atomic orbitals as a list of instances of the class AtomicOrbital,
        extracted from the full set of basis functions defined in the instance of Orbitals. """

        # collect the lists of AOs of the shells and combine them in a unique list
        AOList = []
        for shell in self.AOshells:
            AOList += shell.atomicOrbitals

        # return the list
        return AOList

    ###########################################################################################
    #                            READ MO FROM GAUSSIAN OUTPUT
    ###########################################################################################

    @classmethod
    def readMOfromGaussianOutput(cls, gFileName):

        # open gaussian output file
        with open(gFileName) as f:
            fileText = f.read()

        chkfilename = gFileName.replace(".log", ".chk")
        subprocess.call(["formchk", chkfilename])
        fchkfilename = chkfilename.replace(".chk", ".fchk")
        with open(fchkfilename) as fchkfile:
            fchklines = fchkfile.readlines()

        # define regular expressions to read geometry from gaussian output
        geometryRegExpr = "Input orientation:" + ".*?" + "--*-" + ".*?" + "--*-" + "(.*?)" + "--*-"
        geometry2RegExpr = "Z-Matrix orientation:" + ".*?" + "--*-" + ".*?" + "--*-" + "(.*?)" + "--*-"

        # extract geometry information with regular expression
        try:
            regExprText = re.findall(geometryRegExpr, fileText, re.DOTALL)[0]
        except IndexError:
            regExprText = re.findall(geometry2RegExpr, fileText, re.DOTALL)[0]
        atomCoords = {}
        for line in regExprText.rstrip().lstrip().splitlines():
            atomCoords[int(line.split()[0])] = [float(x) / constants.Bohr2Ang for x in line.split()[3:]]

        # define regular expressions to read type of radial function
        radialRegExpr = "(\([0-9]*D, [0-9]*F\))"

        # reading cartesian or spherical harmonics functions
        radialType = re.findall(radialRegExpr, fileText)[0]
        cartesianD, cartesianF = True, True
        if radialType == "(5D, 7F)":
            cartesianD, cartesianF = False, False
        elif radialType == "(6D, 10F)":
            cartesianD, cartesianF = True, True
        elif radialType == "(6D, 7F)":
            cartesianD, cartesianF = True, False
        else:
            logwrt.fatalerror("the radial definitions of basis set '{0}' is not implemented!\n".format(radialType) +
                              "Please, use either (5D, 7F), (6D, 10F) or (6D, 7F) in Gaussian input file.")

        # define regular expressions to read the basis set definition
        gauBasisRegExpr = "AO basis set in the form of general basis input.*?:" + \
                          "(.*?)" + "\*\*\*\*\n\n"

        # reading basis set definition from input file
        try:
            regExprText = re.findall(gauBasisRegExpr, fileText, re.DOTALL)[0]
        except IndexError:
            logwrt.fatalerror("cannot fine the basis set definition in the Gaussian log!\n" +
                              "Please add the 'gfinput' option in the Gaussian input file.")

        # now loop over all the atomic shells defined in the basis set definition
        shellList = []
        for atomShell in regExprText.split("****"):

            # process text, removing initial and trailing newlines, and then split line-by-line
            shellText = atomShell.rstrip().lstrip()
            shellLines = shellText.splitlines()

            # get the atom ID and from the ID get the atomic center
            atomID = int(shellLines[0].split()[0])
            atomCenter = atomCoords[atomID]

            # define regular expression to find the shell descriptor lines
            shellDescriptorRegExpr = "([S,P,D,F,G]P?  ) *?([0-9]*) *?([0-9]\.[0-9]*)"

            nline = 2
            for shellDescriptor in re.findall(shellDescriptorRegExpr, shellText):

                shellType, nGauss, scaling = shellDescriptor
                shellType = shellType.rstrip().lstrip()
                nGauss = int(nGauss)
                scaling = float(scaling)

                alphas = [float(line.split()[0].replace('D', 'E')) * scaling ** 2
                          for line in shellLines[nline:nline + nGauss]]
                coeffs = [float(line.split()[1].replace('D', 'E')) for line in shellLines[nline:nline + nGauss]]

                if shellType[0] == "D":
                    cartesian = cartesianD
                elif shellType[0] == "F":
                    cartesian = cartesianF
                else:
                    cartesian = True
                shellList.append(Shell(idAtom=atomID, center=atomCenter, angMom=shellType[0], cartesian=cartesian,
                                       alphaList=alphas, coeffList=coeffs))

                if shellType == "SP":
                    coeffsP = [float(line.split()[2].replace('D', 'E')) for line in shellLines[nline:nline + nGauss]]
                    shellList.append(Shell(idAtom=atomID, center=atomCenter, angMom="P", cartesian=True,
                                           alphaList=alphas, coeffList=coeffsP))

                nline += nGauss + 1

        # read MOs from fchk file
        READMOs = False
        molOrb, MOoccupation = [], []
        Nalpha, Nbeta, Nbasis, N = 0, 0, 0, 0

        for i in range(len(fchklines)):
            if "Number of alpha electrons" in fchklines[i]:
                Nalpha = int(fchklines[i].split()[-1])
            if "Number of beta electrons" in fchklines[i]:
                Nbeta = int(fchklines[i].split()[-1])
            if "Number of basis functions"in fchklines[i]:
                Nbasis = int(fchklines[i].split()[-1])
                N_O = max([Nalpha, Nbeta])
                N_V = Nbasis - N_O
                MOoccupation += ["O" for k in range(N_O)]
                MOoccupation += ["V" for k in range(N_V)]
            if "Alpha MO coefficients" in fchklines[i]:
                READMOs = True
                N = int(fchklines[i].split()[-1])
                continue
            if READMOs:
                molOrb += [float(j) for j in fchklines[i].split()]
                if len(molOrb) == N:
                    break
        molOrb = [molOrb[j:j + Nbasis] for j in range(0, N, Nbasis)]

        return cls(AOshells=shellList, MOcoeffs=molOrb, MOoccup=MOoccupation)


def _primitiveOverlap(alpha, beta, A, B, la, lb):
    """Compute the overlap between two primitive gaussian functions, using a recursion formula to
    evaluate the integrals of the cartesian angular components. Note that the integral here is defined
    assuming that the gaussian functions do not contain any normalization factor.
    For all the definition, refer to Ho, Hernandez-Perez
    'Evaluation of Gaussian Molecular Integrals. I. Overlap Integrals'
    https://www.mathematica-journal.com/2012/02/16/evaluation-of-gaussian-molecular-integrals/ """

    # compute the vector between the centers of the gaussian functions
    AB = np.array(A) - np.array(B)
    # compute the vector P, which is the center of the gaussian product of the two gaussians A and B
    P = (alpha * np.array(A) + beta * np.array(B)) / (alpha + beta)

    # compute the term that is independent on the angular components
    EAB = math.exp(-(alpha * beta) * np.linalg.norm(AB) ** 2 / (alpha + beta))

    # define the tables of the integrals of the angular parts, computed with the recursion relation
    # (see below the function recursionTable)
    sx = recursionTable(alpha, beta, P[0], A[0], B[0], la[0] + lb[0])
    sy = recursionTable(alpha, beta, P[1], A[1], B[1], la[1] + lb[1])
    sz = recursionTable(alpha, beta, P[2], A[2], B[2], la[2] + lb[2])

    # now return the overlap integral of the (not-normalized) atomic orbitals
    return math.sqrt((math.pi / (alpha + beta)) ** 3) * EAB * sx[(la[0], lb[0])] * sy[(la[1], lb[1])] * sz[
        (la[2], lb[2])]


def recursionTable(alpha, beta, P, A, B, imax):
    """Compute the integrals of the cartesian angular components up to the elements with na + nb = imax"""

    # define starting elements of the recursion relation
    s = {(0, 0): 1.0, (1, 0): -(A - P), (0, 1): -(B - P)}

    # now construct the rest of the elments (when are needed) with the recursion relation
    for i in range(1, imax):
        # first build the element with labels (i,0)
        s[(i + 1, 0)] = -(A - P) * s[(i, 0)] + float(i) / (2 * (alpha + beta)) * s[(i - 1, 0)]
        # and then construct all the elements corresponding to the same a+b value
        for b in range(i + 1):
            a = i - b
            s[(a, b + 1)] = s[(a + 1, b)] + (A - B) * s[(a, b)]

    # return the whole dictionary of s values
    return s


def normalization(alpha, la):
    """Compute the normalization factor of a primitive atomic orbital, given the exponent alpha and the set of
    three exponents defining the cartesian components of the angular part"""

    denom = math.sqrt(dfactorial(2 * la[0] - 1) * dfactorial(2 * la[1] - 1) * dfactorial(2 * la[2] - 1))
    return (2. * alpha / math.pi) ** (3. / 4) * (4. * alpha) ** (float(la[0] + la[1] + la[2]) / 2) / denom


def dfactorial(i):
    """Semifactorial of a number i defined with its recursion relation within a recursive function"""

    if i == -1 or i == 0:
        return 1
    else:
        return i * dfactorial(i - 2)


if __name__ == '__main__':

    # start COBRAMM with message to log
    logwrt.cobramstart()

    gauOrb = Orbitals.readMOfromGaussianOutput("gaussian-QM.log")

    # for shell in gauOrb.AOshells:
    #     print(shell.prettyLog())

    # for AO in gauOrb.atomicOrbitals:
    #     print(AO.prettyLog())

    atomOrb = gauOrb.atomicOrbitals
    overlapMatrix = np.zeros([len(atomOrb),len(atomOrb)])
    orbMatrix = np.array(gauOrb.MOcoeffs)

    sys.stdout.write("\n\nOVERLAP BETWEEN ATOMIC ORBITALS\n\n")
    for ileft in range(len(atomOrb)):
        for iright in range(len(atomOrb)):
            overlap = AtomicOrbital.AOoverlap(atomOrb[ileft], atomOrb[iright])
            overlapMatrix[ileft, iright] = overlap
        #     sys.stdout.write("{0:6d}{1:6d}{2:20.6f}\n".format(ileft, iright, overlap))
        # sys.stdout.write("\n")

    sys.stdout.write("\n\nOVERLAP BETWEEN MOLECULAR ORBITALS\n\n")
    MOoverlap = np.matmul(orbMatrix, np.matmul(overlapMatrix, orbMatrix.T))
    for ileft in range(len(atomOrb)):
        for iright in range(len(atomOrb)):
            sys.stdout.write("{0:6d}{1:6d}{2:20.6f}\n".format(ileft, iright, MOoverlap[ileft, iright]))
        sys.stdout.write("\n")

    # AOlist = gauOrb.atomicOrbitals
    # pippo = AtomicOrbital.AOoverlap(AOlist[0], AOlist[6])
    # print("overlap between function 0 and 6: {0}\n".format(pippo))
