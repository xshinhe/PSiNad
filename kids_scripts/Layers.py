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

# standard library
import sys  # operating system utilities
import os
# import os                       # system-specific parameters and functions
import copy  # shallow and deep copy operations

# external modules
import numpy as np  # numpy library for scientific computation

# cobramm local modules
import kids_log  # manages log file output + start/end procedures
import constants  # values of physical constants and conversion factors


class Layers:
    """
    Description
        This is the geometry object for the system
        Replaces AXYZ_real (now self.cartesian) and all derived geometries like AXYZ_model, AXYZ_modelH
        Replaces CBF.ReadInitialGeom, CBF.updatereal, amber.makerealcrd, CBF.rebuildrealcrd
        Replaces parallel_numerics.displace (and restore function is not used anymore)
        Partially replaces "splitted" from CBF.split (geometry related fields)
        Replaces "lists" from CBF.ReadInitialLists without command (formerly lists[6])
    How To
        __init__ method takes geometry filename as argument (defaults to 'real_layers.xyz')
        static properties are computed once @__init__ and stored in memory
        sub-geometries are computed on-demand each time, can be accessed as attributes and DON'T contain atomLabel
    """

    # hardcoded values for bond lengths
    dCH = 1.090  # float(command[21])
    dOH = 0.947  # float(command[22])
    dNH = 1.008  # float(command[23])
    dCC = 1.526  # float(command[24])
    dCO = 1.410  # float(command[25])
    dCN = 1.475  # float(command[26])

    def __init__(self, atomLabels: list, atomCoords: list, atomTypes: list, atomLinks: list = None):
        """
        Constructor method, initialize the geometry attribute (self.cartesian) and all the static lists
        :param atomLabels    list of the atomic labels of the system
        :param atomCoords    list of the atomic coordinates, in format  [[x0,y0,z0], [x1,y1,z1], ...]
        :param atomTypes     list of the labels for the H/M/L assignment of the atoms
        :param atomLinks     list of the atomic links between layers, [[a1, a2 ..., aN], [b1, b2, ..., bN]]
        """

        self.PTE = {
                "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36  
        }
        # constructing static lists of layer atoms
        # loop over atoms, IMPORTANT: integer labels of the atoms are i+1 !!!
        _list_HIGH, _list_MEDIUM, _list_LOW, _list_MEDIUM_HIGH = [], [], [], []
        _list_MM, _list_QM = [], []
        _H_M_AL = []
        for i in range(len(atomTypes)):
            # atom: i+1 is of type: atom_type
            atom_id = i+1
            atom_type = atomTypes[i]
            # assign the atom to the appropriate lists
            if atom_type == 'H':
                _list_HIGH.append(atom_id)
                _list_QM.append(atom_id)
            elif atom_type == 'M':
                _list_MEDIUM.append(atom_id)
            else:
                _list_LOW.append(atom_id)
            if atom_type in 'ML':
                _list_MM.append(atom_id)
            if atom_type in 'HM':
                _H_M_AL.append(atom_type)
                _list_MEDIUM_HIGH.append(atom_id)

        # process geometry data from [[x0,y0,z0], ...] top [[x0, x1 ...], [...], [...]] format
        _x, _y, _z = [], [], []
        for coord in atomCoords:
            _x.append(float(coord[0]))
            _y.append(float(coord[1]))
            _z.append(float(coord[2]))

        # init dummy variables
        if atomLinks is None:
            _atom_link = [[], []]
        else:
            _atom_link = atomLinks

        # define instance variables using np.array
        self.cartesian = [np.array(_x), np.array(_y), np.array(_z)]  # x,y,z coordinates
        self.atomType = np.array(atomTypes)  # H/M/L
        self.atomLabel = np.array(atomLabels)  # atom labels
        self.list_HIGH = np.array(_list_HIGH)  # H type atoms index
        self.list_MEDIUM = np.array(_list_MEDIUM)  # M type atoms index
        self.list_LOW = np.array(_list_LOW)  # L type atoms index
        self.list_MM = np.array(_list_MM)  # M/L type atoms index
        self.list_QM = np.array(_list_QM)  # duplicate of list_HIGH
        self.list_MEDIUM_HIGH = np.array(_list_MEDIUM_HIGH)  # HM type atoms index
        self.H_M_AL = np.array(_H_M_AL)  # HM types
        self.atomLink_BA = [[], []]  # atom links with H type atoms
        self.atomLink_DC = [[], []]  # atom links without H type atoms

        # A: HIGH atom @HM/HL interface
        # B: MEDIUM/LOW atom @HM/HL interface
        # C: MEDIUM atom @ML interface
        # D: LOW atom @ML interface
        for i in range(len(_atom_link[0])):
            if _atom_link[0][i] in self.list_QM:
                self.atomLink_BA[1].append(_atom_link[0][i])  # A (QM)
                self.atomLink_BA[0].append(_atom_link[1][i])  # B (?)
            elif _atom_link[0][i] in self.list_MM:
                self.atomLink_DC[1].append(_atom_link[0][i])  # C (MM)
                self.atomLink_DC[0].append(_atom_link[1][i])  # D

        # useful array lengths
        self.atomNum = len(self.cartesian[0])  # total number of atoms
        self.NsubH = len(self.atomLink_BA[0])  # n. of hydrogen substituted atoms
        self.NatomQM = len(self.list_QM)  # n. of H layer atoms
        self.NatomMM = len(self.list_MM)  # n. of M/L layer atoms
        self.NatomH = len(self.list_HIGH)  # n. of H layer atoms
        self.NatomM = len(self.list_MEDIUM)  # n. of M layers atoms
        self.NatomHM = len(self.list_MEDIUM_HIGH)  # n. of H/M layer atoms
        self.NatomL = len(self.list_LOW)  # n. of L layers atoms

        # H/M/L calculation type
        self.calculationType = ''.join(['H' if len(self.list_HIGH) != 0 else '',
                                        'M' if len(self.list_MEDIUM) != 0 else '',
                                        'L' if len(self.list_LOW) != 0 else ''])
        if self.calculationType == 'L':
            Log.fatalError('No atom in the HIGH and MEDIUM levels... check your geometry definition!')

    # =============================================================================================================

    @classmethod
    def from_only_xyz(cls, filetext: (str, list)):
        """
        Alternative constructor for the instances of the class (defined as a class method that returns
        an instance of the class). This method initializes a Layers instance by parsing the text of
        a real_layers.xyz file
        :rtype: Layers
        :param filetext: string containing the text of the "real_layers.xyz" file
        """

        # init dummy variables
        _atom_type, _atom_label, _atom_coords = [], [], []
        _atom_link = [[], []]

        # define the list of lines in the filetext, if it's a list, assume that it has already been split
        try:
            lineslist = filetext.splitlines()
        except AttributeError:
            lineslist = filetext

        # loop over the lines of real_layers.xyz, "i" is the atom id and starts from 1
        for i, line in enumerate(lineslist, start=1):
            if line.strip() == '': continue

            # split line with spaces
            s_line = line.split()
            # store coordinates, labels and layer definition for each atom
            _atom_coords.append([float(s_line[1]), float(s_line[2]), float(s_line[3])])
            _atom_label.append(s_line[0])
            _atom_type.append('H')
        # return instance of Layers
        return cls(_atom_label, _atom_coords, _atom_type, atomLinks=_atom_link)

    @classmethod
    def from_real_layers_xyz(cls, filetext: (str, list)):
        """
        Alternative constructor for the instances of the class (defined as a class method that returns
        an instance of the class). This method initializes a Layers instance by parsing the text of
        a real_layers.xyz file
        :rtype: Layers
        :param filetext: string containing the text of the "real_layers.xyz" file
        """

        # init dummy variables
        _atom_type, _atom_label, _atom_coords = [], [], []
        _atom_link = [[], []]

        # define the list of lines in the filetext, if it's a list, assume that it has already been split
        try:
            lineslist = filetext.splitlines()
        except AttributeError:
            lineslist = filetext

        # loop over the lines of real_layers.xyz, "i" is the atom id and starts from 1
        for i, line in enumerate(lineslist, start=1):
            # split line with spaces
            s_line = line.split()
            # store coordinates, labels and layer definition for each atom
            _atom_coords.append([float(s_line[2]), float(s_line[3]), float(s_line[4])])
            _atom_label.append(s_line[0])
            _atom_type.append(s_line[5])
            # when the line contains more than 6 field, there is an atom link
            if len(s_line) > 6:
                _atom_link[0].append(int(s_line[7]))
                _atom_link[1].append(i)

        # return instance of Layers
        return cls(_atom_label, _atom_coords, _atom_type, atomLinks=_atom_link)

    # =============================================================================================================

    @property
    def reallayertext(self):

        # define line format of the real_layers file (except the final piece with the link definition)
        realLayersFormat = "{0:3s}  0      {1:10.6f}    {2:10.6f}    {3:10.6f}  {4:1s}"

        # construct list of the lines of the real_layers.xyz file
        flines = []
        for lab, lay, x, y, z in zip(self.atomLabel, self.atomType, *self.cartesian):
            flines.append(realLayersFormat.format(lab, x, y, z, lay))

        # add the indications for the atom links
        for b, a in zip(*self.atomLink_BA):
            flines[b] += ' {0} {1}'.format(self.atomType[a-1], a)
        for d, c in zip(*self.atomLink_DC):
            flines[d] += ' {0} {1}'.format(self.atomType[c-1], c)

        # return string with lines separated by line break
        return "\n".join(flines)

    # =============================================================================================================

    def updatereal(self, filename="real.crd"):
        """
        This method updates self.cartesian reading coordinates from input file
        The input file can be the real.crd file (AMBER style) or the QM.in file (SHARC style)
        When the filename is "QMMM.in", tries to read the file assuming the latter option,
        otherwise tries to read it as a standard real.crd file.

        :param filename: input file containing new coordinates
        :return: None
        """

        # read the content of the file, and store the list of lines
        try:
            with open(filename) as f:
                crd = f.readlines()
        except IOError:
            Log.fatalError('{0} not found in current directory'.format(filename))

        if filename == "QMMM.in":  # read the input in the SHARC style
            # check that the number of coordinates is the same of the actual geometry
            N_atoms = int(crd[0].strip())
            if N_atoms != self.atomNum:
                Log.fatalError('The number of atoms in real_layers.xyz'
                                  '({0}) and in {1} ({2}) are different'.format(self.atomNum, filename, N_atoms))

            # read the following lines that contain the coordinates of the atoms
            x, y, z = [], [], []
            for iatom in range(N_atoms):
                row = crd[2 + iatom].split()
                x.append(float(row[1])), y.append(float(row[2])), z.append(float(row[3]))

            # check which distance unit is used for the input coordinates
            bohr = True
            for line in crd[2 + N_atoms:]:
                if line.split() and line.split()[0] == "unit":
                    if line.split()[1] == "bohr":
                        bohr = True
                    elif line.split()[1] == "angstrom":
                        bohr = False

            # if unit is atomic unit, convert to angstrom
            if bohr:
                self.cartesian = [np.array(x) * constants.Bohr2Ang, np.array(y) * constants.Bohr2Ang,
                                  np.array(z) * constants.Bohr2Ang]
            else:
                self.cartesian = [np.array(x), np.array(y), np.array(z)]

        else:  # read the input in the AMBER style
            # check that the number of coordinates is the same of the actual geometry
            N_atoms = int(crd[1].strip().split()[0])
            if N_atoms != self.atomNum:
                Log.fatalError('The number of atoms in real_layers.xyz'
                                  '({0}) and in {1} ({2}) are different'.format(self.atomNum, filename, N_atoms))

            # store the coordinates of the crd file
            x, y, z = [], [], []
            for line in crd[2:]:
                row = line.split()
                try:
                    x.append(float(row[0])), y.append(float(row[1])), z.append(float(row[2]))
                    x.append(float(row[3])), y.append(float(row[4])), z.append(float(row[5]))
                except IndexError:
                    pass
            self.cartesian = [np.array(x), np.array(y), np.array(z)]

    # =============================================================================================================

    def makerealcrd(self, filename='real.crd', filename2='geometry.xyz'):
        """
        This method generate a coordinates file using self.cartesian

        :param filename: coordinates file name, defaults to real.crd (??? unit is angstrom???)
        """
        with open(filename, 'w') as crd:
            crd.write('\n{0}\n'.format(self.atomNum))
            for i in range(self.atomNum):
                crd.write('{0:12.7f}{1:12.7f}{2:12.7f}{3}'.format(
                    self.cartesian[0][i], self.cartesian[1][i], self.cartesian[2][i], '\n' if i % 2 else ''))

        with open(filename2, 'w') as crd:
            crd.write('{0}\n\n'.format(self.NatomHM))
            for i in range(self.atomNum):        
                if i+1 in self.list_MEDIUM_HIGH:
                    crd.write('{0:3s}{1:12.7f}{2:12.7f}{3:12.7f}{4}'.format(
                        self.atomLabel[i], self.cartesian[0][i], self.cartesian[1][i], self.cartesian[2][i], '\n'))

    # =============================================================================================================

    def makemodelHcrd(self, filename):
        """
        This method writes into "filename" the coordinates of modelH geometry
        :param filename: file to write coordinates
        :return: None
        """
        modelH = self.modelH
        with open(filename, 'w') as crd:
            crd.write('\n{0}\n'.format(len(modelH[0])))
            for i in range(len(modelH[0])):
                try:
                    crd.write('{0:12.7f}{1:12.7f}{2:12.7f}{3}'.format(modelH[0][i],
                                                                      modelH[1][i],
                                                                      modelH[2][i],
                                                                      '\n' if i % 2 else ''))
                except IndexError:
                    pass

    def updateHMlayers(self, newgeom):
        """
        This method updates self.cartesian with new HM atom positions after an optimization/MD step
        :param newgeom: new geometry obtained from optimization or MD
        :return:
        """
        # init empty lists
        x, y, z = [], [], []
        x1, y1, z1 = [], [], []
        x2, y2, z2 = [], [], []
        # populate x,y,z with H atoms coordinates and x1,y1,z1 with M atoms coordinates
        for i in range(len(newgeom[0])):
            if self.H_M_AL[i] == 'H':
                x.append(float(newgeom[0][i]))
                y.append(float(newgeom[1][i]))
                z.append(float(newgeom[2][i]))
            elif self.H_M_AL[i] == 'M':
                x1.append(float(newgeom[0][i]))
                y1.append(float(newgeom[1][i]))
                z1.append(float(newgeom[2][i]))
        j = 0
        k = 0
        # populate x2,y2,z2 with all atoms coordinates
        for i in range(self.atomNum):
            # if i+1 in self.list_HIGH
            if np.add.reduce(np.equal(self.list_HIGH, i + 1)) == 1:
                x2.append(x[j])
                y2.append(y[j])
                z2.append(z[j])
                j += 1
            # elif i+1 in self.list_MEDIUM
            elif np.add.reduce(np.equal(self.list_MEDIUM, i + 1)) == 1:
                x2.append(x1[k])
                y2.append(y1[k])
                z2.append(z1[k])
                k += 1
            else:
                x2.append(self.cartesian[0][i])
                y2.append(self.cartesian[1][i])
                z2.append(self.cartesian[2][i])

        # update self.cartesian
        self.cartesian = [np.array(x2), np.array(y2), np.array(z2)]

    # =============================================================================================================

    def updateHlayer(self, newgeom):
        """
        This method updates self.cartesian with new H atom positions
        :param newgeom: new geometry to put into self
        :return:
        """

        # counter to loop over the lists of newgeom
        i_atom_input = 0

        # loop over the indices of the H atoms in the full geometry
        for i_atom_full in self.list_HIGH:
            for idim in range(3):
                self.cartesian[idim][i_atom_full-1] = newgeom[idim][i_atom_input]
            # increment the counter for newgeom
            i_atom_input += 1

    # =============================================================================================================

    def getModel(self, key):
        """
        Method called by other class internal methods, generate the requested set of coordinates
        :param key: type of geometry to generate
        :return: list of np.array containing atoms coordinates for selected sub-geometry
        """
        k_dict = {'model': self.list_QM,
                  'pod': self.list_MM,
                  'HIGH': self.list_HIGH,
                  'MEDIUM': self.list_MEDIUM,
                  'LOW': self.list_LOW,
                  'MEDIUM_HIGH': self.list_MEDIUM_HIGH}
        x, y, z = [], [], []
        for N_atom in k_dict[key]:
            x.append(self.cartesian[0][N_atom - 1])
            y.append(self.cartesian[1][N_atom - 1])
            z.append(self.cartesian[2][N_atom - 1])
        return [np.array(x), np.array(y), np.array(z)]

    def getAtomLabels(self, key):
        """
        Method called by other class internal methods, generate the requested set of atomic labels
        :param key: type of model for which the function generates the list of atomic labels
        :return: list of strings containing the atomic labels for the selected sub-geometry
        """

        # define in a dictionary the options that are already available and do not need an ad-hoc definition
        k_dict = {'model': self.list_QM,
                  'pod': self.list_MM,
                  'HIGH': self.list_HIGH,
                  'MEDIUM': self.list_MEDIUM,
                  'LOW': self.list_LOW,
                  'MEDIUM_HIGH': self.list_MEDIUM_HIGH}

        # for modelH, the list needs to be constructed by taking the QM atoms + a number of H atoms
        if key == "modelH":
            return [self.atomLabel[i - 1] for i in self.list_QM] + ["H"] * self.NsubH
        # the other key values should correspond to atom lists that are already defined
        else:
            return [self.atomLabel[i - 1] for i in k_dict[key]]

    @property
    def model(self):
        """
        Property-method that calls self.getMethod
        :return: a list of np.arrays containing the coordinates of QM atoms
        """
        # return only QM atoms
        return self.getModel('model')

    @property
    def modelH(self):
        """
        modelH adds hydrogen positions to self.model (coordinates of QM atoms),
        H coordinates are appended at the END of self.model geometry.
        The method does NOT return the atomLabel list at current version, may be added later
        :return: a list of np.array cointaining the coordinates of QM atoms saturated including hydrogens
        """

        # init x,y,z from self.model
        tmp_cartesian = self.model
        x = list(tmp_cartesian[0])
        y = list(tmp_cartesian[1])
        z = list(tmp_cartesian[2])

        # if atom links are present in the H layer
        if self.NsubH > 0:
            # loop over self.atomLink_BA
            for i in range(self.NsubH):
                j = self.atomLink_BA[1][i] - 1
                k = self.atomLink_BA[0][i] - 1
                # "a" is the HIGH layer atom
                # "b" is the MEDIUM/LOW layer atom, replaced with hydrogen
                Xa = float(self.cartesian[0][j])
                Ya = float(self.cartesian[1][j])
                Za = float(self.cartesian[2][j])
                Xb = float(self.cartesian[0][k])
                Yb = float(self.cartesian[1][k])
                Zb = float(self.cartesian[2][k])
                # calculate hydrogen coordinates
                if self.atomLabel[j] == 'C':
                    Xc = ((Xb - Xa) * (self.dCH / self.dCC)) + Xa
                    Yc = ((Yb - Ya) * (self.dCH / self.dCC)) + Ya
                    Zc = ((Zb - Za) * (self.dCH / self.dCC)) + Za
                elif self.atomLabel[j] == 'O':
                    Xc = ((Xb - Xa) * (self.dOH / self.dCO)) + Xa
                    Yc = ((Yb - Ya) * (self.dOH / self.dCO)) + Ya
                    Zc = ((Zb - Za) * (self.dOH / self.dCO)) + Za
                elif self.atomLabel[j] == 'N':
                    Xc = ((Xb - Xa) * (self.dNH / self.dCN)) + Xa
                    Yc = ((Yb - Ya) * (self.dNH / self.dCN)) + Ya
                    Zc = ((Zb - Za) * (self.dNH / self.dCN)) + Za
                else:
                    Log.writeLog(
                        '\nCut through a bond involving atoms different than C,O,N in the high layer not implemented\n')
                    sys.exit(1)

                x.append(Xc)
                y.append(Yc)
                z.append(Zc)

        return [np.array(x), np.array(y), np.array(z)]

    @property
    def pod(self):
        """
        Property-method that calls self.getMethod
        :return: a list of np.arrays containing the coordinates of MM atoms
        """
        return self.getModel('pod')

    @property
    def pod_mobileatoms(self):
        """
        Property-method that returns a boolean vector to identify the mobile atoms
        within the MM atom list
        :return: a list of True(=mobile)/False(=fixed), one for each MM atom
        """

        return [(True if natom in self.list_MEDIUM else False) for natom in self.list_MM]

    @property
    def HIGH(self):
        """
        Property-method that calls self.getMethod
        :return: a list of np.arrays containing the coordinates of QM atoms
        """
        return self.getModel('HIGH')

    @property
    def MEDIUM(self):
        """
        Property-method that calls self.getMethod
        :return: a list of np.arrays containing the coordinates of MEDIUM atoms
        """
        return self.getModel('MEDIUM')

    @property
    def LOW(self):
        """
        Property-method that calls self.getMethod
        :return: a list of np.arrays containing the coordinates of LOW atoms
        """
        return self.getModel('LOW')

    @property
    def MEDIUM_HIGH(self):
        """
        Property-method that calls self.getMethod
        :return: a list of np.arrays containing the coordinates of MEDIUM + HIGH atoms
        """
        return self.getModel('MEDIUM_HIGH')

    def distance(self, iatom, jatom):
        """
        Return Euclidean distance between two atoms defined in the current Layers instance.
        :param iatom: integer, index of the first atom
        :param jatom: integer, index of the second atom
        :return: floating point number, cartesian distance between the two atoms
        """
        return np.sqrt((self.cartesian[0][iatom] - self.cartesian[0][jatom]) ** 2 +
                       (self.cartesian[1][iatom] - self.cartesian[1][jatom]) ** 2 +
                       (self.cartesian[2][iatom] - self.cartesian[2][jatom]) ** 2)

    def hbondratio(self, atomlabel): # unused function
        """
        Return the ratio between the average length of the X-H atom bonds and the X-C atom bond.
        :param atomlabel: string with the name of the X atom (currently only C, O and N are implemented)
        :return: a float with d(X-H)/d(X-C), or None if the value of atomlabel is not supported
        """

        # depending on the value of atomlabel, compute the ration between the values stored in the object instance
        if atomlabel == 'C':
            ratio = self.dCH / self.dCC  # in former COBRAMM implementation, this was command[21] / command[24]
        elif atomlabel == 'O':
            ratio = self.dOH / self.dCO  # in former COBRAMM implementation, this was command[22] / command[25]
        elif atomlabel == 'N':
            ratio = self.dNH / self.dCN  # in former COBRAMM implementation, this was command[23] / command[26]
        else:
            ratio = None

        # return the value to caller
        return ratio

    # =============================================================================================================

    def modelHdisplace(self, iAt, iCoord, iDir, disp, Hlinkfollow):
        """
        This routine determines which atom to move, displaces it by user-specified value (default is 0.001 Angstrom)
        It is a new version of displace, that has all the necessary input information explicitly defined in the
        input arguments.
        :param iAt: integer index of the High layer atom to displace
        :param iCoord: components (0 = x, 1 = y, 2 = z) to displace
        :param iDir: direction of displacement (+1 = +, -1 = -)
        :param disp: magnitude of displacement, in the same units used to store the geometry
        :param Hlinkfollow: H atoms of the atom links follow the atoms to which they are bound
        :return: displaced modelH coordinates
        """

        # move the H atoms with QM atoms that they are bound to
        if Hlinkfollow:

            # temporary placeholders for geometry
            _tmp = copy.deepcopy(self.cartesian)

            # displace geometry of the H atom BEFORE generating the modelH
            atomId = self.list_HIGH[iAt]
            # displace the given coordinate of the given atom, with the given direction
            self.cartesian[iCoord][atomId - 1] += iDir * disp

            # create modelH geometry AFTER displacing atom
            modelH = self.modelH

            # restore geometry
            self.cartesian = _tmp

        # otherwise the H atoms are fixed in the position generated at equilibrium
        else:

            # create modelH geometry BEFORE displacing atom
            modelH = self.modelH
            # displace the given coordinate of the given atom, with the given direction
            modelH[iCoord][iAt] += iDir * disp

        # return displaced modelH geometry
        return modelH

    # =============================================================================================================

    def displace(self, command, step): # repeat with above but only used in molcas.py
        """
        This routine determines which atom to move, displaces it by user-specified value (default is 0.001 Angstrom)
        Replaces parallel_numerics.displace
        :param command: list of commands
        :param step: step number
        :return: displaced atom index and coordinates
        """
        numerical_displace = command[16]
        disp = float(command[12])

        # initialize coord
        coord = [None, None, None]

        # define atom and coordinate to displace
        if command[1] in ['optxg', 'mdv', 'irc', 'ci', 'ts'] and command[10] == '0':
            iAt = (step - 1) // 3
            iCoord = (step - 1) % 3
            iDir = 0
        else:
            iAt = (step - 1) // 6
            iCoord = ((step - 1) % 6) // 2
            iDir = ((step - 1) % 6) % 2

        # explicitly compute numerical derivative also for H atoms!
        if numerical_displace == "1":

            # create modelH geometry BEFORE displacing atom
            modelH = self.modelH
            coord = [modelH[0], modelH[1], modelH[2]]

            if not iDir:  # plus displacement (iDir == 0)
                coord[iCoord][iAt] += disp
            elif iDir:  # minus displacement (iDir == 1)
                coord[iCoord][iAt] -= disp

        # do not compute numerical derivative of H atoms, but move them
        # with QM atoms that they are bound to
        elif numerical_displace == "0":

            # temporary placeholders for geometry
            _tmp = copy.deepcopy(self.cartesian)

            # displace geometry of the H atom BEFORE generating the modelH
            atomId = self.list_HIGH[iAt]

            if not iDir:  # plus displacement (iDir == 0)
                self.cartesian[iCoord][atomId - 1] += disp
            elif iDir:  # minus displacement (iDir == 1)
                self.cartesian[iCoord][atomId - 1] -= disp

            # create modelH geometry AFTER displacing atom
            modelH = self.modelH
            coord = [modelH[0], modelH[1], modelH[2]]

            # restore geometry
            self.cartesian = _tmp

        # return displaced modelH geometry with info on atom, coord, direction
        return iAt, iCoord, iDir, coord[0], coord[1], coord[2]

    # =============================================================================================================

    def read_displacement(self, step): # only used by old molcas
        """
        This routine reads the displacement from external file "dqn.xyz" where n==step. The user has to make sure that the displacement module matched the user-specified value command[12] (default is 0.001 Angstrom)
        The routine is only accessed when command[210] > 0
        :param step: step number
        :return: displaced atom index and coordinates
        """
        Log.writeLog('reading external displacement for step '+str(step)+' from external file dq'+str(step)+'.xyz\n')
        iAt = (step - 1) // 6
        iCoord = ((step - 1) % 6) // 2
        iDir = ((step - 1) % 6) % 2
        # initialize coord
        coord = [[], [], []]
        f=open('dq'+str(step)+'.xyz', 'r')
        for line in f:
            if len(line.split())==4:
                coord[0].append(float(line.split()[1]))
                coord[1].append(float(line.split()[2]))
                coord[2].append(float(line.split()[3]))

        # return displaced modelH geometry with info on atom, coord, direction
        return iAt, iCoord, iDir, coord[0], coord[1], coord[2]

    # =============================================================================================================

    def to_string(self, geom_model='cartesian', labels=False, precision=6):
        """
        Method that returns a printable string for selected sub-geometry
        :param geom_model: wanted geometry model (model, modelH, pod,
        HIGH, MEDIUM, LOW, MEDIUM_HIGH)
        :param labels: boolean, print atom labels
        :param precision: int, number of decimals to print
        :return: formatted string with requested sub-geometry
        """
        # retrieve coordinates and store them in _geom
        _geom = eval("self.{0}".format(geom_model))

        # compute atom labels
        if labels:
            # same as getModel's k_dict
            model_dict = {'cartesian': range(1, len(self.cartesian[0]) + 1),
                          'model': self.list_QM,
                          'pod': self.list_MM,
                          'HIGH': self.list_HIGH,
                          'MEDIUM': self.list_MEDIUM,
                          'LOW': self.list_LOW,
                          'MEDIUM_HIGH': self.list_MEDIUM_HIGH}
            # crafting atom labels for modelH
            if geom_model == 'modelH':
                AT = [self.atomLabel[i - 1] for i in self.list_QM] + \
                     ['H' for _ in range(self.NsubH)]
            else:
                AT = [self.atomLabel[i - 1] for i in model_dict[geom_model]]

        else:
            AT = ['' for _ in range(len(_geom[0]))]

        # transpose _geom and store in geom2
        geom2 = np.transpose(_geom)

        # craft formatted string
        final_string = ''
        form_text = '{' + '0:15.{0}f'.format(precision) + '}'
        for at, line in zip(AT, geom2):
            line = [form_text.format(c) for c in line]
            final_string += '{0}\t{1}\n'.format(at, ''.join(line))

        return final_string
