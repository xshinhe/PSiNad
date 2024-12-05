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

# import statements of module from python standard library

import os  # filesystem utilities
import shutil  # filesystem utilities
import subprocess  # run external program as child process
import re  # process output with regular expressions
import math  # import mathematical functions
import shlex       # simple lexical analysis for unix syntax

# imports of local modules

import logwrt  # write messages and output to log (std out or cobramm.log)
import cobrammenv  # general purpose functions in use in COBRAMM
import constants  # value of physical constants of common use

# imports of local objects

from amberCalculator import AmberCalculator, AmberSnapshot, AmberTopology  # classes that constitute the Amber wrapper
from layers import Layers  # class to define the geometry of a QM/MM calculation
from output import Output  # class to parse and write the cobramm.xml file

# math libraries

import numpy as np  # numpy: arrays and math utilities


######################################################################################################################

class CobrammError(Exception):
    """A custom exception used to report errors in use of COBRAMM"""


######################################################################################################################

class CobrammCalculator:
    """ wrapper for the execution of COBRAMM and related utilities """

    _BOND_TOLERANCE = 1.3  # tolerance in defining bonds between atoms
    #                         (1.5 means that atoms at a distance 1.5*reference value are considered bonded)

    # class variable that controls when the check on the environment has been done
    _cobrammCheck = False

    # names of the COBRAMM input/output files
    REALTOP = "real.top"
    MODELHTOP = "model-H.top"
    REALLAYERS = "real_layers.xyz"
    COBRAMMCOMMAND = "cobram.command"
    COBRAMMLOG = "cobramm.log"

    def __init__(self):

        if not CobrammCalculator._cobrammCheck:

            # environmental variable with the COBRAMM path
            try:
                cobramm_path = os.environ['COBRAM_PATH']
            except KeyError:
                raise CobrammError("environment variable $COBRAM_PATH is not defined")

            # check if COBRAMM executables are available
            for exe in ["cobram.py", "cobramext", "freqext"]:
                if not cobrammenv.which(exe):
                    raise CobrammError(exe + " executable is not available")

            # prepare execution of AMBER tools by creating an AmberCalculator instance
            AmberCalculator()

            # write message about the AMBER version currently in use
            logwrt.writelog("Running COBRAMM {0} from {1}\n".format(cobrammenv.getVersion(), cobramm_path))

            # set class variable to true, to skip this check at next initialization
            CobrammCalculator._cobrammCheck = True

    # ================================================================================================================

    @staticmethod
    def run(inputfile, geometry, modeltopology, realtopology, calc_dir="cobrammCalc", n_cores=1, store=False):
        """ run COBRAMM with the given input files and store the results in outputData.

        :rtype: CobrammOutput
        :param inputfile: CobrammInput class instance with the input file for the COBRAMM run
        :param geometry: Layers class instance with the input geometry and layer definition
        :param modeltopology: AmberTopology class instance with the real topology file
        :param realtopology: AmberTopology class instance with the model topology file
        :param calc_dir: string that specify where the calculation is run
        :param n_cores: number of cores used in a parallel run
        :param store: do not remove output files when execution ends
        :return: CobrammOutput class instance with the results of the calculation
        """

        # store starting dir
        start_dir = os.getcwd()

        # when the existing working directory is not re-used, move there and remove it
        if os.path.isdir(calc_dir):
            logwrt.writewarning("overwriting previous " + calc_dir + " directory ")
            shutil.rmtree(calc_dir)

        # if needed create directory and move there
        if not os.path.isdir(calc_dir):
            os.mkdir(calc_dir)

        # when the required number of cores for the QM calculation is more than 1, modify the input file
        if n_cores > 1:
            inputfile["nproc"] = "{0}".format(n_cores)

        # write topologies to file
        with open(os.path.join(calc_dir, CobrammCalculator.MODELHTOP), "w") as f:
            f.write(modeltopology.topotext)
        with open(os.path.join(calc_dir, CobrammCalculator.REALTOP), "w") as f:
            f.write(realtopology.topotext)
        # write geometry and layers to the real layers file
        with open(os.path.join(calc_dir, CobrammCalculator.REALLAYERS), "w") as f:
            f.write(geometry.reallayertext)
        # write cobramm command file
        with open(os.path.join(calc_dir, CobrammCalculator.COBRAMMCOMMAND), "w") as f:
            f.write(inputfile.inputtext)

        # run the COBRAMM calculation
        os.chdir(calc_dir)
        with open(CobrammCalculator.COBRAMMLOG, "w") as fstdout:
            subprocess.call(["cobram.py"], stdout=fstdout, stderr=fstdout)
        os.chdir(start_dir)

        # now create the outputData instance reading from the output of the AMBER calculation
        # AMBER files are stored only if requested with the argument store
        return CobrammOutput(calc_dir, store_files=store)

    # ================================================================================================================

    @staticmethod
    def rattle(geometry: Layers):
        """
        :param geometry:
        :rtype: (str)
        :return:
        """

        # initialize list to store rattle bonds
        rattleBonds = []

        # loop over all the couples of atoms in the M layer
        for j, iAtom2 in enumerate(geometry.list_MEDIUM):
            for iAtom1 in geometry.list_MEDIUM[0:j]:

                # extract the atomic labels
                lab1 = geometry.atomLabel[iAtom1 - 1].upper()
                lab2 = geometry.atomLabel[iAtom2 - 1].upper()

                # consider only the X-H possible bonds
                if lab1 == "H" or lab2 == "H":

                    # get the relevant reference distance for this type of atoms
                    refDistance = constants.sngBond[(lab1, lab2)] * constants.Bohr2Ang
                    # compute the actual distance between the atoms
                    actDistance = geometry.distance(iAtom1 - 1, iAtom2 - 1)

                    # when the distance does not exceed the reference + tolerance, then this is a constrained X-H bond
                    if actDistance < refDistance * CobrammCalculator._BOND_TOLERANCE:
                        rattleBonds.append((iAtom1, iAtom2))
                        logwrt.writelog("Constraining bond {0}-{1} between atoms {2} and {3} at {4} Ang\n".format(
                            lab1, lab2, iAtom1, iAtom2, actDistance
                        ))

        # when at least one bond was identified, print the rattle group
        if len(rattleBonds) > 0:
            rattleStr = "!rattle\n"
            for iAtom1, iAtom2 in rattleBonds:
                rattleStr += "{0} {1}\n".format(iAtom1, iAtom2)
            rattleStr += "?rattle\n"
        else:  # otherwise return None
            rattleStr = None
        return rattleStr

    # ================================================================================================================

    @staticmethod
    def qmmmLayers(realSnapshot: AmberSnapshot, realTopology: AmberTopology, nrResidue: int = 1,
                   calcDir: str = "tmpLayers", mobileThreshold: float = 3.0, reorderRes: bool = True, dry: bool = False):
        """
        :param realSnapshot: AmberSnapshot instance defyining the initial structure
        :param realTopology: AmberTopology instance defyining the initial topology
        :param nrResidue: ordinar number (1,2...) that indicates the residues that will be the QM HIGH layer
        :param calcDir: temporary directory where the scratch files of this function are written
        :param mobileThreshold: distance threshold to divide the MEDIUM from the LOW layer, in Angstrom
        :param reorderRes: logical flag to control the reordering of the residues
        :rtype: (Layers, AmberTopology, AmberTopology)
        :return: the data needed to initialize a COBRAMM calculation: the topologies of the QM and full structure
                 and a Layers instance that store all the geometrical data and the layer assignment
        """
        # store starting dir and move to directory of the calculation
        startDir = os.getcwd()
        # then move to the work d where preprocessing files are stored, if the dir exists, remove first
        if os.path.isdir(calcDir):
            logwrt.writewarning("overwriting previous " + calcDir + " directory ")
            shutil.rmtree(calcDir)
        os.mkdir(calcDir), os.chdir(calcDir)

        # write topology to file
        with open("input00.top", "w") as f:
            f.write(realTopology.topotext)
        # write coordinates and unit cell constants to file
        with open("input00.crd", "w") as f:
            f.write(realSnapshot.crdtext)

        # *********************************************************************************************

        # write input file for cpptraj (create topology for the QM part: molecule without solvent)
        with open("cpptrajscript.01", "w") as f:
            f.write(" parm {0} \n trajin {1} \n reference {1} \n".format("input00.top", "input00.crd"))
            f.write(" strip !(:{0}) outprefix nosolv \n go \n".format(nrResidue))

        # run cpptraj
        with open("cpptrajscript.stdout", "w") as fstdout:
            subprocess.call(["cpptraj", "-i", "cpptrajscript.01"], stdout=fstdout, stderr=fstdout)
        if dry:
            modelTopology = realTopology
        else:
            with open("nosolv.input00.top", "r") as f:
                modelTopology = AmberTopology(f.read())

        # *********************************************************************************************

        # extract atomic symbols, coordinates and residues
        symbolList, coordsMatrix, residueList = AmberCalculator.extractAtomsAndResidues(realTopology, realSnapshot)

        # extract the matrix of the coordinates for the QM residue
        QMResidue = [atomCoords for atomCoords, atomRes in zip(coordsMatrix, residueList) if atomRes == nrResidue]

        # residue map (assignment of each residue to the H,M,L layer)
        residueMap = []
        # loop over the residues of the system
        for iRes in range(max(residueList)):
            # skip QM residue
            if iRes + 1 == nrResidue:
                residueMap.append("H")
                continue
            # initialize the list of the coordinates of the residue iRes+1
            currentMMResidue = []
            # collect coordinates of the residue iRes+1
            for atomCoords, atomRes in zip(coordsMatrix, residueList):
                if atomRes == iRes + 1:
                    currentMMResidue.append(atomCoords)
            # compute distance
            distance = CobrammCalculator.residueDistance(QMResidue, currentMMResidue, "min_min")
            # based on the distance, assign all the atoms of the residue to the M or L layer
            if distance < mobileThreshold:
                residueMap.append("M")
            else:
                residueMap.append("L")
        # now using the residueMap list, define the layer to which each atom belongs
        layerList = [residueMap[resnr - 1] for resnr in residueList]

        # *********************************************************************************************

        # when requested, reorder the residues in such a way that the M residues are concentrated
        # at the beginning of the structure, while the L are at the end
        # this is achieved by exchanging residues with the same structure among themselves, then it
        # works only for a large number of identical residues, as for the case of a solvent
        # data that is processed here is: coordsMatrix and layerList, residueList is not used in the
        # following and is then not reordered
        if reorderRes:

            # define a reordering mask for the residues
            nresidues = max(residueList)
            reorder = [i for i in range(nresidues)]
            # construct a dictionary of lists, with the atoms composing each residue
            residuesAtomsLabels = {i: [] for i in range(nresidues)}
            for lab, res in zip(symbolList, residueList):
                residuesAtomsLabels[res - 1].append(lab)
            # now loop over the residues, from the beginning and from the end at the same time
            for iback in range(nresidues - 1, -1, -1):
                if residueMap[reorder[iback]] == "M":
                    for iforw in range(0, iback, +1):
                        # when we find two identical residues, and the furthest is M and the closest is L, swap them
                        if residueMap[reorder[iforw]] == "L" and \
                                residuesAtomsLabels[reorder[iforw]] == residuesAtomsLabels[reorder[iback]]:
                            reorder[iforw], reorder[iback] = reorder[iback], reorder[iforw]
                            break

            # construct a dictionary of lists, with the indices of the atoms composing each residue
            residuesAtomsIndices = {i: [] for i in range(nresidues)}
            for i, res in enumerate(residueList):
                residuesAtomsIndices[res - 1].append(i)
            # now construct new geometry matrix and new layerList
            # note that by very definition of the reordering, there is no need to reorder that atom labels
            newCoordsMatrix = [[] for _ in range(len(symbolList))]
            newLayerList = [[] for _ in range(len(symbolList))]
            # newResidueList = [[] for _ in range(len(symbolList))]
            # and finally reorder data according to the reorder map
            for srcRes, tarRes in enumerate(reorder):
                for iAtSrc, iAtTar in zip(residuesAtomsIndices[srcRes], residuesAtomsIndices[tarRes]):
                    newCoordsMatrix[iAtTar] = coordsMatrix[iAtSrc]
                    newLayerList[iAtTar] = layerList[iAtSrc]
                    # newResidueList[iAtTar] = residueList[iAtSrc]

            # residueList = newResidueList
            coordsMatrix = newCoordsMatrix
            layerList = newLayerList

        # *********************************************************************************************

        # store the coordinates, the atomic labels and the layer assignment in a Layer instance
        geometry = Layers(symbolList, coordsMatrix, layerList)

        # move back to starting directory and clean up temporary data
        os.chdir(startDir), shutil.rmtree(calcDir)

        # return both topologies (real and modelH) and the layer object
        return geometry, modelTopology, realTopology

    # ================================================================================================================

    @staticmethod
    def extendHlayer(realTopology, geometry, highlayerThreshold):

        # write crd to file
        geometry.makerealcrd()

        # write full topology to files
        with open("real.top", "w") as f:
            f.write(realTopology.topotext)

        if os.path.isfile("model-H.top") == True:
            logwrt.writewarning("overwriting previous model-H.top file")
            os.remove("model-H.top")

        # write input file for parmED
        with open("parmed.script", "w") as f:
            f.write("loadcoordinates real.crd\nstrip !(:1<:{0})\noutparm model-H.top".format(highlayerThreshold))

        # run parmed
        with open("parmedscript.out", "w") as out:
            command = shlex.split("parmed real.top parmed.script")
            subprocess.run(command, stdout=out, stderr=subprocess.STDOUT)

        #modify the real_layers
        with open("model-H.top", "r") as f:
            mh = (f.readlines()[6])
            natoms=int((mh.split())[0])

        standard_text = geometry.reallayertext
        modified_text = ""
        index = 0
        for line in standard_text.splitlines():
            if index in range(natoms):
                if line.split()[5] == "M":
                    modified_text += line.replace(" M", " H") + "\n"
                elif line.split()[5] == "L":
                    modified_text += line.replace(" L", " H") + "\n"
                else:
                    modified_text += line + "\n"
            else:
                modified_text += line + "\n"
            index += 1

        new_layers = Layers.from_real_layers_xyz(modified_text)

        with open("real_layers.xyz", "w") as f:
            f.write(new_layers.reallayertext)

        #return new_layers

    # ================================================================================================================
    @staticmethod
    def residueDistance(QMResidue, MMResidue, distMeasure):

        # compute center of the QM part
        QMCenter = np.array([0.0, 0.0, 0.0])
        for atom in QMResidue:
            QMCenter += np.array(atom)
        QMCenter /= len(QMResidue)

        # compute center of the MM part
        MMCenter = np.array([0.0, 0.0, 0.0])
        for atom in MMResidue:
            MMCenter += np.array(atom)
        MMCenter /= len(MMResidue)

        # distance between geometric centers of the two residues
        if distMeasure == "center_center":

            # return euclidean distance between two two centers
            return np.linalg.norm(QMCenter - MMCenter)

        elif distMeasure == "min_min":

            # find the atom in QMResidue which is closest to MMCenter
            distVec = [np.linalg.norm(np.array(atom) - MMCenter) for atom in QMResidue]
            indexMinQM = distVec.index(min(distVec))

            # find the atom in MMResidue which is closest to QMCenter
            distVec = [np.linalg.norm(QMCenter - np.array(atom)) for atom in MMResidue]
            indexMinMM = distVec.index(min(distVec))

            # return euclidean distance between two two closest atoms
            return np.linalg.norm(np.array(QMResidue[indexMinQM]) - np.array(MMResidue[indexMinMM]))

        else:

            return None


######################################################################################################################

class CobrammInput:
    """This class defines the input keywords for a COBRAMM run. It can be instantiated in two ways,
    either with the standard constructor that sets the values of the keywords based on a series of
    input arguments, or with the readinput method that set the keywords by reading an already
    existing Amber input file. """

    def __init__(self, calc_type="sp", nr_steps=1, qm_charge=0, qm_basis="sto-3g", qm_functional="b3lyp", n_el_states=1, n_root=0, tda=False, pcm=False, implicit_solvent="", charges=False):
        """ Prepare input file by setting reasonable options depending on the input arguments of this function"""

        # the _keywords dictionary contains the entries of the !keyword ... ?keyword section
        self._keywords = {}
        # the _commands dictionary contains all the other sections of the input file
        self._command = {}

        # define general part of the cobram.command
        if calc_type == "sp":
            self._keywords["type"] = "optxg"
            self._keywords["nsteps"] = "sp"
            self._keywords["qm-type"] = "gauss"
            self._keywords["qmem"] = "2000MB"
        elif calc_type == "optxg":
            self._keywords["type"] = "optxg"
            self._keywords["nsteps"] = "{0}".format(nr_steps)
            self._keywords["qm-type"] = "gauss"
            self._keywords["qmem"] = "5000MB"
            self._keywords["geomem"] = "5000MB"
        elif calc_type == "freqxg":
            self._keywords["type"] = "freqxg"
            self._keywords["qm-type"] = "gauss"
            self._keywords["qmem"] = "500MB"
            self._keywords["geomem"] = "500MB"
            self._keywords["savnum"] = "1"
        elif calc_type == "tsh":
            self._keywords["type"] = "mdv"
            self._keywords["qm-type"] = "gauss"
            self._keywords["qmem"] = "1000MB"
            self._keywords["geomem"] = "1000MB"
            self._keywords["nsteps"] = "{0}".format(nr_steps)
            self._keywords["tstep"] = "0.5"
            self._keywords["surhop"] = "persico"
            self._keywords["nacs"] = "tdc"
            self._keywords["tdctype"] = "2"
            self._keywords["hoptogs"] = "0"
            self._keywords["backhop"] = "1"
            self._keywords["velafterhop"] = "1"


        # define sander section
        self._command["sander"] = "comment line\n&cntrl\nimin = 1,\nmaxcyc = 0,\nntb = 0,\nigb = 0,\n" \
                                  "ntr = 0,\nibelly = 1,\ncut = 10\n/"

        # define section for QM with Gaussian
        if tda:
            if n_el_states > 1 and n_root == 0:
                tddirective = "tda(nstates={0}) IOp(9/40=16)".format(n_el_states)
            elif n_el_states > 0 and n_root > 0:
                tddirective = "tda(nstates={0}, root={1}), density=current, iop(6/7=3) IOp(9/40=16)".format(n_el_states, n_root)
            else:
                tddirective = ""
        else:
            if n_el_states > 1 and n_root == 0:
                tddirective = "td(nstates={0})".format(n_el_states)
            elif n_el_states > 0 and n_root > 0:
                tddirective = "td(nstates={0}, root={1}), density=current, iop(6/7=3)".format(n_el_states, n_root)
            else:
                tddirective = ""

        if pcm:
            pcmdirective = "scrf(solvent={})".format(implicit_solvent)
        else:
            pcmdirective = ""

        if charges:
            chargesdirective = "Pop=MK IOp(6/33=2,6/41=10,6/42=17)"
        else:
            chargesdirective = ""

        self._command["gaussian"] = "#p {0} {1} nosym {2} {3} {4}\n\n".format(qm_basis, qm_functional, tddirective, pcmdirective, chargesdirective) + \
                                    "gaussian input generated by COBRAMM \n\n{0} 1".format(qm_charge)

        #iself._command["gaussian"] = "#p {0} {1} nosym {2}\n\n".format(qm_basis, qm_functional, tddirective) + \
                                   # "gaussian input generated by COBRAMM \n\n{0} 1".format(qm_charge)

        # define optimizer options
        if calc_type == "optxg":
            self._command["optxg"] = "#p opt=(Loose, Z-matrix)"

        # define frequency calculation options
        if calc_type == "freqxg":
            self._command["optxg"] = "#p freq(numer,step=10) iop(7/10=8)"

    # ================================================================================

    def __del__(self):
        """Destroy the AmberInput data"""
        del self._keywords
        del self._command

    # ================================================================================

    def __getitem__(self, item):
        """Returns variables stored in the self.input dictionary with the
        "item" key, by using the syntax self[item]. """
        return self._keywords[item]

    def __setitem__(self, key, value):
        """Change the value of a directive of the sander input file: if the key already
        exists, its value will be changed, otherwise, the new key will be added"""
        self._keywords[key] = value

    # ================================================================================

    @property
    def inputtext(self):
        """ Return a string with sander input file constructed with
        the directives stored in the self instance of the AmberInput class

        :return: string with the input file text
        """

        # add the keyword section (which should always be present)
        txt = "!keyword\n"
        for key, value in self._keywords.items():
            txt += "{0}={1}\n".format(key, value)
        txt += "?keyword\n\n"

        # add the other sections
        for command, value in self._command.items():
            txt += "!{0}\n{1}\n?{0}\n\n".format(command, value)

        return txt

    # ================================================================================

    @classmethod
    def readinput(cls, inputtext):
        """Parse the text of a COBRAMM input file, extracting the keywords and their value
        to set up a new instance of the AmberInput class

        :param inputtext: COBRAMM input file text
        :return: new CobrammInput instance with the input file directives """

        # initialize dictionary to save results of the parsing
        cobrammdirectives = {}

        # TODO: implement the parsing of the COBRAMM input file
        if not cobrammdirectives:
            raise CobrammError("the CobrammInput.readinput() method is not yet implemented")

        # create a new instance of CobrammInput and overwrite the dictionary with the keywords
        newinput = cls()
        newinput._keywords = cobrammdirectives
        # return the CobrammInput instance
        return newinput


#####################################################################################################################

class CobrammOutput:

    def __init__(self, out_dir, store_files=False):
        """ Constructor of the CobrammOutput instance. Needs as argument the path of the directory where the
        calculation is stored. The logical variables storeFiles defines whether the COBRAMM files are removed or not
        when cleaning the CobrammOutput instance. It should be set to true when the user is interested in keeping
        the original COBRAMM files for future analysis of the results.

        :param out_dir:  path of the directory where the calculation is stored
        :param store_files: COBRAMM files are removed or not when cleaning the AmberOutput instance
        """

        # initialize storeFiles to True, in this way when the constructor dyes unexpectedly we can save the files
        self.storeFiles = True

        # store the path of the directory where the COBRAMM files are read
        self.outDir = out_dir

        # parse the COBRAMM xml file
        start_dir = os.getcwd()
        os.chdir(out_dir)
        self.xmloutput = Output(parse=True)
        os.chdir(start_dir)

        # when the calculation is a single point, extract the TD electronic state information when available
        if len(self.xmloutput.grep_all('step')) == 1:
            self.elstates = self._states_from_gaussian_tddft_out(os.path.join(out_dir, "QM_data", "qmALL.log"))
        else:
            self.elstates = None

        # extract the charges from the qmALL.log file
        if os.path.exists(os.path.join(out_dir, "QM_data", "qmALL.log")):
            self.qmcharges = self._atomcharges_from_gaussian_out(os.path.join(out_dir, "QM_data", "qmALL.log"))

        # rewrite the geometry.xyz to a 3D vector with the coordinates
        if os.path.exists(os.path.join(out_dir, "QM_data", "qmALL.log")):
            self.coords= self._xyzgeometry_to_vector(os.path.join(out_dir, "geometry.xyz"))

        # read the masses of the atoms from the geometry.log file
        if os.path.exists(os.path.join(out_dir, "geometry.log")):
            self.coord_masses = self._atommasses_from_gaussian_out(os.path.join(out_dir, "geometry.log"))

        # when the normal modes have been computed, extract the raw force matrix from the geometry.chk file
        if self.xmloutput.normal_modes:

            # to use the tools from gaussian to read the chk file, we need to have the gaussian environment set
            cobrammenv.setCobrammProfile()
            # check that the environment for gaussian is actually defined
            env_defined, error_msg = cobrammenv.checkGaussianOptEnv()
            if not env_defined:
                logwrt.fatalerror(error_msg)

            os.chdir(out_dir)
            # format chk file to read the force matrix
            subprocess.call(["gunzip", "-kf", "geometry.chk.gz"])
            subprocess.call(["formchk", "geometry.chk"])
            # open formatted chk file and store its content
            with open("geometry.fchk") as f:
                chk_text = []
                for line in f:
                    chk_text.append(line)

            # process the gaussian chk file to extract the force constant matrix
            self.force_matrix,self.coord_masses = self._process_chk(chk_text)
            os.chdir(start_dir)
        
        # store the value of storeFiles
        self.storeFiles = store_files

    # ================================================================================

    def __del__(self):
        """ destroy the data stored in the class instance, and possibly
        remove the directory of the calculation unless it is requested otherwise """

        if not self.storeFiles:
            shutil.rmtree(self.outDir)

        del self.outDir
        del self.storeFiles
        del self.xmloutput

    # ================================================================================

    def snapshot(self, nr_snap=-1):
        """ extract one of the snapshots stored in the scratch data of the COBRAMM calculation """

        # set the step to extract, when the last step is requested
        if nr_snap <= 0:
            nr_snap = self.xmloutput.steps

        # fetch the wanted step
        step = self.xmloutput.get_step(nr_snap)
        # extract the Layers instance from the chosen step
        return step.geometry

    # ================================================================================

    def eletronicspectrum(self, spectral_grid, width=0.01):
        """return the electronic spectrum computed by convoluting vertical transitions with a gaussian linefunction"""

        # when the electronic transitions have been computed
        if self.elstates:
            # initialize array with the value of the function on the grid
            spec_value = np.zeros(len(spectral_grid))
            # use data from TD-DFT calculation extracted from gaussian to compute spectrum
            for center, intensity in zip(*self.elstates):
                spec_value += intensity * self._line_function(spectral_grid, center, width)
            return spec_value

        # otherwise return none
        else:
            return None

    # ================================================================================

    @staticmethod
    def _states_from_gaussian_tddft_out(file_name, decomposition=False, state=-1):
        """ use regexpr to extract the list of electronic state energies and the corresponding oscillator strengths"""

        # open output file and read content
        with open(file_name, "r") as f:
            g_out_text = f.read()

        if decomposition:
            if state < 10:
                re_string = "Excited State   {}: *.*? *([0-9]*\.[0-9]*) *eV *[0-9]*\.[0-9]* *nm *" \
                    "f=([0-9]*\.[0-9]*) *\<S\*\*2\>=[0-9]*\.[0-9]*".format(state)
            elif state > 9 and state < 99:
                re_string = "Excited State  {}: *.*? *([0-9]*\.[0-9]*) *eV *[0-9]*\.[0-9]* *nm *" \
                    "f=([0-9]*\.[0-9]*) *\<S\*\*2\>=[0-9]*\.[0-9]*".format(state)

        else:
            # excited state string
            re_string = "Excited State *[0-9]*: *.*? *([0-9]*\.[0-9]*) *eV *[0-9]*\.[0-9]* *nm *" \
                    "f=([0-9]*\.[0-9]*) *\<S\*\*2\>=[0-9]*\.[0-9]*"
        # extract results with regular expression
        results = re.findall(re_string, g_out_text)

        # define and return lists of electronic state energies (convert from eV to au) and oscillator strengths
        if results:
            el_energies = [float(i[0]) / constants.Hartree2eV for i in results]
            osci_strength = [float(i[1]) for i in results]
            return el_energies, osci_strength
        else:
            return None

    # ================================================================================

    @staticmethod
    def _line_function(egrid, center, width):
        """gaussian line function, for convoluting vertical lines to create an electronic spectrum"""

        linef = [1. / (width * math.sqrt(2.0 * math.pi)) *
                 np.exp(-(e - center) ** 2 / (2.0 * width ** 2)) for e in egrid]
        return np.array(linef)

    # ================================================================================

    @staticmethod
    def _atommasses_from_gaussian_out(filename):
        """ use regexpr to extract the list of frequencies and atomic weights from GAUSSIAN output file """

        with open(filename, "r") as f:
            file_text = f.read()

        w_string = "AtmWgt=(.*)"
        w_find = re.findall(w_string, file_text)
        atom_masses = [float(m) * constants.amu2au for m in ' '.join(w_find).split()]
        coord_masses = sum([[m, m, m] for m in atom_masses], [])

        return coord_masses

    # ================================================================================

    @staticmethod
    def _atomcharges_from_gaussian_out(filename):
        """ use regexpr to extract the list of charges from GAUSSIAN output file """

        QMcharges = []
        with open(filename) as cf:
            ESPchargesec = False
            for line in cf:
                if "Sum of ESP charges =" in line:
                    ESPchargesec = False
                if ESPchargesec and len(line.split()) == 3:
                    QMcharges.append(float(line.split()[2]))
                if "ESP charges:" in line:
                    ESPchargesec = True



        return QMcharges

    # ================================================================================

    @staticmethod
    def _xyzgeometry_to_vector(filename):
        with open(filename) as f:
            f.readline(), f.readline()
            atomlabels, coords = [], []
            for ln in f:
                if ln.split():
                    atomlabels.append(ln.split()[0])
                    coords.append([float(f) for f in ln.split()[1:]])

        return coords
    # ================================================================================

    @staticmethod
    def _process_chk(chkfile_text, mask=None):

        # extract masses from chk file
        masslines = 0
        masses_list = []
        for line in range(len(chkfile_text)):
            # number of lines to read (assumed fixed format with 5 elements per line)
            if chkfile_text[line].find("Vib-AtMass") != -1:
                masslines = int(math.ceil(float(chkfile_text[line].split()[3])/5))
            elif chkfile_text[line].find("Real atomic weights") != -1:
                masslines = int(math.ceil(float(chkfile_text[line].split()[5])/5))
            if masslines > 0:
                for massline in range(line+1, line+masslines+1):
                    for element in chkfile_text[massline].split():
                        masses_list.extend([float(element)* constants.amu2au, float(element)* constants.amu2au, float(element)* constants.amu2au])
                break

        # extract section of chk file that contains the force constants matrix
        force_list = []
        for line in range(len(chkfile_text)):
            if chkfile_text[line].find("Cartesian Force Constants") != -1:
                # number of lines to read (assumed fixed format with 5 elements per line)
                fclines = int(math.ceil(float(chkfile_text[line].split()[5])/5))
                # number of force matrix elements
                ntriang = int(chkfile_text[line].split()[5])
                nmatrix = int((math.sqrt(8 * ntriang + 1) - 1) / 2)
                # now extract the list of force constants from the text section
                for fcline in range(line+1, line+fclines+1):
                    for element in chkfile_text[fcline].split():
                        force_list.append(float(element))

        # initialize matrix and store with 2 indices the values that have been read from chk file
        mat = np.zeros([nmatrix, nmatrix])
        n = 0
        for i in range(nmatrix):
            for j in range(0, i + 1):
                mat[i, j] = force_list[n]
                if i != j:
                    mat[j, i] = force_list[n]
                n += 1

        # if mask is defined, return the corresponding section of the matrix
        if mask:
            return np.take(np.take(mat, mask, axis=0), mask, axis=1),masses_list
        else:  # otherwise, return the full matrix
            return mat,masses_list


######################################################################################################################

if __name__ == '__main__':
    # test run of the classes contained in this module

    # create an instance of the COBRAMM and Amber interfaces
    AmberCalculator()
    CobrammCalculator()

    # use files from gaussian_036 example calculations
    testdir = os.path.join(os.environ['COBRAM_PATH'], "test", "gaussian_036")
    with open(os.path.join(testdir, CobrammCalculator.MODELHTOP)) as ff:
        testmodelhtop = AmberTopology(ff.read())
    with open(os.path.join(testdir, CobrammCalculator.REALTOP)) as ff:
        testrealtop = AmberTopology(ff.read())
    with open(os.path.join(testdir, CobrammCalculator.REALLAYERS)) as ff:
        testgeometry = Layers.from_real_layers_xyz(ff.read())

    # define cobram.command for a ground state optimization with DFT
    cobcomm = CobrammInput(calc_type="optxg", nr_steps=7)

    # run the cobramm optimization and extract the optimized geometry
    optresult = CobrammCalculator.run(cobcomm, testgeometry, testmodelhtop, testrealtop,
                                      calc_dir="optimization", store=True)
    optgeometry = optresult.snapshot(-1)

    # # define cobram.command for a single point calculation, with TD-DFT and 5 electronic states
    cobcomm = CobrammInput(n_el_states=5)

    # run the cobramm optimization and extract the optimized geometry
    tddftresult = CobrammCalculator.run(cobcomm, optgeometry, testmodelhtop, testrealtop,
                                        calc_dir="tddft", store=True)

    # now process the electronic state results to print the spectrum
    grid = np.linspace(0., 0.5, 200)
    spect = tddftresult.eletronicspectrum(grid)
    with open("spectrum.dat", "w") as ff:
        for en, ints in zip(grid, spect):
            ff.write("{} {}\n".format(en, ints))
