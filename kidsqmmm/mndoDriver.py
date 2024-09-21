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

######################################################################################################################

# import statements of module from python standard library

import os  # filesystem utilities
import shutil  # filesystem utilities
import re  # process output with regular expressions
import math  # import mathematical functions
import pprint
from copy import deepcopy

# imports of local modules

import logwrt  # manages log file output + start/end procedures
import constants  # values of physical constants and conversion factors

# import of local classes
from QMOutput import QMOutput

class MNDOInput:
    """Store the data for the execution of MNDO and process it to check its consistency"""

    def __init__(self, keys, mndoadd, coords, symbols, otheropt):
        """Constructor of the mndoianInput class, requires these arguments:
        - keys: string with the text of the mndoian inp file from the route section to the charge/mult line
        - mndoadd: other options (e.g. the weights of the state averaged CASSCF)
        - coords: coordinates of the atoms
        - symbols: symbols of the atoms
        - otheropt: dictionary with other required options for the QM calculation
                    ("forces", "restart", "suppress_nosymm_error" ...)
        -
        """

        self.keys, self.mndoadd = None, None
        self.coords, self.symbols = None, None
        # dictionary with other options
        self.otheropt = otheropt
        self.calctype = 'OM2-DEBUG'

        # store the sections of the mndoian input files already defined
        self.keys = keys
        self.mndoadd = mndoadd
        # store the information about the molecular geometry
        self.coords = coords
        self.symbols = symbols

        for key in ["iroot", "restart", "mminp", "mmcoup", "jop"]:
            if key in otheropt:
                self.otheropt[key] = otheropt[key]
            else:
                self.otheropt[key] = None

    # =============================================================================================================

    def fileText(self, inpFileName, memory="500MB", nproc=1, gversion="mndo2020"):
        """Prepare the text of mndoian input file named after the inpFileName variable (input file
        should be named <inpFileName>.com, so inpFileName should not include the extension!), using the
        input data defined in the instance of mndoianInput. The memory (var memory) and
        the number of cores (var nproc) to use in the execution and the version of MNDO (gversion)
         are given as input because they are decided when running the MNDO calculation and are not
         considered parameters of the QM calculation"""

        ncharge = 0
        if self.otheropt['charges']:
            chrg = self.otheropt["charges"][1]
            ncharge = len(chrg)

        # initialize the text of the input file
        inputText = ""
        # add the mndoian keys using the variable self.keys up to the first empty line
        end_idx = 0
        keys_copy = deepcopy(self.keys)
        for i in range(len(keys_copy)):
            if self.keys[i].strip()[-1] != '+':
                keys_copy[i] = keys_copy[i].split()
                end_idx = i+1 
                break
            else:
                keys_copy[i] = keys_copy[i].split('+')[0]
                keys_copy[i] = keys_copy[i].split()
        keys_copy = keys_copy[0:end_idx]

        ncigrd = 0
        iroot = 0
        all_keys = []
        for i in range(len(keys_copy)):
            for j in range(len(keys_copy[i])):
                k,v = keys_copy[i][j].split('=')
                all_keys += [k.upper()]
                if k.upper() == 'NCIGRD':
                    ncigrd = int(v)
                if k.upper() == 'IROOT':
                    iroot = int(v)
                if k.upper() == 'LROOT':
                    lroot = int(v)
                if k.upper() == 'NUMATM':
                    keys_copy[i][j] = 'NUMATM=%d'%(ncharge)

        if 'NSAV7' not in all_keys:
            keys_copy[0] += ['NSAV7=7']
        if 'NSAV13' not in all_keys:
            keys_copy[0] += ['NSAV13=2']
        if 'NSAV15' not in all_keys:
            keys_copy[0] += ['NSAV15=4']
        if 'MMINP' not in all_keys:
            keys_copy[1] += ['MMINP=2']
        if 'MMCOUP' not in all_keys:
            keys_copy[1] += ['MMCOUP=2']
        if 'NUMATM' not in all_keys:
            keys_copy[1] += ['NUMATM=%d'%ncharge]
        if 'MPRINT' not in all_keys:
            keys_copy[2] += ['MPRINT=5']
        if 'NPRINT' not in all_keys:
            keys_copy[2] += ['NPRINT=5']
        if 'IUVCD' not in all_keys:
            keys_copy[2] += ['IUVCD=2']

        for i in range(len(keys_copy)):
            keys_copy[i] = ' '.join(keys_copy[i])
        inputText += ' +\n'.join(keys_copy)
        inputText += '\n\n'
        inputText += 'autogenerated for mndo\n'

        
        support_dict = {'C':6, 'H':1, 'O':8, 'N':7}
        znumber = [support_dict[i.upper()] for i in self.symbols]

        # write molecular structure
        for atom in zip(znumber, *self.coords):
            inputText += "{0} {1:16.8f} 0 {2:16.8f} 0 {3:16.8f} 0\n".format(*atom)
        inputText += "{0} {1:16.8f} 0 {2:16.8f} 0 {3:16.8f} 0\n".format(0,0,0,0)
        inputText += ' '.join(['%d'%i for i in range(1,1+ncigrd)]) + '\n'
        # inserting charges
        if self.otheropt['charges']:
            coords = self.otheropt["charges"][0]
            chrg = self.otheropt["charges"][1]
            for atom in zip(chrg, *coords):
                inputText += "{1:14.8f} {2:14.8f} {3:14.8f} {0:16.8f}\n".format(*atom)
            inputText += "\n"

        # return the complete text of the input file
        return inputText

###################################################################################################################

class MNDOOutput(QMOutput):

    def __init__(self, name, calcdir, calctype, SPcalc=False, verbose=False):
        # define options of the base class
        QMOutput.__init__(self)

        self.dataDict = {}

        # FILENAMES, DIRECTORIES, INPUT/OUTPUT
        self.dataDict["name"] = name  # name of the calc, it is the base name of the I/O files (.com, .chk, .log)
        self.dataDict["dir"] = calcdir  # name of the directory where the output files are stored
        self.dataDict["outfile"] = None   # string with the text of the MNDO log file

        # ADDITIONAL PHYSICAL PROPERTIES EXTRACTED FROM THE OUTPUT FILE
        self.dataDict["charges"] = []   # charges / population analysis on QM atoms 
        self.dataDict["elfield"] = {}   # electric field at the position of the point charges (/MM atoms)
        self.dataDict["dipole"] = []    # dipole of system in XYZ and total
        self.dataDict["osc_strength"] = {}   # osc_strength to each electric state
        self.dataDict["natoms"] = None   # number of atoms in the QM calculation
        self.dataDict["nroots"] = 0     # number of roots (/electric states) in the QM calculation
        self.dataDict["optstate"] = 0   # actually == lroot
        
        #
        self.dataDict["energy"] = {} # dictionary of state energy
        self.dataDict["gradient"] = {} # dictionary of state gradients
        self.dataDict["hess"] = {} # dictionary of state hess
        self.dataDict["nac"] = {} # dictionary of NACs
        self.dataDict["fullgradcharge"] = {} # grad on charge
        self.dataDict["fullnaccharge"] = {} # NACs on charge
        self.dataDict["gradcharges_extraterm"] = False # if embedding is performed by QM

        # INFO ON MNDO EXECUTION
        self.dataDict["termination"] = None   # final termination: 0 = Normal termination, 1 = Error termination
        self.dataDict["errormsg"] = []   # in case of Error termination, error message
        self.dataDict["signs"] = [] # list of WF signs from previous step to apply in phase correction
        self.dataDict["scf_mo_maps"] = None
        self.dataDict["psioverlap"] = None # array of the WF overlaps between consecutive steps
        
        self.dataDict["eigenvectors"] = None
        #self.dataDict["SS_root_order"] = None
        self.calctype = calctype # --> can be "OMx-MRCI"

        # store the mndo log removing AO and MO definitons, CI expansions, etc.
        # if SPcalc == False and verbose == False:
        #     self.dataDict["outfile"] = self._filterOutput(os.path.join(calcdir, name + ".log"))
        # else:
        #     with open(os.path.join(calcdir, name + ".log")) as fout:
        #         self.dataDict["outfile"] = fout.read()

        # use full file for reading properties 
        with open(os.path.join(calcdir, name + ".log")) as fout:
            self.dataDict["outfile"] = fout.read()
            # print(self.dataDict["outfile"])

        # loop over the last lines of the log file
        output = self.dataDict["outfile"].splitlines()
        if 'ERROR' in ''.join(output[-11:]).upper():
            self.dataDict["termination"] = 1
            self.dataDict["errormsg"] = ''.join(output[-11:])
            return

        au_2_ang = 5.291772104260590e-01
        au_2_kcal_1mea = 6.275094742071363e+02
        au_2_kcal_1mea_per_ang =  1.185821047928172e+03
        au_2_ev =  2.7211386033e+01
        au_2_wn =  2.1947463147e+05

        # initialize variables
        tempstate = None
        readOvMat = False        
        FLAGCHARGE = False
        ci_currect_state = None
        states_for_NAC = None
        readtimes_for_grad = 0
        full_grad = False # include MM Atoms
        full_nac = False  # include MM Atoms
        nstate = 0 # counter for total states
        i = 0
        while i < len(output):
            # determine the current state
            if "CI CALCULATION FOR STATE:" in output[i]:
                # print('1')

                ci_currect_state = int(output[i].split()[-1]) - 1
                i += 1

            # determine the current stated for analytical NAC
            elif "CI CALCULATION FOR INTERSTATE COUPLING OF STATES:" in output[i]:
                # print('2')

                istate, jstate = map(int, output[i].split()[-2:])
                istate -= 1; jstate -= 1
                states_for_NAC = (istate, jstate)

                i += 1

            # get energy from CI calculation
            elif "SUMMARY OF MULTIPLE CI CALCULATIONS" in output[i]:
                # print('3')

                k = i + 5
                nstate = 0
                while output[k].strip()[0:2] != '--':
                    self.dataDict["energy"][nstate] = float(output[k].split()[2]) / au_2_kcal_1mea
                    nstate += 1
                    k += 1
                i = k + 1

            # get gradient or nac (including MM points)
            elif "GRADIENTS (KCAL/(MOL*ANGSTROM))" in output[i]:
                # print('4')

                k = i + 4
                iatom = 0
                skline = 0
                gradient = [[],[],[]]
                gradientcharge = [[],[],[]]
                while output[k].strip() != '':
                    element = output[k].split()[5:]
                    gradient[0].append(float(element[0]) / au_2_kcal_1mea_per_ang )
                    gradient[1].append(float(element[1]) / au_2_kcal_1mea_per_ang)
                    gradient[2].append(float(element[2]) / au_2_kcal_1mea_per_ang)

                    k += 1

                for _ in range(10):
                    if 'EXTERNAL POINT CHARGES' in output[k]:
                        self.dataDict["gradcharges_extraterm"] = True
                        while "GRADIENTS (KCAL/(MOL*ANGSTROM))" not in output[k]: k += 1
                        k += 4
                        while output[k].strip() != '':
                            element = output[k].split()[5:]
                            gradientcharge[0].append(float(element[0]))
                            gradientcharge[1].append(float(element[1]))
                            gradientcharge[2].append(float(element[2]))
                            k += 1
                        break
                    k += 1

                if readtimes_for_grad == ci_currect_state:
                    self.dataDict["gradient"][ci_currect_state] = gradient      
                    if self.dataDict["gradcharges_extraterm"]:
                        self.dataDict["fullgradcharge"][ci_currect_state] = gradientcharge       
                else:
                    try:
                        self.dataDict["nac"][states_for_NAC[0]][states_for_NAC[1]] = gradient
                    except KeyError:
                        self.dataDict["nac"][states_for_NAC[0]] = {}
                        self.dataDict["nac"][states_for_NAC[0]][states_for_NAC[1]] = gradient
                    try:
                        self.dataDict["nac"][states_for_NAC[1]][states_for_NAC[0]] = [[-x for x in y] for y in gradient]
                    except KeyError:
                        self.dataDict["nac"][states_for_NAC[1]] = {}
                        self.dataDict["nac"][states_for_NAC[1]][states_for_NAC[0]] = [[-x for x in y] for y in gradient]
                    if self.dataDict["gradcharges_extraterm"]:          
                        try:
                            self.dataDict["fullnaccharge"][states_for_NAC[0]][states_for_NAC[1]] = gradientcharge
                        except KeyError:
                            self.dataDict["fullnaccharge"][states_for_NAC[0]] = {}
                            self.dataDict["fullnaccharge"][states_for_NAC[0]][states_for_NAC[1]] = gradientcharge
                        try:
                            self.dataDict["fullnaccharge"][states_for_NAC[1]][states_for_NAC[0]] = [[-x for x in y] for y in gradientcharge]
                        except KeyError:
                            self.dataDict["fullnaccharge"][states_for_NAC[1]] = {}
                            self.dataDict["fullnaccharge"][states_for_NAC[1]][states_for_NAC[0]] = [[-x for x in y] for y in gradientcharge]


                readtimes_for_grad += 1
                i = k + 1
            
            # get (mulliken) charges on QM (HighLayer) atoms (need nprint>=2)
            elif "MULLIKEN POPULATION ANALYSIS." in output[i]:
                # print('5')

                self.dataDict["charges"] = []
                k = i + 4 # skip 4 lines
                while output[k].strip() != '':
                    self.dataDict["charges"].append(float(output[k].split()[2]))
                    k += 1
                i = k + 1

            # get electrostatic field at MM point charges (need mminp = 1, or revised mndo-code with mmcoup>=3)
            elif "ELECTROSTATIC POTENTIAL (VOLT) AND ELECTRIC FIELD COMPONENTS" in output[i]:
                # print('6')

                self.dataDict["elfield"] = {}
                k = i + 4 # skip 4 lines
                elfield = [], [], []
                while output[k].strip() != '':
                    element = output[k].split()
                    Ex, Ey, Ez = float(element[2]), float(element[3]), float(element[4])
                    elfield[0].append(Ex), elfield[1].append(Ey), elfield[2].append(Ez)
                    k += 1
                # MNDO only provide elfield for ground state, here we may copy it to all states
                self.dataDict["elfield"][ci_currect_state] = elfield
                i = k + 1

            if "Properties of transitions   1 -> #" in output[i]:
                self.dataDict["osc_strength"] = {}
                
                k = i + 3
                while output[k].strip() != '':
                    terms = output[k].strip().split()
                    to_state = int(terms[0]) - 1
                    read_pos = -2 # read f_rp only
                    self.dataDict["osc_strength"][to_state] = float(terms[read_pos])
                    k += 1

                i = k + 1

            # get dipole moments
            elif " DIPOLE              X           Y           Z         TOTAL" in output[i]:
                # print('7')

                k = i + 4 # skip 4 lines
                dip = output[k].split()
                self.dataDict["dipole"] = \
                                        float(dip[1]) * constants.Debye2AU, \
                                        float(dip[2]) * constants.Debye2AU, \
                                        float(dip[3]) * constants.Debye2AU, \
                                        float(dip[4]) * constants.Debye2AU

                i = k + 1

            elif "     INPUT GEOMETRY" in output[i]:
                k = i + 6
                xyz_local = []
                while output[k].strip() != '':
                    terms = output[k].strip().replace('*', ' ').split()
                    # print(terms)
                    xyz_local += [float(terms[2]),float(terms[3]), float(terms[4])]
                    k += 1
                self.dataDict['inputXYZ'] = xyz_local
                i = k + 1

            elif "State  1,  Mult. 1,  E-E(1)=  0.00000" in output[i]:
                Elmin = float(output[i].split()[-2]) / au_2_ev
                self.dataDict['Elmin'] = Elmin    
                i += 1

            elif "GRADIENT NORM =" in output[i]:
                self.dataDict["gradNORM"] = float(output[i].split()[3])
                i += 1

            elif 'EIGENVECTORS OF THE MASS-WEIGHTED' in output[i]: # just in AU

                k = i + 3
                colidx = np.array(output[k].split()).astype(np.int32) - 1
                k += 2
                freq[colidx] = np.array(output[k].split()).astype(np.float64) / au_2_wn
                k += 2
                while 'CARTESIAN DISPLACEMENT' not in output[k]: 
                    if output[k].strip() == '': 
                        k += 1
                        continue
                    terms = output[k].split()
                    if len(terms) == len(colidx) + 1:
                        rowidx = int(terms[0]) - 1
                        Tmod[rowidx, colidx] = np.array(terms[1:]).astype(np.float64)
                        k += 1
                    if len(terms) <= len(colidx):
                        colidx = np.array(terms).astype(np.int32) - 1
                        k += 2
                        terms = output[k].split()
                        freq[colidx] = np.array(terms).astype(np.float64)
                        k += 1
                
                self.dataDict['freq'] = freq 
                self.dataDict['Tmod'] = Tmod
                # self.dataDict['hess'] = np.einsum('ik,k,kj->ij', Tmod, freq, Tmod.T)
                i = k + 1

            elif "ERROR" in output[i]:
                self.dataDict['errormsg'] += [output[i]]
                i = i + 1

            elif "TERMINAT" in output[i]: # failed by termination
                self.dataDict['termination'] = 1
                break

            # terminate the reading of the main output file
            elif "COMPUTATION TIME" in output[i]:
                self.dataDict['termination'] = 0
                self.dataDict['errormsg'] = []
                # print('8')
                break

            else:
                i += 1

        # pprint.pprint(self.dataDict)
        # pprint.pprint({'energy':self.dataDict['energy']})
        # pprint.pprint({'gradient':self.dataDict['gradient']})
        # pprint.pprint({'nac':self.dataDict['nac']})
        # pprint.pprint({'hess':self.dataDict['hess']})
        # pprint.pprint({'elfield':self.dataDict['elfield']})
        # pprint.pprint({'charges':self.dataDict['charges']})
        # pprint.pprint({'dipole':self.dataDict['dipole']})
        # pprint.pprint({'osc_strength':self.dataDict['osc_strength']})
        

    # =============================================================================================================

    def __del__(self):
        """Destructor for the mndoianOutput class: not only the memory allocated for the oject attributes
        (the dictionary self.dataDict) needs to be released, also the MNDO I/O files stored on disk
        can be safely removed... """

        # permanently remove directory and mndo files
        # try:
            # if not logwrt.DEBUG_COBRAMM_RUN:
                # if self.dataDict["termination"] == 1:
                    # source=self.dataDict["dir"]+"/mndo-QM.log"
                    # shutil.move(source,"mndo-QM_err.log")
                # shutil.rmtree(self.dataDict["dir"])
        # except FileNotFoundError:
            # pass
        # clean up the log
        del self.log
        # destroy the dictionary that contains the output data
        del self.dataDict    

import sys 

if __name__ == "__main__":
    # check mndoDriver's parser 

    name = sys.argv[1]
    test = MNDOOutput(name, '.', 'OM2-MRCI')
