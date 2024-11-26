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

import os, sys  # filesystem utilities
import shutil   # filesystem utilities
import re  # process output with regular expressions
import math  # import mathematical functions
import pprint
from copy import deepcopy
import numpy as np

# imports of local modules

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from kids_env import which, source
import kids_log  # manages log file output + start/end procedures
import constants  # values of physical constants and conversion factors
from QMOutput import QMOutput


class BdfInput:
    """Store the data for the execution of Bdf and process it to check its consistency"""

    def __init__(self, keys, coords, symbols, otheropt):
        """Constructor of the bdfInput class, requires these arguments:
        - keys: string with the text of the bdf inp file from the route section to the charge/mult line
        - coords: coordinates of the atoms
        - symbols: symbols of the atoms
        - otheropt: dictionary with other required options for the QM calculation
                    ("forces", "restart", "suppress_nosymm_error" ...)
        -
        """

        # store the sections of the bdf input files already defined
        self.keys = keys
        self.coords = coords
        self.symbols = symbols
        self.otheropt = otheropt
        self.calctype = 'TDDFT'

    @staticmethod
    def checkEnv():
        try:
            bdfscript = os.environ['BDF_SCRIPT']  # path of the MOLCAS launch script
            bdfpath = os.environ['BDFHOME']  # path of the MOLCAS code
        except KeyError:
            return False, "environment variables for BDF, $BDF_SCRIPT, $BDFHOME " \
                          "are not defined"

        # bdf script exists and has x permissions
        if not os.path.isfile(bdfscript) or not os.access(bdfscript, os.X_OK):
            return False, f'bdf launch script {bdfscript} not found, please check the env ' \
            'variable BDF_SCRIPT!'
        # bdf directory exists
        if not os.path.isdir(bdfpath):
            return False, 'bdf directory {} not found, please check the env variable BDFHONE!'.format(bdfpath)

        # return True if all the checks were OK
        return True, bdfpath


    def fileText(self, inpFileName, memory="500MB", nproc=1, version="bdfdrv.py"):
        """Prepare the text of bdf input file named after the inpFileName variable (input file
        should be named <inpFileName>.com, so inpFileName should not include the extension!), using the
        input data defined in the instance of bdfInput. The memory (var memory) and
        the number of cores (var nproc) to use in the execution and the version of Bdf (gversion)
         are given as input because they are decided when running the Bdf calculation and are not
         considered parameters of the QM calculation"""

        ncharge = 0
        if 'charges' in self.otheropt and self.otheropt['charges']:
            ncharge = len(self.otheropt["charges"][1])

        # initialize the text of the input file
        inputText = ""
        keys_copy = self.keys.split('\n')
        
        for l in keys_copy:
            if l == '$COORD_XYZ':
                for atom in zip(self.symbols, *self.coords):
                    inputText += '%s   %12.8e  %12.8e  %12.8e\n'%(atom[0], atom[1], atom[2], atom[3])
            elif l == '$EMBED_XYZ' and 'charges' in self.otheropt and self.otheropt['charges']:
                self.embedQ = True
                coords = self.otheropt["charges"][0]
                chrg = self.otheropt["charges"][1]
                for atom in zip(chrg, *coords):
                    inputText += '%12.8e  %12.8e  %12.8e  %12.8e\n'%(atom[0], atom[1], atom[2], atom[3])
            else:
                inputText += l + '\n'
        
        return inputText

###################################################################################################################

class BdfOutput(QMOutput):

    def __init__(self, name, calcdir, SPcalc=False, verbose=False):
        # define options of the base class
        QMOutput.__init__(self)

        self.dataDict = {}

        # FILENAMES, DIRECTORIES, INPUT/OUTPUT
        self.dataDict["inpfile"] = name  # name of the calc, it is the base name of the I/O files (.com, .chk, .log)
        self.dataDict["calcdir"] = calcdir  # name of the directory where the output files are stored
        self.dataDict["outfile"] = None   # string with the text of the Bdf log file

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

        # INFO ON Bdf EXECUTION
        self.dataDict["termination"] = None   # final termination: 0 = Normal termination, 1 = Error termination
        self.dataDict["errormsg"] = []   # in case of Error termination, error message
        self.dataDict["signs"] = [] # list of WF signs from previous step to apply in phase correction
        self.dataDict["scf_mo_maps"] = None
        self.dataDict["psioverlap"] = None # array of the WF overlaps between consecutive steps
        
        self.dataDict["eigenvectors"] = None
        self.dataDict["selfenergy"] = 0.0
        #self.dataDict["SS_root_order"] = None

        # use full file for reading properties 
        with open(os.path.join(calcdir, name)) as fout:
            self.dataDict["outfile"] = fout.read()
            # print(self.dataDict["outfile"])

        # loop over the last lines of the log file
        output = self.dataDict["outfile"].splitlines()
        if 'ERROR' in ''.join(output[-11:]).upper():
            self.dataDict["termination"] = 1
            self.dataDict["errormsg"] = ''.join(output[-11:])
            return
        if 'Total wall time' not in ''.join(output[-11:]):
            self.dataDict["termination"] = 1
            return
        else:
            self.dataDict["termination"] = 0

        au_2_ang = 5.291772104260590e-01
        au_2_kcal_1mea = 6.275094742071363e+02
        au_2_kcal_1mea_per_ang =  1.185821047928172e+03
        au_2_ev =  2.7211386033e+01
        au_2_wn =  2.1947463147e+05

        # initialize variables
        tempstate = None
        readOvMat = False        
        FLAGCHARGE = False
        state_for_GRAD = None
        states_for_NAC = None
        readtimes_for_grad = 0
        full_grad = False # include MM Atoms
        full_nac = False  # include MM Atoms
        has_GS = False
        istate = 0
        nac_istate = 0
        nac_jstate = 0
        natom = None

        i = 0
        while i < len(output):
            # collect error msg
            if 'ERROR' in output[i]:
                self.dataDict['errormsg'] += output[i]
                i += 1 

            # find osc_strength
            elif "State          X           Y           Z          Osc." in output[i]:
                k = i + 1
                self.dataDict["osc_strength"][0] = 0
                cnt = 1
                while output[k].strip() != '':
                    terms = output[k].strip().split()
                    self.dataDict["osc_strength"][cnt] = float(terms[4])
                    k += 1
                    cnt += 1
                i = k + 1
    
            elif "GSGRAD_estate=" in output[i] or "EXGRAD_estate=" in output[i]:
                if 'GSGRAD_estate' in output[i]: has_GS = True

                E = float(output[i].split('=')[1])
                self.dataDict['energy'][istate] = E

                grad_local = [[],[],[]]
                j = i + 2
                terms = output[j].strip().split()
                while terms[0] != 'Sum':
                    grad_local[0] +=  [float(terms[1])]
                    grad_local[1] +=  [float(terms[2])]
                    grad_local[2] +=  [float(terms[3])]
                    j += 1
                    terms = output[j].strip().split()
                self.dataDict["gradient"][istate] = np.array(grad_local)
                istate += 1
                i = j +1

            elif 'FNAC_cpair=' in output[i]:
                terms = output[i].strip().split()
                nac_istate, nac_jstate = int(terms[3]), int(terms[6])
                if not has_GS:
                    nac_istate -= 1
                    nac_jstate -= 1
                i += 1

            elif 'Gradient contribution from Final-NAC(S)-Escaled' in output[i]:
                grad_local = [[],[],[]]
                j = i + 1
                terms = output[j].strip().split()
                while terms[0] != 'Sum':
                    grad_local[0] +=  [float(terms[1])]
                    grad_local[1] +=  [float(terms[2])]
                    grad_local[2] +=  [float(terms[3])]
                    j += 1
                    terms = output[j].strip().split()
                try:
                    self.dataDict["nac"][nac_istate][nac_jstate] = grad_local
                except KeyError:
                    self.dataDict["nac"][nac_istate] = {}
                    self.dataDict["nac"][nac_istate][nac_jstate] = grad_local
                try:
                    self.dataDict["nac"][nac_jstate][nac_istate] = [[-x for x in y] for y in grad_local]
                except KeyError:
                    self.dataDict["nac"][nac_jstate] = {}
                    self.dataDict["nac"][nac_jstate][nac_istate] = [[-x for x in y] for y in grad_local]

                i = j + 1

            elif "Cartesian coordinates (Angstrom)" in output[i]:
                sym = []
                xyz = []
                j = i + 4
                while j < len(output):
                    ll = output[j].strip()
                    if ll == '': continue
                    if ll[0].strip()[0] == '-':
                        break
                    terms = ll.strip().split()
                    sym += [terms[1]] 
                    xyz += [ np.array(terms[3:6]).astype(np.float64) ]
                    j += 1
                
                xyz = np.array(xyz)
                natom = len(xyz)

                self.dataDict['sym'] = sym
                self.dataDict['x0'] = xyz / au_2_ang

                mass = []
                mass = np.zeros(3*len(sym))
                mdict = {'C':12.01, 'H':1.008, 'O': 16.00}
                for k in range(len(mass)):
                    mass[k] = mdict[sym[k//3]] * 1850
                self.dataDict['mass'] = mass
                
                i = j + 1

            elif ' Results of vibrations:' in output[i]:
                hess = np.zeros((N,N))
                Tmod = np.zeros((N,N))
                w = np.zeros((N))
                j = i + 1
                read_col_num = 0
                while j < len(output) and read_col_num < 3*natom:
                    if 'Irreps' in output[j]:
                        # read freq
                        w[read_col_num:read_col_num+3] = np.array(output[j+1].strip().split()[1:4]).astype(np.float64)
                        # read Tmod 
                        for k in range(natom):
                            tget = np.array(output[j+5+k].strip().split()[2:11]).astype(np.float64)
                            Tmod[3*k:3*k+3, read_col_num] = tget[0:3]
                            Tmod[3*k:3*k+3, read_col_num+1] = tget[3:6]
                            Tmod[3*k:3*k+3, read_col_num+2] = tget[6:9]
                        read_col_num += 3
                        j += 5 + natom
                    else:
                        j += 1
                wnew = copy.deepcopy(w)
                Tnew = copy.deepcopy(Tmod)
                wnew[0:6] = w[-6:]
                wnew[6:] = w[:-6]
                Tnew[:,0:6] = Tmod[:,-6:]
                Tnew[:,6:] = Tmod[:,:-6]
                self.dataDict['Tmod'] = Tnew
                self.dataDict['w'] = wnew / au_2_wn

                i = j + 1
            
            # get (mulliken) charges on QM (HighLayer) atoms (need nprint>=2)
            elif "MULLIKEN POPULATION ANALYSIS." in output[i]:
                # print('5')

                self.dataDict["charges"] = []
                k = i + 4 # skip 4 output
                while output[k].strip() != '':
                    self.dataDict["charges"].append(float(output[k].split()[2]))
                    k += 1
                i = k + 1

            # get electrostatic field at MM point charges (need mminp = 1, or revised bdf-code with mmcoup>=3)
            elif "ELECTROSTATIC POTENTIAL (VOLT) AND ELECTRIC FIELD COMPONENTS" in output[i]:
                # print('6')

                self.dataDict["elfield"] = {}
                k = i + 4 # skip 4 output
                elfield = [], [], []
                while output[k].strip() != '':
                    element = output[k].split()
                    Ex, Ey, Ez = float(element[2]), float(element[3]), float(element[4])
                    elfield[0].append(Ex), elfield[1].append(Ey), elfield[2].append(Ez)
                    k += 1
                # Bdf only provide elfield for ground state, here we may copy it to all states
                self.dataDict["elfield"][state_for_GRAD] = elfield
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

                k = i + 4 # skip 4 output
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
        """Destructor for the bdfOutput class: not only the memory allocated for the oject attributes
        (the dictionary self.dataDict) needs to be released, also the Bdf I/O files stored on disk
        can be safely removed... """

        # permanently remove directory and bdf files
        # try:
            # if not kids_log.DEBUG_RUN:
                # if self.dataDict["termination"] == 1:
                    # source=self.dataDict["calcdir"]+"/bdf-QM.log"
                    # shutil.move(source,"bdf-QM_err.log")
                # shutil.rmtree(self.dataDict["calcdir"])
        # except FileNotFoundError:
            # pass
        # clean up the log
        del self.log
        # destroy the dictionary that contains the output data
        del self.dataDict    

if __name__ == "__main__":
    from kids_log import Log
    from pprint import pformat

    Log.startSection(f'[TEST] {__file__}')
    file = 'bdf-QM.log'
    calcdir = 'qm_bdf/qmCalc00001'
    OUT = BdfOutput(file, calcdir, 'TDDFT')
    for k,v in OUT.dataDict.items():
        if k not in ['log', 'outfile']:
            Log.writeLog(f"{k}: " + pformat(v) + '\n')
