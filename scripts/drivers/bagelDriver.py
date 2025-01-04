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
import json
from pprint import pformat, pprint

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from kids_env import which
import kids_log  # manages log file output + start/end procedures
import constants  # values of physical constants and conversion factors
from QMOutput import QMOutput

class BagelInput:
    """Store the data for the execution of Bagel and process it to check its consistency"""

    def __init__(self, keys, coords, symbols, otheropt):
        """Constructor of the bagelianInput class, requires these arguments:
        - keys: string with the text of the bagelian inp file from the route section to the charge/mult line
        - bageladd: other options (e.g. the weights of the state averaged CASSCF)
        - coords: coordinates of the atoms
        - symbols: symbols of the atoms
        - otheropt: dictionary with other required options for the QM calculation
                    ("forces", "restart", "suppress_nosymm_error" ...)
        -
        """

        self.keys, self.bageladd = None, None
        self.coords, self.symbols = None, None
        # dictionary with other options
        self.otheropt = otheropt

        # store the sections of the mndoian input files already defined
        self.keys = keys
        # store the information about the molecular geometry
        self.coords = coords
        self.symbols = symbols

        self.calctype = 'BAGEL-DEBUG'

        tmp = self.keys.replace('$COORD_XYZ', '0')
        tmp = tmp.replace('$EMBED_XYZ', '0')
        self.config = json.loads(tmp)


    @staticmethod
    def checkEnv():
        """ the function checks if the environment for BAGEL execution
           for QM calculations is properly defined and if all the
           necessary executables are available """

        # get the name of the executable and the path from the environmental variable
        try:
            bagelversion = os.environ["BAGEL_EXE_QM"]
            bagelpath = os.environ["BAGEL_DIR_QM"]
        except KeyError:
            return False, "environment variables for BAGEL, $BAGEL_EXE_QM, " \
                          "$BAGEL_DIR_QM are not defined"

        # check if bagel executable is available
        if not which(bagelversion):
            message = bagelversion + " executable is not available"
            return False, message

        # return True if all the checks were OK
        return True, bagelpath


    def fileText(self, inpFileName, memory="500MB", nproc=1, version="bagel2020"):
        """Prepare the text of bagelian input file named after the inpFileName variable (input file
        should be named <inpFileName>.com, so inpFileName should not include the extension!), using the
        input data defined in the instance of bagelianInput. The memory (var memory) and
        the number of cores (var nproc) to use in the execution and the version of Bagel (gversion)
         are given as input because they are decided when running the Bagel calculation and are not
         considered parameters of the QM calculation"""

        self.embedQ = False

        # write molecular structure
        if self.config['bagel'][0]['title'] == 'molecule':
            self.config['bagel'][0]['geometry'] = []
            for atom in zip(self.symbols, *self.coords):
                self.config['bagel'][0]['geometry'] += [
                    {'atom': atom[0], "xyz": [atom[1], atom[2], atom[3]]}
                ]
            # append charges
            if 'charges' in self.otheropt and self.otheropt['charges']:
                self.embedQ = True

                coords = self.otheropt["charges"][0]
                chrg = self.otheropt["charges"][1]
                for atom in zip(chrg, *coords):
                    self.config['bagel'][0]['geometry'] += [
                        {'atom': 'Q', "xyz": [atom[1], atom[2], atom[3]], "charge": atom[0]}
                    ]
        if not os.path.exists('laststep.archive') and self.config['bagel'][1]['title'] == 'load_ref':
            del self.config['bagel'][1]
        inputText = json.dumps(self.config)
        # return the complete text of the input file
        return inputText

###################################################################################################################

class BagelOutput(QMOutput):

    def __init__(self, name, calcdir, calctype=False, verbose=False):

        # define options of the base class
        QMOutput.__init__(self)

        # Add additional bagelian-specific attributes to the data dictionary
        # FILE NAMES, DIRECTORIES, INPUT/OUTPUT
                # FILENAMES, DIRECTORIES, INPUT/OUTPUT
        self.dataDict["inpfile"] = name  # name of the calc, it is the base name of the I/O files (.com, .chk, .log)
        self.dataDict["calcdir"] = calcdir  # name of the directory where the output files are stored
        self.dataDict["outfile"] = None   # string with the text of the MNDO log file

        # ADDITIONAL PHYSICAL PROPERTIES EXTRACTED FROM THE OUTPUT FILE
        self.dataDict["charges"] = []   # charges / population analysis on QM atoms 
        self.dataDict["elfield"] = {}   # electric field at the position of the point charges (/MM atoms)
        self.dataDict["dipole"] = []    # dipole of system in XYZ and total
        self.dataDict["osc_strength"] = {}   # osc_strength to each electric state
        self.dataDict["natoms"] = None   # number of atoms in the QM calculation
        self.dataDict["nroots"] = 0     # number of roots (/electric states) in the QM calculation
        self.dataDict["optstate"] = 0   # target??
        
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
        self.calctype = calctype # --> can be "CASSCF or CASPT2"

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
        else:
            self.dataDict["termination"] = 0

        # initialize variables
        tempstate = None

        # loop over the lines of the log file
        if 'Error' in ''.join(output[-4:]):
            self.dataDict["termination"] = 1
            self.dataDict["errormsg"] = output[-11:-4]
            return

        qmdim = 0
        mmdim = 0
        i = 0
        while i < len(output):
            # determine the current state
            if "*** Geometry ***" in output[i]:
                k = i + 1
                sym = []
                ii = k
                for k in range(i+2, len(output)):
                    ll = output[k].strip()
                    if len(ll)==0 or ll[0] != '{':
                        break
                    s = ll.split(',')[0].split(':')[1].strip().split('"')[1]
                    if s == 'Q':
                        mmdim += 1
                    else:
                        qmdim += 1
                    sym += [s]
                self.dataDict["sym"] = sym
                i = k

            if '"atom" : "Q"' in output[i]:
                self.dataDict["gradcharges_extraterm"] = True        
                i += 1

            if "* NACME Target states:" in output[i]:
                terms = output[i].split()
                istate, jstate = int(terms[-3]), int(terms[-1])
                states_for_NAC = (istate, jstate)

                i += 1

            if "* Oscillator strength for transition between 0" in output[i]:
                terms = output[i].strip().split()
                to_state = int(terms[8])                
                self.dataDict['osc_strength'][0] = 0
                self.dataDict['osc_strength'][to_state] = float(terms[9])

                i += 1

            # how to read charged?
            elif "ERROR" in output[i]:
                self.dataDict['errormsg'] += [output[i]]
                i = i + 1

            elif "TERMINAT" in output[i]: # failed by termination
                self.dataDict['termination'] = 1
                break

            # terminate the reading of the main output file
            elif "COMPUTAT" in output[i]:
                self.dataDict['termination'] = 0
                self.dataDict['errormsg'] = []
                # print('8')
                break

            else:
                i += 1

        nroots = 0
        if os.path.exists(calcdir + '/ENERGY.out'):
            lines = open(calcdir + '/ENERGY.out', 'r').readlines()
            nroots = len(lines)
            for i in range(len(lines)):
                self.dataDict['energy'][i] = float(lines[i])
        self.dataDict['nroot'] = nroots

        for i in range(nroots): # try read grad
            if os.path.exists(calcdir + '/FORCE_%d.out'%i):
                has_emb = False
                lines = open(calcdir + '/FORCE_%d.out'%i, 'r').readlines()
                g = [],[],[]
                gch = [],[],[]
                for k in range(1, len(lines)):
                    terms = lines[k].strip().split()
                    if len(terms) == 0: continue
                    if k <= qmdim: # check it
                        g[0].append(float(terms[1]))
                        g[1].append(float(terms[2]))
                        g[2].append(float(terms[3]))
                    else:
                        has_emb = True
                        gch[0].append(float(terms[1]))
                        gch[1].append(float(terms[2]))
                        gch[2].append(float(terms[3]))
                self.dataDict['gradient'][i] = g
                if has_emb: self.dataDict['fullgradcharge'][i] = gch

        for i in range(nroots): # try read grad
            for k in range(i+1, nroots): # try read grad
                if os.path.exists(calcdir + '/NACME_%d_%d.out'%(i,k)):
                    has_emb = False
                    lines = open(calcdir + '/NACME_%d_%d.out'%(i,k), 'r').readlines()
                    g = [],[],[]
                    gch = [],[],[]
                    for j in range(1, len(lines)):
                        terms = lines[j].strip().split()
                        if len(terms) == 0: continue
                        if j <= qmdim: # check it
                            g[0].append(float(terms[1]))
                            g[1].append(float(terms[2]))
                            g[2].append(float(terms[3]))
                        else:
                            has_emb = True
                            gch[0].append(float(terms[1]))
                            gch[1].append(float(terms[2]))
                            gch[2].append(float(terms[3]))
                    if i not in self.dataDict['nac']:
                        self.dataDict['nac'][i] = {}
                        if has_emb: self.dataDict['fullnaccharge'][i] = {}
                    if k not in self.dataDict['nac']:
                        self.dataDict['nac'][k] = {}
                        if has_emb: self.dataDict['fullnaccharge'][k] = {}
                    self.dataDict['nac'][i][k] = g
                    self.dataDict['nac'][k][i] = [[-x for x in y] for y in g]
                    if has_emb: self.dataDict['fullnaccharge'][i][k] = gch
                    if has_emb: self.dataDict['fullnaccharge'][k][i] = [[-x for x in y] for y in gch]

        # pprint(self.dataDict)
        # pprint({'energy':self.dataDict['energy']})
        # pprint({'gradient':self.dataDict['gradient']})
        # pprint({'nac':self.dataDict['nac']})
        # pprint({'hess':self.dataDict['hess']})
        # pprint({'elfield':self.dataDict['elfield']})
        # pprint({'charges':self.dataDict['charges']})
        # pprint({'dipole':self.dataDict['dipole']})
        # pprint({'osc_strength':self.dataDict['osc_strength']})
        
           

    def __del__(self):
        """Destructor for the bagelianOutput class: not only the memory allocated for the oject attributes
        (the dictionary self.dataDict) needs to be released, also the Bagel I/O files stored on disk
        can be safely removed... """

        # permanently remove directory and bagelian files
        # try:
        #     if not kids_log.DEBUG_RUN:
        #         if self.dataDict["termination"] == 1:
        #             source=self.dataDict["calcdir"]+"/bagel-QM.log"
        #             shutil.move(source,"bagel-QM_err.log")
        #         shutil.rmtree(self.dataDict["calcdir"])
        # except FileNotFoundError:
        #     pass
        # # clean up the log
        del self.log
        # destroy the dictionary that contains the output data
        del self.dataDict



if __name__ == '__main__':
    from kids_log import Log
    from pprint import pformat

    Log.startSection(f'[TEST] {__file__}')
    file = 'bagel-QM.log'
    calcdir = 'qm_bagel/qmCalc00001'
    OUT = BagelOutput(file, calcdir, 'CASSCF')
    for k,v in OUT.dataDict.items():
        if not k in ['log', 'outfile']:
            Log.writeLog(f"{k}: " + pformat(v) + '\n')

