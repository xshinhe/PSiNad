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
# imports of local modules

import logwrt  # manages log file output + start/end procedures
import constants  # values of physical constants and conversion factors

# import of local classes
from QMOutput import QMOutput
import json 
from pprint import pformat, pprint

class BagelInput:
    """Store the data for the execution of Bagel and process it to check its consistency"""

    def __init__(self, keys, bageladd, coords, symbols, otheropt):
        """Constructor of the bagelianInput class, requires these arguments:
        - keys: string with the text of the bagelian inp file from the route section to the charge/mult line
        - bageladd: other options (e.g. the weights of the state averaged CASSCF)
        - coords: coordinates of the atoms
        - symbols: symbols of the atoms
        - otheropt: dictionary with other required options for the QM calculation
                    ("forces", "restart", "suppress_nosymm_error" ...)
        -
        """
        pprint(keys)
        self.keys, self.bageladd = None, None
        self.coords, self.symbols = None, None
        # dictionary with other options
        self.otheropt = otheropt

        # store the sections of the mndoian input files already defined
        self.keys = keys
        self.bageladd = bageladd
        # store the information about the molecular geometry
        self.coords = coords
        self.symbols = symbols

        self.calctype = 'BAGEL-DEBUG'

        print('\n'.join(self.keys))

        self.config = json.loads('\n'.join(self.keys))
        pprint(self.config)

    # =============================================================================================================

    def fileText(self, inpFileName, memory="500MB", nproc=1, gversion="bagel2020"):
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
            if self.otheropt['charges']:
                self.embedQ = True

                coords = self.otheropt["charges"][0]
                chrg = self.otheropt["charges"][1]
                for atom in zip(chrg, *coords):
                    self.config['bagel'][0]['geometry'] += [
                        {'atom': 'Q', "xyz": [atom[1], atom[2], atom[3]], "charge": atom[0]}
                    ]
        inputText = json.dumps(self.config)
        print(inputText)
        # return the complete text of the input file
        return inputText

###################################################################################################################

class BagelOutput(QMOutput):

    def __init__(self, name, calcdir, SPcalc=False, verbose=False):

        # define options of the base class
        QMOutput.__init__(self)

        # Add additional bagelian-specific attributes to the data dictionary
        # FILE NAMES, DIRECTORIES, INPUT/OUTPUT
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
        with open(os.path.join(calcdir, name + ".log")) as fout:
            self.dataDict["outfile"] = fout.read()
            # print(self.dataDict["outfile"])

        # loop over the last lines of the log file
        output = self.dataDict["outfile"].splitlines()
        if 'ERROR' in ''.join(output[-11:]).upper():
            self.dataDict["termination"] = 1
            self.dataDict["errormsg"] = ''.join(output[-11:])
            return

        # initialize variables
        tempstate = None

        # loop over the lines of the log file
        output = outfile.splitlines()
        if output[-4].split()[0].strip() == 'Error':
            self.dataDict["termination"] = 1
            self.dataDict["errormsg"] = output[-11:-4]
            return

        while i < len(output):
            # determine the current state
            if '"atom" : "Q"' in output[i]:
                self.dataDict["gradcharges_extraterm"] = True        
                i += 1

            if "* NACME Target states:" in output[i]:
                terms = output[i].split()
                istate, jstate = int(terms[-3]), int(terms[-1])
                states_for_NAC = (istate, jstate)

                i += 1

            if "* Oscillator strength for transition between 0" in output[i]:
                terms = lines[i].strip().split()
                to_state = int(terms[8])                
                self.dataDict['osc_strength'][0] = 0
                self.dataDict['osc_strength'][to_state] = [float(terms[9])]

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
        if os.path.exists(calcDir + '/ENERGY.out'):
            lines = open(calcDir + '/ENERGY.out', 'r').readlines()
            nroots = len(lines)
            for i in range(len(lines)):
                self.dataDict['energy'][i] = float(lines[i])
        self.dataDict['nroot'] = nroots

        for i in range(nroots): # try read grad
            if os.path.exists(calcDir + '/FORCE_%d.out'%i):
                lines = open(calcDir + '/FORCE_%d.out'%i, 'r').readlines()
                g = [],[],[]
                gch = [],[],[]
                for k in range(1, len(lines)):
                    terms = lines[k].strip.split()
                    if k <= geometry.NatomHDDED: # check it
                        g[0].append(terms[1])
                        g[1].append(terms[2])
                        g[2].append(terms[3])
                    else:
                        gch[0].append(terms[1])
                        gch[1].append(terms[2])
                        gch[2].append(terms[3])
                self.dataDict['gradient'] = g
                self.dataDict['fullgradcharge'] = gch

        for i in range(nroots): # try read grad
            for k in range(i+1, nroots): # try read grad
                if os.path.exists(calcDir + '/NACME_%d_%d.out'%(i,k)):
                    lines = open(calcDir + '/NACME_%d_%d.out'%(i,k), 'r').readlines()
                    g = [],[],[]
                    gch = [],[],[]
                    for k in range(1, len(lines)):
                        terms = lines[k].strip.split()
                        if k <= geometry.NatomHDDED: # check it
                            g[0].append(terms[1])
                            g[1].append(terms[2])
                            g[2].append(terms[3])
                        else:
                            gch[0].append(terms[1])
                            gch[1].append(terms[2])
                            gch[2].append(terms[3])
                    if i not in self.dataDict['nac']:
                        self.dataDict['nac'][i] = {}
                        self.dataDict['fullnaccharge'][i] = {}
                    if k not in self.dataDict['nac']:
                        self.dataDict['nac'][k] = {}
                        self.dataDict['fullnaccharge'][k] = {}
                    self.dataDict['nac'][i][k] = g
                    self.dataDict['nac'][k][i] = [[-x for x in y] for y in g]
                    self.dataDict['fullnaccharge'][i][k] = gch
                    self.dataDict['fullnaccharge'][k][i] = [[-x for x in y] for y in gch]

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
        """Destructor for the bagelianOutput class: not only the memory allocated for the oject attributes
        (the dictionary self.dataDict) needs to be released, also the Bagel I/O files stored on disk
        can be safely removed... """

        # permanently remove directory and bagelian files
        # try:
        #     if not logwrt.DEBUG_COBRAMM_RUN:
        #         if self.dataDict["termination"] == 1:
        #             source=self.dataDict["dir"]+"/bagel-QM.log"
        #             shutil.move(source,"bagel-QM_err.log")
        #         shutil.rmtree(self.dataDict["dir"])
        # except FileNotFoundError:
        #     pass
        # # clean up the log
        del self.log
        # destroy the dictionary that contains the output data
        del self.dataDict

    

    
