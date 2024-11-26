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
import shutil  # filesystem utilities
import re  # process output with regular expressions
import math  # import mathematical functions
import copy # deep copy of lists and arrays
import numpy as np

# imports of local modules

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from kids_env import which, source
from kids_log import Log  # manages log file output + start/end procedures
import constants  # values of physical constants and conversion factors
from QMOutput import QMOutput


######################################################################################################################

def _optionInRoute(keys, options):
    """This function returns True whenpassed option is
    found in the routine section of the keys text, that is when the option is contained
    in one of the lines that preceed the first empty line"""

    # initialize to false
    optionIsInRoute = False

    # loop over the lines
    for line in keys:
        # when one of the options is found in the line, set optionIsInRoute to True (search case-insensitive!)
        optionIsInRoute = optionIsInRoute or any([re.search(s, line, re.IGNORECASE) for s in options])

    # return value
    return optionIsInRoute


#####################################################################################################

class MolcasInput:
    """Store the data for the execution of Molcas and process it to check its consistency"""

    def __init__(self, keys, basisset, sewardkey, coords, symbols, otheropt):
        """Constructor of the molcasInput class, requires these arguments:
        - keys: string with the text of the molcas inp file from the route section to the charge/mult line
        - basisset: string with the text of the basis set definition (used when !basisset section is specified in the input file)
        - sewardkeys: string with the text of the additional seward keywords (used when !seward section is specified in the input file)
        - coords: coordinates of the atoms
        - symbols: symbols of the atoms
        - otheropt: dictionary with other required options for the QM calculation
                    ("forces", "restart", ...)
        - step: calculation step (needed to build differentaited input)
        """

        # initialize MolcasInput instance attributes
        # sections of the inp file given : keys (route section + charge/mult), basisset (basis set), seward (other)
        self.keys, self.basisset, self.sewardkey = None, None, None
        # coordinates and symbols of the molecule
        self.coords, self.symbols = None, None
        # dictionary with other options
        self.otheropt = {}

        # store the sections of the molcas input files already defined
        self.keys = keys
        self.basisset = basisset
        self.sewardkey = sewardkey

        #calculation type (i.e., level of theory, needed for MolcasOutput)
        self.calctype = self.readCalcType()
        Log.writeLog("QM Molcas calculation: requested calculation type is "+self.calctype+"\n")
        #number of roots requested (determined through molcas input section)
        self.nroots = self.readNroots()

        # store the information about the molecular geometry
        self.coords = coords
        self.symbols = symbols
        self.step = 0

        # process the dictionary with the other options, to check and complete it
        for key in ["forces", "restart", "charges", "field", "tdcouplings", "expert", "ricd", "cdth",
                    "nstate", "forcestate","basis", "cartesian", "SP", "command1", "mdv_type", "newBP", "ROST", "highest_NAC_state", "gradOnly", "NAC_for_vel_rescale"]:
            if key in otheropt:
                self.otheropt[key] = otheropt[key]
            else:
                self.otheropt[key] = None

    @staticmethod
    def checkEnv():
        """ the function checks if the environment for MOLCAS or OPENMOLCAS execution
            is properly defined and if all the necessary executables are available """

        # get the name of the executable and the path from the environmental variable
        try:
            molcasscript = os.environ['MOLCAS_SCRIPT']  # path of the MOLCAS launch script
            molcaspath = os.environ['MOLCAS']  # path of the MOLCAS code
        except KeyError:
            return False, "environment variables for Molcas/Openmolcas, MOLCAS_SCRIPT and MOLCAS, are not defined"

        # molcas script exists and has x permissions
        if not os.path.isfile(molcasscript) or not os.access(molcasscript, os.X_OK):
            return False, 'molcas launch script {} not found, please check the env variable MOLCAS_SCRIPT!'.format(
                molcasscript)
        # molcas directory exists
        if not os.path.isdir(molcaspath):
            return False, 'molcas directory {} not found, please check the env variable MOLCAS!'.format(molcaspath)

        # return True if all the checks were OK
        return True, molcaspath


    def readCalcType(self):
        if _optionInRoute(self.keys, ['&caspt2']):
            """F: for now SS, MS, XMS and RMS are the only supported types of CASPT2,
            therefore the routine checks for XMS, RMS and SS, otherwise the calc is assumed
            to be of MS type. Change this setting if more PT2 flavors are to be implemented in Cobramm"""
            if _optionInRoute(self.keys, ['xmul']):
                calc = "XMSPT2"
            elif _optionInRoute(self.keys, ['rmul']):
                calc = "RMSPT2"
            elif _optionInRoute(self.keys, ['nomu']):
                calc = "SSPT2"
            else:
                calc = "MSPT2"
        elif _optionInRoute(self.keys, ['&rasscf']):
            calc = "CAS"
        elif _optionInRoute(self.keys, ['&mbpt2']):
            calc = "MBPT2"
        elif _optionInRoute(self.keys, ['&scf']):
            calc = "SCF"
        else:
            Log.fatalError("Could not determine molcas calculation type: check molcas input section.\n Currently supported routines: &scf, &mbpt2, &rasscf, &caspt2 (SS, MS or XMS)\n")

        return calc

    # =============================================================================================================

    def readNroots(self):
        #returns number of electronic states requested by molcas input
        if _optionInRoute(self.keys, ['&caspt2']):
            """F: for now SS, MS and XMS are the only supported types of CASPT2,
            therefore the routine checks for XMS and SS, otherwise the calc is assumed
            to be of MS type. Change this setting if more PT2 flavors are to be implemented in Cobramm"""
            for i in range(len(self.keys)):
                if re.search('mul', self.keys[i], re.IGNORECASE) and not re.search('nomu', self.keys[i], re.IGNORECASE):
                    try:
                        Nroots = int(self.keys[i + 1].split()[0])
                    except:
                        Nroots = int(self.keys[i].split('=')[1].split()[0])

        elif _optionInRoute(self.keys, ['&rasscf']):
            for i in range(len(self.keys)):
                if re.search('ciro', self.keys[i], re.IGNORECASE):
                    try:
                        Nroots = int(self.keys[i + 1].split()[0])
                    except:
                        Nroots = int(self.keys[i].split('=')[1].split()[0])
        elif _optionInRoute(self.keys, ['&mbpt2']):
            Nroots = 1
        elif _optionInRoute(self.keys, ['&scf']):
            Nroots = 1
        else:
            Log.fatalError("Could not determine molcas calculation type: check molcas input section.\n Currently supported routines: &scf, &mbpt2, &rasscf, &caspt2 (SS, MS or XMS)\n")

        return Nroots

    # =============================================================================================================

    def restartCheck(self):
        '''Check compatibility between used-defined RASSCF keyword (LUMORB/CIRESTART)
        and identified restart file (RasOrb/JobIph). Return a fatal error if they are
        not compatible'''
        if _optionInRoute(self.keys, ['LUMO']):
            if self.otheropt["restart"] == "INPORB" or  self.otheropt["restart"] == "molcas.RasOrb":
                return True
            else:
                Log.writewarning("LUMORB keyword was used in &RASSCF routine but no orbital file was found! Molcas will use guess orbitals")
                return True
        elif _optionInRoute(self.keys, ['CIRE']):
            if self.otheropt["restart"] == "molcas.JobIph":
                return True
            else:
                Log.fatalError("CIRESTART keyword requires molcas.JobIph restart file. Please provide restart file and try again.")
        elif _optionInRoute(self.keys, ['FILE']):
            Log.fatalError("FILEorb keyword is not supported by COBRAMM. Please use LUMOrb keyword instead and provide molcas.RasOrb or INPORB file.")
        else:
            Log.writewarning("No LUMORB or CIRESTART keywords found in Molcas input! Molcas will use scratch orbitals despite orbital file...")


    # =============================================================================================================

    def getRoutineKeys(self, routine):
        '''returns the keywords of a specific molcas routine'''
        # initialize indexes for first and last keyword
        if not _optionInRoute(self.keys,[routine]):
            return []
        start, end = None, None
        for i in range(len(self.keys)):
            if self.keys[i].lower() == "&" + routine.lower():
                start = i+1
            elif "&" in self.keys[i].lower() and start != None:
                #end of routine is identified by the first appearance of '&' after the desired routine name
                end = i
                break

        return [i.lower() for i in self.keys[start:end] if (i != "" and ">" not in i)]

    # =============================================================================================================

    def getMolcasRoutines(self):
        '''returns the molcas routine names present in the cobram.command file'''
        routines = []
        for i in self.keys:
            if "&" in i:
                routines.append(i[1:].lower())

        return routines

    # =============================================================================================================

    def getUserDefinedRassiKeys(self, routine):
        '''detects and returns the keywords of any additional user-defined RASSI routine
        present in the cobram.command file following a selected routine'''
        routines = self.getMolcasRoutines()
        if routine in routines:
            routine_index = routines.index(routine)
            if routine_index != len(routines) - 1  and routines[routine_index + 1] == 'rassi':
                #The presence of the passed routine activates READ flag (to avoid reading the wrong RASSI routine)
                #the flag is switched off after reading the first following RASSI
                READ = False
                start, end = None, None
                for i in range(len(self.keys)):
                    if self.keys[i].lower() == "&" + routine.lower():
                        READ = True
                    if self.keys[i].lower() == "&rassi" and READ:
                        start = i+1
                    elif "&" in self.keys[i].lower() and start != None and READ:
                        READ = False
                        #end of routine is identified by the first appearance of '&' after the desired routine name
                        end = i
                        break
                return [i.lower() for i in self.keys[start:end] if i != ""]

        #in all other cases return False
        return False


    # =============================================================================================================

    def fileText(self, inpFileName, memory="500MB", nproc=1):
        """Prepare the text of molcas input file named after the inpFileName variable (input file
        should be named <inpFileName>.inp, so inpFileName should not include the extension!), using the
        input data defined in the instance of molcasInput. The memory (var memory) and
        the number of cores (var nproc) to use in the execution and the version of molcas (mversion)
         are given as input because they are decided when running the Molcas calculation and are not
         considered parameters of the QM calculation"""

        # initialize the text of the input file
        inputText = ""

        # add the number of cores to use and the memory of the calculation
        inputText += '>>>EXPORT MOLCAS_MEM={0}\n\n'.format(memory)
        
        #SEWARD section is automatically written by Cobramm
        inputText += "&SEWARD  &END\nTitle\nCOBRAMM\n"
        if self.otheropt["ricd"]: inputText += "RICD\nCDTH\n{}\nDoana\n".format(str(self.otheropt["cdth"]))
        # PT2 analytical gradients require RICD keyword: return fatal error if it is not active
        if self.calctype in ("SSPT2", "MSPT2", "XMSPT2", "RMSPT2") and not self.otheropt["SP"] and not self.otheropt["ricd"]:
            Log.fatalError("PT2 analytical gradients require Cholesky decomposition. Please activate command 196=1 (keyword ricd=1)\n")

        #additional seward keys (user-defined through optional !seward section)
        if len(self.sewardkey) != 0:
            for line in self.sewardkey:
                inputText += line + '\n'
            inputText += '\n'

        basisfile = self.readbassisset()
        contractions = self.extract_contraction_scheme(basisfile)
        user_routines = self.getMolcasRoutines()

        element = 0
        for atom in zip(self.symbols, *self.coords):
            inputText += "Basis set\n"
            inputText += str(contractions[element])
            inputText += "{0}".format(atom[0])+"{0}".format(element+1)+"{0:16.8f} {1:16.8f} {2:16.8f}\n".format(atom[1]/ constants.Bohr2Ang,atom[2]/ constants.Bohr2Ang,atom[3]/ constants.Bohr2Ang)
            if self.otheropt["cartesian"]: inputText +="Cartesian all\n"
            else: inputText += "Spherical all\n"
            inputText += "End of Basis\n"
            element +=1
        inputText += "\n"

        for line in self.sewardkey:
            inputText += line + "\n"

        # add external charges and points for electric field computation (if this is a QMMM calculation)
        if self.otheropt["charges"]:
            inputText += "XField\n{0}\n".format(len(self.otheropt["charges"][1]))
            for ch, x, y, z in zip(self.otheropt["charges"][1], *self.otheropt["charges"][0]):
                inputText += "{1:14.8f}{2:14.8f}{3:14.8f}{0:16.8f} 0.0 0.0 0.0\n".format(
                    ch, x / constants.Bohr2Ang, y / constants.Bohr2Ang, z / constants.Bohr2Ang)
        if self.otheropt["field"]:
            inputText += "EFLD\n{0}\n".format(len(self.otheropt["field"][1]))
            for x, y, z in zip(*self.otheropt["field"][0]):
                inputText += "{0:14.8f}{1:14.8f}{2:14.8f}\n".format(
                    x / constants.Bohr2Ang, y / constants.Bohr2Ang, z / constants.Bohr2Ang)

        #in first step, check consistency of restart keyword and restart file
        if self.step == 0 and self.calctype in ("CAS", "SSPT2", "MSPT2", "XMSPT2", "RMSPT2"):
                self.restartCheck() 
        
        # add the molcas keys using the variable self.keys
        if self.otheropt["SP"]:
            for line in self.keys:
                inputText += line + '\n'
        else:
            '''if calculation is not of SP type, check the presence of Molcas keywords to specify the state of interest
            and die if they are found (state of interest must be set through command[13] to avoid messy code)'''
            if _optionInRoute(self.getRoutineKeys('alaska'), ['root']) or _optionInRoute(self.keys, ['rlxr']):
                Log.fatalError('ROOT/RLXR keywords detected! They will be automatically added by COBRAMM. \nPlease remove them and set initial state of interest through command 13 (keyword numrlx)')
            # SCF of SCF + MBPT2
            if self.calctype in ['SCF', 'MBPT2']:
                inputText += '\n&SCF\n'
                for keyword in self.getRoutineKeys('scf'):
                    inputText += keyword + '\n'
                if self.calctype == 'MBPT2':
                    inputText += '\n&MBPT2\n'
                    inputText += 'grdt\n'
                for keyword in self.getRoutineKeys('mbpt2'):
                    if keyword != 'grdt':
                        inputText += keyword + '\n'
                inputText += '\n&ALASKA\n'
                inputText += 'show\n'
                for keyword in self.getRoutineKeys('alaska'):
                    if keyword != 'show':
                        inputText += keyword + '\n'
            # RASSCF
            elif self.calctype == 'CAS':
                ''' In case of mdv with THS scheme (for step > 0), before running RASSCF copy the JOBIPH file used for restart 
                    (i.e., WF at previous step) as QMPrev to be used for overlap (later) 
                    This should be avoided to save time if this is the newstate gradient calculation after a hop (gradOnly)'''
                if self.otheropt['command1'] == 'mdv' and int(self.otheropt['mdv_type']) == 1 and self.step > 0 and not self.otheropt['gradOnly']:
                    inputText += "\n>>> COPY JOBIPH QMPrev.JobIph\n"
                inputText += "\n&RASSCF\n"
                # for step 0: simply add the user's RASSCF keys
                if self.step == 0:
                    for keyword in self.getRoutineKeys('rasscf'):inputText += keyword + '\n'
                # for step > 0, restart from JOBIPH and add the user's keys (except fpr LUMO or CIRE)
                else:
                    inputText += "CIREstart\n"
                    for keyword in self.getRoutineKeys('rasscf'):
                        if keyword[:4] != 'lumo' and keyword[:4] != 'cire':
                            inputText += keyword + '\n'
                inputText += 'RLXROOT\n'
                inputText += str(self.otheropt["nstate"]+1)+"\n"
                #check for any additional RASSI on top of RASSCF WF
                if self.getUserDefinedRassiKeys('rasscf') and not self.otheropt['gradOnly']:
                    Log.writewarning("Additional user-defined RASSI routine identified after &RASSCF: check that this does not conflict with COBRAMM workflow!")
                    inputText += '\n&RASSI\n'
                    for keyword in self.getUserDefinedRassiKeys('rasscf'):
                        inputText += keyword + '\n'
                #gradient of reference state
                inputText += "\n&ALASKA\n"
                inputText += "show\n"
                # additional alaska keys (if present in cobram.command)
                for keyword in self.getRoutineKeys('alaska'):
                    if keyword != 'show':
                        inputText += keyword + '\n'

                if self.otheropt['command1'] == 'ci':
                    #gradient of lower state
                    inputText += "\n&ALASKA\n"
                    inputText += "ROOT="+str(self.otheropt["nstate"])+"\n"
                    inputText += "show\n"
                    # additional alaska keys (if present in cobram.command)
                    for keyword in self.getRoutineKeys('alaska'):
                        if keyword != 'show':
                            inputText += keyword + '\n'
                    # add ALASKA routine for NAC (only needed for CI search with traditional branching plane)
                    if not self.otheropt["newBP"]:
                        inputText += "\n&ALASKA\n"
                        inputText += "NAC = {0} {1}\n".format(self.otheropt["nstate"]+1, self.otheropt["nstate"])
                        inputText += "show\n"
                
                if self.otheropt['command1'] == 'mdv' and not self.otheropt['gradOnly']:
                    # dynamics with TDC (mdv_type == 1)
                    # step 0 of TDC dyn requires only RASSCF and ALASKA. step > 0 require RASSI with WF of previous step
                    if int(self.otheropt['mdv_type']) == 1 and self.step > 0:
                        inputText += "\n>>>COPY QMPrev.JobIph JOB001\n"
                        inputText += ">>>COPY JOBIPH JOB002\n"
                        inputText += "\n&RASSI\n"
                        inputText += "NROF = 2 "+str(self.nroots)+" "+str(self.nroots)+"\n"
                        for i in range(2):
                            for j in range(1,self.nroots+1):
                                inputText += str(j)+" "
                            inputText +="\n"
                        inputText += "ONEL\n"
                    # dynamics with NAC (mdv_type == 2)
                    elif int(self.otheropt['mdv_type']) == 2:
                        #add NAC calculation for all couples of states
                        if not self.otheropt['highest_NAC_state']:
                            if self.step == 0 :
                                Log.writewarning("Adding ALASKA to compute spatial NACs between all couples of states.\nThis might increase the computational cost: consider using THS hopping scheme (command 14=1, keyword nacs=tdc)\n or excluding some roots from amplitudes propagation (command 80, keyword dynroots)\n")
                            highest_root = self.nroots
                        else:
                            highest_root = self.otheropt['highest_NAC_state']
                        for i in range(1,highest_root+1):
                            for j in range(i+1,highest_root+1):
                                inputText += "\n&ALASKA\n"
                                inputText += "NAC\n"
                                inputText += str(i)+" "+str(j)+"\n"
                if self.otheropt['NAC_for_vel_rescale']:
                #add NAC calculation for old and new state after a hop states
                        Log.writeLog("Computing NAC between old and new state for velocity rescaling...\n")
                        inputText += "\n&ALASKA\n"
                        inputText += "NAC\n"
                        inputText += str(self.otheropt['NAC_for_vel_rescale'][0] +1 )+" "+str(self.otheropt['NAC_for_vel_rescale'][1] + 1)+"\n"
                else:
                    inputText += "\n"

            # MULTISTATE (NORMAL, EXTENDED OR ROTATED) CASPT2
            elif self.calctype in ['MSPT2', 'XMSPT2', 'RMSPT2']:
                inputText += "\n&RASSCF\n"
                # for step 0: simply add the user's RASSCF keys
                if self.step == 0 and not self.otheropt["ROST"]:
                    for keyword in self.getRoutineKeys('rasscf'):inputText += keyword + '\n'
                # for step > 0, restart from JOBIPH and add the user's keys (except fpr LUMO or CIRE)
                else:
                    inputText += "CIREstart\n"
                    for keyword in self.getRoutineKeys('rasscf'):
                        if keyword[:4] != 'lumo' and keyword[:4] != 'cire':
                            inputText += keyword + '\n'
                inputText += 'RLXROOT\n'
                inputText += str(self.otheropt["nstate"]+1)+"\n"
                # if this is the second run for QM/MM with MS, XMS or RMS PT2: apply state rotation (ROST) and retrieve electric field through SS-CASPT2
                if self.otheropt["ROST"]:
                    #add ROST keyword to CASSCF
                    inputText += "ROST\n"
                    #run SS-CASPT2 on the rotates states only for root of interest (to obtain electric field)
                    inputText += "\n&CASPT2\n"
                    caspt2keys = self.getRoutineKeys('caspt2')
                    for key in range(len(caspt2keys)):
                        # add IPEA and IMAG SHIFT if present in cobram.command
                        # remember that MOLCAS keywords can be specified in-line (with =) or over multiple lines
                        if 'imag' in caspt2keys[key]:
                            inputText += caspt2keys[key] + "\n"
                            if '=' not in caspt2keys[key]:
                                inputText += caspt2keys[key+1]  + "\n"
                        if 'ipea' in caspt2keys[key]:
                            inputText += caspt2keys[key]  + "\n"
                            if '=' not in caspt2keys[key]:
                                inputText += caspt2keys[key+1] + "\n"
                    # in case of CI opt we need both gradients (and both electric fields)
                    if self.otheropt['command1'] == 'ci':
                        inputText += 'MULT= 2 '+str(self.otheropt["nstate"])+" "+str(self.otheropt["nstate"]+1)+"\n"
                    # otherwise compute only root of interest
                    else:
                        inputText += 'MULT= 1 '+str(self.otheropt["nstate"]+1)+"\n"
                    inputText += "NOMULT\n"
                    inputText += "PROP\n"
                # normal (first run) PT2 calculation, following the user's input
                else:
                    #check for any additional RASSI on top of RASSCF WF
                    if self.getUserDefinedRassiKeys('rasscf') and not self.otheropt['gradOnly']:
                        Log.writewarning("Additional user-defined RASSI routine identified after &RASSCF: check that this does not conflict with COBRAMM workflow!")
                        inputText += '\n&RASSI\n'
                        for keyword in self.getUserDefinedRassiKeys('rasscf'):
                            inputText += keyword + '\n'
                    inputText += "\n&CASPT2\n"
                    for keyword in self.getRoutineKeys('caspt2'):
                        if keyword[:4] != 'prop':
                            inputText += keyword + '\n'
                    #add GRDT keyword if not provided
                    if not _optionInRoute(self.getRoutineKeys('caspt2'), ['grdt']):
                        inputText += "GRDT\n"
                    inputText += 'RLXROOT\n'
                    inputText += str(self.otheropt["nstate"]+1)+"\n"
                    inputText += "NOPROP\n"
                    #check for any additional RASSI on top of CASPT2 WF
                    if self.getUserDefinedRassiKeys('caspt2') and not self.otheropt['gradOnly']:
                        Log.writewarning("Additional user-defined RASSI routine identified after &CASPT2: check that this does not conflict with COBRAMM workflow!")
                        inputText += '\n>>>COPY molcas.JobMix JOB001\n'
                        inputText += '\n&RASSI\n'
                        for keyword in self.getUserDefinedRassiKeys('caspt2'):
                            inputText += keyword + '\n'
                    #gradient of upper state
                    inputText += "\n&ALASKA\n"
                    inputText += "show\n"
                    for keyword in self.getRoutineKeys('alaska'):
                        if keyword != 'show':
                            inputText += keyword + '\n'
                    ''' CI opt requires both gradients, but unfortunately CASPT2 gradients require a run of CASPT2 each (with GRDT and RLXR keys)
                        so we need to run CASPT2 twice (for upper and lower state).
                        In case of CI opt with traditional branching plane, NAC is also required and therefore one more CASPT2 routine with NAC key'''
                    if self.otheropt['command1'] == 'ci':
                        # second CASPT2 run (lower state)
                        inputText += "\n&CASPT2\n"
                        for keyword in self.getRoutineKeys('caspt2'):
                            if keyword[:4] != 'prop':
                                inputText += keyword + '\n'
                        # add GRDT keyword if not provided
                        if not _optionInRoute(self.getRoutineKeys('caspt2'), ['grdt']):
                            inputText += "GRDT\n"
                        inputText += 'RLXROOT\n'
                        inputText += str(self.otheropt["nstate"])+"\n"
                        inputText += "NOPROP\n"
                        # gradient of lower state
                        inputText += "\n&ALASKA\n"
                        inputText += "show\n"
                        for keyword in self.getRoutineKeys('alaska'):
                            if keyword != 'show':
                                inputText += keyword + '\n'
                        # third CASPT2 run (NAC; only needed for CI search with traditional branching plane)
                        if not self.otheropt["newBP"]:
                            inputText += "\n&CASPT2\n"
                            for keyword in self.getRoutineKeys('caspt2'):
                                if keyword[:4] != 'prop':
                                    inputText += keyword + '\n'
                            # add GRDT keyword if not provided
                            if not _optionInRoute(self.getRoutineKeys('caspt2'), ['grdt']):
                                inputText += "GRDT\n"
                            inputText += "NAC="+str(self.otheropt["nstate"])+" "+str(self.otheropt["nstate"]+1)+"\n"
                            inputText += "\n&ALASKA\n"
                            inputText += "NAC="+str(self.otheropt["nstate"])+" "+str(self.otheropt["nstate"]+1)+"\n"
                    if self.otheropt['command1'] == 'mdv' and not self.otheropt['gradOnly']:
                        # dynamics with TDC (mdv_type == 1)
                        # step 0 of TDC dyn requires only RASSCF and ALASKA. step > 0 require RASSI with WF of previous step
                        if int(self.otheropt['mdv_type']) == 1 and self.step > 0:
                            inputText += "\n>>>COPY Prev.JobMix JOB001\n"
                            inputText += ">>>COPY molcas.JobMix JOB002\n"
                            inputText += "\n&RASSI\n"
                            inputText += "NROF = 2 "+str(self.nroots)+" "+str(self.nroots)+"\n"
                            for i in range(2):
                                for j in range(1,self.nroots+1):
                                    inputText += str(j)+" "
                                inputText +="\n"
                        # dynamics with NAC (mdv_type == 2)
                        # NAC requires one more CASPT2 routine with NAC key'''
                        elif int(self.otheropt['mdv_type']) == 2:
                            if not self.otheropt['highest_NAC_state']:
                                Log.writewarning("Adding an additional CASPT2 routine for each couple of roots (spatial NACs).\nThis might increase significantly the computational cost...\n...consider using THS hopping scheme for CASPT2 dynamics (command 14=1, keyword nacs=tdc)\n or excluding some roots from amplitudes propagation (command 80, keyword dynroots)\n")
                                highest_root = self.nroots
                            else:
                                highest_root = self.otheropt['highest_NAC_state']
                            for i in range(1,highest_root+1):
                                for j in range(i+1,highest_root+1):
                                    inputText += "\n&CASPT2\n"
                                    for keyword in self.getRoutineKeys('caspt2'):
                                        if keyword[:4] != 'prop':
                                            inputText += keyword + '\n'
                                    # add GRDT keyword if not provided
                                    if not _optionInRoute(self.getRoutineKeys('caspt2'), ['grdt']):
                                        inputText += "GRDT\n"
                                    inputText += "NAC="+str(i)+" "+str(j)+"\n"
                                    inputText += "\n&ALASKA\n"
                                    inputText += "NAC="+str(i)+" "+str(j)+"\n"
                    if self.otheropt['NAC_for_vel_rescale']:
                    #add NAC calculation for old and new state after a hop states
                        Log.writeLog("Computing NAC between old and new state for velocity rescaling...\n")
                        inputText += "\n&CASPT2\n"
                        for keyword in self.getRoutineKeys('caspt2'):
                            if keyword[:4] != 'prop':
                                inputText += keyword + '\n'
                        # add GRDT keyword if not provided
                        if not _optionInRoute(self.getRoutineKeys('caspt2'), ['grdt']):
                            inputText += "GRDT\n"
                        inputText += "NAC=" + str(self.otheropt['NAC_for_vel_rescale'][0] +1 ) + " " + str(self.otheropt['NAC_for_vel_rescale'][1] + 1) + "\n"
                        inputText += "\n&ALASKA\n"
                        inputText += "NAC=" + str(self.otheropt['NAC_for_vel_rescale'][0] +1 ) + " " + str(self.otheropt['NAC_for_vel_rescale'][1] + 1) + "\n"
            elif self.calctype == 'SSPT2':
                inputText += "\n&RASSCF\n"
                # for step 0: simply add the user's RASSCF keys
                if self.step == 0 and not self.otheropt["ROST"]:
                    for keyword in self.getRoutineKeys('rasscf'):inputText += keyword + '\n'
                # for step > 0, restart from JOBIPH and add the user's keys (except fpr LUMO or CIRE)
                else:
                    inputText += "CIREstart\n"
                    for keyword in self.getRoutineKeys('rasscf'):
                        if keyword[:4] != 'lumo' and keyword[:4] != 'cire':
                            inputText += keyword + '\n'
                if self.otheropt["ROST"]: 
                    inputText += 'RLXROOT\n'
                    inputText += str(self.otheropt["nstate"]+1)+"\n"
                    inputText += "ROST\n"
                else:
                    #any additional RASSI will be performed in 1st run
                    #check for any additional RASSI on top of RASSCF WF
                    if self.getUserDefinedRassiKeys('rasscf') and not self.otheropt['gradOnly']:
                        Log.writewarning("Additional user-defined RASSI routine identified after &RASSCF: check that this does not conflict with COBRAMM workflow!\n NOTE: RASSI will be performed on RASSCF states BEFORE rotation (order might be different from SS-CASPT2)")
                        inputText += '\n&RASSI\n'
                        for keyword in self.getUserDefinedRassiKeys('rasscf'):
                            inputText += keyword + '\n'
                inputText += "\n&CASPT2\n"
                for keyword in self.getRoutineKeys('caspt2'):
                    if keyword[:4] != 'prop' and keyword[:4] != 'grdt':
                        inputText += keyword + '\n'
                if self.otheropt["ROST"]:
                    inputText += "GRDT\n"
                    inputText += 'RLXROOT\n'
                    inputText += str(self.otheropt["nstate"]+1)+"\n"
                    inputText += "PROP\n"
                    #gradient of upper state
                    inputText += "\n&ALASKA\n"
                    inputText += "show\n"
                    for keyword in self.getRoutineKeys('alaska'):
                        if keyword != 'show':
                            inputText += keyword + '\n'
                ''' CI opt requires both gradients, but unfortunately CASPT2 gradients require a run of CASPT2 each (with GRDT and RLXR keys)
                    so we need to run CASPT2 twice (for upper and lower state). --> this is done ONLY in the second run
                    In case of CI opt with traditional branching plane, NAC is also required and therefore one more CASPT2 routine with NAC key'''
                if self.otheropt['command1'] == 'ci' and self.otheropt["ROST"]:
                    # second CASPT2 run (lower state)
                    inputText += "\n&CASPT2\n"
                    for keyword in self.getRoutineKeys('caspt2'):
                        if keyword[:4] != 'prop':
                            inputText += keyword + '\n'
                    # add GRDT keyword if not provided
                    if not _optionInRoute(self.getRoutineKeys('caspt2'), ['grdt']):
                        inputText += "GRDT\n"
                    inputText += 'RLXROOT\n'
                    inputText += str(self.otheropt["nstate"])+"\n"
                    inputText += "NOPROP\n"
                    # gradient of lower state
                    inputText += "\n&ALASKA\n"
                    inputText += "show\n"
                    for keyword in self.getRoutineKeys('alaska'):
                        if keyword != 'show':
                            inputText += keyword + '\n'
                    # third CASPT2 run (NAC; only needed for CI search with traditional branching plane)
                    if not self.otheropt["newBP"]:
                        inputText += "\n&CASPT2\n"
                        for keyword in self.getRoutineKeys('caspt2'):
                            if keyword[:4] != 'prop':
                                inputText += keyword + '\n'
                        # add GRDT keyword if not provided
                        if not _optionInRoute(self.getRoutineKeys('caspt2'), ['grdt']):
                            inputText += "GRDT\n"
                        inputText += "NAC="+str(self.otheropt["nstate"])+" "+str(self.otheropt["nstate"]+1)+"\n"
                        inputText += 'NOPROP'
                        inputText += "\n&ALASKA\n"
                        inputText += "NAC="+str(self.otheropt["nstate"])+" "+str(self.otheropt["nstate"]+1)+"\n"
                if self.otheropt['command1'] == 'mdv' and self.otheropt["ROST"] and not self.otheropt['gradOnly']:
                    # dynamics with TDC (mdv_type == 1)
                    # step 0 of TDC dyn requires only RASSCF and ALASKA. step > 0 require RASSI with WF of previous step
                    if int(self.otheropt['mdv_type']) == 1 and self.step > 0:
                        inputText += "\n>>>COPY Prev.JobIph JOB001\n"
                        inputText += ">>>COPY JOBIPH JOB002\n"
                        inputText += "\n&RASSI\n"
                        inputText += "NROF = 2 "+str(self.nroots)+" "+str(self.nroots)+"\n"
                        for i in range(2):
                            for j in range(1,self.nroots+1):
                                inputText += str(j)+" "
                            inputText +="\n"
                    # dynamics with NAC (mdv_type == 2)
                    # NAC requires one more CASPT2 routine with NAC key'''
                    elif int(self.otheropt['mdv_type']) == 2:
                        if not self.otheropt['highest_NAC_state']:
                            Log.writewarning("Adding an additional CASPT2 routine for each couple of roots (spatial NACs).\nThis might increase significantly the computational cost...\n...consider using THS hopping scheme for CASPT2 dynamics (command 14=1, keyword nacs=tdc)\n or excluding some roots from amplitudes propagation (command 80, keyword dynroots)\n")
                            highest_root = self.nroots
                        else:
                            highest_root = self.otheropt['highest_NAC_state']
                        for i in range(1,highest_root+1):
                            for j in range(i+1,highest_root+1):
                                inputText += "\n&CASPT2\n"
                                for keyword in self.getRoutineKeys('caspt2'):
                                    if keyword[:4] != 'prop':
                                        inputText += keyword + '\n'
                                # add GRDT keyword if not provided
                                if not _optionInRoute(self.getRoutineKeys('caspt2'), ['grdt']):
                                    inputText += "GRDT\n"
                                inputText += "NAC="+str(i)+" "+str(j)+"\n"
                                inputText += "\n&ALASKA\n"
                                inputText += "NAC="+str(i)+" "+str(j)+"\n"
                if self.otheropt["ROST"] and self.otheropt['NAC_for_vel_rescale']:
                #add NAC calculation for old and new state after a hop states
                    Log.writeLog("Computing NAC between old and new state for velocity rescaling...\n")
                    inputText += "\n&CASPT2\n"
                    for keyword in self.getRoutineKeys('caspt2'):
                        if keyword[:4] != 'prop':
                            inputText += keyword + '\n'
                    # add GRDT keyword if not provided
                    if not _optionInRoute(self.getRoutineKeys('caspt2'), ['grdt']):
                        inputText += "GRDT\n"
                    inputText += "NAC=" + str(self.otheropt['NAC_for_vel_rescale'][0] +1 ) + " " + str(self.otheropt['NAC_for_vel_rescale'][1] + 1) + "\n"
                    inputText += "\n&ALASKA\n"
                    inputText += "NAC=" + str(self.otheropt['NAC_for_vel_rescale'][0] +1 ) + " " + str(self.otheropt['NAC_for_vel_rescale'][1] + 1) + "\n"

        return inputText

    # =============================================================================================================
    def readbassisset(self):

        natom = len(self.symbols)
        basissetkeys = self.basisset
        if len(basissetkeys) == 0:
            basisset_info = [[0 for cols in range(2)] for rows in range(natom)]  # a 2D list of length natom with
            # first column basis set label and second column the basis set contraction; basis set contraction is
            # useful for generally contracted basis sets used by Molcas, it is not used with Molpro
            for i in range(natom):
                basisset_info[i][0] = self.otheropt["basis"]
        else:
            basisset_info = [[0 for cols in range(2)] for rows in range(natom)]
            for i in range(len(basissetkeys)):
                tmp_elements = []
                try:
                    tmp_elements = basissetkeys[i].split()[0].strip().split(',')
                    for k in range(len(tmp_elements)):
                        tmp2_elements = []
                        tmp2_elements = tmp_elements[k].split('-')
                        if len(tmp2_elements) == 1:
                            basisset_info[int(tmp2_elements[0]) - 1][0] = basissetkeys[i].split()[1].strip()
                            if basissetkeys[i].split()[1].strip() not in ['6-31G', '6-31G*', '6-31Gp',
                                                                                  '6-31G**','6-31Gpp', 'STO-3G', 
                                                                                  'CC-PVDZ', 'CC-PVTZ', 
                                                                                  'DEF2-SVP', 'DEF2-TZVP', 'DEF2-TZVPP', 
                                                                                  'PC-0', 'PC-1', 'PC-2', 'PC-3', 'PC-4']:
                                basisset_info[int(tmp2_elements[0]) - 1][1] = basissetkeys[i].split()[2].strip()
                            else:
                                basisset_info[int(tmp2_elements[0]) - 1][1] = 'dummy'
                        else:
                            for l in range(int(tmp2_elements[0]) - 1, int(tmp2_elements[1])):
                                basisset_info[l][0] = basissetkeys[i].split()[1].strip()
                                if basissetkeys[i].split()[1].strip() not in ['6-31G', '6-31G*', '6-31Gp',
                                                                                  '6-31G**','6-31Gpp', 'STO-3G', 
                                                                                  'CC-PVDZ', 'CC-PVTZ', 
                                                                                  'DEF2-SVP', 'DEF2-TZVP', 'DEF2-TZVPP', 
                                                                                  'PC-0', 'PC-1', 'PC-2', 'PC-3', 'PC-4']:
                                    basisset_info[l][1] = basissetkeys[i].split()[2].strip()
                                else:
                                    basisset_info[l][1] = 'dummy'
                except:
                    pass
            for i in range(natom):
                if not basisset_info[i][0] or not basisset_info[i][1]:
                    Log.fatalError("Basis set definition for atom " + str(i+1) + " is missing")

        return basisset_info

    #############################################################################################################

    def extract_contraction_scheme(self, basis_set_info):

            contractions = []
            for atom in range(len(self.symbols)):
            # in the filenames for the Pople's polarized basis set, *'s are substituted with p's
                if basis_set_info[atom][0] == '6-31G**':
                    basisName = '6-31Gpp'
                elif basis_set_info[atom][0] == '6-31G*':
                    basisName = '6-31Gp'
                else:
                    basisName = basis_set_info[atom][0]
                try:
                    molcaspath = os.getenv('MOLCAS')
                    with open(molcaspath + '/basis_library/' + basisName) as basisfile:
                        output = basisfile.readlines()
                except IOError:
                    Log.writeLog(molcaspath + '/basis_library/' + basisName+"\n")
                    Log.fatalError("The given basis set {0} has not been found (path: "+ molcaspath + '/basis_library/' + basisName + "), please change basis definition".format(
                        basis_set_info[atom][0]))

                        # now extract contraction scheme from the basis set file
                for j in range(len(output)):
                    # in the files for the Pople's polarized basis set, names with *'s should be used
                    if basis_set_info[atom][0] == '6-31Gpp':
                        basisName = '6-31G**'
                    elif basis_set_info[atom][0] == '6-31Gp':
                        basisName = '6-31G*'
                    else:
                        basisName = basis_set_info[atom][0]
                    lookfor = output[j].lower().find('/' + self.symbols[atom].lower() + '.' + str(basisName.lower()))
                    # replace default contraction with the user-defined unless the basis set has only one contraction
                    if ((basis_set_info[atom][0] not in ['6-31G', '6-31G*', '6-31Gp',
                                                            '6-31G**','6-31Gpp', 'STO-3G', 
                                                            'CC-PVDZ', 'CC-PVTZ', 
                                                            'DEF2-SVP', 'DEF2-TZVP', 'DEF2-TZVPP', 
                                                            'PC-0', 'PC-1', 'PC-2', 'PC-3', 'PC-4'])
                            and (lookfor != -1)):
                        tmp_label = output[j][1:].split('.')
                        if basis_set_info[atom][1] != 0:
                            tmp_label[4] = str(basis_set_info[atom][1])  # set user-defined basis set contraction
                        contraction_string = '.'.join(tmp_label) + '\n'
                    elif lookfor != -1:
                        contraction_string = output[j][1:] + '\n'

                contractions.append(contraction_string)

            return contractions


###################################################################################################################

class MolcasOutput(QMOutput): #MolcasOutput is subclass of QMOutput

    def __init__(self, name, calcdir, calctype, SPcalc=False): #--> F: addded calctype, to be set according to input! If not possible, add routine to determine calctype

        # define options of the base class
        QMOutput.__init__(self)

        # Add additional gaussian-specific attributes to the data dictionary -->F: check if some can be deleted and/or Molcas needs additional ones
        # FILE NAMES, DIRECTORIES, INPUT/OUTPUT
        self.dataDict["inpfile"] = name  # name of the calc, it is the base name of the I/O files (.com, .chk, .log)
        self.dataDict["calcdir"] = calcdir  # name of the directory where the output files are stored
        self.dataDict["outfile"] = None   # string with the text of the Molcas log file, initialized to None --> F: what's the difference with the log attribute of QMOutput?
        # ADDITIONAL PHYSICAL PROPERTIES EXTRACTED FROM THE OUTPUT FILE
        self.dataDict["elfield"] = {}   # electric field at the position of the point charges -->F: check for MOLCAS
        self.dataDict["natoms"] = None   # number of atoms in the QM calculation
        self.dataDict["nroots"] = 0   # number of roots in the QM calculation
        # INFO ON MOLCAS EXECUTION
        self.dataDict["termination"] = None   # final termination: 0 = Normal termination, 1 = Error termination
        self.dataDict["errormsg"] = None   # in case of Error termination, error message
        self.dataDict["signs"] = [] # list of WF signs from previous step to apply in phase correction
        self.dataDict["psioverlap"] = None # array of the WF overlaps between consecutive steps
        self.dataDict["gradient"] = {} # dictionary of state gradients
        self.dataDict["nac"] = {} # dictionary of NACs
        self.dataDict["eigenvectors"] = None
        #self.dataDict["SS_root_order"] = None
        self.calctype = calctype # --> F: calculation type can be "SCF", "MBPT2", "CAS", "SSPT2", "MSPT2", "XMSPT2" or ""RMSPT2"

        # store the Gaussian log removing AO and MO definitons, CI expansions, etc. -->F: is filter needed for Molcas?
        if SPcalc == False:
            self.dataDict["outfile"] = self._filterOutput(os.path.join(calcdir, name + ".log")) # -->F: for now I leave the filter as it will have no effect on Molcas out
        else:
            with open(os.path.join(calcdir, name + ".log")) as fout:
                self.dataDict["outfile"] = fout.read()

        # use full file for reading properties 
        with open(os.path.join(calcdir, name + ".log")) as fout:
            outfile = fout.read()

        # initialize variables
        tempstate = None
        Efieldstate = None

        # loop over the lines of the log file
        output = outfile.splitlines()
        #if output[-4].split()[0].strip() == 'Error': -->F: old "if"
        if "Happy landing!" not in output[-4]: #-->F: "Happy landing!" is printed in the 4th line from bottom whenever calculation is successfull
            self.dataDict["termination"] = 1
            self.dataDict["errormsg"] = output[-30:-4]
            return
        else:
            self.dataDict["termination"] = 10

        ''' F: created flag variable for reading electric field: it is switched on in case of SCF or MBPT2 calculatrion (only one state for electric fied)
            in case of CAS or PT2 calculation, it turns on only when electric field is printed in the routine of intrest (MCLR reprints electric potential 
            and field but we do not want to it to overwrite the CAS one)'''
        readEfield = True
        # F: created flag variable for reading overlap matrix: it is usually switched off, but it is turned on when "RASSI" section is encountered
        readOvMat = False
        
        #initialize flag for reading QM charges 
        FLAGCHARGE = False

        #initialize state for gradient (GS by default, edited later automatically for CASSCF and CASPT2)
        state_for_grad = 0

        for i in range(len(output)):

            # extract SCF energy from log file, when available
            if self.calctype == "SCF" and "Total SCF energy" in output[i]:
                self.dataDict["energy"][0] = float(output[i].split()[-1].strip())
                self.dataDict["optstate"] = 0
                self.log += "SCF energy is {0:12.8f} Hartree\n".format(self.dataDict["energy"][0])
            
            # extract MP2 energy from log file, when available
            if self.calctype == "MBPT2" and "Total MBPT2 energy" in output[i]:
                self.dataDict["energy"][0] = float(output[i].split()[-1].strip())
                self.dataDict["optstate"] = 0
                self.log += "MP2 energy is {0:12.8f} Hartree\n".format(self.dataDict["energy"][0])

            # extract CAS energies from log file, when available
            if self.calctype == "CAS" and "RASSCF root number" in output[i]:
                rootnr = int(output[i].split()[4])
                energy = float(output[i].split()[-1].strip())
                self.dataDict["energy"][rootnr - 1] = energy
                self.log += 'CASSCF energy for state {0:3d} is {1:12.8f} Hartree\n'.format(rootnr, energy)

            # extract SSPT2 energies from log file, when available
            """ in case of SS-PT2 states need to be reordered and mapped to CASSCF states. 
            This is done later in this code (after finishing reading the output line by line).
            The mapping information is stored for subsequent overlap matrix reordering (if present)"""
            if self.calctype == "SSPT2" and "::    CASPT2 Root" in output[i]:
                rootnr = int(output[i].split()[3])
                energy = float(output[i].split()[-1].strip())
                self.dataDict["energy"][rootnr - 1] = energy

            # extract MSPT2 energies from log file, when available
            if self.calctype == "MSPT2" and "::    MS-CASPT2 Root" in output[i]:
                rootnr = int(output[i].split()[3])
                energy = float(output[i].split()[-1].strip())
                #avoid rewriting and overcounting of nroots when there are multiple CASPT2 routines
                if rootnr -1 not in self.dataDict["energy"]:
                    self.dataDict["energy"][rootnr - 1] = energy
                    self.log += 'MS-PT2 energy for state {0:3d} is {1:12.8f} Hartree\n'.format(rootnr, energy)
                    self.dataDict["nroots"] += 1

            # extract XMSPT2 energies from log file, when available
            if self.calctype == "XMSPT2" and "::    XMS-CASPT2 Root" in output[i]:
                rootnr = int(output[i].split()[3])
                energy = float(output[i].split()[-1].strip())
                #avoid rewriting and overcounting of nroots when there are multiple CASPT2 routines
                if rootnr -1 not in self.dataDict["energy"]:
                    self.dataDict["energy"][rootnr - 1] = energy
                    self.log += 'XMS-PT2 energy for state {0:3d} is {1:12.8f} Hartree\n'.format(rootnr, energy)
                    self.dataDict["nroots"] += 1

            # extract RMSPT2 energies from log file, when available
            if self.calctype == "RMSPT2" and "::    RMS-CASPT2 Root" in output[i]:
                rootnr = int(output[i].split()[3])
                energy = float(output[i].split()[-1].strip())
                #avoid rewriting and overcounting of nroots when there are multiple CASPT2 routines
                if rootnr -1 not in self.dataDict["energy"]:
                    self.dataDict["energy"][rootnr - 1] = energy
                    self.log += 'RMS-PT2 energy for state {0:3d} is {1:12.8f} Hartree\n'.format(rootnr, energy)
                    self.dataDict["nroots"] += 1

            # extract eigenvectors from log file (to be used for ROST in case of QM/MM X/R/MSPT2 calculation)
            if self.calctype in ['MSPT2'] and ("Eigenvectors:" in output[i]) and (self.dataDict["eigenvectors"] is None):
                nroots = self.dataDict["nroots"]
                nblocks = nroots // 5
                if nroots % 5 != 0:
                    nblocks += 1
                self.dataDict["eigenvectors"] = self.readEigenvectors(output[i: i+nroots*nblocks+nblocks+2], nroots, "Eigenvectors")

            #in case of XMS or RMS we extract the eigenvectors in ters of the original CASSCF states
            elif self.calctype in ['XMSPT2', 'RMSPT2'] and ("In terms of the input states:" in output[i]) and (self.dataDict["eigenvectors"] is None):
                nroots = self.dataDict["nroots"]
                nblocks = nroots // 5
                if nroots % 5 != 0:
                    nblocks += 1
                self.dataDict["eigenvectors"] = self.readEigenvectors(output[i: i+nroots*nblocks+nblocks+2], nroots, "In terms of the input states:")


            # determine the state of interest for optimization --> F: check that this is ok for cobramm...
            if "Root passed to geometry opt" in output[i]:
                 self.dataDict["optstate"] = int(output[i].split()[-1]) - 1 

            # set root for analytical gradient (from MCLR output in case of RASSCF or CASPT2 gradient)
            if "Lagrangian multipliers are calculated for state no" in output[i]:
                state_for_grad = int(output[i].split()[-1]) - 1

            # get analytical gradient (from ALASKA output)
            if "Molecular gradients" in output[i]:
                # initialize the dictionary entry for the gradient
                gradient = [[], [], []]
                # loop over following lines
                j = 0
                while True:
                    # try to read the block of numbers with the gradient, if an exception is raised finish reading
                    try:
                        element = output[i + j + 8].split()
                        gradient[0].append(float(element[1]))
                        gradient[1].append(float(element[2]))
                        gradient[2].append(float(element[3]))
                    except (IndexError, ValueError):
                        break
                    # move to following line
                    j += 1
                    # when the gradient is too large, abort calculation
                    if abs(gradient[0][-1]) > 1.0 or abs(gradient[1][-1]) > 1.0 or abs(gradient[2][-1]) > 1.0:
                        Log.fatalError('Large gradient during analytical computation. Check last QM '
                                          'output {0}'.format(os.path.join(calcdir, name + ".log")))
                # now store the gradient extracted from the molcas log
                #self.dataDict["gradient"] = {self.dataDict["optstate"]: gradient}
                self.dataDict["gradient"][state_for_grad] = gradient

            # set roots for analytical NAC (from MCLR output in case of RASSCF or CASPT2 gradient)
            if "Lagrangian multipliers are calculated for states no" in output[i]:
                states_for_NAC = (int(output[i].split()[-2].strip("/")) - 1, int(output[i].split()[-1]) - 1)

            # get analytical NACs (from ALASKA output)
            if "Total derivative coupling" in output[i]:
                # initialize the dictionary entry for the NAC
                NAC = [[], [], []]
                # loop over following lines
                j = 0
                while True:
                    # try to read the block of numbers with the gradient, if an exception is raised finish reading
                    try:
                        element = output[i + j + 8].split()
                        NAC[0].append(float(element[1]))
                        NAC[1].append(float(element[2]))
                        NAC[2].append(float(element[3]))
                    except (IndexError, ValueError):
                        break
                    # move to following line
                    j += 1
                # now store the NAC extracted from the molcas log
                try:
                    self.dataDict["nac"][states_for_NAC[0]][states_for_NAC[1]] = NAC
                except KeyError:
                    self.dataDict["nac"][states_for_NAC[0]] = {}
                    self.dataDict["nac"][states_for_NAC[0]][states_for_NAC[1]] = NAC
                try:
                    self.dataDict["nac"][states_for_NAC[1]][states_for_NAC[0]] = [[-x for x in y] for y in NAC]
                except KeyError:
                    self.dataDict["nac"][states_for_NAC[1]] = {}
                    self.dataDict["nac"][states_for_NAC[1]][states_for_NAC[0]] = [[-x for x in y] for y in NAC]


            # extract High layer charges --> F: for the future: change Mulliken to ESP
            if self.calctype == 'SCF' or self.calctype == 'MBPT2':
                if "Mulliken charges per centre and basis function type" in output[i]:
                    FLAGCHARGE = True
                    self.dataDict["charges"] = []
            else:
                if self.calctype == 'CAS' and "Mulliken population analysis for root number" in output[i] and int(output[i].split()[-1]) == self.dataDict["optstate"]+1:
                    FLAGCHARGE = True
                    self.dataDict["charges"] = []
                elif self.calctype in ['SSPT2', 'MSPT2', 'XMSPT2', 'RMSPT2'] and "Compute H0 matrices for state" in output[i] and int(output[i].split()[-1]) == self.dataDict["optstate"]+1:
                    FLAGCHARGE = True
                    self.dataDict["charges"] = []
            if "Molecular properties" in output[i]:
                FLAGCHARGE = False
            if "N-E" in output[i] and FLAGCHARGE:
                for k in range(1,len(output[i].split())):
                    self.dataDict["charges"].append(float(output[i].split()[k]))

            #Turn on electric field reading for state of interest: for CAS, there is explicit printing of state number before the properties
            # for SSPT2 one needs to refer to the PT2 correction string "Compute H0 matrices for state"
            if self.calctype == 'CAS' and "Expectation values of various properties for root number" in output[i] and readEfield:
                Efieldstate = int(output[i].split()[-1]) -1
            elif (self.calctype in ['SSPT2', 'MSPT2', 'XMSPT2', 'RMSPT2']) and ("Compute H0 matrices for state" in output[i]) and readEfield:
                Efieldstate = int(output[i].split()[-1]) -1
            elif (self.calctype in ['SCF', 'MBPT2']) and ("Molecular properties" in output[i]):
                Efieldstate = 0

            # get electrostatic field at MM point charges 
            if "Electric field:" in output[i] and readEfield and (Efieldstate is not None):
                # loop over following lines
                j = 0
                elfield = [], [], []
                while True:
                    # try to read the block of numbers with the electric field, if an exception is raised finish reading
                    try:
                        element = output[i + j + 2].split()
                        Ex, Ey, Ez = float(element[1]), float(element[2]), float(element[3])
                        elfield[0].append(Ex), elfield[1].append(Ey), elfield[2].append(Ez)
                    except (IndexError, ValueError):
                        break
                    # move to following line
                    j += 1
                if elfield[0]:
                    #self.dataDict["elfield"] = {self.dataDict["optstate"]: elfield}
                    self.dataDict["elfield"][Efieldstate] = elfield
                else:
                    self.dataDict["elfield"] = None
            if "Start Module: mclr" in output[i]:
                #turn off electric field reading in MCLR output
                readEfield = False

            # get dipole moments --> F: same as elfield...which one to read?
            if "Dipole Moment (Debye)" in output[i] and readEfield and Efieldstate == self.dataDict["optstate"]:
                dip = output[i + 2].split()
                self.dataDict["dipole"] = float(dip[1]) * constants.Debye2AU, float(dip[3]) * constants.Debye2AU, \
                    float(dip[5]) * constants.Debye2AU, float(dip[7]) * constants.Debye2AU
                #self.dataDict["originaldipole"] = float(dip[1]), float(dip[3]), \
                #    float(dip[5]), float(dip[7])

            #get overlap matrix if RASSI routine is present
            if "RASSI" in output[i]:
                readOvMat = True

            if readOvMat and "Nr of states:" in output[i]:
                n_job1 = 0
                n_job2 = 0
                read_nstates = True
                j = 0
                while read_nstates:
                    j += 1
                    if "JobIph:" in output[i+j]:
                        n_job1 += output[i+j].split().count('1')
                        n_job2 += output[i+j].split().count('2')
                    if "OVERLAP MATRIX FOR THE ORIGINAL STATES" in output[i+j]:
                        read_nstates = False
                if n_job2 != 0:
                    #create a shorter output to pass to readOverlap() function (RASSI is usually at the end of the calculation)
                    tmp_out = output[i:]
                    nstates = int(int(output[i].split()[-1])/2)
                    # F: o is the "raw" overlap matrix read from molcas output, which needs to be reordered in case of SS-PT2 calculation
                    self.dataDict["psioverlap"] = self.readOverlap(tmp_out, n_job1, n_job2)

            # terminate the reading of the main output file
            if "Happy landing!" in output[i]:
                nmain = i
                break

        #set attribute number of states according to length of energy list
        self.dataDict["nroots"] = len(self.dataDict["energy"])
        
        # reordering the SSPT2 energies and corresponding overlap matrix (if present) and gradients
        # correct order of RASSCF states is saved in dataDict for subsequent gradient run
        if self.calctype == 'SSPT2':
            eigenvectors = np.zeros((self.dataDict["nroots"], self.dataDict["nroots"]))
            unsorted_energies = list(self.dataDict["energy"].values())
            sorted_energies= sorted(unsorted_energies)
            for i in range(len(sorted_energies)):
                for j in range(len(unsorted_energies)):
                    if sorted_energies[i] == unsorted_energies[j]:
                        eigenvectors[i][j] = 1
            self.dataDict["eigenvectors"] = np.array(eigenvectors)
            #print energies only for second run (identified by the presence of the gradient)
            if self.dataDict["gradient"] != {}:
                for i in self.dataDict["energy"]:
                    self.log += 'Sorted SS-PT2 energy for state {0:3d} is {1:12.8f} Hartree\n'.format(i+1, self.dataDict["energy"][i])


        # F: I do not know exactly what is the following "for" cycle needed for...maybe to read output of SP calculations on top of dynamics step (command 200)?...commented for now...
        #addEne = {}
        #for i in range(nmain, len(output)):

        #    # extract SCF energy from log file, when available
        #    try:
        #        if output[i].split()[0].strip() == 'SCF' and output[i].split()[1].strip() == 'Done:':
        #            scf = float(output[i].split()[4].strip())
        #            self.log += "SCF energy for secondary output is {0:12.8f} Hartree\n".format(scf)
        #            addEne[0] = scf
        #    except (IndexError, ValueError):
        #        pass

        #    # extract MP2 energy from log file, when available
        #    try:
        #        for mp2str1 in output[i].split('\\'):
        #            if mp2str1.find('EUMP2') != -1:
        #                mp2 = float(mp2str1.split()[5])
        #                self.log += "MP2 energy for secondary output is {0:12.8f} Hartree\n".format(mp2)
        #                addEne[0] = mp2
        #    except (IndexError, ValueError):
        #        pass

        #nstart = 100
        #if addEne:
        #    for i, e in addEne.items():
        #        if i == 0:
        #            pass
        #        else:
        #            if i in self.dataDict["energy"]:
        #                self.dataDict["energy"][nstart+i] = e
        #            else:
        #                self.dataDict["energy"][i] = e

        self.log += "\n"

        # finish collecting data from .log file
        # now collect some info on the job execution

        # check if Molcas termination is normal or error, and store excerpt of log file with error message -->F: this stuff has already beed done before! Commented...
        #if output[-1].split()[0].strip() == 'Normal':
        #    self.dataDict["termination"] = 0
        #elif output[-4].split()[0].strip() == 'Error':
        #    self.dataDict["termination"] = 1
        #    self.dataDict["errormsg"] = output[-11:-4]

    # =============================================================================================================

    def __del__(self):
        """Destructor for the molcasOutput class: not only the memory allocated for the oject attributes
        (the dictionary self.dataDict) needs to be released, also the Gaussian I/O files stored on disk
        can be safely removed... """

        # permanently remove directory and gaussian files
        try:
            if not kids_log.DEBUG_RUN:
                if self.dataDict["termination"] == 1:
                    source=self.dataDict["calcdir"]+"/molcas.log"
                    shutil.move(source,"molcas_err.log")
                    source=self.dataDict["calcdir"]+"/molcas.RasOrb"
                    shutil.move(source,"molcas_err.RasOrb")
                    source=self.dataDict["calcdir"]+"/molcas.rasscf.molden"
                    shutil.move(source,"molcas_err.rasscf.molden")
                shutil.rmtree(self.dataDict["calcdir"])
        except FileNotFoundError:
            pass
        # clean up the log
        del self.log
        # destroy the dictionary that contains the output data
        del self.dataDict

    # =============================================================================================================

    def _filterOutput(self,log):
        """Remove AO and MO definitions, CI expansions from output """
        cleanlog = []

        with open(log) as fout:
            log = fout.read()
        fulllog = log.splitlines()

        toPrint = 1
        for iLine in range(len(fulllog)):
            if fulllog[iLine].find("++    Molecular orbitals") != -1:
                toPrint = 0
            elif fulllog[iLine] == "--":
                toPrint = 1
            elif fulllog[iLine].find("++    Molecular charges") != -1:
                toPrint = 0
            elif fulllog[iLine].find("++    Molecular properties") != -1:
                toPrint = 0
            elif fulllog[iLine].find("Von Neumann Entropy") != -1:
                toPrint = 0
            elif fulllog[iLine].find("Mulliken population analysis for root number") != -1:
                toPrint = 0
            elif fulllog[iLine].find("Expectation values of various properties for root number") != -1:
                toPrint = 0


            if toPrint == 1: 
                cleanlog.append(fulllog[iLine])

        return "\n".join(cleanlog)

    # =============================================================================================================

    def readOverlap(self, log, dim1, dim2):
        o_tmp=np.zeros((dim2,dim1))
        for j in range(len(log)):
            if log[j].find('OVERLAP MATRIX FOR THE ORIGINAL STATES') != -1 and log[j + 2].find('Diagonal, with elements') == -1:
                j+=1
                #routine for reading the lower triangular overlap matrix with arbitrary size
                #<1|1>              # upper triangular matrix corresponds to overlaps at the reference <Psi_i(r)|Psi_j(r)>
                #<2|1> <2|2>        # in the case of SS-PT2 the WFs are not orthogonal, os this should introduce a small correction to SS-PT2 NACs
                #<3|1> <3|2> <3|3>  # values are stored in o_ref
                #<1'|1> <1'|2> <1'|3> <1'|1'>          # lower square matrix (nroots x nroots) corresponds to overlap <Psi_i(r+dr)|Psi_j(r)>
                #<2'|1> <2'|2> <2'|3> <2'|1'> <2'|1'>  # values are stored in o
                #<3'|1> <3'|2> <3'|3> <3'|1'> <3'|2'>  #
                #<3'|3'>
                #note, that only 5 elements are printed per line
                #m runs over 2 x nroots (.eq raws in the overlap matrix)
                for m in range(dim1+dim2):
                    #reading next line of the output
                    j+=1
                    #the first nroot raws correspond to overlaps at the reference <Psi_i(r)|Psi_j(r)>
                    if m < dim1:
                        #la is an auxiliary iterator which runs from 0 to 4 (i.e. over the elements in one line)
                        la = 0
                        #l runs over the 0:m+1, i.e. reading only a triangluar matrix (including diagonal elements, therefore m+1)
                        for l in range(m+1):
                            #if the 5th elements is reached, CR to 0 and continue reading next line from the first element
                            if la >= 5:
                                j += 1
                                la = 0 
                            la += 1
                    ##the second nroot raws correspond to overlaps at the reference <Psi_i(r+dr)|Psi_j(r)>
                    elif m >= dim1:
                        #as m runs over 2 x nroots, k maps back to a particular root
                        k = m-dim1
                        la = 0
                        #l runs over the 0:nroots+k+1, i.e. reading only a triangluar matrix (including diagonal elements, therefore m+1)
                        for l in range(m+1):
                            if la >= 5:
                                j += 1
                                la = 0
                            if l < dim1:
                                o_tmp[k][l] = float(log[j].split()[la])
                            la += 1
        #to avoid problems with step 0 (in case an additional rassi routine is present),
        #return "o" only if we have effectively read an overlap matrix, else return None
        if np.any(o_tmp):
            return o_tmp
        else:
            return None

    # =============================================================================================================

    def readEigenvectors(self,log, nroots, match):
        """ read eigenvector matrix from &CASPT2 output (RMS, XMS, MS)"""
        #flag to activate matrix reading
        READ = False 
        # if there are more that 5 stsates, we have multiple blocks
        block = 0
        #j keeps track of block line
        j=0
        #initialize empty matrix
        eigenvectors = np.zeros((nroots,nroots))

        for line in log:
            if match in line:
                READ = True
                continue
            elif "--" in line and READ:
                READ = False
                break
            if READ:
                values = line.strip().split()
                if values:
                    for i in range(len(values)):
                        eigenvectors[i + 5*block][j] = float(values[i])
                    j += 1
    
                if j == nroots:
                    #we are at the end of a block: increase block value and reset j
                    j = 0
                    block += 1
    
        return eigenvectors
    # =============================================================================================================

    def orbfiles(self):
        """Check their existance and then return a list with the names of the files to be
         saved to store the orbitals to disk"""

        files_to_save = []
        files_string = ""

        #in order to save disk space, only files relevant for each type of calculation will be saved
        if self.calctype == 'SCF' or self.calctype == 'MBPT2':
            extensions = (".ScfOrb", ".scf.molden")
        elif self.calctype == 'CAS' or self.calctype == 'SSPT2':
            extensions = (".RasOrb", ".JobIph", ".rasscf.molden", ".rasscf.h5", ".rassi.h5")
        elif self.calctype == 'MSPT2' or self.calctype == 'XMSPT2' or self.calctype == 'RMSPT2':
            extensions = (".RasOrb", ".JobMix", ".rasscf.molden", ".rasscf.h5", ".rassi.h5")

        # name of the orbital/wavefunction files that needs to be saved
        for ext in extensions:
            filename = os.path.join(self.dataDict["calcdir"], self.dataDict["inpfile"]+ext)

            # if the file exists, add its path to the list
            if os.path.exists(filename):
                files_to_save.append(filename)
                files_string += self.dataDict["inpfile"]+ext+ " "
        if len(files_to_save) != 0:
            Log.writeLog("Molcas orbital/wavefunction files to be saved: "+files_string+"\n")
        else:
            Log.writewarning("Molcas orbital/wavefunction not found! It will not be saved\n")

        return files_to_save

    # =============================================================================================================

    def restartfile(self):
        """Check their existance and then return a the names of the files to be
         saved to store the orbitals to disk for restart purpose"""

        # name of the chk file that needs to be saved for restart
        restname1 = os.path.join(self.dataDict["calcdir"], "JOBIPH")
        restname2 = os.path.join(self.dataDict["calcdir"], self.dataDict["inpfile"]+".JobIph")
        restname3 = os.path.join(self.dataDict["calcdir"], self.dataDict["inpfile"]+".RasOrb")

        # return restart file name with priority to JobIph (RasOrb is used otherwise)
        if os.path.exists(restname1):
            shutil.move(restname1, os.path.join(self.dataDict["calcdir"], self.dataDict["inpfile"]+".JobIph"))
            return restname2
        elif os.path.exists(restname2):
                return restname2
        elif os.path.exists(restname3):
                return restname3
        else:
            # if no restart file is present, return None
            return None

    # =============================================================================================================

    def SS_CASPT2_reorder(self):
        """Reorder energies, gradients, electric fields, psioverlaps and nacs according
        to SS-CASPT2 order"""
        # NOW THIS FUNCTION IS NOT USED...however, it took me some time to write it, so I leave it here in case it is needed in the future
        cas_order = []
        unsorted_energies = list(self.dataDict["energy"].values())
        sorted_energies= sorted(unsorted_energies)
        #if the states are already in the correct order, save time and don't do nothing
        if sorted_energies == unsorted_energies:
            self.dataDict["SS_root_order"] = list(range(self.dataDict["nroots"]))
            return
        #in all other cases, re-order all necessary quantities
        #ENERGIES
        for i in range(len(sorted_energies)):
            self.dataDict["energy"][i] = sorted_energies[i]
            #when reordering the energies, build the cas_order vector
            for j in range(len(unsorted_energies)):
                if sorted_energies[i] == unsorted_energies[j]:
                    cas_order.append(j)
        self.dataDict["SS_root_order"] = cas_order
        #PSI OVERLAPS and NACS
        for a in ("psioverlap", "nac"):
            if self.dataDict[a]:
                o_tmp = copy.deepcopy(self.dataDict[a])
                for i in range(len(o_tmp)):
                    for j in range(len(o_tmp)):
                        self.dataDict[a][i][j] = o_tmp[cas_order[i]][cas_order[j]]
        #GRADIENTS and ELECTRIC FIELDS
        for a in ("gradient", "elfield"):
            if self.dataDict[a]:
                a_tmp = {}
                for i in range(len(cas_order)):
                    try:
                        a_tmp[i] = self.dataDict[a][cas_order[i]]
                    except KeyError:
                        pass
                self.dataDict[a] = a_tmp
        #CHANGE ALSO OPTSTATE ACCORDINGLY (to avoid mismach in cobram.py workflow)
        optstate = None
        for i in range(len(cas_order)):
            if cas_order[i]  == self.dataDict["optstate"]:
                optstate = i 
        self.dataDict["optstate"] = optstate









