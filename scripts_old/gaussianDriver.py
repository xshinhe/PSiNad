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


######################################################################################################################

def _optionInRoute(keys, options):
    """This function returns True when one of the options of the list options is
    found in the route section of the keys text, that is when the option is contained
    in one of the lines that preceed the first empty line"""

    # initialize to false
    optionIsInRoute = False

    # loop over the lines
    for line in keys:
        # when one of the options is found in the line, set optionIsInRoute to True (search case-insensitive!)
        optionIsInRoute = optionIsInRoute or any([re.search(s, line, re.IGNORECASE) for s in options])
        # if the line is empty, break loop
        if line.strip() == '': break

    # return value
    return optionIsInRoute


#####################################################################################################

def _energiesFromGaussianOptOut(fileText):
    """ use regexpr to extract the list of total energies from the text of a GAUSSIAN output file """

    eString = 'Energy= *(-?[0-9]*.[0-9]*) * NIter='
    fString = 'Maximum Force *([0-9]*.[0-9]*)'
    dString = 'Maximum Displacement *([0-9]*.[0-9]*)'

    eFind = re.findall(eString, fileText)
    fFind = re.findall(fString, fileText)
    dFind = re.findall(dString, fileText)

    enList, forceList, displList = [], [], []
    for en, f, d in zip(eFind, fFind, dFind):
        enList.append(float(en))
        forceList.append(float(f))
        displList.append(float(d))

    return enList, forceList, displList


#####################################################################################################

class GaussianInput:
    """Store the data for the execution of Gaussian and process it to check its consistency"""

    def __init__(self, keys, gen, gaussadd, gaussweights, coords, symbols, otheropt):
        """Constructor of the gaussianInput class, requires these arguments:
        - keys: string with the text of the gaussian inp file from the route section to the charge/mult line
        - gen: string with the text of the basis set definition (used when Gen is specified in the input file)
        - gaussadd: other options (e.g. the weights of the state averaged CASSCF)
        - coords: coordinates of the atoms
        - symbols: symbols of the atoms
        - otheropt: dictionary with other required options for the QM calculation
                    ("forces", "restart", "suppress_nosymm_error" ...)
        -
        """

        # initialize gaussianInput instance attributes
        # sections of the inp file given : keys (route section + charge/mult), gen (basis set), gaussadd (other)
        self.keys, self.gen, self.gaussadd, self.gaussweights = None, None, None, None
        # coordinates and symbols of the molecule
        self.coords, self.symbols = None, None
        # dictionary with other options
        self.otheropt = {}

        # store the sections of the gaussian input files already defined
        self.keys = keys
        self.gen = gen
        self.gaussadd = gaussadd
        self.gaussweights = gaussweights

        # store the information about the molecular geometry
        self.coords = coords
        self.symbols = symbols

        # process the dictionary with the other options, to check and complete it
        for key in ["forces", "SP", "verbose", "restart", "suppress_nosymm_error", "suppress_TDDFT", "charges", "field", "tdcouplings",
                    "nstate", "forcestate", "CIopt", "upper_state"]:
            if key in otheropt:
                self.otheropt[key] = otheropt[key]
            else:
                self.otheropt[key] = None

        # check if the gaussian command line contains the NOSYM option
        symmetry = not _optionInRoute(self.keys, ['NoSymmetry', 'NoSymm', 'NoSym',
                                                  'Symmetry=None', 'Symmetry=\(.*?None.*?\)'])

        # if symmetry is turned on and the calculation is of QM/MM type, give error/warning depend on command 192
        if symmetry and self.otheropt["charges"]:
            if otheropt["suppress_nosymm_error"]:
                logwrt.writewarning("""Allowing Gaussian calculation with symmetry... 
        Use this option with care: Gaussian expects all charges and geometries \
        in the standard orientation!!""")
            else:
                logwrt.fatalerror("""The NoSymmetry keyword was not found in the !gaussian input section!
        If you know what you are doing, you can disable this error with command 192, 
        otherwise please add NoSymmetry/NoSymm to the !gaussian ... ?gaussian input section. """)

        # check if the gaussian command line contains the GEN option
        genBasis = _optionInRoute(self.keys, ['GEN'])
        if genBasis and len(gen) == 0:
            logwrt.fatalerror("""User-specified input of the basis set is required with the GEN keyword, 
                    but no !gen ... ?gen group was found in the COBRAMM command file!
                    Please specify a basis set in the Gaussian route section or provide a 
                    !gen ... ?gen group in the command file. """)
        
        #check if keywords are not double written because are present in the cobram.command and they should not be
        ischarge = _optionInRoute(self.keys, ["charge"])
        if ischarge:
            logwrt.fatalerror("""keyword 'charge' defined in the cobram.command. Please, remove it.
            COBRAMM will take care of it!""")

    # =============================================================================================================

    def fileText(self, inpFileName, memory="500MB", nproc=1, gversion="g16"):
        """Prepare the text of gaussian input file named after the inpFileName variable (input file
        should be named <inpFileName>.com, so inpFileName should not include the extension!), using the
        input data defined in the instance of gaussianInput. The memory (var memory) and
        the number of cores (var nproc) to use in the execution and the version of Gaussian (gversion)
         are given as input because they are decided when running the Gaussian calculation and are not
         considered parameters of the QM calculation"""

        # initialize the text of the input file
        inputText = ""

        # add the section with the names of the auxiliary files
        inputText += '%chk={0}\n%rwf={0}\n%int={0}\n%d2e={0}\n'.format(inpFileName)
        # add the number of cores to use and the memory of the calculation
        inputText += '%Nproc={0}\n%Mem={1}\n'.format(nproc, memory)

        # add the gaussian keys using the variable self.keys up to the first empty line
        for line in self.keys:
            inputText += line.lower()
            if line.strip() == '':
                break
            else:
                inputText += "\n"

        # if forces are required in the output
        if self.otheropt["forces"]: inputText += " force"
        # if the calculation has background charges
        if self.otheropt["charges"]: inputText += " charge"
        # if computing electric field is required
        if self.otheropt["field"]: inputText += " prop=(field,read)"
        # check if a command for producing charges is there, otherwise add the option POP=(FULL,CHELPG)
        if not _optionInRoute(self.keys, ['CHELP', 'MK']): inputText += " pop=(full,chelpg)"
        # if a chk file is provided for restarting the orbitals
        if self.otheropt["restart"]: inputText += " guess=(read)"
        # options that are always needed
        if not _optionInRoute(self.keys, ['gfinput']): inputText += " gfinput"
        inputText += " scf=tight"

        # when it is requested to compute td couplings, require tamm-damkov and force printing of all excitations
        if self.otheropt["tdcouplings"]:
            if _optionInRoute(self.keys, ['td']):  # when the TD option is active, substitute it with TDA (Tamm-Damkov)
                logwrt.writewarning("Calculation of time derivative couplings is only possible with TDA...replacing TD keyword in Gaussian input\n")
                inputText = inputText.replace(" td ", " tda ")  # without subsequent options
                inputText = inputText.replace(" td(", " tda(")  # with subsequent options (syntax 1)
                inputText = inputText.replace(" td=(", " tda(")  # with subsequent options (syntax 2)
            else:
                inputText += " tda"
            inputText += " IOp(9/40=16)"  # add option to print the whole set of excitation coefficients
        inputText += "\n"

        # when it is requested to change the electronic state number
        addTDSection = None
        if self.otheropt["nstate"] is not None:
            if _optionInRoute(self.keys, ['td', 'tda']):  # when the calculation is TD or TDA
                if self.otheropt["nstate"] == 0:
                    findtdgroup = re.findall("(tda?\(.*?\))", inputText)
                    if findtdgroup:
                        addTDSection = findtdgroup[0]
                        inputText = inputText.replace(findtdgroup[0], "")
                    else:  # in this case there is an isolated td / tda keyword
                        addTDSection = "tda"
                        inputText = inputText.replace(" td ", " ")
                        inputText = inputText.replace(" tda ", " ")
                else:
                    findroot = re.findall("(root *= *[0-9]*)", inputText)  # look for an expression of type "root = N"
                    if findroot:
                        if self.otheropt["nstate"] != int(findroot[0].split("=")[1]) and \
                                not self.otheropt["forcestate"]:
                            logwrt.fatalerror("Mismatch of the electronic state between command file and QM input")
                        inputText = inputText.replace(findroot[0], "root={0}".format(self.otheropt["nstate"]))
                    else:
                        findtdgroup = re.findall("(tda?=?\(.*?\))", inputText)  # look for a group tda( ) or td( )
                        if findtdgroup:
                            inputText = inputText.replace(findtdgroup[0], findtdgroup[0][:-1] +
                                                          ",root={0})".format(self.otheropt["nstate"]))
                        else:  # in this case there is an isolated td / tda keyword
                            inputText = inputText.replace(" td ", " td(root={0}) ".format(self.otheropt["nstate"]))
                            inputText = inputText.replace(" tda ", " tda(root={0}) ".format(self.otheropt["nstate"]))
            else:  # this means that this is not a td/tda calculation, the "nstate" option is not implemented
                logwrt.fatalerror("cannot change the electronic state number for these Gaussian options")

        # add a blank line, the title, another blank line and then charge and multeplicity
        inputText += "\nQM single-point calculation for COBRAMM\n\n"
        inputText += self.keys[-1] + "\n"

        # write molecular structure
        for atom in zip(self.symbols, *self.coords):
            inputText += "{0} {1:16.8f} {2:16.8f} {3:16.8f}\n".format(*atom)
        inputText += "\n"

        # inserting charges
        if self.otheropt["charges"]:
            coords = self.otheropt["charges"][0]
            chrg = self.otheropt["charges"][1]
            for atom in zip(chrg, *coords):
                inputText += "{1:14.8f} {2:14.8f} {3:14.8f} {0:16.8f}\n".format(*atom)
            inputText += "\n"

        # insert basis set information only when the GEN keyword is found in the given route options
        if _optionInRoute(self.keys, ['GEN']):
            for line in self.gen:
                inputText += line + "\n"
            inputText += "\n"

        # in g16 the gaussweights section (weights of the CASSCF) should be before the EF
        if gversion == 'g16' and len(self.gaussweights) != 0:
            for line in self.gaussweights:
                inputText += line + '\n'
            inputText += '\n'

        # inserting EF for gradient computation
        if self.otheropt["field"]:
            for coords in zip(*self.otheropt["field"][0]):
                inputText += "{0:14.8f} {1:14.8f} {2:14.8f}\n".format(*coords)
            inputText += '\n'

        # in g16 the gaussweights section (weights of the CASSCF) should be after the EF
        if gversion != 'g16' and len(self.gaussadd) != 0:
            for line in self.gaussweights:
                inputText += line + '\n'
            inputText += '\n'

        # add the remaining part of the input, stored in the gaussadd variable
        for line in self.gaussadd:
            inputText += line + '\n'
        inputText += '\n'

        # in case of a TD calculation with ground state gradient, append a second input at the end
        if not self.otheropt["suppress_TDDFT"] and addTDSection is not None:
            # add the "link" separation
            inputText += "--Link1--\n"
            # add the section with the names of the auxiliary files
            inputText += '%chk={0}\n%rwf={0}\n%int={0}\n%d2e={0}\n'.format(inpFileName)
            # add the number of cores to use and the memory of the calculation
            inputText += '%Nproc={0}\n%Mem={1}\n'.format(nproc, memory)
            # add the route line
            inputText2 = "" 
            for line in self.keys:
                inputText2 += line.lower()
                if line.strip() == '':
                    break
                else:
                    inputText2 += "\n"
            inputText2 = inputText2.replace(" td ", " tda ")  # without subsequent options
            inputText2 = inputText2.replace(" td(", " tda(")  # with subsequent options
            inputText2 = inputText2.replace(" force ", " ")   # remove force  
            inputText2 = inputText2.replace(",eqsolv", "")   # ES not in equilibrium when dynamics in GS
            inputText2 = inputText2.replace("eqsolv,", "")   # ES not in equilibrium when dynamics in GS
            inputText2 = inputText2.replace("eqsolv", "")   # ES not in equilibrium when dynamics in GS
            # if the calculation has background charges
            if self.otheropt["charges"]: inputText2 += " charge"
            inputText2 += " gfinput"
            inputText2 += " scf=tight"
            inputText += inputText2
            # add restart options
            #inputText += " guess=(read) geom=(check) " + addTDSection + " IOp(9/40=16) \n"
            inputText += " guess=(read) geom=(check) IOp(9/40=16) \n"
            # add a blank line, the title, another blank line and then charge and multeplicity
            inputText += "\nQM single-point calculation for COBRAMM\n\n"
            inputText += self.keys[-1] + "\n\n"

            # inserting charges also for second calculation
            if self.otheropt["charges"]:
                coords = self.otheropt["charges"][0]
                chrg = self.otheropt["charges"][1]
                for atom in zip(chrg, *coords):
                    inputText += "{1:14.8f} {2:14.8f} {3:14.8f} {0:16.8f}\n".format(*atom)
                inputText += "\n"

        # in case of a CI optimization (with gmean branching plane) we also need gradient of lower state
        if self.otheropt["CIopt"]:
            # set lower state and determine if it is GS (DFT only) or not (TD required)
            lower_state = self.otheropt["upper_state"] - 1
            if lower_state == 1:
                GS = True
            else:
                GS = False
            # add the "link" separation
            inputText += "--Link1--\n"
            # add the section with the names of the auxiliary files
            inputText += '%chk={0}\n%rwf={0}\n%int={0}\n%d2e={0}\n'.format(inpFileName)
            # add the number of cores to use and the memory of the calculation
            inputText += '%Nproc={0}\n%Mem={1}\n'.format(nproc, memory)
            # add the route line
            inputText2 = "" 
            for line in self.keys:
                inputText2 += line.lower()
                if line.strip() == '':
                    break
                else:
                    inputText2 += "\n"
            if GS:
                findtdgroup = re.findall("(tda?=?\(.*?\))", inputText2)
                if findtdgroup:
                    inputText2 = inputText2.replace(findtdgroup[0], "")
                else:  # in this case there is an isolated td / tda keyword
                    inputText2 = inputText2.replace(" td ", " ")
                    inputText2 = inputText2.replace(" tda ", " ")  
            else:
                findroot = re.findall("(root *= *[0-9]*)", inputText2)  # look for an expression of type "root = N"
                if findroot:
                    inputText2 = inputText2.replace(findroot[0], "root={0}".format(lower_state - 1))
                else:
                    findtdgroup = re.findall("(tda?\(.*?\))", inputText2)  # look for a group tda( ) or td( )
                    if findtdgroup:
                        inputText2 = inputText2.replace(findtdgroup[0], findtdgroup[0][:-1] +
                                                      ",root={0})".format(lower_state - 1))
                    else:  # in this case there is an isolated td / tda keyword
                        inputText2 = inputText2.replace(" td ", " td(root={0}) ".format(lower_state - 1))
                        inputText2 = inputText2.replace(" tda ", " tda(root={0}) ".format(lower_state - 1))
            if self.otheropt["forces"]: inputText2 += " force"
            # if the calculation has background charges
            if self.otheropt["charges"]: inputText2 += " charge"
            # if computing electric field is required
            if self.otheropt["field"]: inputText2 += " prop=(field,read)"
            # check if a command for producing charges is there, otherwise add the option POP=(FULL,CHELPG)
            if not _optionInRoute(self.keys, ['CHELP', 'MK']): inputText2 += " pop=(full,chelpg)"
            inputText2 += " gfinput"
            inputText2 += " scf=tight"
            inputText += inputText2
            # add restart options
            inputText += " guess=(read) geom=(check)\n"
            # add a blank line, the title, another blank line and then charge and multeplicity
            inputText += "\nQM single-point calculation for COBRAMM\n\n"
            inputText += self.keys[-1] + "\n\n"

            # inserting charges also for second calculation
            if self.otheropt["charges"]:
                coords = self.otheropt["charges"][0]
                chrg = self.otheropt["charges"][1]
                for atom in zip(chrg, *coords):
                    inputText += "{1:14.8f} {2:14.8f} {3:14.8f} {0:16.8f}\n".format(*atom)
                inputText += "\n"

            # insert basis set information only when the GEN keyword is found in the given route options
            if _optionInRoute(self.keys, ['GEN']):
                for line in self.gen:
                    inputText += line + "\n"
                inputText += "\n"
    
            # in g16 the gaussweights section (weights of the CASSCF) should be before the EF
            if gversion == 'g16' and len(self.gaussweights) != 0:
                for line in self.gaussweights:
                    inputText += line + '\n'
                inputText += '\n'
    
            # inserting EF for gradient computation
            if self.otheropt["field"]:
                for coords in zip(*self.otheropt["field"][0]):
                    inputText += "{0:14.8f} {1:14.8f} {2:14.8f}\n".format(*coords)
                inputText += '\n'
    
            # in g16 the gaussweights section (weights of the CASSCF) should be after the EF
            if gversion != 'g16' and len(self.gaussadd) != 0:
                for line in self.gaussweights:
                    inputText += line + '\n'
                inputText += '\n'
    
            # add the remaining part of the input, stored in the gaussadd variable
            for line in self.gaussadd:
                inputText += line + '\n'
            inputText += '\n'


        # return the complete text of the input file
        return inputText


###################################################################################################################

class GaussianOutput(QMOutput):

    def __init__(self, name, calcdir, SPcalc=False, verbose=False):

        # define options of the base class
        QMOutput.__init__(self)

        # Add additional gaussian-specific attributes to the data dictionary
        # FILE NAMES, DIRECTORIES, INPUT/OUTPUT
        self.dataDict["name"] = name  # name of the calc, it is the base name of the I/O files (.com, .chk, .log)
        self.dataDict["dir"] = calcdir  # name of the directory where the output files are stored
        self.dataDict["outfile"] = None   # string with the text of the Gaussian log file, initialized to None
        # ADDITIONAL PHYSICAL PROPERTIES EXTRACTED FROM THE OUTPUT FILE
        self.dataDict["elfield"] = None   # electric field at the position of hte point charges
        self.dataDict["natoms"] = None   # number of atoms in the QM calculation
        # INFO ON GAUSSIAN EXECUTION
        self.dataDict["termination"] = None   # final termination: 0 = Normal termination, 1 = Error termination
        self.dataDict["errormsg"] = None   # in case of Error termination, error message
        self.dataDict["signs"] = [] # list of WF signs from previous step to apply in phase correction
        self.dataDict["psioverlap"] = None # array of the WF overlaps between consecutive steps
        self.dataDict["nroots"] = 1 # number of states (GS + ESs): initialized to 1 and then increased in case of TDDFT

        # store the gaussian log removing AO and MO definitons, CI expansions, etc.
        if SPcalc == False and verbose == False:
            self.dataDict["outfile"] = self._filterOutput(os.path.join(calcdir, name + ".log"))
        else:
            with open(os.path.join(calcdir, name + ".log")) as fout:
                self.dataDict["outfile"] = fout.read()

        # use full file for reading properties 
        with open(os.path.join(calcdir, name + ".log")) as fout:
            outfile = fout.read()

        # initialize variables
        tempstate = None

        # loop over the lines of the log file
        output = outfile.splitlines()
        if output[-4].split()[0].strip() == 'Error':
            self.dataDict["termination"] = 1
            self.dataDict["errormsg"] = output[-11:-4]
            return
        for i in range(len(output)):

            #extract the number of atoms of the QM calculation
            try:
                if output[i].split()[0].strip() == 'NAtoms=':
                    self.dataDict["natoms"] = int(output[i].split()[1].strip())
            except (IndexError, ValueError):
                pass

            # extract the number of occupied orbitals 
            try:
                if output[i].split()[1].strip() == 'alpha' and output[i].split()[2].strip() == 'electrons':
                    alpha = int(output[i].split()[0].strip())
                    if int(output[i].split()[3].strip()) != alpha:
                        beta = int(output[i].split()[3].strip())
                        self.log += 'Unrestricted GS calculation with {0:d} alpha and {1:d} beta electrons\n'.format(alpha, beta)
                        nocc = alpha + 1  
                    else:
                        nocc = alpha
            except (IndexError, ValueError):
                pass

            # extract SCF energy from log file, when available
            try:
                if output[i].split()[0].strip() == 'SCF' and output[i].split()[1].strip() == 'Done:':
                    self.dataDict["energy"][0] = float(output[i].split()[4].strip())
                    self.dataDict["optstate"] = 0
                    self.log += "SCF energy is {0:12.8f} Hartree\n".format(self.dataDict["energy"][0])
            except (IndexError, ValueError):
                pass

            # extract MP2 energy from log file, when available
            try:
                for mp2str1 in output[i].split('\\'):
                    if mp2str1.find('EUMP2') != -1:
                        self.dataDict["energy"][0] = float(mp2str1.split()[5])
                        self.dataDict["optstate"] = 0
                        self.log += "MP2 energy is {0:12.8f} Hartree\n".format(self.dataDict["energy"][0])
            except (IndexError, ValueError):
                pass

            # get the TD-DFT excited state energies, when available
            try:
                if output[i].split()[0].strip() == 'Convergence' and output[i].split()[1].strip() == 'achieved':
                    TDstates = int(output[i - 1].split()[1])
                    count = TDstates
                    for j in range(TDstates):
                        self.dataDict["energy"][j + 1] = float(output[i - count].split()[3]) / constants.Hartree2eV + \
                                   self.dataDict["energy"][0]
                        count = count - 1
                    self.dataDict["nroots"] = TDstates + 1
            except (IndexError, ValueError):
                pass

            # TDDFT occasionally calculates many more startes, therefore we need to find out the real number of requested states
            try:
                if output[i].split()[5].strip() == 'velocity' and output[i].split()[6].strip() == 'dipole':
                    nstates = int(output[i - 1].split()[0])
                    self.log += 'Number of excited states: {0:3d}\n'.format(nstates)
                    for j in range(nstates+1,TDstates+1):
                        del self.dataDict["energy"][j] 
                    for j in range(nstates):
                        self.log += 'DFT energy for state {0:3d} (ES {1:2d}) is {2:12.8f} Hartree\n'.format(
                            j + 2, j + 1, self.dataDict["energy"][j+1])
                    self.dataDict["nroots"] = nstates + 1
            except (IndexError, ValueError):
                pass

            # determine the state of interest in TD-DFT calculation
            try:
                if output[i].split()[0].strip() == 'Excited' and output[i].split()[1].strip() == 'State':
                    tempstate = output[i].split()[2]
                elif output[i].split()[0].strip() == 'This' and output[i].split()[1].strip() == 'state':
                    self.dataDict["optstate"] = int(tempstate.split(':')[0])
            except (IndexError, ValueError):
                pass

            # get the CASSCF energies, when available
            try:
                if output[i].split()[2].strip() == 'EIGENVALUE' and output[i].split()[0].strip() == '(':
                    energy = float(output[i].split()[3].strip())
                    istate = int(output[i].split()[1].strip(" ()"))
                    self.log += 'CASSCF energy for state {0:3d} is {1:12.8f} Hartree\n'.format(istate, energy)
                    self.dataDict["energy"][istate - 1] = energy
                    self.dataDict["optstate"] = istate - 1
            except (IndexError, ValueError):
                pass

            # # get energy from solvation
            # try:
            #     if output[i].find('with all non electrostatic terms') != -1:
            #         # Reset energy
            #         self.log += "Energy after re-setting {0:12.8f}\n".format(self.dataDict["energy"][0])
            # except (IndexError, ValueError):
            #     pass

            # get self-energy of the background charges
            if output[i].split('=')[0].strip() == 'Self energy of the charges':
                self.dataDict["selfenergy"] = float(output[i].split()[6].strip())
                self.log += "Gaussian energy includes electrostatic interaction between external point charges\n" + \
                            "For consistency this self-energy {0:12.8f} Hartree will be substracted.\n\n".format(
                                    self.dataDict["selfenergy"])

            # get gradient
            if output[i].strip() == 'Center     Atomic                   Forces (Hartrees/Bohr)':
                # initialize the dictionary entry for the gradient
                gradient = [[], [], []]
                # loop over following lines
                j = 0
                while True:
                    # try to read the block of numbers with the gradient, if an exception is raised finish reading
                    try:
                        element = output[i + j + 3].split()
                        gradient[0].append(-float(element[2]))
                        gradient[1].append(-float(element[3]))
                        gradient[2].append(-float(element[4]))
                    except (IndexError, ValueError):
                        break
                    # move to following line
                    j += 1
                    # when the gradient is too large, abort calculation
                    if abs(gradient[0][-1]) > 1.0 or abs(gradient[1][-1]) > 1.0 or abs(gradient[2][-1]) > 1.0:
                        logwrt.fatalerror('Large gradient during analytical computation. Check last QM '
                                          'output {0}'.format(os.path.join(calcdir, name + ".log")))
                # now store the gradient extracted from the gaussian log
                self.dataDict["gradient"] = {self.dataDict["optstate"]: gradient}

            # extract ESP charges
            if output[i].strip() == 'Fitting point charges to electrostatic potential':
                # initialize the dictionary entry for the charges
                self.dataDict["charges"] = []
                # loop over following lines
                j = 0
                while True:
                    # try to read the block of numbers with the gradient, if an exception is raised finish reading
                    try:
                        element = output[i + j + 4].split()
                        self.dataDict["charges"].append(float(element[2]))
                    except (IndexError, ValueError):
                        break
                    # move to following line
                    j += 1

            # get electrostatic field at MM point charges
            if output[i].strip() == 'Electrostatic Properties (Atomic Units)':
                # loop over following lines
                j = 0
                elfield = [], [], []
                while True:
                    # try to read the block of numbers with the gradient, if an exception is raised finish reading
                    try:
                        element = output[i + j + 6 + self.dataDict["natoms"]].split()
                        Ex, Ey, Ez = float(element[2]), float(element[3]), float(element[4])
                        elfield[0].append(Ex), elfield[1].append(Ey), elfield[2].append(Ez)
                    except (IndexError, ValueError):
                        break
                    # move to following line
                    j += 1
                if elfield[0]:
                    self.dataDict["elfield"] = {self.dataDict["optstate"]: elfield}
                else:
                    self.dataDict["elfield"] = None

            # get dipole moments
            if output[i].strip() == 'Dipole moment (field-independent basis, Debye):':
                dip = output[i + 1].split()
                self.dataDict["dipole"] = float(dip[1]) * constants.Debye2AU, float(dip[3]) * constants.Debye2AU, \
                    float(dip[5]) * constants.Debye2AU, float(dip[7]) * constants.Debye2AU

            # terminate the reading of the main output file
            if "Normal termination" in output[i]:
                nmain = i
                break

        addEne = {}
        for i in range(nmain, len(output)):

            # extract SCF energy from log file, when available
            try:
                if output[i].split()[0].strip() == 'SCF' and output[i].split()[1].strip() == 'Done:':
                    scf = float(output[i].split()[4].strip())
                    self.log += "SCF energy for secondary output is {0:12.8f} Hartree\n".format(scf)
                    self.dataDict["lowerstate"] = 0
                    addEne[0] = scf
            except (IndexError, ValueError):
                pass

            # extract MP2 energy from log file, when available
            try:
                for mp2str1 in output[i].split('\\'):
                    if mp2str1.find('EUMP2') != -1:
                        mp2 = float(mp2str1.split()[5])
                        self.log += "MP2 energy for secondary output is {0:12.8f} Hartree\n".format(mp2)
                        self.dataDict["lowerstate"] = 0
                        addEne[0] = mp2
            except (IndexError, ValueError):
                pass


            # get the TD-DFT excited state energies, when available
            try:
                if output[i].split()[0].strip() == 'Convergence' and output[i].split()[1].strip() == 'achieved':
                    TDstates = int(output[i - 1].split()[1])
                    count = TDstates
                    for j in range(TDstates):
                        addEne[j+1] = float(output[i - count].split()[3]) / constants.Hartree2eV + scf
                        count = count - 1
            except (IndexError, ValueError):
                pass

            # TDDFT occasionally calculates many more startes, therefore we need to find out the real number of requested states
            try:
                if output[i].split()[5].strip() == 'velocity' and output[i].split()[6].strip() == 'dipole':
                    nstates = int(output[i - 1].split()[0])
                    self.log += 'Number of excited states in secondary output: {0:3d}\n'.format(nstates)
                    for j in range(nstates+1,TDstates+1):
                        del addEne[j]
                    for j in range(nstates):
                        self.log += 'DFT energy for state {0:3d} (ES {1:2d}) is {2:12.8f} Hartree\n'.format(
                            j + 2, j + 1, addEne[j+1])
            except (IndexError, ValueError):
                pass

            # determine the state of interest in TD-DFT calculation
            try:
                if output[i].split()[0].strip() == 'Excited' and output[i].split()[1].strip() == 'State':
                    tempstate = output[i].split()[2]
                elif output[i].split()[0].strip() == 'This' and output[i].split()[1].strip() == 'state':
                    self.dataDict["lowerstate"] = int(tempstate.split(':')[0])
            except (IndexError, ValueError):
                pass

            # get the CASSCF energies, when available
            try:
                if output[i].split()[2].strip() == 'EIGENVALUE' and output[i].split()[0].strip() == '(':
                    if not casscf: casscf = []
                    self.log += 'CASSCF energies from secondary output\n'
                    energy = float(output[i].split()[3].strip())
                    istate = int(output[i].split()[1].strip(" ()"))
                    self.log += 'CASSCF energy for state {0:3d} is {1:12.8f} Hartree\n'.format(istate, energy)
                    addEne[istate - 1] = energy
            except (IndexError, ValueError):
                pass

            # get gradient
            if output[i].strip() == 'Center     Atomic                   Forces (Hartrees/Bohr)':
                # initialize the dictionary entry for the gradient
                gradient = [[], [], []]
                # loop over following lines
                j = 0
                while True:
                    # try to read the block of numbers with the gradient, if an exception is raised finish reading
                    try:
                        element = output[i + j + 3].split()
                        gradient[0].append(-float(element[2]))
                        gradient[1].append(-float(element[3]))
                        gradient[2].append(-float(element[4]))
                    except (IndexError, ValueError):
                        break
                    # move to following line
                    j += 1
                    # when the gradient is too large, abort calculation
                    if abs(gradient[0][-1]) > 1.0 or abs(gradient[1][-1]) > 1.0 or abs(gradient[2][-1]) > 1.0:
                        logwrt.fatalerror('Large gradient during analytical computation. Check last QM '
                                          'output {0}'.format(os.path.join(calcdir, name + ".log")))
                # now store the gradient extracted from the gaussian log
                self.dataDict["gradient"][self.dataDict["lowerstate"]] = gradient

            # get electrostatic field at MM point charges
            if output[i].strip() == 'Electrostatic Properties (Atomic Units)':
                # loop over following lines
                j = 0
                elfield = [], [], []
                while True:
                    # try to read the block of numbers with the gradient, if an exception is raised finish reading
                    try:
                        element = output[i + j + 6 + self.dataDict["natoms"]].split()
                        Ex, Ey, Ez = float(element[2]), float(element[3]), float(element[4])
                        elfield[0].append(Ex), elfield[1].append(Ey), elfield[2].append(Ez)
                    except (IndexError, ValueError):
                        break
                    # move to following line
                    j += 1
                if elfield[0]:
                    self.dataDict["elfield"][self.dataDict["lowerstate"]] = elfield
                else:
                    self.dataDict["elfield"] = None

        nstart = 100
        if addEne:
            for i, e in addEne.items():
                if i == 0:
                    pass
                else:
                    #old version was saving a double copy of existing ESs from item 100 onwards
                    #if i in self.dataDict["energy"]:
                    #    self.dataDict["energy"][nstart+i] = e
                    #else:
                    #    self.dataDict["energy"][i] = e
                    # F: I think we do not need this double copy (which is messing up the energy printout in cobramm.log), so I will save only states which are not present in primary output
                    # restore to old version if this creates troubles 
                    if i not in self.dataDict["energy"]:
                        self.dataDict["energy"][i] = e

        # extracting CIS/TD excitation coefficients
        if "Excited states from <AA,BB:AA,BB> singles matrix:" in outfile:
            self.log += '\nExtracting CIS/TD(A) coefficients from gaussian output\n'
            # convert "cis_coeffs" to from None to empty dictionary:
            self.dataDict["cis_coeffs"] = {}
            # define regular expressions to isolate the sections of each excited state and to list excit coeff.s
            stateRegExp = "Excited State *([0-9]*): *Singlet\-\?Sym *(-?[0-9]*\.[0-9]*) eV *(-?[0-9]*\.[0-9]*) nm *f=(-?[0-9]*\.[0-9]*) *<S\*\*2>=([0-9]*\.[0-9]*)\n([\s\S]*?)\n ?\n" 
            excitationRegExp = "([0-9]*) *-> *([0-9]*) *(-?[0-9]*\.[0-9]*)"
            # loop over states and then over each excitation line
            for nstate, eV, nm, f, S, text in re.findall(stateRegExp, outfile):
                cisExpansion = {}
                for occorb, virorb, coeff in re.findall(excitationRegExp, text):
                    cisExpansion[(int(occorb), int(virorb))] = float(coeff) * math.sqrt(2)
                self.dataDict["cis_coeffs"][int(nstate)] = cisExpansion
                cisExpansion = self._truncateExpansion(cisExpansion, 0.9)
                self.log += " Excited State   {0:d} (ES {1:d}): {2:.2f} eV ({3:d} nm) f={4:.3f}\n".format(int(nstate)+1, int(nstate), float(eV), int(float(nm)), float(f))
                count=0
                for ex, coeff in cisExpansion.items():
                    occ = 'HOMO'
                    if ex[0] != nocc:
                        occ += '-{0:d}'.format(nocc-ex[0])
                    virt = 'LUMO'
                    if ex[1] != nocc+1:
                        virt += '+{0:d}'.format(ex[1]-nocc-1)
                    self.log += '{0:d} -> {1:d} ({2} -> {3}): {4:5.2f} ({5:5.2f})\n'.format(ex[0], ex[1], occ, virt, coeff, coeff**2)
                    count += 1
                    if count == 5:
                        break
                self.log += "\n"

        self.log += "\n"

        # finish collecting data from .log file
        # now collect some info on the job execution

        # check if Gaussian termination is normal or error, and store excerpt of log file with error message
        if output[-1].split()[0].strip() == 'Normal':
            self.dataDict["termination"] = 0
        elif output[-4].split()[0].strip() == 'Error':
            self.dataDict["termination"] = 1
            self.dataDict["errormsg"] = output[-11:-4]

    # =============================================================================================================

    def __del__(self):
        """Destructor for the gaussianOutput class: not only the memory allocated for the oject attributes
        (the dictionary self.dataDict) needs to be released, also the Gaussian I/O files stored on disk
        can be safely removed... """

        # permanently remove directory and gaussian files
        try:
            if not logwrt.DEBUG_COBRAMM_RUN:
                if self.dataDict["termination"] == 1:
                    source=self.dataDict["dir"]+"/gaussian-QM.log"
                    shutil.move(source,"gaussian-QM_err.log")
                shutil.rmtree(self.dataDict["dir"])
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
            if fulllog[iLine].find("Background charge distribution read from input stream") != -1:
                toPrint = 0
            elif fulllog[iLine].find("Pt Chg Charge") != -1:
                toPrint = 1
            if fulllog[iLine].find("AO basis set in the form of general basis input") != -1:
                toPrint = 0
            elif fulllog[iLine].find("primitive gaussians") != -1:
                toPrint = 1
            elif fulllog[iLine].find("Excitation energies and oscillator strengths") != -1:
                cleanlog.append(self.extractCIcoeff(log))
                toPrint = 0
            elif fulllog[iLine].find(" SavETr:  write IOETrn=") != -1 or fulllog[iLine].find("Leave Link  914") != -1:
                toPrint = 1
            elif fulllog[iLine].find("Molecular Orbital Coefficients") != -1:
                toPrint = 0
            elif fulllog[iLine].find("Full Mulliken population analysis") != -1 or fulllog[iLine].find("Leave Link  601") != -1:
                toPrint = 1

            if toPrint == 1: 
                cleanlog.append(fulllog[iLine])

        return "\n".join(cleanlog)

    # =============================================================================================================

    def extractCIcoeff(self,outfile):
        """Print leading CI coefficients"""

        cleanlog = [" Excitation energies and oscillator strengths:"]
        stateRegExp = "Excited State *([0-9]*): *Singlet\-\?Sym *(-?[0-9]*\.[0-9]*) eV *(-?[0-9]*\.[0-9]*) nm *f=(-?[0-9]*\.[0-9]*) *<S\*\*2>=([0-9]*\.[0-9]*)\n([\s\S]*?)\n ?\n"
        excitationRegExp = "([0-9]*) *-> *([0-9]*) *(-?[0-9]*\.[0-9]*)"
        # loop over states and then over each excitation line
        for nstate, eV, nm, f, S, text in re.findall(stateRegExp, outfile):
            cisExpansion = {}
            for occorb, virorb, coeff in re.findall(excitationRegExp, text):
                cisExpansion[(int(occorb), int(virorb))] = float(coeff) * math.sqrt(2)
            cisExpansion = self._truncateExpansion(cisExpansion, 0.95)
            cleanlog.append(" Excited State   {0:d}:      Singlet-?Sym    {1:.4f} eV  {2:.2f} nm  f={3:.4f}  <S**2>={4:.3f}".format(int(nstate), float(eV), float(nm), float(f), float(S)))
            count=0
            for ex, coeff in cisExpansion.items():
                cleanlog.append("   {0:d} -> {1:d}        {2:9.5f}".format(ex[0], ex[1], coeff))
                count += 1
                if count == 10:
                    break
            cleanlog.append("\n")

        return "\n".join(cleanlog) 

    # =============================================================================================================

    def orbfiles(self):
        """Check their existance and then return a list with the names of the files to be
         saved to store the orbitals to disk"""

        # name of the chk file that needs to be saved
        # the chk is considered as the orbital file because it can be used to produce cube files of the orbitals
        chkname = os.path.join(self.dataDict["dir"], self.dataDict["name"]+".chk")

        # if the file is not there, print warning and return an empty list
        if not os.path.exists(chkname):
            logwrt.writewarning("gaussian file {0} cannot be found: it will not be stored".format(chkname))
            return []
        # otherwise return a list with the chk file alone
        else:
            return [chkname]

    # =============================================================================================================

    def restartfile(self):
        """Check their existance and then return a the names of the files to be
         saved to store the orbitals to disk for restart purpose"""

        # name of the chk file that needs to be saved for restart
        chkname = os.path.join(self.dataDict["dir"], self.dataDict["name"]+".chk")

        # if the file is not there, print warning and return an empty list
        if not os.path.exists(chkname):
            return None
        # otherwise return a list with the chk file alone
        else:
            return chkname

    # =============================================================================================================

    @staticmethod
    def _truncateExpansion(expansion: dict, threshold: float) -> dict:
        """
        Truncate the CI expansion so that the total sum_i |CI_i|^2 > threshold
        """

        # sort CI coeff. in decreasing order
        expansion_sort = sorted(expansion.items(), key=lambda coeff: abs(coeff[1]),reverse=True)

        # determine the length of the truncated expansion
        cumulative_weight = 0
        # default length (determines the length of the CI vector in case no mdv is run)
        maxel = len(expansion_sort)
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

