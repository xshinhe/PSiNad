#!/usr/bin/env python3
#   Coding=utf-8

#   KIDS SCRIPTS
#   Author: xshinhe
#   
#   Copyright (c) 2024 PeKing Univ. - GNUv3 License

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
import os
import shutil
import gzip
import shelve
from typing import Union, List, Dict, Optional, Any, TypedDict
import numpy.typing as npt
from pprint import pprint

# imports of local modules
from Layers import Layers
from constants import element_list
import amberDriver
import molcasDriver
import molproDriver

from kids_log import Timing, Log
from kids_config import Config
import kids_arg

INPUTSTORAGEDIR = "input_files"

class TomlDict(dict):
    def get_nested(self, key_path, default=None):
        keys = key_path.split('.')
        try:
            return reduce(getitem, keys, self)
        except (KeyError, TypeError):
            return default

def readFile(filename: str) -> List[str]:
    """
    Reads the content of a file and returns it as a list of strings.

    Parameters:
    filename (str): The name of the file to be read.

    Returns:
    List[str]: A list where each element is a line from the file.

    Raises:
    SystemExit: Exits the program if the file does not exist.
    """
    try:
        with open(filename, 'r') as f:
            out = f.readlines()
    except IOError:
        pprint.pprint(f'File {filename} does not exist!')
        traceback.print_exc()
        sys.exit(-1)
    return out

def writeFile(filename: str, content: Union[str, List[str]]):
    """
    Writes content to a file. Content can be a string or a list of strings.

    Parameters:
    filename (str): The name of the file to write to.
    content (Union[str, list[str]]): The content to write to the file.

    Raises:
    SystemExit: Exits the program if writing to the file fails.
    """
    try:
        with open(filename, 'w') as f:
            if isinstance(content, List):
                for line in content:
                    f.write(line)
            elif isinstance(content, str):
                f.write(content)
            else:
                print(f'Content {content} cannot be written to file!')
    except IOError:
        print(f'Could not write to file {filename}!')
        traceback.print_exc()
        sys.exit(-1)

class QMout(TypedDict):
    natom: int
    energy: npt.ArrayLike
    gradient: npt.ArrayLike 
    nacv: npt.ArrayLike


def saveinputs():
    if INPUTSTORAGEDIR in os.listdir('.'):
        shutil.rmtree(INPUTSTORAGEDIR)
    os.mkdir(INPUTSTORAGEDIR)
    tocopy = ['real_layers.xyz', 'real.top', 'model-H.top', 'real.crd',
              'geometry.chk', 'gaussian-QM.chk', 'velocity.dat', 'INPORB', 'turboinp']
    # save file when present
    for f in tocopy:
        if f in os.listdir('.'): shutil.copy(f, INPUTSTORAGEDIR)

def GZIP(file):
    with open(file, 'rb') as r_file:
        w_file = gzip.GzipFile(file + '.gz', 'w', 9)
        w_file.write(r_file.read())
        w_file.flush()
        w_file.close()
    os.remove(file)

def GUNZIP(file):
    r_file = gzip.GzipFile(file, 'r')
    with open(os.path.splitext(file)[0], 'w') as w_file:
        w_file.write(r_file.read().decode())
    r_file.close()
    os.unlink(file)

def garbager(geometry, ks_data):
    # the calculation has a QM part, clean up files from QM run
    if "H" in geometry.calculationType:
        if ks_data['QM']['exec'] == 'MOLCAS':  # QM calulation by Molcas
            molcasDriver.clean()
        if ks_data['QM']['exec'] == 'MOLPRO':  # QM calulation by MOLPRO
            molproDriver.clean()

    # the calculation has a MM part, clean up files from MM run
    if "M" in geometry.calculationType or "L" in geometry.calculationType:
        amberDriver.clean()

    # remove optimization files and zip chk file
    toremove2 = ["geometry.com", 'geometry.dat', 'geometry.rwf', 'geometry.d2e', 'geometry.int']
    for f in toremove2:
        if os.path.isfile(f): os.remove(f)
    if os.path.isfile('geometry.chk'): GZIP('geometry.chk')

    # rm geomXXX directories remaining from a parallel calculation
    for i in range(999):
        if os.path.isdir("geom" + str(i)): shutil.rmtree("geom" + str(i))

    toremove = ["velocity_0.dat", "cobram-sef", 'frozenmmatoms.pdb', 'basis', 'auxbasis', 'mos',
                'statistics', 'fort.99', 'test.out', 'infile', 'newfile', 'oldfile']
    for f in toremove:
        if os.path.isfile(f): os.remove(f)

    toremove = ["velocity.dat", "EKIN", "ETOT", "velocityOLD.dat", "TSTEP", "states"]
    for f in toremove:
        if os.path.isfile(f): os.remove(f)

def getCalcType( flag:str ):
   """ get a string that describes the the type of calculation that is requested"""
   if flag == "optxg":
      calcType = "An optimization run"
   elif flag == "ts":
      calcType = "A transition state (TS) search"
   elif flag == "irc":
      calcType = "An intrinsic reaction coordinate (IRC) calculation"
   elif flag == "ci":
      calcType = "A conical intersection (CI) search"
   elif flag == "freqxg":
      calcType = "A harmonic frequency calculation"
   elif flag == "mdv":
      calcType = "A molecular dynamics run"
   elif flag == "nad":
      calcType = "An interface for nad"
   elif flag.isdigit():
      calcType = "An interface for level%d"%int(flag) 
   else:
      fatalError( "calculation type (ks_data 1) '{0}' is not defined. ".format(flag) )

   return calcType

#================================================================================

def getLayers( geometry ):
   """ get a string that describes the layers that are defined in the calculation"""

   if geometry.calculationType == "H":
      layers = "HIGH layer only"
   elif geometry.calculationType == "M":
      layers = "MEDIUM layer only"
   elif geometry.calculationType == "HM":
      layers = "HIGH and MEDIUM layers"
   elif geometry.calculationType == "HL":
      layers = "HIGH and LOW layers"
   elif geometry.calculationType == "ML":
      layers = "MEDIUM and LOW layers"
   elif geometry.calculationType == "HML":
      layers = "HIGH, MEDIUM and LOW layers"
   else:
      fatalError( "layer scheme '{0}' is not defined. ".format(geometry.calculationType) )

   return layers

#================================================================================

def getOptionSummary( ks_data ):
    """ returns a string with a summary of the options defined in the input ks_data list """

    summary = ""
    summary += "*** SYSTEM OPTIONS ***\n"
    try:
        summary += "* allocated memory for Gaussian optimizer: {0}\n".format(ks_data['config']['opt_mem'] )
        summary += "* allocated memory for QM computation: {0}\n".format(ks_data['config']['QM']['mem'] )
        summary += "* nr. of cores to use for QM third-party software: {0}\n".format(ks_data['config']['QM']['ncore'] )
        summary += "* nr. of cores to use for COBRAMM internal parallel numerics routine: {0}\n".format(ks_data['config']['ncore'] )
        if ks_data['config']['MM']['AMBER']['type'] == "serial":
            summary += "* use AMBER serial executable sander\n"
        elif ks_data['config']['MM']['AMBER']['type'] == "gpu":
            summary += "* use AMBER GPU code pmemd.cuda\n"
        elif ks_data['config']['MM']['AMBER']['type'] == "mpi":
            summary += "* use AMBER parallel executable sander.MPI\n"
    except KeyError:
        pass 

    summary += "\n"
    summary += "*** GENERAL OPTIONS ***\n"
    # if int(ks_data[200]) > 0:
    #     summary += "* perform a single point computation every {0} steps (for Molpro and Molcas, requires !singlepoint  block)\n".format(ks_data[200])
    # elif ks_data[200] == "-1":
    #     summary += "* perform a single point computation at every converged IRC geometry (for Molpro and Molcas, requires !singlepoint block)\n"

    try:
        summary += "* number of optimization steps / microsteps in each IRC step / MD steps: {0}\n".format(ks_data[60])
        if ks_data['config']['QM']['store_log_step']  == "0":
            summary += "* do not store QM log files \n"
        else:
            summary += "* save QM log files every {0} steps\n".format(ks_data['config']['QM']['store_log_step'] )
        if ks_data['config']['QM']['store_wfn_step']  == "0":
            summary += "* do not store QM wavefunction files \n"
        else:
            summary += "* save QM wavefunction files every {0} steps\n".format(ks_data['config']['QM']['store_wfn_step'] )
    except KeyError:
        pass

    try:
        summary += "\n"
        summary += "*** MOLECULAR MODEL OPTIONS ***\n"
        if ks_data['config']['MM']['border_atom_flag']  == "1":
            summary += "* border atoms are frozen\n"
        if ks_data['config']['QMMM']['border_atom_force']  == "1":
            summary += "* project the force of the link atoms on the border atoms\n"

        if ks_data['args']['type']  in ["optxg", "ts", "irc", "ci"] and int(ks_data[60]) > 1:
            summary += "\n"
            summary += "*** OPTIMIZATION OPTIONS ***\n"
            summary += "* Gaussian optimization convergence threshold: RMS Force = {0}\n".format(float(ks_data[67])*0.0001)
            summary += "* Gaussian optimization convergence threshold: Max Force = {0}\n".format(float(ks_data[67])*0.0001*1.5)
            summary += "* Gaussian optimization convergence threshold: RMS Displ = {0}\n".format(float(ks_data[67])*0.0001*4.)
            summary += "* Gaussian optimization convergence threshold: Max Displ = {0}\n".format(float(ks_data[67])*0.0001*6.)
            summary += "* maximum step size of an optimization step: {0} bohr\n".format(float(ks_data[68])*0.01)

        if ks_data['args']['type']  == "irc":
            summary += "\n"
            summary += "*** IRC OPTIONS ***\n"
            summary += "* number of intermediate geometries along the MEP computed with IRC: {0} \n".format(ks_data[65])
            summary += "* step size along the MEP: {0} bohr\n".format(float(ks_data[66])*0.01)

        if ks_data['args']['type']  == "mdv":
            summary += "\n"
            summary += "*** MOLECULAR DYNAMICS OPTIONS ***\n"
            summary += "* regular time-step of the molecular dynamics: {0} fs\n".format(ks_data[83])
            if int(ks_data[85]) > 0  and ks_data[14] == '0': 
                summary += "* short time-step when Tully or Tully-Persico is active: {0} fs\n".format(ks_data[84])
                summary += "* energy difference between active state and near-by states for activating short time step: {0} kcal/mol\n".format(ks_data[86])
            if ks_data[85] == "0":
                summary += "* perform adiabatic dynamics\n"
            elif ks_data[85] == "1":
                summary += "* dynamics with hoppings according to Tully's FSSH, with Persico's decoherence correction\n"
            elif ks_data[85] == "2":
                summary += "* dynamics with hoppings according on only energy difference\n"
                summary += "* energy threshold (kcal/mol) to hop: {0}\n".format(ks_data[207])
            if ks_data[85] == '1':
                summary += "* energy threshold (kcal/mol) to activate smaller timestep and hoppings: {0} kcal/mol\n".format(ks_data[86])
            if ks_data[90] == "0":
                summary += "* turn off hopping to the GS\n"
            elif ks_data[90] == "0":
                summary += "* compute hopping probability to the GS in the same way as between ES\n"
            elif ks_data[90] == "2":
                summary += "* force hop to the ground state when the GS-ES energy gap is < {0} kcal/mol\n".format(ks_data[91])
            if ks_data[87] == "0":
                summary += "* switch off back-hopping during Tully's FSSH\n"
            elif ks_data[87] == "1":
                summary += "* switch on back-hopping during Tully's FSSH\n"
            if ks_data[14] == "0":
                summary += "* compute (spatial) non-adiabatic couplings (NACs)\n"
            elif ks_data[14] == "1":
                summary += "* compute time-derivative couplings (TDC)\n"
            if ks_data[206] == "0":
                summary += "* rescale velocity after hop along the non-adiabatic coupling (NAC) vector\n"
            elif ks_data[206] == "2":
                summary += "* rescale velocity after hop along the gradient difference (GD) vector\n"
            elif ks_data[206] == "1":
                summary += "* simple rescaling of the velocity after hop\n"
            if ks_data[80] != "0":
                summary += "* use defined number in Tully FSSH, for testing purpose\n"
            if ks_data[193] == "1":
                summary += "* scale initial velocities for isotopes to match initial moment of normal atom (only if !isotopes is present)\n"
            if ks_data[199] == "1":
                summary += "* activate single state calculation after the hopping region\n"
            if ks_data[200] == "-1":
                summary += "* activate single state calculation after the hopping region\n"

        if ks_data['args']['type']  == "freqxg" or ks_data[8] == 1:
            summary += "\n"
            summary += "*** NUMERICAL DERIVATIVES ***\n"
            summary += "* displacement for finite difference computation of derivatives: {0} Ang\n".format(ks_data[12])
            if ks_data[10] == "0":
                summary += "* displace only in + direction for finite difference derivatives\n"
            elif ks_data[10] == "1":
                summary += "* displace in +/- direction for finite difference derivatives\n"
            elif ks_data[10] == "2":
                summary += "* displace in +/- direction only when hopping is active, in + direction otherwise\n"
            if ks_data[11] == "0":
                summary += "* optimize only CI coefficients during numerical computation of forces\n"
            elif ks_data[11] == "1":
                summary += "* optimize both MO & CI coefficients during numerical computation of forces\n"
            if ks_data[16] == "0":
                summary += "* do not perform numerical differentiation for link atoms\n"
            elif ks_data[16] == "1":
                summary += "* perform numerical differentiation for link atoms\n"

    # GAUSSIAN : no option
        #if ks_data['config']['QM']['exec'] == "1":
            #summary += "\n"
            #summary += "*** GAUSSIAN OPTIONS ***\n"

    # MOLPRO : 92, 93, 100, 194, 198
        if ks_data['config']['QM']['exec'] == "7":
            summary += "\n"
            summary += "*** MOLPRO OPTIONS ***\n"
            if ks_data[92] == "0":
                summary += "* do not perform smooth switching off of state-averaging\n"
            else:
                summary += "* smooth switching off of state-averaging at dE = {0} kcal/mol (only for 2 states)\n".format(ks_data[92])
            if ks_data[93] == "0":
                summary += "* do not perform abrupt switching off of state-averaging\n"
            else:
                summary += "* abrupt switching off of state-averaging at dE = {0} kcal/mol\n".format(ks_data[93])
            if ks_data[194] == "0":
                summary += "* use spherical harmonics to construct the basis set\n"
            elif ks_data[194] == "1":
                summary += "* use Cartesian functions to construct the basis set\n"
            if ks_data[198] == "0":
                summary += "* do not bomb the job when molpro WF does not converge\n"
            elif ks_data[198] == "1":
                summary += "* bomb the job when molpro WF does not converge\n"

    # MOLCAS : 93, 101, 102, 194, 195, 196, 197
        if ks_data['config']['QM']['exec'] == "6":
            summary += "\n"
            summary += "*** MOLCAS OPTIONS ***\n"
            if ks_data[93] == "0":
                summary += "* do not perform abrupt switching off of state-average\n"
            else:
                summary += "* switch off state-averaging at dE = {0} kcal/mol\n".format(ks_data[93])
            if ks_data[101] == "1":
                summary += "* switch on internal Molcas numerics\n"
            if ks_data[102] == "0":
                summary += "* allow only down-hop in CI vector rotation protocol (85 = 1)\n"
            else:
                summary += "* forbid back hopping for {0} number of steps\n".format(ks_data[102])
            if ks_data[194] == "0":
                summary += "* use spherical harmonics to construct the basis set\n"
            elif ks_data[194] == "1":
                summary += "* use Cartesian functions to construct the basis set\n"
            summary += "* threshold for discarding integrals when using Cholesky: {0}\n".format(ks_data[195])
            if ks_data[196] == "0":
                summary += "* do not use Cholesky decomposition\n"
            elif ks_data[196] == "1":
                summary += "* use Cholesky decomposition to speed up integrals\n"
            summary += "* basis set definition (possible overwritten by !basisset): {0}\n".format(ks_data[197])

        # TURBOMOLE : 190, 191
        if ks_data['config']['QM']['exec'] == "5":
            summary += "\n"
            summary += "*** TURBOMOLE OPTIONS ***\n"
            summary += "* the gradient will be computed for the {0} excited state\n".format(ks_data[190])
            if ks_data[191] == "0":
                summary += "* use grad and dscf for running TURBOMOLE\n"
            elif ks_data[191] == "1":
                summary += "* use RI and rigrad for running TURBOMOLE\n"
            elif ks_data[191] == "2":
                summary += "* use ricc2 for running TURBOMOLE\n"

        # SHARC INTERFACE
        if ks_data['config']['QM']['exec'] == "11":
            summary += "\n"
            summary += "*** SHARC INTERFACE OPTIONS ***\n"
            summary += "* using Sharc interface to compute QM using {0}\n".format(ks_data[110])
    except KeyError:
        pass 

    summary += "\n"

    return summary

#####################################################################################################
#            FUNCTIONS TO MANAGE COBRAMM COMMAND FILE
#####################################################################################################

def getCobrammCommand(filename):
    # write in the cobram.log the cobram.ks_data file
    filein = open(filename)
    x = filein.readlines()
    filein.close()
    return x

def makehard():
    # create the list of cobram ks_datas
    # if you set here a value it will be used as default in all calculations
    ks_data = []
    # initialize to 0 all value
    for i in range(300):
        ks_data.append(str(0))
    # insert default ks_data in cobram
    ks_data['args']['type']  = str('optxg')  # calc type (optxg, optxgp, ts, tsp, irc, ircp, ci, cip, freqxg, freqxgp, mdv, mdvp)
    ks_data[2] = str(0)  # verbosity of cobram.out and cobram.log (0=non verbose, 1=verbose, 2=very verbose)
    ks_data[3] = str(0)  # 0=do not freeze; 1=freeze all MEDIUM part; 2=freeze all HIGH part;
    # 3= active MEDIUM (for old G03 optimizer)
    ks_data['config']['MM']['border_atom_flag']  = str(0)  # 0=do not freeze; 1=freeze border atoms
    ks_data['config']['opt_mem']  = str('500MB')  # allocated memory for geometry optimization
    ks_data[6] = str('auto')  # if read or not  QM wavefunc (gaussian) from check;
    ks_data['config']['QM']['ncore']  = str(1)  # Nproc for parallel runs
    ks_data[8] = 0  # run with parallel finite difference: 0 = standard serial calculation, 1 = parallel finite diff
    ks_data['config']['ncore']  = str(1)  # Nproc for parallel frequency calculations
    ks_data[10] = str(1)  # 0: displace only in + direction for numerical computations (3*N);
    # 1: dispace in +/- direction for numerical computations (6*N);
    # 2: displace only in + direction when hopping is not active, dispace in +/- direction when hopping is active
    ks_data[11] = str(1)  # 0: optimize only CI coefficients during numerical compuation of the gradient and NACs;
    # 1: optimize MO and CI coeffcients (works only with MOLCAS so far)
    ks_data[12] = str(0.001)  # displacement (in Angstrom) for numerical computations (accuracy up to 0.0001)
    ks_data[13] = str(1)  # state to relax in numerical optxg or follow in numerical mdv (GS = 1, S1 = 2, etc.)
    ks_data[14] = str(0)  # type of derivative coupling 0: compute DC through spatial NACs;
    # 1: compute DC numerically through Tully-Hammes-Schaefer from WF at t+dt and t-dt (scaling along GD);
    ks_data[15] = str(0)  # 0: do not save single points during numerics
    # 1: save every single point during numerics (for debugging) to molcasALLnum.log
    ks_data[16] = str(1)  # 0: do not perform numerical differentiation for link atoms in numerics;
    # 1: perform numerical differentiation for link atoms in numerics
    ks_data[18] = str(0)  # 0: do not perform QM calculatrions for displacements of M-layer atoms in a frequency calculation (do MM calculation only for such steps) - 1: perform all QM calculations of the freq run (H+M atoms)
    ks_data[19] = str(0)  # branching plane definition in CI optimizations 0: use NACs 1: use gradient mean (according to eq.6 of J. Chem. Theory Comput., Vol. 6, No. 5, 2010)
    ks_data['config']['MM']['AMBER']['type'] = str(0)  # compute forces using the amber velocities (for big systems this is often
    # faster than using the implemented forcedump routine)
    # 1 = use pmemd GPU code with a dielektrikum (carefully compare your results!)
    # 2 = use sander.MPI parallel code
    # 3 = use serial amber but compute velocities
    ks_data[21] = str(1.090)  # CH link distance
    ks_data[22] = str(0.947)  # OH link distance
    ks_data[23] = str(1.008)  # NH link distance
    ks_data[24] = str(1.526)  # CC equillibrium distance
    ks_data[25] = str(1.410)  # CO equillibrium distance
    ks_data[26] = str(1.475)  # CN equillibrium distance
    ks_data[40] = str(0)  # frequency below which the normal modes are used to construct the sub-space the project the gradient; 
    ks_data[41] = str(999) #cut off value for MM calculations
    #used for Wigner sampling through QM/MM dynamics for low frequency modes 
    ks_data['config']['QM']['exec'] = str(1)  # QM calculation type: 0=none, 1=Gaussian, 2=orca,3=DFTB, 4=DFTB+, 5=TURBOMOLE, 6=molcas,
    # 7=molpro, 8=DFT-MRCI, 11=SHARC QM interfaces
    ks_data['config']['QM']['mem']  = str('700MB')  # Memory Definition for Gaussian OPT and QM part
    ks_data[60] = str(100)  # number of optimization cicles (in case of IRC equal to "number of points along the
    # reaction path" x "number of cycles required for each step"
    ks_data['config']['QMMM']['border_atom_force']  = str(1)  # correct the gradient of Q1 (1 = yes)
    ks_data[65] = str(10)  # number of points along the reaction path for IRC computations
    ks_data[66] = str(10)  # step size along the reaction path; .gt 0: in units of 0.01 Bohr;
    # .lt 0: in units of 0.01 amu1/2 Bohr (beware that the convergence threshold also depends on the stepsize:
    # threshold == 30 x step; this can be modified with ks_data[67])
    ks_data[67] = str(300)  # convergence threshold of the first derivative: RMS=1.5x, displ=4x, max disp=6x
    ks_data[68] = str(30)  # maximum size for an optimization step in units of 0.01 Bohr (default is 30, i.e. 0.30 Bohr)
    ks_data[71] = str(0)  # MM calculation type: 0=amber, no other option at the moment
    ks_data[80] = str(0)  # not so random number (for testing); 0: use random number, .ne 0: user defined random number
    ks_data[81] = str(0)  # in case of hopping scheme based on NACs, highest state to include in amplitude propagation for dynamics; 0: include all states in propagation, .ne 0: states above this number are spectators (may be needed for stability reasons but they do not participate in pop.dynamics, numbering: S1 = 2, etc.))
    ks_data[83] = str(1.0)  # first time step in the MD
    ks_data[84] = str(0.25)  # second timestep
    ks_data[85] = str(0)  # activate surface hopping 0= no, 1=Tully surface hopping with Persico decoherence, 2=hop based on energy difference
    ks_data[86] = str(1000.0)  # energy treshold to activate the surface hopping
    ks_data[87] = str(1)  # active the back surface hopping 0= no 1= yes
    ks_data[88] = str(0)  # activate surface hopping 0= no 1= charge amount
    ks_data[89] = str(15.0)  # tresh. for the S0/S1 energy gap to deactivate the charge surface hopping
    ks_data[90] = str(0)  # type of ES-GS hop with TDDFT (0 = no hop, 1 = use GS-ES tdc, 2 = hop when deltaE<threshold)
    ks_data[91] = str(2.0)  # threshold for ES->GS hop when ks_data[90]=="2" in kcal/mol
    ks_data[92] = str(0)  # OW switch off weighting smoothly.  n=0: off, n>0 start to switch off at dE=n kcal/mol
    ks_data[93] = str(0)  # OW switch off weighting abruptly.  n=0: off, n>0 switch off at dE=n kcal/mol
    ks_data[94] = str(0)  # OW intruder state detection  n=0: off, n>0 mix intruder
    # state from dE=40kcal/mol using sigmoid function, scaled by n
    ks_data[95] = ''  # OW scale sigmoid exponent by n (defaul=-0.4)
    ks_data[96] = ''  # OW scale half height position of sigmoid (kcal/mol, default=18)
    ks_data[97] = str(0)  # type of TD coupling for TD-DFT: Izmaylov (0) Hammes-Schiffer Tully (1) or MO unitary transformation (2)  
    ks_data[98] = str(1.0) # threshold for truncating WF expansion in WF overlap calculation in HST 
    ks_data['config']['QM']['store_log_step']  = str(1)  # save log file every X steps for all the QM calculation types
    ks_data['config']['QM']['store_wfn_step']  = str(-1)  # save wavefunction files every X steps for all the QM calculation types
    ks_data[101] = str(0)  # 0: switch off internal Molcas numerics (CASSCF, CASPT2);
    # 1: switch on internal Molcas numerics (CASSCF, CASPT2)
    ks_data[102] = str(0)  # 0: allow only down-hop in CI vector rotation protocol (85 = 1);
    ks_data[103] = str(1)  # 0: compute only GS after hopping to the GS (DFT); 1: always compute TDDFT after hopping to the GS
    # ne. 0: forbid back hopping for n steps
    ks_data[110] = ''  # name of the QM code used through the SHARC interface
    ks_data[120] = str(0)  # OW Use molcas GHOST atoms instead of XField when "1" (Use this for molcas RI!)
    ks_data[122] = str('0')  # criterion stopping trajectory computation, when the system is in the GS for X fs
    ks_data[130] = str(0)  # whether to compute gradient of HighLayer: 0 yes (normal run), 1 freeze HighLayer (compute only MediumLayer gradient based on electric field)
    ks_data[190] = str(0)  # OW select excited state gradient (egrad) for Turbomole (default '0', use grad)
    ks_data[191] = str(0)  # Use RI and rigrad for Turbomole if '1', or ricc2 if '2' (default '0', use grad and dscf)
    ks_data[192] = str(0)  # Disable NOSYM warning and continue
    ks_data[193] = str(0)  # scale initial vel for isotopes to match initial moment of normal atom (only for isotopes!)
    ks_data[194] = str(0)  # OW Use Cartesian functions in molcas if '1', default is Spherical ('0')
    ks_data[195] = str('1.0E-4')  # OW Threshold for RICD in molcas (CDTH) Defaults to 1.0E-4
    ks_data[196] = str(0)  # OW use RI approximation in molcas when '1'
    ks_data[197] = str('6-31Gp')  # OW Basis set used for molcas
    ks_data[198] = str(1)  # OW 0: Do not bomb the job when molpro WF does not converge (Be careful!);
    # 1: bomb the job when molpro WF does not converge
    ks_data[199] = str(0)  # 0: do nothing; .ne 0: use !molcasSS in the GS when the 0-1 gap is above threshold
    ks_data[200] = str(0)  # use single point on top of a SS traj every n steps only for molpro/molcas
    ks_data[201] = str(0)  # 0=Massey, Tully or Tully corrected with NAC scheme 1,
    # 1=Massey, Tully or Tully corrected with NAC scheme 2
    ks_data[203] = str(0)  # 0: do nothing; .ne 0: switch to MP2 in GS when the 0-1 gap is above threshold
    ks_data[204] = str(2.5)  # threshold (kcal/mol) for including states beyond the photoactive state during
    # the computation of the numerical gradients at SS-PT2 and PM-CASPT2 level to account for possible state swapping
    ks_data[205] = str(0)  # 0: fixed number of states in computation of the numerical gradients at
    # SS-PT2 and PM-CASPT2 level; 1: include only state below threshold 204
    ks_data[206] = str(0)  # correct velocity after hop: (0) along the DC (NAC, default for all types of dyn);
    # (1) rescale velocity; (2) along GD vector
    ks_data[207] = str(0) # threshold for ES->ES hop when ks_data[85]=="2" in kcal/mol
    ks_data.pop(0)
    return ks_data


def ReadCobramCommand(cobcom, getpart, filename2):
    # check the input files presence in the working directory
    # if not take it from cobram.commnd
    if os.path.isfile(filename2):
        # get file from working directory
        filein = open(filename2)
        x = filein.readlines()
        filein.close()
    else:
        # read a part of the cobram.ks_data file
        # get only the part defined in the variable getpart
        x = []
        start, end = 0, 0
        for i in range(len(cobcom)):
            row = cobcom[i].strip()
            if (row == '!' + getpart) or (row == '!' + getpart + 'p'):
                start = i
            if (row == '?' + getpart) or (row == '?' + getpart + 'p'):
                end = i
        try:
            for i in range(start + 1, end):
                x.append(cobcom[i].strip())
        except:
            pass
    return x


def key2hard(key_soft, ks_datahard):
    # parse the !keywords section and then merge it with
    # the ks_data_hard (the default keyword in of the code)
    # beware that !ks_datas will always overwrite !keywords
    #
    # first create a dictionary of keywords and IOPs
    key_dict = {'type': 1, 'verbose': 2, 'freeze': 3, 'border': 4, 'geomem': 5, 'nproc': 7, 'numproc': 9, 'distype': 10, 
                'cicopt': 11, 'displace': 12, 'numrlx': 13, 'DC': 14, 'savnum': 15, 'difflink': 16, 'skipQM' : 18, 'BranchPlane': 19,'amber': 20, 'cutoff': 41, 
                'savchk': 50, 'qm-type': 51,
                'qmem': 53, 'nsteps': 60, 'q1cor': 61, 'ircpoints': 65, 'ircstepsize': 66, 'conver': 67, 'stepsize': 68,
                'mm-type': 71, 'rand': 80, 'dynroots': 81, 'tstep': 83, 'tsshort': 84, 'surhop': 85, 'ediff': 86, 'backhop': 87,
                'hoptogs': 90, 'gsesth': 91,
                'smooth': 92, 'brute': 93, 'intruder': 94, 'sigm': 95, 'height': 96, 'tdctype': 97, 'thresOv': 98, 'savQMlog': 99,
                'savwfu': 100, 'molcasnum': 101,
                'cirotbackhop': 102, 'tdafterhop': 103, 'sharc-qm': 110, 'ghost': 120, 'enefail': 121, 'stop': 122, 'egrad': 190,
                'rigrad': 191, 'scale-iso': 193,
                'basisfunc': 194, 'cdth': 195, 'ricd': 196, 'basis': 197, 'bomb': 198, 'single': 199, 'sstep': 200,
                'velafterhop': 206,
                'delta_en': 207, 'ta' : 208
                }

    # now special dict's for key-options
    key_freeze = {'no': 0, 'medium': 1, 'high': 2, 'active': 3}
    key_distype = {'plus': 0, 'both': 1, 'phop': 2}
    key_cicopt = {'ci': 0, 'cimo': 1}
    key_DC = {'nac': 0, 'tdc': 1}
    key_qmtype = {'none': 0, 'gauss': 1, 'orca': 2, 'dftb': 3, 'dftb+': 4, 'turbo': 5, 'molcas': 6, 'molpro': 7,
                  'dft-mrci': 8, 'dalton': 9,
                  'sharc': 11,
                  'mndo':10086, 'bagel':10087} #
    key_surhop = {'none': 0, 'tully': 1, 'deltaE': 2}
    key_amber = {'serial': 0, 'cuda': 1, 'mpi': 2, 'velo': 3}
    key_mmtype = {'amber': 0, 'gromac':1, 'openmm':2, 'lammps':3} #
    key_basisfunc = {'spher': 0, 'cart': 1}
    key_BranchPlane = {'nac': 0, 'gmean': 1}
    key_velafterhop = {'nac': 0, 'uniform': 1, 'GD': 2}
    key_other = {'yes': 1, 'true': 1, 'on': 1, 'active': 1, 'no': 0, 'false': 0, 'off': 0, 'inactive': 0}
    row = []
    comm = []
    key = []
    value = []
    for i in range(len(key_soft)):
        tmp = key_soft[i].split()
        for j in range(len(tmp)):
            tmp2 = tmp[j].split('=')
            # print tmp2
            key.append(tmp2[0])
            # print "appending key:",tmp2[0]
            value.append(tmp2[1])
            # print "with value",tmp2[1]
    # now resolve the keywords into iops
    for i in range(len(key)):
        try:
            row.append(key_dict[key[i]])
        except:
            fatalError("Key {} not recognized\n".format(key[i]))

        if key[i] == 'freeze':
            comm.append(str(key_freeze[value[i]]))
        elif key[i] == 'distype':
            comm.append(str(key_distype[value[i]]))
        elif key[i] == 'cicopt':
            comm.append(str(key_cicopt[value[i]]))
        elif key[i] == 'DC':
            comm.append(str(key_DC[value[i]]))
        elif key[i] == 'qm-type':
            comm.append(str(key_qmtype[value[i]]))
        elif key[i] == 'surhop':
            comm.append(str(key_surhop[value[i]]))
        elif key[i] == 'amber':
            comm.append(str(key_amber[value[i]]))
        elif key[i] == 'mm-type':
            comm.append(str(key_mmtype[value[i]]))
        elif key[i] == 'basisfunc':
            comm.append(str(key_basisfunc[value[i]]))
        elif key[i] == 'BranchPlane':
            comm.append(str(key_BranchPlane[value[i]]))
        elif key[i] == 'velafterhop':
            comm.append(str(key_velafterhop[value[i]]))
        elif value[i] in key_other:
            comm.append(str(key_other[value[i]]))
        else:
            comm.append(str(value[i]))
    ks_data_soft = [row, comm]
    # print "ks_data_soft=",ks_data_soft
    # lists[6].insert(0,'0')
    for i in range(1, len(ks_datahard)):
        # print "checking element", lists[6][i]
        # compare and substitute hard and soft keywords
        for j in range(len(ks_data_soft[0])):
            if i == ks_data_soft[0][j]:
                # if (lists[6][i].strip()!=ks_data_soft[1][j].strip()):
                #    x.append('default ks_data value:\t'+str(i)+' '+lists[6][i].strip()+
                #             '\t new value:\t'+str(i)+' '+ks_data_soft[1][j]+'\n')
                ks_datahard[i] = ks_data_soft[1][j]
    return ks_datahard


def soft2hard(ks_data_soft, ks_datahard):
    # merge the commad_soft(users cobram keywords) and
    # the commad_hard (the default keyword in of the code)
    # initialization of the lists
    row = []
    comm = []
    for i in range(len(ks_data_soft)):
        # clean up the cobram_soft from user comments
        element = ks_data_soft[i].split()
        row.append(int(element[0]))
        comm.append(element[1])
    ks_data_soft = [row, comm]
    # print "ks_data_soft in soft2hard",ks_data_soft
    # lists[6].insert(0,'0')
    for i in range(1, len(ks_datahard)):
        # compare and substitute hard and soft keywords
        for j in range(len(ks_data_soft[0])):
            if i == ks_data_soft[0][j]:
                # if (lists[6][i].strip()!=ks_data_soft[1][j].strip()):
                #    x.append('default ks_data value:\t'+str(i)+' '+lists[6][i].strip()+
                #             '\t new value:\t'+str(i)+' '+ks_data_soft[1][j]+'\n')
                ks_datahard[i] = ks_data_soft[1][j]
    return ks_datahard


#####################################################################################################
#            OLD INTERFACES TO QM AND MM
#####################################################################################################

# "par_num" controls the behavior of Cobram in the QM routine
# 0: sequential computation
# 1: computation at the reference geometry during parallel numerics (behaves like a sequential computation)
# 2: collecting data for computing FREQs, GRADs and NACs during parallel numerics
# 3: perform surface hopping
def QM(ks_data, cobcom, charges, geometry, step, par_num):
    # inizialize QM_results with zeros
    QM_results = 0.0, [], [], 0.0, [0.0, 0.0, 0.0, 0.0]

    if geometry.calculationType in ['M', 'ML']:
        logwrt.writeLog('No QM calculation will be done, it is a ' + geometry.calculationType + ' calculation\n')
        return QM_results

    if par_num == 0 or step == 0: logwrt.writeLog('Starting QM calculation ')

    if step != 0 or int(ks_data[211]) > 0:
        sef = shelve.open("cobram-sef", 'r')
        calc_coupl = sef['calc_coupl']
        sef.close()
    else:
        calc_coupl = []

    # ---------------------------------
    if ks_data['config']['QM']['exec'] == '6':  # QM calculation by MOLCAS

        if par_num == 0 or par_num == 1:
            molcas.prepare(ks_data, step)
            molcas.makeinp(charges, geometry, ks_data, step, cobcom, par_num)

        # if par_num == 1 for a freq calc, and this is not the first step, only prepare input but do not run Molcas
        # parallel Molcas instances will be run later in the routine QM2
        if ks_data['args']['type']  == 'freqxg' and par_num == 1 and step != 0:
            QM_results = []

        # if par_num == 1 for another type of parallel calc, only prepare input but do not run Molcas
        # parallel Molcas instances will be run later in the routine QM2
        elif ks_data['args']['type']  in ['optxgp', 'mdvp', 'ircp', 'cip', 'tsp'] and par_num == 1:
            QM_results = []

        # if par_num == 2 and this is a frequency calc, collect all the information necessary to
        # compute numerical frequencies (via Gaussian): ENERGY, DIPOLE and GRAD
        elif ks_data['args']['type']  == 'freqxg' and par_num == 2:
            # extract results from calculations
            QM_results = molcas.molcasEneGradCrg(ks_data, geometry, charges, step, par_num)

        # if par_num == 2 and this is another type of calculationa,
        # collect all the information necessary to compute numerical GRADs: ENERGY
        elif ks_data['args']['type']  in ['optxg', 'mdv', 'irc', 'ci', 'ts'] and par_num == 2:
            # compute gradient
            molcasGradient = molcas.molcasEne(ks_data, geometry)
            # compute NAC when needed
            molcasNAC = []
            if (ks_data['args']['type']  == 'mdv' and ks_data[85] > '3' and calc_coupl != [] and ks_data[14] == '0') or \
                    ks_data['args']['type']  == 'ci':
                molcasNAC = molcas.molcasNAC(ks_data, geometry)
            # QM_Results returned at this point is not complete as it contains only GRAD (and NACs), while all the
            # remaining data has already been collected at the reference geometry... therefore, the QM_Results is
            # returned to parallel_numerics.run as a temporary tmp_Results, serving to update the GRAD and NACs
            QM_results = [molcasGradient, molcasNAC]

        # the molcas.molcasMDV routine is called once the GRADs and NACs have been computed
        elif ks_data['args']['type']  == 'mdv' and par_num == 3:
            with Timer("molcasMDV1"):
                molcasResults = molcas.molcasMDV(ks_data, geometry, charges, cobcom, step, par_num)
            # the QM_Results is returned to parallel_numerics.run as a temporary tmp_Results
            # in case of a hop the QM_Results are updated with data for the new state (done in molcas.molcasMDV)
            # if there is no hop the returned QM_results is empty and QM_Results is not overwritten
            QM_results = molcasResults

        else:
            if os.path.exists(os.getcwd() + '/molcas.log') and step == 0 and ks_data['args']['type']  in ['freqxgp', 'freqxg']:
                pass
            else:
                # launch MOLCAS calculation
                molcas.launch(ks_data, step)
            # extract results
            molcasResults = molcas.molcasEneGradCrg(ks_data, geometry, charges, step, par_num)
            molcas.singlepoint(ks_data, step, cobcom)
            sef = shelve.open("cobram-sef", 'r')
            state = sef['state']
            newstate = sef['newstate']
            sef.close()
            # if MCLR fails we skip reading all infomation from this step and go directly to Vverlet
            # in case of a HOP molcasMDV will be repeated only when couplings are to be computed with Molpro (14 2)
            if ks_data['args']['type']  == 'mdv' and par_num == 0 and (
                    newstate == state or (newstate != state and ks_data[14] == '2')):
                with Timer("molcasMDV1"):
                    tmp_Results = molcas.molcasMDV(ks_data, geometry, charges, cobcom, step, par_num)
                # the QM_Results is returned as a temporary tmp_Results
                # in case of a hop the QM_Results are updated with data for the new state (done in molcas.molcasMDV)
                # if there is no hop the returned QM_results is empty and QM_Results is not overwritten
                if tmp_Results != 0:
                    molcasResults = tmp_Results
            QM_results = molcasResults

    # ---------------------------------
    elif ks_data['config']['QM']['exec'] == '7':  # QM calculation by MOLPRO

        molpro.prepare(ks_data, step)
        chMED = molpro.makeinp(ks_data, geometry, charges.CRG_emb, step)
        molpro.launch(ks_data, step)
        QM_results = molpro.molproEneGradCrg(ks_data, geometry, charges.CRG_emb, step, chMED, cobcom)

    if par_num == 0 or par_num == 2 or step == 0: logwrt.writeLog("\n")

    return QM_results


# =============================================================================================================


def save_step_temp(step, geometry, ks_data, copyLog, copyOrb):
    """ at the end of each iteration in the main loop, collect and save files of the MM and QM calculations """

    # the calculation has a QM part, collect files with orbitals
    if "H" in geometry.calculationType:

        if ks_data['config']['QM']['exec'] == '6':  # QM calulation by Molcas
            molcas.save_QM_molcas_step(step, ks_data, copyLog, copyOrb)
        if ks_data['config']['QM']['exec'] == '7':  # QM calulation by MOLPRO
            molpro.save_QM_molpro_step(step, ks_data, copyLog, copyOrb)


# =============================================================================================================


def MM(step, geometry, ks_data, cobcom):
    """ wrapper function for the MM part """

    # initialize MMresults with a series of dummy arguments
    MMresults = [], 0.0, 0.0, [], [], 0.0

    # when the scheme actually contains a MM part
    if "M" in geometry.calculationType or "L" in geometry.calculationType:

        logwrt.writeLog('Starting MM calculation ')

        # AMBER calculation (only option implemented now)
        if ks_data[71] == '0':
            logwrt.writeLog('using AMBER ... \n')
            MMresults = amber.MMcalculations(step, geometry, ks_data, cobcom)
        else:
            fatalError("MM software (ks_data 71) '{0}' is not defined. ".format(ks_data[71]))
        logwrt.writeLog('\n')

    return MMresults

if __name__ == '__main__':
    args = kids_arg.parser.parse_args()
    if args.type == kids_arg.CalculationType.NAD.value:
        ks_data = KSInput(args.input)
        pprint(ks_data)
