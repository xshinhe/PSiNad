#! /usr/bin/env python3
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

# imports of local modules

import logwrt

#####################################################################################################
# get pretty output strings of the command options
#####################################################################################################

def getCalcType( command ):
   """ get a string that describes the the type of calculation that is requested"""

   if command[1] == "optxg" and command[60] == '1':
      calcType = "A single point calculation"
   elif command[1] == "optxg" and command[60] != '1':
      calcType = "An optimization run"
   elif command[1] == "ts":
      calcType = "A transition state (TS) search"
   elif command[1] == "irc":
      calcType = "An intrinsic reaction coordinate (IRC) calculation"
   elif command[1] == "ci":
      calcType = "A conical intersection (CI) search"
   elif command[1] == "freqxg":
      calcType = "A harmonic frequency calculation"
   elif command[1] == "mdv":
      calcType = "A molecular dynamics run"
   elif command[1] == "nad":
      calcType = "An interface for nad"
   else:
      logwrt.fatalerror( "calculation type (command 1) '{0}' is not defined. ".format(command[1]) )

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
      logwrt.fatalerror( "layer scheme '{0}' is not defined. ".format(geometry.calculationType) )

   return layers

#================================================================================

def getQMCode( command ):
   """ get a string that describes the QM code that is used for the calculation"""

   if command[51] == "0":
      QMsoftware = None
   elif command[51] == "1" or command[51]=='10':
      QMsoftware = "Gaussian"
   elif command[51] == "2":
      QMsoftware = "ORCA"
   elif command[51] == "3":
      QMsoftware = "DFTB"
   elif command[51] == "4":
      QMsoftware = "DFTB+"
   elif command[51] == "5":
      QMsoftware = "TURBOMOLE"
   elif command[51] == "6":
      QMsoftware = "MOLCAS / OPENMOLCAS"
   elif command[51] == "7":
      QMsoftware = "MOLPRO"
   elif command[51] == "8":
      QMsoftware = "TURBOMOLE-DFT-MRCI"
   elif command[51] == "11":
       QMsoftware = "SHARC QM interfaces"
   elif command[51] == "10086":
       QMsoftware = "MNDO QM interfaces"
   elif command[51] == "10087":
       QMsoftware = "Bagel QM interfaces"
   else:
      logwrt.fatalerror( "QM software (command 51) '{0}' is not defined. ".format(command[51]) )

   return QMsoftware

#================================================================================

def getMMCode( command ):
   """ get a string that describes the QM code that is used for the calculation"""

   if command[71] == "0":
      MMsoftware = "Amber"
   else:
      logwrt.fatalerror( "MM software (command 71) '{0}' is not defined. ".format(command[71]) )

   return MMsoftware

#================================================================================

def getOptionSummary( command ):
    """ returns a string with a summary of the options defined in the input command list """

    summary = ""
    summary += "*** SYSTEM OPTIONS ***\n"
    summary += "* allocated memory for Gaussian optimizer: {0}\n".format(command[5])
    summary += "* allocated memory for QM computation: {0}\n".format(command[53])
    summary += "* nr. of cores to use for QM third-party software: {0}\n".format(command[7])
    summary += "* nr. of cores to use for COBRAMM internal parallel numerics routine: {0}\n".format(command[9])
    if command[20] == "0":
        summary += "* use AMBER serial executable sander\n"
    elif command[20] == "1":
        summary += "* use AMBER GPU code pmemd.cuda\n"
    elif command[20] == "2":
        summary += "* use AMBER parallel executable sander.MPI\n"

    summary += "\n"
    summary += "*** GENERAL OPTIONS ***\n"
    summary += "* number of optimization steps / microsteps in each IRC step / MD steps: {0}\n".format(command[60])
    if int(command[200]) > 0:
        summary += "* perform a single point computation every {0} steps (for Molpro and Molcas, requires !singlepoint  block)\n".format(command[200])
    elif command[200] == "-1":
        summary += "* perform a single point computation at every converged IRC geometry (for Molpro and Molcas, requires !singlepoint block)\n"
    if command[99] == "0":
        summary += "* do not store QM log files \n"
    else:
        summary += "* save QM log files every {0} steps\n".format(command[99])
    if command[100] == "0":
        summary += "* do not store QM wavefunction files \n"
    else:
        summary += "* save QM wavefunction files every {0} steps\n".format(command[100])

    summary += "\n"
    summary += "*** MOLECULAR MODEL OPTIONS ***\n"
    if command[4] == "1":
        summary += "* border atoms are frozen\n"
    if command[61] == "1":
        summary += "* project the force of the link atoms on the border atoms\n"

    if command[1] in ["optxg", "ts", "irc", "ci"] and int(command[60]) > 1:
        summary += "\n"
        summary += "*** OPTIMIZATION OPTIONS ***\n"
        summary += "* Gaussian optimization convergence threshold: RMS Force = {0}\n".format(float(command[67])*0.0001)
        summary += "* Gaussian optimization convergence threshold: Max Force = {0}\n".format(float(command[67])*0.0001*1.5)
        summary += "* Gaussian optimization convergence threshold: RMS Displ = {0}\n".format(float(command[67])*0.0001*4.)
        summary += "* Gaussian optimization convergence threshold: Max Displ = {0}\n".format(float(command[67])*0.0001*6.)
        summary += "* maximum step size of an optimization step: {0} bohr\n".format(float(command[68])*0.01)

    if command[1] == "irc":
        summary += "\n"
        summary += "*** IRC OPTIONS ***\n"
        summary += "* number of intermediate geometries along the MEP computed with IRC: {0} \n".format(command[65])
        summary += "* step size along the MEP: {0} bohr\n".format(float(command[66])*0.01)

    if command[1] == "mdv":
        summary += "\n"
        summary += "*** MOLECULAR DYNAMICS OPTIONS ***\n"
        summary += "* regular time-step of the molecular dynamics: {0} fs\n".format(command[83])
        if int(command[85]) > 0  and command[14] == '0': 
            summary += "* short time-step when Tully or Tully-Persico is active: {0} fs\n".format(command[84])
            summary += "* energy difference between active state and near-by states for activating short time step: {0} kcal/mol\n".format(command[86])
        if command[85] == "0":
            summary += "* perform adiabatic dynamics\n"
        elif command[85] == "1":
            summary += "* dynamics with hoppings according to Tully's FSSH, with Persico's decoherence correction\n"
        elif command[85] == "2":
            summary += "* dynamics with hoppings according on only energy difference\n"
            summary += "* energy threshold (kcal/mol) to hop: {0}\n".format(command[207])
        if command[85] == '1':
            summary += "* energy threshold (kcal/mol) to activate smaller timestep and hoppings: {0} kcal/mol\n".format(command[86])
        if command[90] == "0":
            summary += "* turn off hopping to the GS\n"
        elif command[90] == "0":
            summary += "* compute hopping probability to the GS in the same way as between ES\n"
        elif command[90] == "2":
            summary += "* force hop to the ground state when the GS-ES energy gap is < {0} kcal/mol\n".format(command[91])
        if command[87] == "0":
            summary += "* switch off back-hopping during Tully's FSSH\n"
        elif command[87] == "1":
            summary += "* switch on back-hopping during Tully's FSSH\n"
        if command[14] == "0":
            summary += "* compute (spatial) non-adiabatic couplings (NACs)\n"
        elif command[14] == "1":
            summary += "* compute time-derivative couplings (TDC)\n"
        if command[206] == "0":
            summary += "* rescale velocity after hop along the non-adiabatic coupling (NAC) vector\n"
        elif command[206] == "2":
            summary += "* rescale velocity after hop along the gradient difference (GD) vector\n"
        elif command[206] == "1":
            summary += "* simple rescaling of the velocity after hop\n"
        if command[80] != "0":
            summary += "* use defined number in Tully FSSH, for testing purpose\n"
        if command[193] == "1":
            summary += "* scale initial velocities for isotopes to match initial moment of normal atom (only if !isotopes is present)\n"
        if command[199] == "1":
            summary += "* activate single state calculation after the hopping region\n"
        if command[200] == "-1":
            summary += "* activate single state calculation after the hopping region\n"

    if command[1] == "freqxg" or command[8] == 1:
        summary += "\n"
        summary += "*** NUMERICAL DERIVATIVES ***\n"
        summary += "* displacement for finite difference computation of derivatives: {0} Ang\n".format(command[12])
        if command[10] == "0":
            summary += "* displace only in + direction for finite difference derivatives\n"
        elif command[10] == "1":
            summary += "* displace in +/- direction for finite difference derivatives\n"
        elif command[10] == "2":
            summary += "* displace in +/- direction only when hopping is active, in + direction otherwise\n"
        if command[11] == "0":
            summary += "* optimize only CI coefficients during numerical computation of forces\n"
        elif command[11] == "1":
            summary += "* optimize both MO & CI coefficients during numerical computation of forces\n"
        if command[16] == "0":
            summary += "* do not perform numerical differentiation for link atoms\n"
        elif command[16] == "1":
            summary += "* perform numerical differentiation for link atoms\n"

# GAUSSIAN : no option
    #if command[51] == "1":
        #summary += "\n"
        #summary += "*** GAUSSIAN OPTIONS ***\n"

# MOLPRO : 92, 93, 100, 194, 198
    if command[51] == "7":
        summary += "\n"
        summary += "*** MOLPRO OPTIONS ***\n"
        if command[92] == "0":
            summary += "* do not perform smooth switching off of state-averaging\n"
        else:
            summary += "* smooth switching off of state-averaging at dE = {0} kcal/mol (only for 2 states)\n".format(command[92])
        if command[93] == "0":
            summary += "* do not perform abrupt switching off of state-averaging\n"
        else:
            summary += "* abrupt switching off of state-averaging at dE = {0} kcal/mol\n".format(command[93])
        if command[194] == "0":
            summary += "* use spherical harmonics to construct the basis set\n"
        elif command[194] == "1":
            summary += "* use Cartesian functions to construct the basis set\n"
        if command[198] == "0":
            summary += "* do not bomb the job when molpro WF does not converge\n"
        elif command[198] == "1":
            summary += "* bomb the job when molpro WF does not converge\n"

# MOLCAS : 93, 101, 102, 194, 195, 196, 197
    if command[51] == "6":
        summary += "\n"
        summary += "*** MOLCAS OPTIONS ***\n"
        if command[93] == "0":
            summary += "* do not perform abrupt switching off of state-average\n"
        else:
            summary += "* switch off state-averaging at dE = {0} kcal/mol\n".format(command[93])
        if command[101] == "1":
            summary += "* switch on internal Molcas numerics\n"
        if command[102] == "0":
            summary += "* allow only down-hop in CI vector rotation protocol (85 = 1)\n"
        else:
            summary += "* forbid back hopping for {0} number of steps\n".format(command[102])
        if command[194] == "0":
            summary += "* use spherical harmonics to construct the basis set\n"
        elif command[194] == "1":
            summary += "* use Cartesian functions to construct the basis set\n"
        summary += "* threshold for discarding integrals when using Cholesky: {0}\n".format(command[195])
        if command[196] == "0":
            summary += "* do not use Cholesky decomposition\n"
        elif command[196] == "1":
            summary += "* use Cholesky decomposition to speed up integrals\n"
        summary += "* basis set definition (possible overwritten by !basisset): {0}\n".format(command[197])

    # TURBOMOLE : 190, 191
    if command[51] == "5":
        summary += "\n"
        summary += "*** TURBOMOLE OPTIONS ***\n"
        summary += "* the gradient will be computed for the {0} excited state\n".format(command[190])
        if command[191] == "0":
            summary += "* use grad and dscf for running TURBOMOLE\n"
        elif command[191] == "1":
            summary += "* use RI and rigrad for running TURBOMOLE\n"
        elif command[191] == "2":
            summary += "* use ricc2 for running TURBOMOLE\n"

    # SHARC INTERFACE
    if command[51] == "11":
        summary += "\n"
        summary += "*** SHARC INTERFACE OPTIONS ***\n"
        summary += "* using Sharc interface to compute QM using {0}\n".format(command[110])

    summary += "\n"

    return summary

#####################################################################################################
