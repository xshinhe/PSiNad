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

from dataclasses import dataclass, field
from typing import Dict, List, Any, Tuple

@dataclass
class MMOutput:
    charge_real: List[float]
    charge_modelH: List[float]

    energy_real: float
    energy_real_nocharge: float
    energy_modelH: float

    grad_real: List[Any]
    grad_real_nocharge: List[Any]
    grad_modelH: List[Any]


@dataclass
class QMOutput0:
    inpfile: str = None
    calcdir: str = None
    outfile: str = None

    status: int = None
    errormsg: List[str] = field(default_factory=list)
    
    natoms: int = None
    nroots: int = None

    iroot0: int = None
    iroot1: int = None

    charges: List[float] = field(default_factory=list)
    elfield: Dict[int, Any] = field(default_factory=dict)
    
    dipole: List[float] = field(default_factory=list)
    osc_strength: Dict[str, Any] = field(default_factory=dict)
    
    energy: Dict[int, float] = field(default_factory=dict)
    gradient: Dict[int, Any] = field(default_factory=dict)
    hess: Dict[int, Any] = field(default_factory=dict)
    nac: Dict[Tuple[int, int], Any] = field(default_factory=dict)
    selfenergy: float = None
    
    gradcharges_extraterm: bool = False
    fullgradcharge: Dict[int, Any] = field(default_factory=dict)
    fullnaccharge: Dict[Tuple[int, int], Any] = field(default_factory=dict)
    
    signs: List[Any] = field(default_factory=list)
    scf_mo_maps: Any = None
    psioverlap: Any = None
    eigenvectors: Any = None

class QMOutput:
    def __init__(self):

        # Setup a string to store the log of the QM output processing
        self.log = ""

        # Setup a dictionary to store the info of the output file
        self.dataDict = {
            # results of the calculation
            "energy": {},  # dictionary of QM energies for the state considered in the Gaussian calculation
            "gradient": {},  # dictionary of the gradients computed with QM

            "optstate": None,  # optimized electronic state
            "charges": [],  # charges of the modelH calculation, to be used to redefine the modelH force field
            "dipole": None,  # dipole moment of the modelH calculation, printed as relevant output information

            "selfenergy": 0.0,  # self energy of the background charges
            "gradcharges": None,  # electric field at the position of the point charges
            # logical variable that triggers the removal of the M-ML charge interaction from the external charge force
            "gradcharges_extraterm": False,

            "cis_coeffs": None,  # dictionary to store the CIS/TD excitation coefficients, initialized to None (will remain None for Molcas)
            "psi_overlap": None # F: overlap matrix between the wavefunctions, initialized to None (will remain None for Gaussian)
        }

    # =============================================================================================================

    def set(self, label, value):
        """ store value in the dictionary with index given by label"""
        self.dataDict[label] = value

    def get(self, label):
        """ get value in the dictionary with index given by label"""
        try:
            return self.dataDict[label]
        except KeyError:
            return None

    # =============================================================================================================

    def __del__(self):
        """Destructor for the QMOutput class: deallocate memory for the object attribute"""

        # clean up the log
        del self.log
        # destroy the dictionary that contains the output data
        del self.dataDict

    # =============================================================================================================

    def orbfiles(self):
        """Return a list with the names of the orbital files. In the case of this QMOutput class that
        represent a generic QM output, return an empty list. """

        return []

    # =============================================================================================================

    def restartfile(self):
        """Return a the names of the file needed for orbital restart. In the case of this
        dummy QMOutput class that represent a generic QM output, return None. """

        return None

###################################################################################################################
