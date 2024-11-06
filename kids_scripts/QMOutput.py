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
