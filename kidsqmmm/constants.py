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

# Constants Definition of mathematical / physical constants and units conversion factors
#
# Use MATH python module \n  
# Sources of the constants: \n  
#       1) The NIST Reference on Constants, Units and Uncertainty ( http://physics.nist.gov/cuu/index.html ) \n 
#       2) The NIST Atomic Weights and Isotopic Compositions ( http://www.nist.gov/pml/data/comp.cfm ) \n

import math

# PHYSICAL CONSTANTS

# Boltzmann's constant in eV/K
kb = 8.617332478e-5
# Atomic Mass Unit (AMU) in kg (from NIST reference)
uma = 1.660538921e-27
# Electron mass in kg (from NIST reference)
mel = 9.10938291e-31
#  Avogadro's number (from NIST reference)
NAvo = 6.02214129E23
#  Thermochemical Calorie ( from NIST reference )
ThermoCal = 4.184
#  Atomic unit of time in fs ( from NIST reference )
AUofTime = 2.418884326502E-2
# Speed of light in the vacuum, in m/s ( from NIST reference )
SpeedOfLight = 299792458

# CONVERSION FACTORS

# Hartree - eV conversion factor (from NIST reference)
Hartree2eV = 27.21138505
# Conversion factor from Hartree to Joule (from NIST reference)
Hatree2joule = 4.35974434e-18
# Conversion factor from Hartree to kJ/mol
Hatree2kjmol = Hatree2joule * NAvo / 1000.
# Conversion factor from Hartree to kcal/mol
Hatree2kcalmol = Hatree2kjmol / ThermoCal

# Decimal degrees - radiants conversion factor
Deg2Rad = math.pi / 180.0

# Bohr - Ang conversion factor (from NIST reference)
Bohr2Ang = 0.52917721092

# UMA - Atomic Units conversion factor
amu2au = uma / mel

# Conversion factor from picosecond to Atomic Units
ps2au = 1000.0 / AUofTime
# Conversion factor from femtosecond to Atomic Units
fs2au = 1.0 / AUofTime
# Conversion factor from time unit of AMBER to Atomic Units
MDtime2au = 1000. * fs2au / 20.455

# Conversion from cm-1 to fs-1 (2PI for angular freq conversion)
wavnr2Hz = 2.0 * math.pi * SpeedOfLight * 1.0E-13
# Conversion from cm-1 to Atomic Units (freq=energy)
wavnr2au = wavnr2Hz / fs2au

# Conversion from K to au
K2AU = kb / Hartree2eV
# Boltzmann's constant in Hartree/K
kbAU = kb / Hartree2eV
# convert dipole moment from Debye to au
Debye2AU = 0.39342215569939517797

# Conversion factor from charge unit of AMBER to Atomic Units
MDcharge2au = 1.0 / math.sqrt(332.0)

# ELEMENTS MASSES ( http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some )
# isotopes with highest abundance are defined
# only for Chlorine, ATOMIC WEIGTH using the 0.76-0.24 standard composition of 35+37

masses = { "h":  1.0078250321 * amu2au,  "d": 2.01410178 * amu2au, 
          "li":  7.0160034366 * amu2au, "be":  9.0121822 * amu2au,  "b": 11.0093054 * amu2au, "c": 12.0000000 * amu2au,  "n": 14.0030740044 * amu2au, "o": 15.99491462 * amu2au,  "f": 18.99840316 * amu2au, 
          "na": 22.9897692820 * amu2au, "mg": 23.9850417 * amu2au, "al": 26.9815386 * amu2au, "p": 30.9737620 * amu2au, "si": 27.9769265325 * amu2au, "s": 31.97207100 * amu2au, "cl": 34.96885268 * amu2au}


# now define function to return mass values, with two functions:
# - identify the atom with case-insensitive labels
# - return None if the atom is not in the dictionary
def atommass(atomname):
    if atomname.lower() in masses:
        return masses[atomname.lower()]
    else:
        return None

# BOND LENGTH

# dictionary for standard X-H single bond length
_bonds = {
    ("H" , "H"):  0.740 / Bohr2Ang,
    ("C" , "H"):  1.090 / Bohr2Ang,
    ("N" , "H"):  1.008 / Bohr2Ang,
    ("O" , "H"):  0.947 / Bohr2Ang,
    ("F" , "H"):  0.920 / Bohr2Ang,
    ("SI", "H"):  1.480 / Bohr2Ang,
    ("P" , "H"):  1.420 / Bohr2Ang,
    ("S" , "H"):  1.340 / Bohr2Ang,
    ("CL", "H"):  1.270 / Bohr2Ang,
}
# make the dictionary symmetric with respect to the key couples
sngBond = {}
for k, v in _bonds.items():
    sngBond[k] = v
    if (k[1], k[0]) not in sngBond:
        sngBond[(k[1], k[0])] = v
