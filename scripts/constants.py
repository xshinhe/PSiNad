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


au_2_ang = 5.291772104260590e-01
au_2_kcal_1mea = 6.275094742071363e+02
au_2_kcal_1mea_per_ang =  1.185821047928172e+03
au_2_ev =  2.7211386033e+01
au_2_wn =  2.1947463147e+05

periodic_table = '''
  -----                                                               -----
1 | H |                                                               |He |
  |---+----                                       --------------------+---|
2 |Li |Be |                                       | B | C | N | O | F |Ne |
  |---+---|                                       |---+---+---+---+---+---|
3 |Na |Mg |3B  4B  5B  6B  7B |    8B     |1B  2B |Al |Si | P | S |Cl |Ar |
  |---+---+---------------------------------------+---+---+---+---+---+---|
4 | K |Ca |Sc |Ti | V |Cr |Mn |Fe |Co |Ni |Cu |Zn |Ga |Ge |As |Se |Br |Kr |
  |---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---|
5 |Rb |Sr | Y |Zr |Nb |Mo |Tc |Ru |Rh |Pd |Ag |Cd |In |Sn |Sb |Te | I |Xe |
  |---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---|
6 |Cs |Ba |LAN|Hf |Ta | W |Re |Os |Ir |Pt |Au |Hg |Tl |Pb |Bi |Po |At |Rn |
  |---+---+---+------------------------------------------------------------
7 |Fr |Ra |ACT|
  -------------
              -------------------------------------------------------------
   Lanthanide |La |Ce |Pr |Nd |Pm |Sm |Eu |Gd |Tb |Dy |Ho |Er |Tm |Yb |Lu |
              |---+---+---+---+---+---+---+---+---+---+---+---+---+---+---|
   Actinide   |Ac |Th |Pa | U |Np |Pu |Am |Cm |Bk |Cf |Es |Fm |Md |No |Lw |
              -------------------------------------------------------------
'''

element_list = {
    "NU": ["NU",  0, -1, 0.0          ],
    "H" : ["H" ,  1, -1, 1.007947     ],
    "HE": ["HE",  2, -1, 4.0026022    ],
    "LI": ["LI",  3, -1, 6.9412       ],
    "BE": ["BE",  4, -1, 9.0121823    ],
    "B" : ["B" ,  5, -1, 10.8117      ],
    "C" : ["C" ,  6, -1, 12.01078     ],
    "N" : ["N" ,  7, -1, 14.00672     ],
    "O" : ["O" ,  8, -1, 15.99943     ],
    "F" : ["F" ,  9, -1, 18.99840325  ],
    "NE": ["NE", 10, -1, 20.17976     ],
    "NA": ["NA", 11, -1, 22.989769282 ],
    "MG": ["MG", 12, -1, 24.30506     ],
    "AL": ["AL", 13, -1, 26.98153868  ],
    "SI": ["SI", 14, -1, 28.08553     ],
    "P" : ["P" , 15, -1, 30.9737622   ],
    "S" : ["S" , 16, -1, 32.0655      ],
    "CL": ["CL", 17, -1, 35.4532      ],
    "AR": ["AR", 18, -1, 39.9481      ],
    "K" : ["K" , 19, -1, 39.09831     ],
    "CA": ["CA", 20, -1, 40.0784      ],
    "SC": ["SC", 21, -1, 44.9559126   ],
    "TI": ["TI", 22, -1, 47.8671      ],
    "V" : ["V" , 23, -1, 50.94151     ],
    "CR": ["CR", 24, -1, 51.99616     ],
    "MN": ["MN", 25, -1, 54.9380455   ],
    "FE": ["FE", 26, -1, 55.8452      ],
    "CO": ["CO", 27, -1, 58.9331955   ],
    "NI": ["NI", 28, -1, 58.69342     ],
    "CU": ["CU", 29, -1, 63.5463      ],
    "ZN": ["ZN", 30, -1, 65.4094      ],
    "GA": ["GA", 31, -1, 69.7231      ],
    "GE": ["GE", 32, -1, 72.641       ],
    "AS": ["AS", 33, -1, 74.921602    ],
    "SE": ["SE", 34, -1, 78.963       ],
    "BR": ["BR", 35, -1, 79.9041      ],
    "KR": ["KR", 36, -1, 83.7982      ],
    "RB": ["RB", 37, -1, 85.46783     ],
    "SR": ["SR", 38, -1, 87.621       ],
    "Y" : ["Y" , 39, -1, 88.905852    ],
    "ZR": ["ZR", 40, -1, 91.2242      ],
    "NB": ["NB", 41, -1, 92.906382    ],
    "MO": ["MO", 42, -1, 95.942       ],
    "TC": ["TC", 43, -1, 98.0         ],
    "RU": ["RU", 44, -1, 101.072      ],
    "RH": ["RH", 45, -1, 102.905502   ],
    "PD": ["PD", 46, -1, 106.421      ],
    "AG": ["AG", 47, -1, 107.86822    ],
    "CD": ["CD", 48, -1, 112.4118     ],
    "IN": ["IN", 49, -1, 114.8183     ],
    "SN": ["SN", 50, -1, 118.7107     ],
    "SB": ["SB", 51, -1, 121.7601     ],
    "TE": ["TE", 52, -1, 127.603      ],
    "I" : ["I" , 53, -1, 126.904473   ],
    "XE": ["XE", 54, -1, 131.2936     ],
    "CS": ["CS", 55, -1, 132.90545192 ],
    "BA": ["BA", 56, -1, 137.3277     ],
    "LA": ["LA", 57, -1, 138.905477   ],
    "CE": ["CE", 58, -1, 140.1161     ],
    "PR": ["PR", 59, -1, 140.907652   ],
    "ND": ["ND", 60, -1, 144.2423     ],
    "PM": ["PM", 61, -1, 145.0        ],
    "SM": ["SM", 62, -1, 150.362      ],
    "EU": ["EU", 63, -1, 151.9641     ],
    "GD": ["GD", 64, -1, 157.253      ],
    "TB": ["TB", 65, -1, 158.925352   ],
    "DY": ["DY", 66, -1, 162.5001     ],
    "HO": ["HO", 67, -1, 164.930322   ],
    "ER": ["ER", 68, -1, 167.2593     ],
    "TM": ["TM", 69, -1, 168.934212   ],
    "YB": ["YB", 70, -1, 173.043      ],
    "LU": ["LU", 71, -1, 174.9671     ],
    "HF": ["HF", 72, -1, 178.492      ],
    "TA": ["TA", 73, -1, 180.947882   ],
    "W" : ["W" , 74, -1, 183.841      ],
    "RE": ["RE", 75, -1, 186.2071     ],
    "OS": ["OS", 76, -1, 190.233      ],
    "IR": ["IR", 77, -1, 192.2173     ],
    "PT": ["PT", 78, -1, 195.0849     ],
    "AU": ["AU", 79, -1, 196.9665694  ],
    "HG": ["HG", 80, -1, 200.592      ],
    "TL": ["TL", 81, -1, 204.38332    ],
    "PB": ["PB", 82, -1, 207.21       ],
    "BI": ["BI", 83, -1, 208.980401   ],
    "PO": ["PO", 84, -1, 210.0        ],
    "AT": ["AT", 85, -1, 210.0        ],
    "RN": ["RN", 86, -1, 220.0        ],
    "FR": ["FR", 87, -1, 223.0        ],
    "RA": ["RA", 88, -1, 226.0        ],
    "AC": ["AC", 89, -1, 227.0        ],
    "TH": ["TH", 90, -1, 232.038062   ],
    "PA": ["PA", 91, -1, 231.035882   ],
    "U" : ["U" , 92, -1, 238.028913   ],
    "NP": ["NP", 93, -1, 237.0        ],
    "PU": ["PU", 94, -1, 244.0        ],
    "AM": ["AM", 95, -1, 243.0        ],
    "CM": ["CM", 96, -1, 247.0        ],
    "BK": ["BK", 97, -1, 247.0        ],
    "CF": ["CF", 98, -1, 251.0        ],
    "ES": ["ES", 99, -1, 252.0        ],
    "FM": ["FM",100, -1, 257.0        ],
    "MD": ["MD",101, -1, 258.0        ],
    "NO": ["NO",102, -1, 259.0        ],
    "LR": ["LR",103, -1, 262.0        ],
    "RF": ["RF",104, -1, 261.0        ],
    "DB": ["DB",105, -1, 262.0        ],
    "SG": ["SG",106, -1, 266.0        ],
    "BH": ["BH",107, -1, 264.0        ],
    "HS": ["HS",108, -1, 277.0        ],
    "MT": ["MT",109, -1, 268.0        ],
    "DS": ["DS",110, -1, 271.0        ],
    "RG": ["RG",111, -1, 272.0        ],
    "CN": ["CN",112, -1, 285.0        ],
    "NH": ["NH",113, -1, 284.0        ],
    "FL": ["FL",114, -1, 289.0        ],
    "MC": ["MC",115, -1, 288.0        ],
    "LV": ["LV",116, -1, 292.0        ],
    "TS": ["TS",117, -1, 291.0        ],
    "OG": ["OG",118, -1, 294.0        ],
    #
    0   : ["NU",  0, -1, 0.0          ],
    1   : ["H" ,  1, -1, 1.007947     ],
    2   : ["HE",  2, -1, 4.0026022    ],
    3   : ["LI",  3, -1, 6.9412       ],
    4   : ["BE",  4, -1, 9.0121823    ],
    5   : ["B" ,  5, -1, 10.8117      ],
    6   : ["C" ,  6, -1, 12.01078     ],
    7   : ["N" ,  7, -1, 14.00672     ],
    8   : ["O" ,  8, -1, 15.99943     ],
    9   : ["F" ,  9, -1, 18.99840325  ],
    10  : ["NE", 10, -1, 20.17976     ],
    11  : ["NA", 11, -1, 22.989769282 ],
    12  : ["MG", 12, -1, 24.30506     ],
    13  : ["AL", 13, -1, 26.98153868  ],
    14  : ["SI", 14, -1, 28.08553     ],
    15  : ["P" , 15, -1, 30.9737622   ],
    16  : ["S" , 16, -1, 32.0655      ],
    17  : ["CL", 17, -1, 35.4532      ],
    18  : ["AR", 18, -1, 39.9481      ],
    19  : ["K" , 19, -1, 39.09831     ],
    20  : ["CA", 20, -1, 40.0784      ],
    21  : ["SC", 21, -1, 44.9559126   ],
    22  : ["TI", 22, -1, 47.8671      ],
    23  : ["V" , 23, -1, 50.94151     ],
    24  : ["CR", 24, -1, 51.99616     ],
    25  : ["MN", 25, -1, 54.9380455   ],
    26  : ["FE", 26, -1, 55.8452      ],
    27  : ["CO", 27, -1, 58.9331955   ],
    28  : ["NI", 28, -1, 58.69342     ],
    29  : ["CU", 29, -1, 63.5463      ],
    30  : ["ZN", 30, -1, 65.4094      ],
    31  : ["GA", 31, -1, 69.7231      ],
    32  : ["GE", 32, -1, 72.641       ],
    33  : ["AS", 33, -1, 74.921602    ],
    34  : ["SE", 34, -1, 78.963       ],
    35  : ["BR", 35, -1, 79.9041      ],
    36  : ["KR", 36, -1, 83.7982      ],
    37  : ["RB", 37, -1, 85.46783     ],
    38  : ["SR", 38, -1, 87.621       ],
    39  : ["Y" , 39, -1, 88.905852    ],
    40  : ["ZR", 40, -1, 91.2242      ],
    41  : ["NB", 41, -1, 92.906382    ],
    42  : ["MO", 42, -1, 95.942       ],
    43  : ["TC", 43, -1, 98.0         ],
    44  : ["RU", 44, -1, 101.072      ],
    45  : ["RH", 45, -1, 102.905502   ],
    46  : ["PD", 46, -1, 106.421      ],
    47  : ["AG", 47, -1, 107.86822    ],
    48  : ["CD", 48, -1, 112.4118     ],
    49  : ["IN", 49, -1, 114.8183     ],
    50  : ["SN", 50, -1, 118.7107     ],
    51  : ["SB", 51, -1, 121.7601     ],
    52  : ["TE", 52, -1, 127.603      ],
    53  : ["I" , 53, -1, 126.904473   ],
    54  : ["XE", 54, -1, 131.2936     ],
    55  : ["CS", 55, -1, 132.90545192 ],
    56  : ["BA", 56, -1, 137.3277     ],
    57  : ["LA", 57, -1, 138.905477   ],
    58  : ["CE", 58, -1, 140.1161     ],
    59  : ["PR", 59, -1, 140.907652   ],
    60  : ["ND", 60, -1, 144.2423     ],
    61  : ["PM", 61, -1, 145.0        ],
    62  : ["SM", 62, -1, 150.362      ],
    63  : ["EU", 63, -1, 151.9641     ],
    64  : ["GD", 64, -1, 157.253      ],
    65  : ["TB", 65, -1, 158.925352   ],
    66  : ["DY", 66, -1, 162.5001     ],
    67  : ["HO", 67, -1, 164.930322   ],
    68  : ["ER", 68, -1, 167.2593     ],
    69  : ["TM", 69, -1, 168.934212   ],
    70  : ["YB", 70, -1, 173.043      ],
    71  : ["LU", 71, -1, 174.9671     ],
    72  : ["HF", 72, -1, 178.492      ],
    73  : ["TA", 73, -1, 180.947882   ],
    74  : ["W" , 74, -1, 183.841      ],
    75  : ["RE", 75, -1, 186.2071     ],
    76  : ["OS", 76, -1, 190.233      ],
    77  : ["IR", 77, -1, 192.2173     ],
    78  : ["PT", 78, -1, 195.0849     ],
    79  : ["AU", 79, -1, 196.9665694  ],
    80  : ["HG", 80, -1, 200.592      ],
    81  : ["TL", 81, -1, 204.38332    ],
    82  : ["PB", 82, -1, 207.21       ],
    83  : ["BI", 83, -1, 208.980401   ],
    84  : ["PO", 84, -1, 210.0        ],
    85  : ["AT", 85, -1, 210.0        ],
    86  : ["RN", 86, -1, 220.0        ],
    87  : ["FR", 87, -1, 223.0        ],
    88  : ["RA", 88, -1, 226.0        ],
    89  : ["AC", 89, -1, 227.0        ],
    90  : ["TH", 90, -1, 232.038062   ],
    91  : ["PA", 91, -1, 231.035882   ],
    92  : ["U" , 92, -1, 238.028913   ],
    93  : ["NP", 93, -1, 237.0        ],
    94  : ["PU", 94, -1, 244.0        ],
    95  : ["AM", 95, -1, 243.0        ],
    96  : ["CM", 96, -1, 247.0        ],
    97  : ["BK", 97, -1, 247.0        ],
    98  : ["CF", 98, -1, 251.0        ],
    99  : ["ES", 99, -1, 252.0        ],
    100 : ["FM",100, -1, 257.0        ],
    101 : ["MD",101, -1, 258.0        ],
    102 : ["NO",102, -1, 259.0        ],
    103 : ["LR",103, -1, 262.0        ],
    104 : ["RF",104, -1, 261.0        ],
    105 : ["DB",105, -1, 262.0        ],
    106 : ["SG",106, -1, 266.0        ],
    107 : ["BH",107, -1, 264.0        ],
    108 : ["HS",108, -1, 277.0        ],
    109 : ["MT",109, -1, 268.0        ],
    110 : ["DS",110, -1, 271.0        ],
    111 : ["RG",111, -1, 272.0        ],
    112 : ["CN",112, -1, 285.0        ],
    113 : ["NH",113, -1, 284.0        ],
    114 : ["FL",114, -1, 289.0        ],
    115 : ["MC",115, -1, 288.0        ],
    116 : ["LV",116, -1, 292.0        ],
    117 : ["TS",117, -1, 291.0        ],
    118 : ["OG",118, -1, 294.0        ],
}
