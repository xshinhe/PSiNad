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

# import statments of module from python standard library

import argparse    # command line arguments parser
import sys         # System-specific parameters and functions
import os          # filesystem utilities
import shutil      # filesystem utilities
import numpy as np
import math

# hijack $PYTHONPATH (sys.path) to find cobramm local modules
try:
    sys.path.append(os.path.join(os.environ['COBRAM_PATH'], 'cobramm'))
except KeyError:
    raise RuntimeError('cannot import cobramm module.\n'
                       'Please check if $COBRAM_PATH is defined in your environment, if not run:\n'
                       'export COBRAM_PATH=/path/to/COBRAMM')

# imports of local modules

import logwrt  # write messages and output to log (std out or cobramm.log)
import cobrammenv  # management of the enviroment for running COBRAMM (and utils)
import constants

from cobrammCalculator import CobrammOutput
#####################################################################################################
def main():

    # this new class defines a style for the parser message that join the two predefined styles
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    def dir_path(string):
        if os.path.isdir(string):
            return string
        else:
            raise NotADirectoryError(f"directory {string} doesn't exist")


    # Define command line arguments with argparser
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="COBRAMM utility script: convolute the UV spectrum",
        epilog="The script convolutes the UV spectrum with a Gaussian function, from the resutls obtained to an ensamble of vertical excitations, each of them convoluted with a Gaussian function", formatter_class=CustomFormatter)

    parser.add_argument("-dir", "--directories", nargs="+", required=True, metavar='DIR', dest="directories_list",
                        type=dir_path, help="path of the directories containing the qmALL.dat file")
    parser.add_argument("-emin", "--emin", default=0.1, metavar='EMIN', dest="emin",
                        type=float, help="minimum energy in the spectrum")
    parser.add_argument("-emax", "--emax", default=10.0, metavar='EMAX', dest="emax",
                        type=float, help="maximum energy in the spectrum")
    parser.add_argument("-ew", "--ewid", default=0.15,  metavar='EWID', dest="ew",
                        type=float, help="broadening in energy domain (eV)")
    parser.add_argument("-g", "--grdid", default=200,  metavar='GRID', dest="grid",
                        type=int, help="number of grid points of the specrum")
    parser.add_argument("-s", "--state", default = -1,metavar='STATE', dest="state",
                        type=int, help="convolute the spectrum targeting one specific state")
    parser.add_argument('--version', action='version',
                        version='%(prog)s is part of COBRAMM {}'.format(cobrammenv.getVersion()))
    # TODO ADD PARSER OUTPUT

    args = parser.parse_args()

    collected_energies = []
    collected_intensities = []
    grid, spectrum = [], []
    

    #get the values
    for folder in args.directories_list:
        if args.state > 0:
            el_energies, osci_strength = CobrammOutput._states_from_gaussian_tddft_out(os.path.join(folder, "QM_data", "qmALL.log"),decomposition=True, state = args.state)
        else:
            el_energies, osci_strength = CobrammOutput._states_from_gaussian_tddft_out(os.path.join(folder, "QM_data", "qmALL.log"))
        for i in range(len(el_energies)):
            collected_energies.append(el_energies[i]*constants.Hartree2eV)
            collected_intensities.append(osci_strength[i])

    #define grid
    spectral_grid = np.linspace(args.emin, args.emax, args.grid) 
    spec_value = np.zeros(len(spectral_grid))

    #convolute spectrum
    for center, intensity in zip(collected_energies,collected_intensities):
        linef = np.array([1. / (args.ew * math.sqrt(2.0 * math.pi)) * np.exp(-(e - center) ** 2 / (2.0 * args.ew ** 2))
                            for e in spectral_grid])
        to_add = intensity*linef
        #print(to_add)
        spec_value += to_add

    #store in list format
    for en, ints_e in zip(spectral_grid, spec_value):
        grid.append(en)
        spectrum.append(ints_e) #/ wavelength**2)

    #write in list format
    comment = "#\n" \
            "#      energy(eV)    wavelength(nm)    intensity(a.u.)\n" \
            "#\n"

    if args.state > 0:
        out = "spectrum_S_{}.dat".format(args.state)
    else: 
        out = "spectrum.dat"
    with open(out, "w") as f:
        f.write(comment)
        for ev, ints in zip(grid, spectrum):
            f.write("{0:16.4f} {1:16.4f} {2:16.4f}\n".format(ev, 1239.8/ev, ints))



#####################################################################################################

if __name__ == '__main__':
    main()
