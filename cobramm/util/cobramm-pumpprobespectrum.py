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
from cobrammCalculator import CobrammCalculator, CobrammInput, CobrammOutput
from transientCalculator import PumpProbeCalculator, MultiwfnCalculator

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
        description="COBRAMM utility script: extract population during a TSH run from quantum amplitudes",
        epilog="The script calculates the population for an ensemble of TSH trajectories from the sum absolute squared of "
               "the quantum amplitudes and the fraction of population for each excited states\n\n",
        formatter_class=CustomFormatter)

    parser.add_argument("-dir", "--directories", nargs="+", required=True, metavar='DIR', dest="directories_list",
                        type=dir_path, help="path of the directories containing the trajectories")
    parser.add_argument("-m", "--method", required=True, metavar="METHOD", dest="method",
                        type=str, help="QM method of choice")
    parser.add_argument("-ns", "--nstates", required=True, metavar='NSTATE', dest="nstates",
                        type=int, help="number of electronic states to consider")
    parser.add_argument("-emin", "--emin", default=0.1, metavar='EMIN', dest="emin",
                        type=float, help="minimum energy in the spectrum")
    parser.add_argument("-emax", "--emax", default=10.0, metavar='EMAX', dest="emax",
                        type=float, help="maximum energy in the spectrum")
    parser.add_argument("-ew", "--ewid", default=0.15,  metavar='EWID', dest="ew",
                        type=float, help="broadening in energy domain (eV)")
    parser.add_argument("-tw", "--twid", default=25,  metavar='TWID', dest="tw",
                        type=int, help="broadening in time domain (fs)")
    parser.add_argument("-g", "--grdid", default=400,  metavar='GRID', dest="grid",
                        type=int, help="number of grid points of the specrum")
    parser.add_argument("-st", "--simtime", required=True, metavar='SIMTIME', dest="simulation_time",
                        type=int, help="simulation time of the dynamics (fs)")
    parser.add_argument("-ts", "--timestep", required=True, metavar='TIMESTEP', dest="time_step",
                        help="time step performed (fs)")
    parser.add_argument("-pol", "--polarization",  default="none", dest="polarization",
                        help="polarization scheme required")
    parser.add_argument("-dec", "--decomposition", dest="decomposition",
                        type=int, help="state from which transitions will be considered")
    parser.add_argument("-Dt", "--deltat", required=True, metavar="DELTATIME", dest="delta_t",
                        type=int, help="delta of time vertical excitations are performed (fs)")
    parser.add_argument('--version', action='version',
                        version='%(prog)s is part of COBRAMM {}'.format(cobrammenv.getVersion()))
    # TODO ADD PARSER OUTPUT

    args = parser.parse_args()

    if args.decomposition:

        selected = args.decomposition
        for target in range(2, 12):

            pp = PumpProbeCalculator(method=args.method, simulation_time=args.simulation_time, delta_t=args.delta_t,
                                     nstates=args.nstates, trjdir_list=args.directories_list,en_min=args.emin,
                                     en_max=args.emax, en_width=args.ew, t_width=args.tw)

            print("convolution ESA from S_{0} to S_{1}\n".format(selected,target))
            for step in range(0, pp.nsteps+1, pp.delta):
                t = step*pp.time_step

                print("Convolution of the spectrum for time {}\n".format(t))

                collected_values = pp.collect_decomposed_values_single_time(nstep=step,selected_state=selected,target_state=target)
                pp.get_spectrum_single_time_all_traj(values=collected_values, currtime=t)

            print("Convolution of the final pump-probe spectrum")

            pp.time_convolution(output="spectrum_from_S{0}_to_S{1}".format(selected, target))

    else:
        #initialize class
        pp = PumpProbeCalculator(method=args.method, simulation_time=args.simulation_time, delta_t=args.delta_t, nstates=args.nstates,
                                 trjdir_list=args.directories_list, en_min=args.emin, en_max=args.emax,
                                 en_width=args.ew, t_width=args.tw)

        print("Convolution of the spectra for each time\n")

        all_init_tdm = pp.extract_init_tdm()

        for step in range(0, pp.nsteps + 1, pp.delta):
            t = step * pp.time_step

            print("Convolution of the spectrum for time {}\n".format(t))

            collected_values = pp.collect_values_single_time(nstep=step, polarization=args.polarization,all_init_tdm=all_init_tdm)

            #print(collected_values)
            pp.get_spectrum_single_time_all_traj(values=collected_values, currtime=t)
        print("Convoluting the final pump-probe spectrum")

        pp.time_convolution(output="spectrumPP.txt")
        print("Spectrum written in 'spectrumPP.txt'\n")
        print("Writing gnuplot script as 'plot_matrix.gp'")
        pp.write_gnuplot(toplot="spectrumPP.txt")

#####################################################################################################

if __name__ == '__main__':
    main()
