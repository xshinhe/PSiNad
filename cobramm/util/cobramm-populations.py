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
    sys.path.insert(0, os.path.join(os.getenv('COBRAM_PATH'),'cobramm'))
except KeyError:
    raise RuntimeError('cannot import cobramm module.\n'
                       'Please check if $COBRAM_PATH is defined in your environment, if not run:\n'
                       'export COBRAM_PATH=/path/to/COBRAMM')

# imports of local modules

import logwrt  # write messages and output to log (std out or cobramm.log)
import cobrammenv  # management of the enviroment for running COBRAMM (and utils)

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
                        type=dir_path, help="path of the directories containing the AmplitudesALL.dat and occupation.dat files")
    parser.add_argument("-ns", "--nstates", required=True, metavar='NSTATE', dest="nstates",
                        type=int, help="number of electronic states to consider")
    parser.add_argument("-st", "--simtime", required=True, metavar='SIMTIME', dest="simulation_time",
                        type=int, help="simulation time of the dynamics (fs)")
    parser.add_argument("-ts", "--timestep", required=True, metavar='TIMESTEP', dest="time_step",
                        help="time step performed (fs)")
    parser.add_argument('--version', action='version',
                        version='%(prog)s is part of COBRAMM {}'.format(cobrammenv.getVersion()))


    args = parser.parse_args()



    try:
        ts = float(args.time_step)
    except:
        try:
            ts = int(args.time_step)
        except:
            raise

    nsteps = int(args.simulation_time/ts+1)

    time = np.arange(0, args.simulation_time, ts).reshape(nsteps-1, 1)

    #AMPLITUDES

    temp = np.zeros((nsteps-1, args.nstates))

    for directory in args.directories_list:
        ampl = np.loadtxt("{}/AmplitudesALL.dat".format(directory))
        sqr_parts = np.array(np.square(ampl[1:, 1:]))
        final_ampl = np.array(([np.add(sqr_parts[0, state], sqr_parts[0, state+1])
                                for state in range(0, args.nstates*2, 2)]))
        
        for step in range(1, nsteps-1):
            summed_parts = np.array(([np.add(sqr_parts[step, state], sqr_parts[step, state+1])
                                      for state in range(0, args.nstates*2, 2)]))
            final_ampl = np.row_stack((final_ampl, summed_parts))
        temp = np.add(temp, final_ampl)

    final_temp = np.divide(temp, len(args.directories_list))
    pop = np.hstack((time, final_temp))

    #ADD PARSER OUTPUT

    pop_amplitudes = np.mat(pop)

    with open("populations.dat", "w") as out:
        comment_string=f'{"## t": <19}'
        for state in range(0,args.nstates):
            comment_string+="S_{singlet: <15}".format(singlet=state)
        comment_string += "\n##\n##\n"
        out.write(comment_string)
        for row in pop_amplitudes:
            np.savetxt(out, row,fmt='%16.12f')

    #FRACTION OF TRAJECTORIES

    temp = np.zeros((nsteps-1, args.nstates))

    for directory in args.directories_list:
        act_states = np.loadtxt("{}/occupation.dat".format(directory), skiprows=1, max_rows=nsteps-1, usecols=(np.arange(1, args.nstates+1)))
        temp = np.add(temp, act_states)

    final_temp = np.divide(temp, len(args.directories_list))
    fraction = np.mat(np.hstack((time, final_temp)))

    with open("fractions_trj.dat", "w") as out:
        comment_string = f'{"## t": <19}'
        for state in range(0, args.nstates):
            comment_string += "S_{singlet: <15}".format(singlet=state)
        comment_string +="\n##\n##\n"
        out.write(comment_string)
        for row in fraction:
            np.savetxt(out, row, fmt='%16.12f')


#####################################################################################################

if __name__ == '__main__':
    main()
