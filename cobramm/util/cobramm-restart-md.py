#! /usr/bin/env python3
# -*- coding: utf-8 -*-

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

import sys
import os
import argparse
import numpy as np

# hijack $PYTHONPATH (sys.path) to find cobramm local modules
sys.path.insert(0, os.path.join(os.getenv('COBRAM_PATH'),'cobramm'))
# now we can import output module
try:
    from output import Output
except ImportError:
    print( 'Error: cannot import cobramm output.py module.\n'
       'Please check if $COBRAM_PATH is defined in your environment, if not run:\n'
       'export COBRAM_PATH=/path/to/COBRAMM\n' )
    sys.exit(1)

r_velocity = 'restart/velocity.dat'
r_crd = 'restart/real.crd'

# open output file
try:
    OUT = Output(parse=True)
except IOError:
    print('ERROR: output file not found in current directory\n'
          'Aborting...\n')
    sys.exit(1)

# check for MD calculation
if not OUT.optimization_type.startswith('mdv'):
    print('ERROR: calculation is not a molecular dynamics run\n'
          'Aborting...\n')
    sys.exit(1)

# print header and steps info
print('\n{0}\n'
      'COBRAMM MD step taker\n'.format('='*40))
if OUT.steps == 0:
    print('Single step found in output file\n')
elif OUT.steps == -1:
    print('ERROR: no steps found\n'
          'Please check if COBRAMM calculation terminated correctly\n')
    sys.exit(1)
else:
    print('{0} steps found in output file [0-{1}]\n'.format(OUT.steps + 1,
                                                            OUT.steps))

# argparse input
parser = argparse.ArgumentParser(description='COBRAMM MD restart utility')
parser.add_argument('-n', '--number', default=None, help='Wanted step number')

args = parser.parse_args()

if args.number is not None:
    try:
        wanted_step = int(args.number)
    except ValueError:
        print('ERROR: An integer is required!\n'
              'Aborting...\n')
        sys.exit(1)
    try:
        wanted_step = range(OUT.steps + 1)[wanted_step]
    except IndexError:
        print('ERROR: Step out of range!\n'
              'Aborting...\n')
        sys.exit(1)
    ask_user_input = False
else:
    ask_user_input = True

# ask user for step number
while ask_user_input:
    try:
        in_wanted_step = input('What step do you need (-1 for last one)? [-1] ')
        if in_wanted_step != '':
            wanted_step = int(in_wanted_step)
        else:
            wanted_step = -1

    except (NameError, ValueError):
        print('Please provide an integer number!\n')
        continue

    except SyntaxError:
        wanted_step = -1

    try:
        wanted_step = range(OUT.steps + 1)[wanted_step]
    except IndexError:
        print('The step you want does not exist:'
              'the number you entered is not in range[0-{0}]\n'.format(OUT.steps))
        continue

    print('The step you want is in range[0-{0}]'.format(OUT.steps))
    print(' OK: writing data from steps [0-{0}]\n'.format(wanted_step))
    break

# create directory
if not os.path.isdir('restart'):
    os.system('mkdir restart')
else:
    print('Directory "restart" already exists, overwriting data...\n')

# fetch Step object
step = OUT.get_step(wanted_step)

# copying files into restart directory
os.system('cp -r real_layers.xyz cobram.command WORK.wfu work.wfu WFU/WORK.wfu WFU/work.wfu wchk real.top model-H.top cobram-sef imomap.dat fort.11 restart 2>/dev/null')

# truncate and copy cobram output into restart directory
with open(OUT.filename,'rb') as f:
    new_text = b''
    for line in f.readlines():
        new_text += line
        if line.strip().decode() == '</STEP {0}>'.format(wanted_step):
            break
    with open('restart/{0}'.format(OUT.filename), 'wb') as fw:
        fw.write(new_text)
os.system('echo " RESTART ----------" >> restart/{0}'.format(OUT.filename))

# write velocity.dat into restart directory
with open(r_velocity, 'w') as fw:
    for vx, vy, vz in zip(*step.velocity):
       fw.write("{0:17.8e}{1:17.8e}{2:17.8e}\n".format(vx, vy, vz))


# write real.crd into restart directory
if os.path.isfile(r_crd):
    os.system('rm {0}'.format(r_crd))
step.geometry.makerealcrd(r_crd)

##################################################


