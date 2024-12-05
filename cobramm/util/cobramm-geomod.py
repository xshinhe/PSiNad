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

# standard library
import sys
import os
import argparse

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

# open output file
try:
    OUT = Output(parse=True)
except IOError:
    print('ERROR: output file not found in current directory\n'
          'Aborting...\n')
    sys.exit(1)

# argparse input
parser = argparse.ArgumentParser(description='COBRAMM geometry edit utility')
parser.add_argument('-n', '--number', default=None, help='Wanted step number')

args = parser.parse_args()


# if edit-me.xyz is found, create new_real.crd
dir_elemets = os.listdir('.')
if 'edit-me.xyz' in dir_elemets and args.number is None:
    print('\nedit-me.xyz found in current directory, new_real.crd will be created')
    with open('edit-me.xyz') as f:
        data = f.readlines()[2:]

    geom = OUT.geometry
    X,Y,Z = geom.cartesian[:]
    HM_list = geom.list_MEDIUM_HIGH
    for i,line in enumerate(data):
        lsplit = line.split()
        atom_index = HM_list[i] - 1
        X[atom_index] = lsplit[1]
        Y[atom_index] = lsplit[2]
        Z[atom_index] = lsplit[3]

    geom.cartesian = [X,Y,Z]

    if 'new_real.crd' in dir_elemets:
        os.remove('new_real.crd')
    geom.makerealcrd('new_real.crd')
    os.remove('edit-me.xyz')
    print('\nnew_real.crd created in current directory\n')
    sys.exit(0)


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

    break

step = OUT.get_step(wanted_step)
geom_string = step.geometry.to_string(geom_model='MEDIUM_HIGH',
                                      labels=True,
                                      precision=7)
atom_num = step.geometry.NatomHM
with open('edit-me.xyz','w') as fw:
    fw.write('{0}\n\n{1}'.format(atom_num, geom_string))
print('\nedit-me.xyz created in current directory\n'
      'Please edit it, then re-launch this script\n')
