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

# check if cobramm input file exists
if not os.path.isfile('cobram.command'):
    print('ERROR: cobram.command not found in current directory\n'
          'Aborting...\n')
    sys.exit(1)

# open output file
try:
    OUT = Output(parse=True)
except IOError:
    print('ERROR: output file not found in current directory\n'
          'Aborting...\n')
    sys.exit(1)

# print header and steps info
print('\n{0}\n'
      'COBRAMM step taker\n'.format('='*40))
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
parser = argparse.ArgumentParser(description='Cobramm step taker')
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
    print(' OK: writing data from step {0}\n'.format(wanted_step))
    break

# create directory
if wanted_step == OUT.steps:
    namedir = 'inputs_from_last_step'
else:
    namedir = 'inputs_from_step_{0}'.format(wanted_step)
if not os.path.isdir(namedir):
    os.system('mkdir {0}'.format(namedir))
else:
    print('Directory {0} already exists, overwriting data...\n'.format(namedir))

# copy input files in target directory
os.system('cp cobram.command {0}/'.format(namedir))
if os.path.isfile('real_layers.xyz'):
    os.system('cp real_layers.xyz {0}/'.format(namedir))
else:
    print('WARNING: real_layers.xyz not found in current directory.\n'
          'Automatic build is not implemented, please copy real_layers.xyz in {0}\n'.format(namedir))

# fetch Step object
step = OUT.get_step(wanted_step)

if OUT.calculation_type in ['HML', 'HM', 'HL']:
    with open('{0}/model-H.top'.format(namedir), 'w') as fw:
        fw.write(OUT.modelH_top)

if OUT.calculation_type in ['HML', 'HM', 'HL', 'ML', 'M']:
    with open('{0}/real.top'.format(namedir), 'w') as fw:
        fw.write(OUT.real_top)

# delete real.crd file if already exists
crd_file = '{0}/real.crd'.format(namedir)
if os.path.isfile(crd_file):
    os.remove(crd_file)
# make real.crd with geometry from current step
step.geometry.makerealcrd(crd_file)

# QM checkpoint data
try:
    QM_dir = os.listdir('QM_data')
    chk_file_names = {'gaussian-QM_{0}.chk.gz': 'gaussian-QM.chk',
                      'molcas_{0}.ScfOrb.gz':   'INPORB',
                      'molcas_{0}.RasOrb.gz':   'INPORB',
                      'work{0}.wfu.gz':         'INP.wfu'}
    # loop over dict keys
    for f in chk_file_names.keys():
        fname = f.format(wanted_step)
        if fname in QM_dir:
            os.system('cp -f QM_data/{0} {1}'.format(fname, namedir))
            os.system('gzip -df {0}/{1}'.format(namedir, fname))
            os.system('mv {0}/{1} {0}/{2}'.format(namedir, fname[:-3], chk_file_names[f]))
except OSError:
    pass

# geometry checkpoint data
list_dir= os.listdir('.')
if 'geometry.chk.gz' in list_dir:
    os.system('/bin/cp -f geometry.chk.gz {0}'.format(namedir))
    os.system('gzip -df {0}/geometry.chk.gz'.format(namedir))
elif 'geometry.chk' in list_dir:
    os.system('/bin/cp -f geometry.chk {0}'.format(namedir))
else:
    print('no geometry checkpoint file in your directory\n')

# print footer
print('COBRAMM step taker ended\n'
      '{0}\n'.format('='*40))
