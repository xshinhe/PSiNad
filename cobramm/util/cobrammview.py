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

import os
import sys
import subprocess
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

# set default print config for numpy
try:
    np.set_printoptions(threshold=sys.maxsize,      # try to avoid truncation for long arrays
                        suppress=True,              # try to avoid scientific notation
                        sign=' ',                   # whitespace in front of positive values
                        floatmode='fixed')          # print all decimals for selected precision (default 8)
except TypeError:
    np.set_printoptions(threshold=sys.maxsize, suppress=True)

# argparse init
parser = argparse.ArgumentParser(description='Cobramm viewer utility')
parser.add_argument('-xyz', '--xyz', help='xyz trajectories ', default=None, action='store_true')
parser.add_argument('-pdb', '--pdb', help='pdb file for selected step (default: 0)', default=None, action='store_true')
parser.add_argument('-mdcrd', '--mdcrd', help='mdcrd trajectories', default=None, action='store_true')
parser.add_argument('-molden', '--molden', help='molden trajectories and forces', default=None, action='store_true')
parser.add_argument('-n', '--numstep', help='select step number (only for xyz and pdb)', default=None)
parser.add_argument('-m', '--model', help='select layers model (H/M/HM/HML/ML)', default=None)

args = parser.parse_args()

# check if at least one mandatory keyword has been supplied by user
if not np.any([args.xyz is not None,
               args.pdb is not None,
               args.mdcrd is not None,
               args.molden is not None]):
    # print help
    subprocess.call('cobrammview.py -h'.split())
    sys.exit(0)

# read cobramm.xml output
OUT = Output(parse=True)
nsteps = OUT.steps
if args.numstep is not None:
    try:
        wanted_step = range(nsteps + 1)[int(args.numstep)]
    except IndexError:
        print('step {0} out of range [0-{1}]\n'
              'Aborting...'.format(args.numstep, nsteps))
        sys.exit(1)
else:
    wanted_step = None

# define steps_range, skip step 0
if OUT.optimization_type.startswith('irc'):
    tmp_range = range(1, nsteps + 2)
    substeps = OUT.grep_all('nIRCstep')
    steps_range = []
    for i in tmp_range:
        # select only converged points for IRC calculations
        if i <= nsteps:
            substep = substeps[i]
            # skipping first step
            if i == 1:
                continue
            elif i < nsteps:
                # skipping non-converged points up to second-last step
                if substep != 1:
                    continue
                else:
                    steps_range.append(i - 1)
            # last step case
            elif substep == 1:
                steps_range.append(i - 1)
            else:
                continue
        else:
            steps_range.append(i - 1)
else:
    steps_range = range(1, nsteps + 1)


if args.xyz:
    models = {'H': 'modelH',
              'M': 'MEDIUM',
              'HM': 'MEDIUM_HIGH',
              'HML': 'cartesian',
              'ML': 'pod',
              'HL': 'not implemented',
              'L': 'not implemented'}
    model = models[args.model] if args.model else models['HM']
    if model == 'not implemented':
        print('HL and L layers filter not implemented... using full geometry instead')
        model = 'cartesian'

    if args.numstep is not None:
        # fetch single step
        xyz_geoms = [OUT.get_step(wanted_step).geometry.to_string(geom_model=model,
                                                                  labels=True)]
    else:
        # get all geometries except first one written in INIT_DATA
        xyz_geoms = [g.to_string(geom_model=model, labels=True) for g in OUT.grep_all('geometry')][1:]
        # skip step 0 and select only converged points for IRC calculations
        if len(xyz_geoms) > 1:
           xyz_geoms = [xyz_geoms[i] for i in steps_range]

    with open('cobramm.xyz', 'w') as xyzfile:
        natoms = len(xyz_geoms[0].strip().split('\n'))
        for geom in xyz_geoms:
            xyzfile.write('{0}\n\n'.format(natoms))
            xyzfile.write(geom)
    snaptraj = 'trajectory' if wanted_step is None else 'snapshot'
    print('cobramm.xyz {0} saved in current directory'.format(snaptraj))


if args.pdb:
    pdbok = True
    if wanted_step is None:
        print('-n argument not supplied, using default value for pdb (0)')
        wanted_step = 0
    if OUT.calculation_type == 'H':
        print('H layer calculation (QM only), skipping pdb...')
        pdbok = False

    if args.model == 'HML' or args.model is None:
        pdb_geom = OUT.get_step(wanted_step).geometry.makerealcrd('tmp.crd')
        topfile = OUT.real_top
    elif args.model == 'H':
        pdb_geom = OUT.get_step(wanted_step).geometry.makemodelHcrd('tmp.crd')
        topfile = OUT.modelH_top
    else:
        pdbok = False
        print('Invalid geometry model ({0}) for pdb,'
              ' only HML and H are supported, skipping...'.format(args.model))

    if pdbok:
        # create pdb file using ambpdb
        with open('tmp.top', 'w') as tmptop:
            tmptop.write(topfile)
        ambpdb_cmd = 'ambpdb -p tmp.top -c tmp.crd'.split()
        with open('cobramm.pdb', 'w') as pdb:
            subprocess.call(ambpdb_cmd, stdout=pdb)
        print('cobramm.pdb saved in current directory')

    # clean tmp files
    try:
        os.remove('tmp.crd')
        os.remove('tmp.top')
    except OSError:
        pass

# abort if user requested mdcrd or molden for a single point
if OUT.steps == 0 and (args.mdcrd or args.molden):
    print('mdcrd and molden trajectory not available for single point. Aborting...')
    sys.exit(1)

if args.mdcrd:
    crdok = True

    if args.model == 'HML' or args.model is None:
        crd_geoms = [g.cartesian for g in OUT.grep_all('geometry')][1:]
    elif args.model == 'H':
        crd_geoms = [g.modelH for g in OUT.grep_all('geometry')][1:]
    else:
        crdok = False
        print('Invalid geometry model ({0}) for mdcrd,'
              ' only HML and H are supported, skipping...'.format(args.model))

    if crdok:
        text = 'Created with cobrammview'
        for crd_geom in crd_geoms:
            text += '\n' if text[-1] != '\n' else ''
            lingeom = np.reshape(crd_geom, (len(crd_geom[0]) * 3, ), order='F')
            for i,crd in enumerate(lingeom, start=1):
                text += '{0:>8.3f}{1}'.format(crd, '' if i%10 else '\n')
        text += '\n' if text[-1] != '\n' else ''
        with open('cobramm.mdcrd', 'w') as mdcrd:
            mdcrd.write(text)
        print('cobramm.mdcrd trajectory saved in current directory')


if args.molden:
    if args.model is not None and args.model != 'HM':
        print('Molden trajectory support only HM layers ({} requested).'
              'Proceeding with HM...')
    # geometry data
    atomnum = OUT.geometry.NatomHM
    atomlinks = OUT.geometry.NsubH
    HM_indexes = OUT.geometry.list_MEDIUM_HIGH

    # get geometries
    geometries = [g.to_string(geom_model='MEDIUM_HIGH', labels=True)
                  for g in OUT.grep_all('geometry')][1:]
    geometries = [geometries[i] for i in steps_range]
    # get charges
    CRG_reals = [c for c in OUT.grep_all('chargeMM')]
    if len(CRG_reals) > nsteps + 1:
        CRG_reals = CRG_reals[1:]
    CRG_reals = [CRG_reals[i] for i in steps_range]

    geom_final = []
    for g,c in zip(geometries, CRG_reals):
        geom_lines = g.strip().split('\n')
        step_geom_text = '{0}\n\n'.format(atomnum)
        for j in range(atomnum):
            line = geom_lines[j]
            ch = c[HM_indexes[j] - 1]
            spacer = '' if ch < 0 else ' '
            step_geom_text += line + '\t{0}{1:.8f}\n'.format(spacer, ch)
        geom_final.append(step_geom_text)

    # energies
    energies = [e for e in OUT.grep_all('E_QMMM')]
    energies = [energies[i] for i in steps_range]
    state_key = 'state' if OUT.optimization_type.startswith('mdv') else 'optstate'
    states = [s for s in OUT.grep_all(state_key)]
    states = [states[i] for i in steps_range]
    if type(energies[0]) == list:
        energies = [e[s] for e,s in zip(energies, states)]
    energies = ['{0:.8f}\n'.format(e) for e in energies]

    # Fmax
    Fmaxs = ['{0:.8f}\n'.format(f) for f in OUT.grep_all('Fmax')]
    diff = len(steps_range) - len(Fmaxs)
    if diff:
        Fmaxs = ['{0:.8f}\n'.format(0.0)] * diff + Fmaxs
    Fmaxs = [Fmaxs[i-1] for i in steps_range]

    # Frms
    Frmss = ['{0:.8f}\n'.format(f if f else 0.0) for f in OUT.grep_all('Frms')]
    diff = len(steps_range) - len(Frmss)
    if diff:
        Frmss = ['{0:.8f}\n'.format(0.0)] * diff + Frmss
    Frmss = [Frmss[i-1] for i in steps_range]

    # Dmax
    Dmaxs = ['{0:.8f}\n'.format(d if d else 0.0) for d in OUT.grep_all('Dmax')]
    diff = len(steps_range) - len(Dmaxs)
    if diff:
        Dmaxs = ['{0:.8f}\n'.format(0.0)] * diff + Dmaxs
    Dmaxs = [Dmaxs[i-1] for i in steps_range]

    # Drms
    Drmss = ['{0:.8f}\n'.format(d if d else 0.0) for d in OUT.grep_all('Drms')]
    diff = len(steps_range) - len(Drmss)
    if diff:
        Drmss = ['{0:.8f}\n'.format(0.0)] * diff + Drmss
    Drmss = [Drmss[i-1] for i in steps_range]

    # QMMM forces
    forces_QMMMM_HM = [f for f in OUT.grep_all('forces_QMMM_HM')]
    diff = len(steps_range) - len(forces_QMMMM_HM)
    if diff > 0:
        forces_QMMMM_HM = [np.zeros(np.shape(forces_QMMMM_HM[0]))] * diff + forces_QMMMM_HM
    elif diff < 0:
        forces_QMMMM_HM = forces_QMMMM_HM[1:]
    forces_QMMMM_HM = [forces_QMMMM_HM[i-1] for i in steps_range]
    forces_QMMMM_HM_final = []
    for i,force in enumerate(forces_QMMMM_HM):
        tforce = np.transpose(force)
        point_text = 'point {0}\n{1}\n'.format(steps_range[i], atomnum)
        for j in range(len(tforce) - atomlinks):
            f_list = ['{0}{1:.8f}'.format(' ' if f>=0.0 else '', f) for f in tforce[j]]
            point_text += '\t'.join(f_list) + '\n'
        forces_QMMMM_HM_final.append(point_text)


    # prepare molden text
    molden = ''

    # geometries
    molden += ' [MOLDEN FORMAT]\n [GEOMETRIES] (XYZ)\n'
    molden += ''.join(geom_final)

    # energies, forces, displacements
    molden += ' [GEOCONV]\n energy\n'
    molden += ''.join(energies)

    molden += ' max-force\n'
    molden += ''.join(Fmaxs)

    molden += ' rms-force\n'
    molden += ''.join(Frmss)

    molden += ' max-step\n'
    molden += ''.join(Dmaxs)

    molden += ' rms-step\n'
    molden += ''.join(Drmss)

    # QMMM forces
    molden += ' [FORCES]\n'
    molden += ''.join(forces_QMMMM_HM_final)

    # write molden text to cobramm.molden
    with open('cobramm.molden', 'w') as molden_file:
        molden_file.write(molden)
    print('cobramm.molden trajectory saved in current directory')
