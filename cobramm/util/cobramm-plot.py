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

import sys
import os
import subprocess
import argparse
import glob
import shutil
import matplotlib.pyplot as plt
import numpy as np

# hijack $PYTHONPATH (sys.path) to find cobramm local modules
sys.path.insert(0, os.path.join(os.getenv('COBRAM_PATH'),'cobramm'))
# now we can import output, constants and cobrammenv modules
try:
    from output import Output
    import constants
    import cobrammenv
except ImportError:
    print( 'Error: cannot import cobramm module.\n'
           'Please check if $COBRAM_PATH is defined in your environment, if not run:\n'
           'export COBRAM_PATH=/path/to/COBRAMM' )
    sys.exit(1)


def sync_remote():
    """
    Retrieve cobramm.xml file from a remote-cluster mounted with autofs/automount
    This method doesn't work out of the box, users need to change variables
    to fit the cluster-storage configuration.
    """
    AUTOFS_CONF = '/etc/auto.nfs'
    SCRATCH_BASE_DIR = '/sge_scratch'
    STORAGE_BASE_DIR = '/storage2'
    JOB_FORMAT = '*ini*e*'
    user = os.getenv('USER')
    hostname = cobrammenv.getHostname()
    workDir = os.getcwd()
    jobDir = os.path.join(STORAGE_BASE_DIR,
                          user,
                          hostname,
                          workDir.split(user)[1][1:])
    jobID = glob.glob(os.path.join(jobDir, JOB_FORMAT))[0][-6:]
    print("Fetching remote results for jobid {0}".format(jobID))

    try:
        os.mkdir('remote_out')
    except OSError:
        pass

    targetDir = os.path.join(workDir, 'remote_out')
    targetFile = os.path.join(targetDir, 'cobramm.xml')

    with open(AUTOFS_CONF) as cnf:
        nodes = [line.split()[0] for line in cnf.readlines()[:-1]]
    # trigger auto-mount and try to fetch output file
    for node in nodes:
        try:
            os.chdir(os.path.join(SCRATCH_BASE_DIR, node))
            filename = os.path.join('{0}_cobram_{1}'.format(user, jobID),
                                    'cobramm.xml')
            try:
                shutil.copy(filename, targetFile)
                break
            except IOError:
                pass
        except OSError:
            pass

    os.chdir(targetDir)


parser = argparse.ArgumentParser(description='Cobramm data extractor and plot utility')
parser.add_argument('-e', '--energymd', help='Total, potential and kinetic energies in eV', action='store_true')
parser.add_argument('-eq', '--energyqm', help='QM potential energies for all states in A.U.', action='store_true')
parser.add_argument('-em', '--energyqmmm', help='QM-MM energies for all states in A.U.', action='store_true')
parser.add_argument('-sh', '--surfacehopping', help='QM-MM energies (A.U.) along a SH trajectory ', action='store_true')
parser.add_argument('-o', '--occupation', help='trajectory states occupations', action='store_true')
parser.add_argument('-opt', '--optimization', help='energies, forces, displacements in optimizations', action='store_true')
parser.add_argument('-all', '--all', help='extract all fields', action='store_true')
parser.add_argument('-p', '--plot', help='automatic plot using matplotlib', action='store_true')
parser.add_argument('-r', '--remote', help='fetch data from remote cluster', action='store_true')

args = parser.parse_args()

if args.remote:
    # fetch remote-cluster output file
    sync_remote()

# read output file
OUT = Output(parse=True)

# init req variable to None
req = None

if args.all:
    if OUT.optimization_type.startswith('mdv'):
        choices = ['-e', '-eq', '-em', '-o']
    elif OUT.steps == 0 or OUT.optimization_type.startswith('freq'):
        choices = ['-eq', '-em']
    else:
        choices = ['-eq', '-em', '-opt']
    for r in choices:
        if args.plot:
            _code = subprocess.call(['cobramm-plot.py', '-p', r])
        else:
            _code = subprocess.call(['cobramm-plot.py', r])
        if _code:
            print('Process terminated with exit code '
                  '{0} for "{1}"'.format(_code, r))
    sys.exit(0)
else:
    choices = [k for k, v in vars(args).items() if v and k not in ['plot', 'remote']]
    # more than one argument
    if len(choices) > 1:
        for r in choices:
            if args.plot:
                _code = subprocess.call(['cobramm-plot.py', '-p', '--{0}'.format(r)])
            else:
                _code = subprocess.call(['cobramm-plot.py', '--{0}'.format(r)])
            if _code:
                print('Process terminated with exit code '
                      '{0} for "{1}"'.format(_code, r))
        sys.exit(0)
    elif len(choices) == 0 and args.remote:
        print('cobramm.xml saved in remote_out directory')
        sys.exit(0)
    elif len(choices) == 0:
        subprocess.call(['cobramm-plot.py', '-h'])
        sys.exit(0)
    else:
        req = choices[0]

graph = plt.subplot('111')
plt_list = []
plt_data = []

Epot0 = 0.0

# compute chosen property
# Etot, Epot, Ekin
if args.energymd:
    if not OUT.optimization_type.startswith('mdv'):
        print('Calculation is not a MD, aborting...')
        sys.exit(1)

    _Etot = [e * constants.Hartree2eV for e in OUT.grep_all('MD_Etot')]
    _Epot = [e * constants.Hartree2eV for e in OUT.grep_all('MD_Epot')]
    Ekin = [e * constants.Hartree2eV for e in OUT.grep_all('MD_Ekin')]
    times = [t / constants.fs2au for t in OUT.grep_all('time')]

    # sanity checks
    diff1 = (OUT.steps + 1) - len(_Etot)
    diff2 = (OUT.steps + 1) - len(_Epot)
    diff3 = (OUT.steps + 1) - len(Ekin)
    diff4 = (OUT.steps + 1) - len(times)
    if diff1:
        _Etot = [0.0] * diff1 + _Etot
    if diff2:
        _Epot = [0.0] * diff2 + _Epot
    if diff3:
        Ekin = [0.0] * diff3 + Ekin
    if diff4:
        times = [0.0] * diff4 + times

    Epot0 = _Epot[0]
    text = '#  time  ' + 'Etot'.rjust(16) + 'Epot'.rjust(16) + 'Ekin'.rjust(16) + '\n'
    Epot = [e - Epot0 for e in _Epot]
    Etot = [e - Epot0 for e in _Etot]
    plt_data.append(times)
    plt_data.append(Etot)
    plt_data.append(Epot)
    plt_data.append(Ekin)

    for i in range(len(times)):
        text += '{0:>7.3f}\t\t{1:>12.8f}\t{2:>12.8f}\t{3:>12.8f}\n'.format(times[i], Etot[i], Epot[i], Ekin[i])

    if args.plot:
        plt_1 = graph.plot(plt_data[0], plt_data[1], label='Etot')
        plt_2 = graph.plot(plt_data[0], plt_data[2], label='Epot')
        plt_3 = graph.plot(plt_data[0], plt_data[3], label='Ekin')

        graph.legend(plt_1 + plt_2 + plt_3, [p.get_label() for p in plt_1 + plt_2 + plt_3])
        graph.set_xlabel('time / fs')
        graph.set_ylabel('energy / eV')
        graph.set_title('ETot, EPot, EKin')


# QM/QM-MM energy
elif args.energyqm or args.energyqmmm:
    MD = OUT.optimization_type.startswith('mdv')
    # store energies in _EN variable
    _EN = OUT.grep_all('E_QM') if args.energyqm else OUT.grep_all('E_QMMM')
    if not type(_EN[0]) == list:
        _EN = [[e] for e in _EN]

    if MD:
        times = [0.0] + [t / constants.fs2au for t in OUT.grep_all('time')]
        t_text = ['{0:>7.3f}'.format(t) for t in times]
    else:
        times = OUT.grep_all('step')
        t_text = ['{0:>5d}'.format(s) for s in times]

    # init text
    text = '#  time  ' if MD else '#  step  '
    for i in range(1, len(_EN[0])+1):
        text += 'state_{0}'.format(i).rjust(16)
    text += '\n'

    # populate plt_data and text
    plt_data.append(times)
    for i in range(len(_EN[0])):
        plt_data.append([e[i] for e in _EN])

    for i in range(len(times)):
        _s = ''.join(['{0:>16.8f}'.format(s) for s in _EN[i]])
        text += '{0}\t\t{1}\n'.format(t_text[i], _s)

    # plot data
    if args.plot:
        if len(plt_data) > 2:
            label = 'state_{0}'
        else:
            label = ''
        for i in range(1, len(plt_data)):
            g = graph.plot(plt_data[0], plt_data[i], label='state_{0}'.format(i))
            plt_list.append(g)
        legend_list = plt_list[0]
        for i in range(1, len(plt_list)):
            legend_list += plt_list[i]
        graph.legend(legend_list, [p.get_label() for p in legend_list])
        graph.set_xlabel('time / fs' if MD else 'optimization step')
        graph.set_ylabel('energy / A.U.')
        graph.set_title('{0} energy'.format('QM' if args.energyqm else 'QM-MM'))


# surface hopping trajectory
elif args.surfacehopping:
    if not OUT.optimization_type.startswith('mdv'):
        print('Calculation is not a MD, aborting...')
        sys.exit(1)

    # store energies in _EN variable
    _EN = OUT.grep_all('E_QMMM')
    if not type(_EN[0]) == list:
        _EN = [[e] for e in _EN]
    # extract the energy of the current state of the trajectory
    states = OUT.grep_all('state')
    _Estate = [_EN[i][os] for i, os in enumerate(states)]

    times = [0.0] + [t / constants.fs2au for t in OUT.grep_all('time')]
    t_text = ['{0:>7.3f}'.format(t) for t in times]

    # init text
    text = '#  time  '
    text += 'current_state'.rjust(16)
    for i in range(1, len(_EN[0])+1):
        text += 'state_{0}'.format(i).rjust(16)
    text += '\n'

    # populate plt_data and text
    plt_data.append(times)
    plt_data.append(_Estate)
    for i in range(len(_EN[0])):
        plt_data.append([e[i] for e in _EN])

    for i in range(len(times)):
        _s = ''.join(['{0:>16.8f}'.format(s) for s in _EN[i]])
        text += '{0}\t\t{1:>16.8f}{2}\n'.format(t_text[i], _Estate[i], _s)

    # plot data
    if args.plot:
        g = graph.plot(plt_data[0], plt_data[1], label='current state', lw=3.0, ls="--")
        plt_list.append(g)
        for i in range(2, len(plt_data)):
            g = graph.plot(plt_data[0], plt_data[i], label='state_{0}'.format(i-1))
            plt_list.append(g)
        legend_list = plt_list[0]
        for i in range(1, len(plt_list)):
            legend_list += plt_list[i]
        graph.legend(legend_list, [p.get_label() for p in legend_list])
        graph.set_xlabel('time / fs')
        graph.set_ylabel('energy / A.U.')
        graph.set_title('QM/MM energy')


# states occupation
elif args.occupation:
    if not OUT.optimization_type.startswith('mdv'):
        print('Calculation is not a MD, aborting...')
        sys.exit(1)

    states = OUT.grep_all('state')
    nstates = len(OUT.get_step(0).E_QM)
    times = [0.0] + [t / constants.fs2au for t in OUT.grep_all('time')]
    _occ = []
    for s in states:
        _occ.append([1.0 if i == s else 0.0 for i in range(nstates)])

    occ = np.transpose(_occ)

    # init text
    text = '#  time        '
    for i in range(nstates):
        text += 'state_{0}'.format(i+1).rjust(12)
    text += '\n'

    # populate plt_data and text
    plt_data.append(times)
    t_text = ['{0:>7.3f}'.format(t) for t in times]

    for e in occ:
        plt_data.append(e)

    for i in range(len(_occ)):
        _s = ''.join(['{0:>13.4f}'.format(s) for s in _occ[i]])
        text += '{0}\t\t{1}\n'.format(t_text[i], _s)

    # plot data
    if args.plot:
        for i in range(1, len(plt_data)):
            g = graph.plot(plt_data[0], plt_data[i], label='state_{0}'.format(i))
            plt_list.append(g)
        legend_list = plt_list[0]
        for i in range(1, len(plt_list)):
            legend_list += plt_list[i]
        graph.legend(legend_list, [p.get_label() for p in legend_list])
        graph.set_xlabel('time / fs' if OUT.calculation_type.startswith('mdv') else 'optimization step')
        graph.set_ylabel('states occupation')
        graph.set_title('Electronic states occupation')


# optimization data
elif args.optimization:
    if OUT.optimization_type.startswith('mdv') or OUT.optimization_type.startswith('freq'):
        print('Calculation is not an optimization, aborting...')
        sys.exit(1)
    elif OUT.steps == 0:
        print('Single step calculation, optimization data not available, aborting...')
        sys.exit(1)

    k_list = ['E_QMMM', 'Fmax', 'Frms', 'Dmax', 'Drms']

    _E = OUT.grep_all('E_QMMM')[1:]
    optstates = OUT.grep_all('optstate')[1:]
    steps = OUT.grep_all('step')[1:]
    if type(_E[0]) == list:
        try:
           _E = [_E[i][os] for i,os in enumerate(optstates)]
        except IndexError:
           _E = [_E[i][os-1] for i,os in enumerate(optstates)]
    elif type(_E[0]) == float:
        pass

    k_data = { 'E_QMMM': _E,
               'Fmax': OUT.grep_all('Fmax'),
               'Frms': OUT.grep_all('Frms'),
               'Dmax': OUT.grep_all('Dmax'),
               'Drms': OUT.grep_all('Drms')}

    # sanity checks
    for k in k_list[1:]:
        diff = len(_E) - len(k_data[k])
        if diff:
            k_data[k] = [0.0] * diff + k_data[k]

    # init text
    text = '#  step  '
    for s in k_list:
        text += '{0}'.format(s).rjust(16)
    text += '\n'

    # populate plt_data and text
    plt_data.append(steps)
    for k in k_list:
        plt_data.append(k_data[k])

    DF = np.asarray([k_data[k] for k in k_list])
    DF = np.transpose(DF)

    for i in range(len(DF)):
        _s = '   '.join(['{0:>13.8f}'.format(v) for v in DF[i]])
        text += '{0:>5d}\t\t{1}\n'.format(steps[i], _s)

    # plot data
    if args.plot:
        g = graph.plot(plt_data[0], plt_data[1], label='E_QMMM_state_{0}'.format(optstates[0]))
        plt_list.append(g)
        graph2 = graph.twinx()
        for i in range(2, len(plt_data)):
            g = graph2.plot(plt_data[0], plt_data[i], label=k_list[i-1], ls=':')
            plt_list.append(g)

        legend_list = plt_list[0]
        for i in range(1, len(plt_list)):
            legend_list += plt_list[i]
        graph.legend(legend_list, [p.get_label() for p in legend_list])
        graph.set_xlabel('optimization step')
        graph.set_ylabel('energy / A.U.')
        graph2.set_yscale('log')
        graph2.set_ylabel('')
        graph.set_title('Optimization data')

# write dat file
with open('{0}.dat'.format(req), 'w') as f:
    f.write(text)

# automatic plot
if args.plot:
    graph.set_xticks(plt_data[0])
    plt.show()
