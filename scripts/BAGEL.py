import argparse
import datetime
import math
import os
import re
import shutil
import subprocess
import sys
import time
import toml
from copy import deepcopy
from multiprocessing import Pool
from pprint import pprint
from traceback import format_exc
from socket import gethostname
import numpy as np
import QMutils
import typing

parser = argparse.ArgumentParser(description='Execute BAGEL Calculation')
parser.add_argument('-d', '--directory', dest='directory', nargs='?', default='BAGEL', type=str,
    help='work directory')
parser.add_argument('-i', '--input', dest='input', nargs='?', default='QM.in.BAGEL', type=str,
    help='input file')
parser.add_argument('-o', '--output', dest='output', nargs='?', default='QM.out.BAGEL', type=str,
    help='output file')
# args = parser.parse_args()

def qm_job(qm_data, args):
    qm_config = qm_data["qm_config"]
    mndo_config = qm_config["QM"]["BAGEL"]
    natom = qm_data["natom"]
    znumber = qm_data["znumber"]
    xyz = qm_data["geom_xyz"]
    
    try:
        F = int(qm_config['QM']['BAGEL']['F'])
        N = int(qm_config['QM']['BAGEL']['N'])
        nciref = F # int(qm_config['QM']['MNDO']['keywords']['nciref'])
    except KeyError:
        print(format_exc())

    job_str = '{"bagel": [\n'

    # molecule part
    job_str += '{\n'
    job_str += '"title": "molecule",\n'
    job_str += '"basis": "%s",\n'%(qm_config['QM']['BAGEL']['basis'])
    job_str += '"df_basis": "%s",\n'%(qm_config['QM']['BAGEL']['df_basis'])
    job_str += '"angstrom": true,\n'
    job_str += '"geometry": [\n'
    for i in range(natom):
        end = ',' if i!=natom-1 else ''
        job_str += '{"atom": "%s", "xyz": [%12.8e, %12.8e, %12.8e]}%s\n'%(
                QMutils.element_list[znumber[i]][0], xyz[i][1], xyz[i][2], xyz[i][3], end)
    job_str += ']\n' # end geometry
    job_str += '},\n' # end molecule part

    # force part
    job_str += '{\n'
    job_str += '"title": "forces",\n'
    job_str += '"dipole": true,\n'
    job_str += '"export": true,\n'
    job_str += '"nproc": 1,\n'
    job_str += '"grads": [\n'
    for i in range(F):
        end = ',' if i!=F-1 else ''
        job_str += '{"title": "force", "target": %d},\n'%i
    for i in range(F):
        for k in range(i+1,F):
            end = '' if i==k-1 and k==F-1 else ','
            job_str += '{"title": "nacme", "target": %d, "target2": %d, "nacmtype": "full"}%s\n'%(i,k,end)
    job_str += '],\n' # end grads
    # force methods
    job_str += '"method": ['

    if qm_config['QM']['BAGEL']['method'] == 'casscf' or qm_config['QM']['BAGEL']['method'] == 'caspt2':
        job_str += '{\n'
        job_str += '"title": "casscf",\n'
        job_str += '"nstate": %d,\n'%(qm_config['QM']['BAGEL']['nstate'])
        job_str += '"nact": %d,\n'%(qm_config['QM']['BAGEL']['nact'])
        # job_str += '"nopen": %d,\n'%(qm_config['QM']['BAGEL']['nopen'])
        job_str += '"nclosed": %d,\n'%(qm_config['QM']['BAGEL']['nclosed'])
        if 'active' in qm_config['QM']['BAGEL']:
            job_str += '['
            cnt = 0
            for istate in qm_config['QM']['BAGEL']['active']:
                if cnt !=0: job_str += ','
                job_str += '%d'%istate
                cnt += 1
            job_str += '],\n'
        job_str += '"charge": %d,\n'%(qm_config['QM']['BAGEL']['charge'])
        job_str += '"nspin": %d\n'%(qm_config['QM']['BAGEL']['nspin'])
        job_str += '}\n'
    
    if qm_config['QM']['BAGEL']['method'] == 'caspt2':
        job_str += ',{\n'
        job_str += '"title": "smith",\n'
        job_str += '"method": "caspt2",\n'
        if 'ms' in qm_config['QM']['BAGEL']:
            job_str += '"ms": "%s",\n'%(qm_config['QM']['BAGEL']['ms'])
        if 'xms' in qm_config['QM']['BAGEL']:
            job_str += '"xms": "%s",\n'%(qm_config['QM']['BAGEL']['xms'])
        if 'sssr' in qm_config['QM']['BAGEL']:
            job_str += '"sssr": "%s",\n'%(qm_config['QM']['BAGEL']['sssr'])
        job_str += '"shift": %.3f\n'%(qm_config['QM']['BAGEL']['shift'])
        job_str += '}\n'

    job_str += ']\n' # end method
    job_str += '}\n' # end force part

    job_str += ']}' # end of bagel parameter

    qm_config["QM"]["env"] = {
        "input_is_ready":True,
        "generated": "QM.run.BAGEL.json",
        "directory": args.directory,
        "output": args.output,
    }

    directory = qm_config['QM']['env']['directory']
    if not os.path.exists(directory): os.makedirs(directory)

    f = open(directory + '/' + qm_config['QM']['env']['generated'], 'w')
    f.write(job_str)
    f.flush()
    f.close()

    exe_str = '%s  %s > %s'%(
        qm_config['QM']['BAGEL']['path'],
        directory + '/' + qm_config['QM']['env']['generated'],
        directory + '/' + qm_config['QM']['env']['output']
        )
    print(exe_str)
    os.system(exe_str)

    parse_result(
        qm_data,
        directory + '/' + qm_config['QM']['env']['output']
    )

def parse_result(qm_data, log_file):
    try:
        natom = qm_data['natom']
        qm_config = qm_data['qm_config']
        F = int(qm_config['QM']['BAGEL']['F'])
        N = int(qm_config['QM']['BAGEL']['N'])
        nciref = F # int(qm_config['QM']['BAGEL']['keywords']['nciref'])
    except KeyError:
        print(format_exc())

    stat = -1
    istate_force = 0
    istate_force_meet = 0
    ERROR_MSG = ""
    f_r = np.zeros((nciref))
    f_p = np.zeros((nciref))
    f_rp = np.zeros((nciref))
    eig = np.zeros((F))
    dE  = np.zeros((F, N))
    nac = np.zeros((F, F, N))
    
    with open(log_file, 'r') as ifs:
        for eachline in ifs:
            if "ERROR" in eachline or "UNABLE" in eachline:
                ERROR_MSG += eachline

            if "Properties of transitions   1 -> #" in eachline:
                ifs.readline()  # blankline
                ifs.readline()  # headline
                for i in range(1, nciref):
                    f_r[i], f_p[i], f_rp[i] = map(float, ifs.readline().split()[-3:])
                stat = 0

            elif "SUMMARY OF MULTIPLE CI CALCULATIONS" in eachline:
                for _ in range(4):
                    ifs.readline()
                for i in range(F):
                    eig[i] = float(ifs.readline().split()[2])
                stat = 0

            elif "CI CALCULATION FOR STATE:" in eachline:
                istate_force = int(eachline.split()[-1]) - 1

            elif "GRADIENTS (KCAL/(MOL*ANGSTROM))" in eachline and istate_force_meet == istate_force:
                istate_force_meet += 1
                for _ in range(3):
                    ifs.readline()
                for i in range(natom):
                    dE[istate_force,3*i:3*i+3] = ifs.readline().split()[5:]

            elif "CI CALCULATION FOR INTERSTATE COUPLING OF STATES:" in eachline:
                istate, jstate = map(int, eachline.split()[-2:])
                istate -= 1; jstate -= 1
                if istate < F and jstate < F and istate != jstate:
                    while True:
                        if "GRADIENTS (KCAL/(MOL*ANGSTROM))" in ifs.readline():
                            for _ in range(3):
                                ifs.readline()
                            for i in range(natom):
                                nac[istate,jstate,3*i:3*i+3] = ifs.readline().split()[5:]
                            nac[jstate,istate,:] = -nac[istate,jstate,:]
                            break
                    stat = 2

    if stat != 2:
        print(f"fail in calling BAGEL! {ERROR_MSG}")

    # convert to au unit
    eig /= QMutils.au_2_kcal_1mea
    dE /= QMutils.au_2_kcal_1mea_per_ang
    nac /= 1.0e0 / QMutils.au_2_ang

    qmout = QMutils.QMout(natom=natom, 
        energy=eig,
        gradient=dE.T,
        nac=np.einsum('ikj->jik', nac)
    )

    f = open(qm_config['QM']['env']['directory'] + '/energy.dat', 'w')
    for i in range(F):
        f.write('{: 12.8e}\n'.format(eig[i]))
    f.close()

    f = open(qm_config['QM']['env']['directory'] + '/gradient.dat', 'w')
    for j in range(N):
        for i in range(F):
            f.write('{: 12.8e} '.format(dE[i,j]))
        f.write('\n')
    f.close()

    f = open(qm_config['QM']['env']['directory'] + '/nacv.dat', 'w')
    for j in range(N):
        for i in range(F):
            for k in range(F):
                f.write('{: 12.8e} '.format(nac[i,k,j]))
        f.write('\n')
    f.close()
    
    #pprint(qmout)
    return qmout

def main():
    return

if __name__ == '__main__':
    args = parser.parse_args()
    pprint(args)
    qm_data_in = QMutils.parseQMinput(args.input)
    qm_job(qm_data_in, args)

    main()
