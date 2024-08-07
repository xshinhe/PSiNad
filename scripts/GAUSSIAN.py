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

parser = argparse.ArgumentParser(description='Execute GAUSSIAN Calculation')
parser.add_argument('-d', '--directory', dest='directory', nargs='?', default='GAUSSIAN', type=str,
    help='work directory')
parser.add_argument('-i', '--input', dest='input', nargs='?', default='QM.in.GAUSSIAN', type=str,
    help='input file')
parser.add_argument('-o', '--output', dest='output', nargs='?', default='QM.out.GAUSSIAN', type=str,
    help='output file')
# args = parser.parse_args()

def qm_job(qm_data, args):
    qm_config = qm_data["qm_config"]
    mndo_config = qm_config["QM"]["GAUSSIAN"]
    natom = qm_data["natom"]
    znumber = qm_data["znumber"]
    xyz = qm_data["geom_xyz"]
    
    try:
        F = int(qm_config['QM']['GAUSSIAN']['F'])
        N = int(qm_config['QM']['GAUSSIAN']['N'])
        nciref = F # int(qm_config['QM']['MNDO']['keywords']['nciref'])
    except KeyError:
        print(format_exc())

    job_str = '{"bagel": [\n'

    # molecule part
    job_str += '{\n'
    job_str += '"title": "molecule",\n'
    job_str += '"basis": "%s",\n'%(qm_config['QM']['GAUSSIAN']['basis'])
    job_str += '"df_basis": "%s",\n'%(qm_config['QM']['GAUSSIAN']['df_basis'])
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

    if qm_config['QM']['GAUSSIAN']['method'] == 'casscf' or qm_config['QM']['GAUSSIAN']['method'] == 'caspt2':
        job_str += '{\n'
        job_str += '"title": "casscf",\n'
        job_str += '"nstate": %d,\n'%(qm_config['QM']['GAUSSIAN']['nstate']) #should >= F (always ==)
        job_str += '"nact": %d,\n'%(qm_config['QM']['GAUSSIAN']['nact'])
        # job_str += '"nopen": %d,\n'%(qm_config['QM']['GAUSSIAN']['nopen'])
        job_str += '"nclosed": %d,\n'%(qm_config['QM']['GAUSSIAN']['nclosed'])
        if 'active' in qm_config['QM']['GAUSSIAN']:
            job_str += '['
            cnt = 0
            for istate in qm_config['QM']['GAUSSIAN']['active']:
                if cnt !=0: job_str += ','
                job_str += '%d'%istate
                cnt += 1
            job_str += '],\n'
        if qm_config['QM']['GAUSSIAN']['method'] == 'caspt2':
            job_str += '"simth": {\n'
            job_str += '"method": "caspt2",\n'
            if 'ms' in qm_config['QM']['GAUSSIAN']:
                job_str += '"ms": "%s",\n'%(qm_config['QM']['GAUSSIAN']['ms'])
            if 'xms' in qm_config['QM']['GAUSSIAN']:
                job_str += '"xms": "%s",\n'%(qm_config['QM']['GAUSSIAN']['xms'])
            if 'sssr' in qm_config['QM']['GAUSSIAN']:
                job_str += '"sssr": "%s",\n'%(qm_config['QM']['GAUSSIAN']['sssr'])
            job_str += '"shift": %.3f\n'%(qm_config['QM']['GAUSSIAN']['shift'])
            job_str += '},\n'

        job_str += '"charge": %d,\n'%(qm_config['QM']['GAUSSIAN']['charge'])
        job_str += '"nspin": %d\n'%(qm_config['QM']['GAUSSIAN']['nspin'])
        job_str += '}\n'
    else:
        raise ValueError("unsupport method")

    job_str += ']\n' # end method
    job_str += '}\n' # end force part

    job_str += ']}' # end of bagel parameter

    qm_config["QM"]["env"] = {
        "input_is_ready":True,
        "generated": "QM.run.GAUSSIAN.json",
        "directory": args.directory,
        "output": args.output,
    }

    directory = qm_config['QM']['env']['directory']
    if not os.path.exists(directory): os.makedirs(directory)

    f = open(directory + '/' + qm_config['QM']['env']['generated'], 'w')
    f.write(job_str)
    f.flush()
    f.close()

    exe_str = 'cd %s && %s  %s > %s && cd -'%(
        qm_config['QM']['env']['directory'],
        qm_config['QM']['GAUSSIAN']['path'],
        qm_config['QM']['env']['generated'],
        qm_config['QM']['env']['output']
        )
    print(exe_str)
    os.system(exe_str)

    parse_result(
        qm_data,
        directory
    )

def parse_result(qm_data, log_path):
    try:
        natom = qm_data['natom']
        qm_config = qm_data['qm_config']
        F = int(qm_config['QM']['GAUSSIAN']['F'])
        N = int(qm_config['QM']['GAUSSIAN']['N'])
        nciref = F # int(qm_config['QM']['GAUSSIAN']['keywords']['nciref'])
    except KeyError:
        print(format_exc())

    ERROR_MSG = ""
    eig = np.zeros((F))
    dE  = np.zeros((F, N))
    nac = np.zeros((F, F, N))
    
    with open(log_path+'/ENERGY.out', 'r') as ifs:
        for i in range(F):
            eig[i] = float(ifs.readline())

    for i in range(F):
        with open(log_path+'/FORCE_%d.out'%i, 'r') as ifs:
            ifs.readline()
            for a in range(natom):
                dE[i,3*a:3*a+3] = ifs.readline().split()[1:]

    for i in range(F):
        for k in range(i+1, F):
            with open(log_path+'/NACME_%d_%d.out'%(i,k), 'r') as ifs:
                ifs.readline()
                for a in range(natom):
                    nac[i,k,3*a:3*a+3] = ifs.readline().split()[1:]
            nac[k,i,:] = -nac[i,k,:]
    
    # convert to au unit(no need)
    
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
