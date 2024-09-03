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
from multiprocessing import Pool, Process
from pprint import pprint
from traceback import format_exc
from socket import gethostname
import numpy as np
import QMutils
import typing

parser = argparse.ArgumentParser(description='Execute BAGEL Calculation')
parser.add_argument('-d', '--directory', dest='directory', nargs='?', default='.', type=str,
    help='work directory')
parser.add_argument('-i', '--input', dest='input', nargs='?', default='QM.in.BAGEL', type=str,
    help='input file')
parser.add_argument('-t', '--task', dest='task', nargs='?', default='nad', type=str,
    help='task type')
parser.add_argument('-o', '--output', dest='output', nargs='?', default='QM.out.BAGEL', type=str,
    help='output file')
# args = parser.parse_args()

def _get_errormsg_from_lines(lines: list[str]):
    find = False
    ERROR_MSG = ''
    for eachline in lines:
        if 'ERROR' in eachline or 'STOP' in eachline: # or warnings
            find = True
            ERROR_MSG += eachline
    return find, ERROR_MSG

def _get_unable_from_lines(lines: list[str]):
    find = False
    SOME_MSG = ''
    for eachline in lines:
        if "UNABLE" in eachline:
            find = True
            SOME_MSG += eachline
        if "ACHIEVED" in eachline:
            find = False
            SOME_MSG += eachline
    return find, SOME_MSG

def _get_oscstrength_from_lines(lines: list[str]):
    find = False
    f_strength = [0,]
    for i in range(len(lines)):
        if "* Oscillator strength for transition between 0" in lines[i]:
            terms = lines[i].strip().split()
            f_strength += [float(terms[9])]
            find = True
    return find, np.array(f_strength)

def _get_energy(log_path:str, F:int):
    find = False
    energies = []
    if os.path.exists(log_path + '/ENERGY.out'):
        find = True
        with open(log_path+'/ENERGY.out', 'r') as ifs:
            for i in range(F):
                energies += [float(ifs.readline())]
        os.rename(log_path + '/ENERGY.out', log_path + '/ENERGY-OLD.out')

    return find, np.array(energies)

def _get_gradient(log_path:str, N:int, F:int):
    find = {}
    natom = N //3 
    gradient = {}
    for i in range(F):
        find[i] = False
        grad_local = []
        if os.path.exists(log_path+'/FORCE_%d.out'%i):
            find[i] = True
            with open(log_path+'/FORCE_%d.out'%i, 'r') as ifs:
                ifs.readline()
                for a in range(natom):
                    grad_local += [ np.array(ifs.readline().split()[1:]).astype(np.float64) ]
            gradient['%d'%i] = np.array(grad_local)
            os.rename(log_path+'/FORCE_%d.out'%i, log_path+'/FORCE_%d-OLD.out'%i)
    allfind = True
    for i in range(F):
        if find[i] == False:
            allfind = False
            break
    return allfind, gradient

def _get_nac(log_path:str, N:int, F:int):
    find = {}
    natom = N //3 
    nac = {}
    for i in range(F):
        for k in range(i+1, F):
            find['%d-%d'%(i,k)] = False
            nac_local = []
            if os.path.exists(log_path+'/NACME_%d_%d.out'%(i,k)):
                find['%d-%d'%(i,k)] = True
                with open(log_path+'/NACME_%d_%d.out'%(i,k), 'r') as ifs:
                    ifs.readline()
                    for a in range(natom):
                        nac_local += [ np.array(ifs.readline().split()[1:]).astype(np.float64) ]
                nac['%d-%d'%(i,k)] = np.array(nac_local)
                nac['%d-%d'%(k,i)] =-np.array(nac_local)
                os.rename(log_path+'/NACME_%d_%d.out'%(i,k), log_path+'/NACME_%d_%d-OLD.out'%(i,k))
    allfind = True
    for i in range(F):
        for k in range(i+1,F):
            if find["%d-%d"%(i,k)] == False:
                allfind = False
                break
    return allfind, nac


def _get_hessian_from_lines(lines: list[str]):
    find = {}
    hess_dict = {}
    geom_begin_idx = -1
    resforce_begin_idx = -1
    symwhess_begin_idx = -1
    whess_eig_begin_idx = -1
    whess_vec_begin_idx = -1
    for i in range(len(lines)):
        if "*** Geometry ***" in lines[i] and geom_begin_idx == -1:
            geom_begin_idx = i
            find['x0'] = True
        if '++ Symmetrized Mass Weighted Hessian ++' in lines[i] and symwhess_begin_idx == -1:
            symwhess_begin_idx = i
            find['hess'] = True
        if 'Mass Weighted Hessian Eigenvalues' in lines[i] and whess_eig_begin_idx == -1:
            find['w'] = True # BAD!!!
            whess_eig_begin_idx = i
        if '* Vibrational frequencies, IR intensities, and corresponding cartesian eigenvectors' in lines[i] and whess_vec_begin_idx == -1:
            find['Tmod'] = True
            whess_vec_begin_idx = i

    if 'x0' in find:
        sym = []
        xyz = []
        for i in range(geom_begin_idx+1, len(lines)):
            ll = lines[i].strip()
            if ll == '': continue
            if ll[0] != '{' and 'first_geom' not in find:
                find['first_geom'] = True
                continue
            if ll[0] != '{' and 'first_geom' in find:
                break
            sym += [ll.split(',')[0].split(':')[1].strip().split('"')[1]] 
            xyz += [ np.array(ll.split('[')[1].split(']')[0].split(',')).astype(np.float64) ]
        xyz = np.array(xyz)
        hess_dict['sym'] = sym
        hess_dict['x0'] = xyz

        mass = []
        mass = np.zeros(3*len(sym))
        mdict = {'C':12.01, 'H':1.008, 'O': 16.00}
        for k in range(len(mass)):
            mass[k] = mdict[sym[k//3]] * 1850
        hess_dict['mass'] = mass

        
    if 'hess' in find:
        natom = len(xyz)
        N = 3*natom
        hess = np.zeros((N,N))
        i = symwhess_begin_idx+1
        read_col_num = 0
        while i < len(lines) and read_col_num < N:
            ll = lines[i].strip()
            if ll == '': 
                i += 1
                continue
            col_idx = np.array(ll.split()).astype(np.int32)
            read_col_num += len(col_idx)
            for k in range(1,1+N):
                data = np.array(lines[i+k].strip().split())
                row_idx = int(data[0])
                data_v = np.array(data[1:]).astype(np.float64)
                hess[row_idx,col_idx] = data_v    
            i += 1+N
        hess_dict['hess'] = hess

    if 'w' in find and False:
        w = np.array(lines[whess_eig_begin_idx+1].split()).astype(np.float64)
        sign = np.sign(w)
        w = np.sqrt(np.abs(w)) * sign
        hess_dict['w'] = w

    if 'Tmod' in find:
        find['w'] = True
        Tmod = np.zeros((N,N))
        wmod = np.zeros((N))
        i = whess_vec_begin_idx+1
        read_col_num = 0
        while i < len(lines) and read_col_num < N:
            ll = lines[i].strip()
            if ll == '': 
                i += 1
                continue
            col_idx = np.array(ll.split()).astype(np.int32)
            read_col_num += len(col_idx)
            wmod[col_idx] = np.array(lines[i+1].strip().split()[2:])
            for k in range(6,6+N):
                data = np.array(lines[i+k].strip().split())
                row_idx = int(data[0])
                data_v = np.array(data[1:]).astype(np.float64)
                Tmod[row_idx,col_idx] = data_v
            i += 6+N
        hess_dict['Tmod'] = Tmod
        hess_dict['w'] = wmod / QMutils.au_2_wn

    return find != {}, hess_dict


def qm_job(qm_data, args):
    qm_config = qm_data["qm_config"]
    mndo_config = qm_config["QM"]["BAGEL"]
    natom = qm_data["natom"]
    znumber = qm_data["znumber"]
    xyz = qm_data["geom_xyz"]
    
    try:
        F = int(qm_config['QM']['BAGEL']['F'])
        N = int(qm_config['QM']['BAGEL']['N'])
    except KeyError:
        print(format_exc())
    nthread = 1
    if 'nthread' in qm_config['QM']:
        nthread = qm_config['QM']['nthread']

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

    level = 0
    try:
        print('current task = ', args.task)
        if args.task.isnumeric():
            print('current task = ', args.task)
            level = int(args.task)
            if level > 0:
                print("failed and try by level = ", level)
    except KeyError:
        print(format_exc())
    if os.path.exists(args.directory + '/laststep.archive') and level ==0:
        job_str += '{\n'
        job_str += '"title": "load_ref",\n'
        job_str += '"file": "laststep",\n'
        job_str += '"continue_geom": false\n'
        job_str += '},\n'

    # force part
    job_str += '{\n'
    job_str += '"title": "forces",\n'
    job_str += '"dipole": true,\n'
    job_str += '"export": true,\n'
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
        job_str += '"nstate": %d,\n'%(qm_config['QM']['BAGEL']['nstate']) #should >= F (always ==)
        job_str += '"nact": %d,\n'%(qm_config['QM']['BAGEL']['nact'])
        # job_str += '"nopen": %d,\n'%(qm_config['QM']['BAGEL']['nopen'])
        job_str += '"nclosed": %d,\n'%(qm_config['QM']['BAGEL']['nclosed'])
        if 'active' in qm_config['QM']['BAGEL']:
            job_str += '"active" :'
            job_str += '['
            cnt = 0
            for istate in qm_config['QM']['BAGEL']['active']:
                if cnt !=0: job_str += ','
                job_str += '%d'%istate
                cnt += 1
            job_str += '],\n'
        if qm_config['QM']['BAGEL']['method'] == 'caspt2':
            job_str += '"simth": {\n'
            job_str += '"method": "caspt2",\n'
            if 'ms' in qm_config['QM']['BAGEL']:
                job_str += '"ms": "%s",\n'%(qm_config['QM']['BAGEL']['ms'])
            if 'xms' in qm_config['QM']['BAGEL']:
                job_str += '"xms": "%s",\n'%(qm_config['QM']['BAGEL']['xms'])
            if 'sssr' in qm_config['QM']['BAGEL']:
                job_str += '"sssr": "%s",\n'%(qm_config['QM']['BAGEL']['sssr'])
            job_str += '"shift": %.3f\n'%(qm_config['QM']['BAGEL']['shift'])
            job_str += '},\n'

        job_str += '"charge": %d,\n'%(qm_config['QM']['BAGEL']['charge'])
        job_str += '"nspin": %d\n'%(qm_config['QM']['BAGEL']['nspin'])
        job_str += '}\n'
    else:
        raise ValueError("unsupport method")

    job_str += ']\n' # end method
    job_str += '},\n' # end force part

    # save part
    job_str += '{\n'
    job_str += '"title": "save_ref",\n'
    job_str += '"file": "laststep"\n'
    job_str += '}\n'

    job_str += ']}' # end of bagel parameter

    qm_config["QM"]["env"] = {
        "input_is_ready":True,
        "generated": "QM.run.BAGEL.json",
        "directory": args.directory,
        "output": args.output,
    }
    #print('config managed by kidsqmmm module:')
    #pprint(qm_config)

    directory = qm_config['QM']['env']['directory']
    if not os.path.exists(directory): os.makedirs(directory)

    f = open(directory + '/' + qm_config['QM']['env']['generated'], 'w')
    f.write(job_str)
    f.flush()
    f.close()

    #print("try to do BAGEL job with parameter as:")
    #print(job_str)

    current_path = os.path.abspath(os.getcwd())
    exe_str = 'export BAGEL_NUM_THREADS=%d && cd %s && %s  %s > %s && cd %s'%(
        nthread,
        qm_config['QM']['env']['directory'],
        qm_config['QM']['BAGEL']['path'],
        qm_config['QM']['env']['generated'],
        qm_config['QM']['env']['output'],
        current_path
    )
    stat = os.system(exe_str)
    
    #def run_command():
    #    os.system(exe_str)
    #p = Process(target=run_command)
    #p.start()
    #stat = p.join()

    #print('os executation status: ', stat)
    #print('now try to parse somthing')
    parse_result(
        qm_data,
        directory
    )

def parse_result(qm_data, log_path):
    try:
        natom = qm_data['natom']
        qm_config = qm_data['qm_config']
        F = int(qm_config['QM']['BAGEL']['F'])
        N = int(qm_config['QM']['BAGEL']['N'])
        output = qm_config['QM']['env']['output']
    except KeyError:
        print(format_exc())

    ERROR_MSG = ""
    eig = np.zeros((F))
    dE  = np.zeros((N, F))
    nac = np.zeros((N, F, F))

    fstrength = np.zeros((F))
    
    lines = open(log_path +'/'+ output, 'r').readlines()
    #pprint(lines)
    finde, ERROR_MSG = _get_errormsg_from_lines(lines)
    findu, UNABLE_MSG = _get_unable_from_lines(lines)
    findf, fstrength = _get_oscstrength_from_lines(lines)

    stat = 0
    if finde: stat = 1
    try:
        find0, energy = _get_energy(log_path, F)
        find1, gradient = _get_gradient(log_path, N, F)
        findc, nacvector = _get_nac(log_path, N, F)
        eig = energy # already in au
        for i in range(F):
            dE[:,i] = gradient['%d'%i].flatten() # already in au
        for i in range(F):
            for k in range(i+1,F):
                nac[:,i,k] = nacvector['%d-%d'%(i,k)].flatten() # already in au
                nac[:,k,i] = nacvector['%d-%d'%(k,i)].flatten() # already in au
    except (IOError, KeyError):
        print(format_exc())
        print('CANNOT PARSE FILE. try to print ERROR MSG')
        print(ERROR_MSG)
        stat = 1

    qmout = QMutils.QMout(natom=natom, 
        energy=eig,
        gradient=dE,
        nac=nac
    )

    # anyway we alway write the interface.ds
    try:
        f = open(qm_config['QM']['env']['directory'] + '/interface.ds', 'w')
        f.write('interface.stat\n')
        f.write(f'kids_int {1}\n')
        f.write(f'{stat}\n\n')

        if find0:
            f.write('interface.eig\n')
            f.write(f'kids_real {F}\n')
            for i in range(F):
                f.write('{: 12.8e}\n'.format(eig[i]))
            f.write('\n')

        if find1:
            f.write('interface.dE\n')
            f.write(f'kids_real {N*F}\n')
            for j in range(N):
                for i in range(F):
                    f.write('{: 12.8e} '.format(dE[j,i]))
                f.write('\n')
            f.write('\n')

        if findc:
            f.write('interface.nac\n')
            f.write(f'kids_real {N*F*F}\n')
            for j in range(N):
                for i in range(F):
                    for k in range(F):
                        f.write('{: 12.8e} '.format(nac[j,i,k]))
                f.write('\n')
            f.write('\n')

        if findf:
            f.write('interface.strength\n')
            f.write(f'kids_real {F}\n')
            for i in range(len(fstrength)):
                f.write('{: 12.8e} '.format(fstrength[i]))
            f.write('\n')
            f.close()
    except (IOError, KeyError):
        print('CANNOT WRITE THE FILE')
        print(format_exc())
    #pprint(qmout)
    if stat != 0:
        current_path = os.path.abspath(os.getcwd())
        tmp = time.strftime('%m-%d-%H-%M-%S', time.localtime())
        print('backup QM.log to QM.log.' + tmp)
        exe_str = 'cd %s && cp QM.log QM.log.%s && cd %s'%(
            qm_config['QM']['env']['directory'],
            tmp,
            current_path
        )
        os.system(exe_str)

    return qmout

def main():
    return

if __name__ == '__main__':
    args = parser.parse_args()
    # pprint(args)
    if args.task == 'nad':
        qm_data_in = QMutils.parseQMinput(args.input)
        qm_job(qm_data_in, args)
    elif args.task == 'hess_ds':
        hess_log = args.input
        lines = open(hess_log, 'r').readlines()
        find, hess = _get_hessian_from_lines(lines)

        fo = open(args.output, 'w')
        if 'vpes' in hess:
            fo.write('model.vpes\n');
            fo.write(f'kids_real 1\n');
            fo.write('{: 12.8e} \n\n'.format(hess['vpes']))
        if 'hess' in hess:
            fo.write('model.hess\n');
            val = hess['hess'].flatten()
            fo.write(f'kids_real {len(val)}\n');
            for i in range(len(val)):
                fo.write('{: 12.8e} '.format(val[i]))
            fo.write('\n\n')
        if 'Tmod' in hess:
            fo.write('model.Tmod\n');
            val = hess['Tmod'].flatten()
            fo.write(f'kids_real {len(val)}\n');
            for i in range(len(val)):
                fo.write('{: 12.8e} '.format(val[i]))
            fo.write('\n\n')
        if 'w' in hess:
            fo.write('model.w\n');
            val = hess['w'].flatten() #?
            fo.write(f'kids_real {len(val)}\n');
            for i in range(len(val)):
                fo.write('{: 12.8e} '.format(val[i]))
            fo.write('\n\n')
        if 'x0' in hess:
            fo.write('model.x0\n');
            val = hess['x0'].flatten() #?
            fo.write(f'kids_real {len(val)}\n');
            for i in range(len(val)):
                fo.write('{: 12.8e} '.format(val[i]))
            fo.write('\n\n')

            fo.write('model.p0\n');
            fo.write(f'kids_real {len(val)}\n');
            for i in range(len(val)):
                fo.write('{: 12.8e} '.format(0))
            fo.write('\n\n')
        fo.close()
        for i in range(len(hess['w'])):
            print(hess['x0'])

            w = np.array(hess['w']).flatten()
            x0 = np.array(hess['x0']).flatten() * QMutils.au_2_ang
            Tmod = np.array(hess['Tmod'])
            mass = np.array(hess['mass']).flatten()
            sym = hess['sym']
            print(sym, mass)
            Natom = len(w)//3
            ft = open('nma-%d.xyz'%i, 'w')
            for t in np.linspace(0, np.pi*10/w[i], 100):
                ft.write('%d\n\n'%Natom)
                dx = x0*0 
                dx[i] = 1*np.sin(w[i]*t)
                xt = x0 + 10 / np.sqrt(mass) * np.dot(Tmod, dx)

                for ia in range(Natom):
                    ft.write('%s %12.8f %12.8f %12.8f\n'%(sym[ia//3], xt[3*ia], xt[3*ia+1], xt[3*ia+2]))
            ft.close()

        # pprint(hess)
