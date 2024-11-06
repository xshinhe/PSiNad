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
import copy 


parser = argparse.ArgumentParser(description='Execute BDF Calculation')
parser.add_argument('-d', '--directory', dest='directory', nargs='?', default='.', type=str,
    help='work directory')
parser.add_argument('-i', '--input', dest='input', nargs='?', default='QM.in.BDF', type=str,
    help='input file')
parser.add_argument('-t', '--task', dest='task', nargs='?', default='nad', type=str,
    help='task type')
parser.add_argument('-o', '--output', dest='output', nargs='?', default='QM.out.BDF', type=str,
    help='output file')
# args = parser.parse_args()

def _get_stat(lines: list[str]):
    finde = True
    for i in range(len(lines)):
        if 'Total wall time' in lines[i]:
            finde = False 
            break 
    return finde

def _get_oscstrength_from_lines(lines: list[str]):
    find = False
    f_strength = [0,]
    for i in range(len(lines)):
        if "State          X           Y           Z          Osc." in lines[i]:
            k = i + 1
            while lines[k].strip() != '':
                terms = lines[k].strip().split()
                f_strength += [float(terms[4])]
                k += 1
            find = True
            break
    return find, np.array(f_strength)

def _get_energy_grad_nac_from_lines(lines: list[str]):
    find = False
    has_GS = False
    energies = []
    gradient = {}
    nac = {}
    istate = 0
    nac_istate = 0
    nac_jstate = 0
    for i in range(len(lines)):
        if "GSGRAD_estate=" in lines[i] or "EXGRAD_estate=" in lines[i]:
            find = True
            if 'GSGRAD_estate' in lines[i]: has_GS = True

            E = float(lines[i].split('=')[1])
            energies += [E]

            grad_local = []
            j = i + 2
            terms = lines[j].strip().split()
            while terms[0] != 'Sum':
                grad_local += [ np.array(terms[1:4]).astype(np.float64) ]
                j += 1
                terms = lines[j].strip().split()            
            gradient['%d'%istate] = deepcopy(np.array(grad_local))

            istate += 1

        if 'FNAC_cpair=' in lines[i]:
            terms = lines[i].strip().split()
            nac_istate, nac_jstate = int(terms[3]), int(terms[6])
            if not has_GS:
                nac_istate -= 1
                nac_jstate -= 1

        if 'Gradient contribution from Final-NAC(S)-Escaled' in lines[i]:
            grad_local = []
            j = i + 1
            terms = lines[j].strip().split()
            while terms[0] != 'Sum':
                grad_local += [ np.array(terms[1:4]).astype(np.float64) ]
                j += 1
                terms = lines[j].strip().split() 
            nac['%d-%d'%(nac_istate, nac_jstate)] = np.array(grad_local)
            nac['%d-%d'%(nac_jstate, nac_istate)] =-np.array(grad_local)  

    return find, np.array(energies), gradient, nac

def _get_hessian_from_lines(lines: list[str]):
    find = {}
    hess_dict = {}
    geom_begin_idx = -1
    resforce_begin_idx = -1
    hess_begin_idx = -1
    for i in range(len(lines)):
        if "Cartesian coordinates (Angstrom)" in lines[i] and geom_begin_idx == -1:
            geom_begin_idx = i
            find['x0'] = True
        if ' Results of vibrations:' in lines[i] and hess_begin_idx == -1:
            hess_begin_idx = i
            find['w'] = True
            find['Tmod'] = True

    if 'x0' in find:
        sym = []
        xyz = []
        for i in range(geom_begin_idx+4, len(lines)):
            ll = lines[i].strip()
            if ll == '': continue
            if ll[0].strip()[0] == '-':
                break
            print(ll)
            terms = ll.strip().split()
            sym += [terms[1]] 
            xyz += [ np.array(terms[3:6]).astype(np.float64) ]
        xyz = np.array(xyz)
        hess_dict['sym'] = sym
        hess_dict['x0'] = xyz / QMutils.au_2_ang
        print(sym)

        mass = []
        mass = np.zeros(3*len(sym))
        mdict = {'C':12.01, 'H':1.008, 'O': 16.00}
        for k in range(len(mass)):
            mass[k] = mdict[sym[k//3]] * 1850
        hess_dict['mass'] = mass
        
    if 'Tmod' in find:
        natom = len(xyz)
        N = 3*natom
        hess = np.zeros((N,N))
        Tmod = np.zeros((N,N))
        w = np.zeros((N))
        i = hess_begin_idx+1
        read_col_num = 0
        while i < len(lines) and read_col_num < N:
            if 'Irreps' in lines[i]:
                # read freq
                w[read_col_num:read_col_num+3] = np.array(lines[i+1].strip().split()[1:4]).astype(np.float64)
                # read Tmod 
                for j in range(natom):
                    tget = np.array(lines[i+5+j].strip().split()[2:11]).astype(np.float64)
                    Tmod[3*j:3*j+3, read_col_num] = tget[0:3]
                    Tmod[3*j:3*j+3, read_col_num+1] = tget[3:6]
                    Tmod[3*j:3*j+3, read_col_num+2] = tget[6:9]
                read_col_num += 3
                i += 5 + natom
            else:
                i += 1
        wnew = copy.deepcopy(w)
        Tnew = copy.deepcopy(Tmod)
        wnew[0:6] = w[-6:]
        wnew[6:] = w[:-6]
        Tnew[:,0:6] = Tmod[:,-6:]
        Tnew[:,6:] = Tmod[:,:-6]
        hess_dict['Tmod'] = Tnew
        hess_dict['w'] = wnew / QMutils.au_2_wn

    return find != {}, hess_dict

def qm_job(qm_data, args):
    qm_config = qm_data["qm_config"]
    bdf_config = qm_config["QM"]["BDF"]
    natom = qm_data["natom"]
    znumber = qm_data["znumber"]
    xyz = qm_data["geom_xyz"]
    
    try:
        F = int(qm_config['QM']['BDF']['F'])
        N = int(qm_config['QM']['BDF']['N'])
    except KeyError:
        print(format_exc())
    nthread = 1
    if 'nthread' in qm_config['QM']:
        nthread = qm_config['QM']['nthread']

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

    job_tpl = qm_config['QM']['BDF']['level%d'%level]
    lines = job_tpl.split('\n')
    job_str = ''
    for l in lines:
        if l == '$COORD_TPL':
            for i in range(natom):
                job_str += '%s   %12.8e  %12.8e  %12.8e\n'%(
                QMutils.element_list[znumber[i]][0], xyz[i][1], xyz[i][2], xyz[i][3])
        else:
            job_str += l + '\n'
    # print(job_str)
    qm_config["QM"]["env"] = {
        "input_is_ready":True,
        "generated": "QM.run.BDF.inp",
        "directory": args.directory,
        "output": args.output,
    }
    #print('config managed by kidsqmmm module:')
    #pprint(qm_config)

    directory = qm_config['QM']['env']['directory']
    if not os.path.exists(directory): os.makedirs(directory)

    print(qm_config['QM']['env']['generated'])
    f = open(directory + '/' + qm_config['QM']['env']['generated'], 'w')
    f.write(job_str)
    f.flush()
    f.close()

    #print("try to do BDF job with parameter as:")
    #print(job_str)

    current_path = os.path.abspath(os.getcwd())
    exe_str = 'cd %s && %s  %s > %s && cd %s'%(
        qm_config['QM']['env']['directory'],
        qm_config['QM']['BDF']['path'],
        qm_config['QM']['env']['generated'],
        qm_config['QM']['env']['output'],
        current_path
    )
    print(exe_str)
    stat = os.system(exe_str)
    parse_result(
        qm_data,
        directory
    )

def parse_result(qm_data, log_path):
    try:
        natom = qm_data['natom']
        qm_config = qm_data['qm_config']
        F = int(qm_config['QM']['BDF']['F'])
        N = int(qm_config['QM']['BDF']['N'])
        output = qm_config['QM']['env']['output']
    except KeyError:
        print(format_exc())

    eig = np.zeros((F))
    dE  = np.zeros((N, F))
    nac = np.zeros((N, F, F))

    fstrength = np.zeros((F))
    
    lines = open(log_path +'/'+ output, 'r').readlines()
    finde = _get_stat(lines)
    findf, fstrength = _get_oscstrength_from_lines(lines)

    stat = 0
    if finde: stat = 1
    try:
        find0, energy, gradient, nacvector = _get_energy_grad_nac_from_lines(lines)
        eig = energy # already in au
        for i in range(F):
            dE[:,i] = gradient['%d'%i].flatten() # already in au
        for i in range(F):
            for k in range(i+1,F):
                try:
                    nac[:,i,k] = nacvector['%d-%d'%(i,k)].flatten() # already in au
                    nac[:,k,i] = nacvector['%d-%d'%(k,i)].flatten() # already in au
                except:
                    pass

    except (IOError, KeyError):
        print(format_exc())
        print('CANNOT PARSE FILE. try to print ERROR MSG')
        # print(ERROR_MSG)
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

            f.write('interface.dE\n')
            f.write(f'kids_real {N*F}\n')
            for j in range(N):
                for i in range(F):
                    f.write('{: 12.8e} '.format(dE[j,i]))
                f.write('\n')
            f.write('\n')

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
    elif args.task == 'test':
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
            print(sym, mass, 'aaaa')
            Natom = len(w)//3
            ft = open('nma-%d.xyz'%i, 'w')
            for t in np.linspace(0, np.pi*10/w[i], 100):
                ft.write('%d\n\n'%Natom)
                dx = x0*0 
                dx[i] = 1*np.sin(w[i]*t)
                xt = x0 + 10 / np.sqrt(mass) * np.dot(Tmod, dx)

                for ia in range(Natom):
                    ft.write('%s %12.8f %12.8f %12.8f\n'%(sym[ia], xt[3*ia], xt[3*ia+1], xt[3*ia+2]))
            ft.close()

        # pprint(hess)
