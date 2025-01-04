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

parser = argparse.ArgumentParser(description='Execute MNDO Calculation')
parser.add_argument('-d', '--directory', dest='directory', nargs='?', default='.', type=str,
    help='work directory')
parser.add_argument('-i', '--input', dest='input', nargs='?', default='QM.in.MNDO', type=str,
    help='input file')
parser.add_argument('-t', '--task', dest='task', nargs='?', default='0', type=str,
    help='task type')
parser.add_argument('-o', '--output', dest='output', nargs='?', default='QM.log', type=str,
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
        if "Properties of transitions   1 -> #" in lines[i]:
            # i+1 is blankline
            # i+2 is headline
            k = 3
            while lines[i+k].strip() != '':
                read_pos = -2 # read f_rp only
                f_strength +=  [ float(lines[i+k].split()[read_pos]) ]
                k += 1
            find = True
            break
    return find, np.array(f_strength)

def _get_energy_from_lines(lines: list[str]):
    find = False
    Emin = 0
    energies = []
    for i in range(len(lines)):
        if "State  1,  Mult. 1,  E-E(1)=  0.000000" in lines[i]:
            Emin = float(lines[i].split()[-2])
        elif "SUMMARY OF MULTIPLE CI CALCULATIONS" in lines[i]:
            k = 5
            while lines[i+k].strip()[0:2] != '--':
                energies += [float(lines[i+k].split()[2])]
                k += 1 
            find = True
    return find, np.array(energies)

def _get_gradient_from_lines(lines: list[str]):
    find = False
    gradient = {}
    istate_force = 0
    i = 0;
    while i < len(lines):
        if "CI CALCULATION FOR STATE:" in lines[i]:
            istate_force = int(lines[i].split()[-1]) - 1
            k = i + 1
            while "GRADIENTS (KCAL/(MOL*ANGSTROM))" not in lines[k]: k += 1
            k += 4
            grad_local = []
            while lines[k].strip() != '':
                grad_local += [ np.array(lines[k].split()[5:]).astype(np.float64) ]
                k += 1
            for _ in range(10):
                k += 1
                if 'EXTERNAL POINT CHARGES' in lines[k]:
                    while "GRADIENTS (KCAL/(MOL*ANGSTROM))" not in lines[k]: k += 1
                    k += 4
                    while lines[k].strip() != '':
                        grad_local += [ np.array(lines[k].split()[5:]).astype(np.float64) ]
                        k += 1
                    break
            gradient['%d'%istate_force] = deepcopy(np.array(grad_local))
            find = True 
            i = k + 1
        else:
            i += 1
    return find, gradient

def _get_nac_from_lines(lines: list[str]):
    find = False
    nac = {}
    i = 0
    while i < len(lines):
        if "CI CALCULATION FOR INTERSTATE COUPLING OF STATES:" in lines[i]:
            istate, jstate = map(int, lines[i].split()[-2:])
            istate -= 1; jstate -= 1
            k = i + 1
            while "GRADIENTS (KCAL/(MOL*ANGSTROM))" not in lines[k]: k += 1
            k += 4
            grad_local = []
            while lines[k].strip() != '':
                grad_local += [ np.array(lines[k].split()[5:]).astype(np.float64) ]
                k += 1
            for _ in range(10):
                k += 1
                if 'EXTERNAL POINT CHARGES' in lines[k]:
                    while "GRADIENTS (KCAL/(MOL*ANGSTROM))" not in lines[k]: k += 1
                    k += 4
                    while lines[k].strip() != '':
                        grad_local += [ np.array(lines[k].split()[5:]).astype(np.float64) ]
                        k += 1
                    break
            nac['%d-%d'%(istate, jstate)] = np.array(grad_local)
            nac['%d-%d'%(jstate, istate)] =-np.array(grad_local)
            find = True
            i = k + 1
        else:
            i += 1
    return find, nac

def _get_hessian_from_lines(lines: list[str]):
    find = False
    hess = {}
    
    ncalc = 0
    for li in lines:
        if 'CURRENT CARTESIAN GRADIENT' in li:
            ncalc += 1
    N = (ncalc - 1)//2

    Emin = 0
    freq = np.zeros(N)
    Tmod = np.zeros((N,N))
    xyz = np.zeros((N))
    i = 0 
    while i < len(lines):
        if "     INPUT GEOMETRY" in lines[i]:
            k = i + 6
            xyz_local = []
            sym = []
            sym_dict = {'1':'H', '6':'C', '7':'N', '8':'O'}
            while lines[k].strip() != '':
                terms = lines[k].split()
                if terms[3] == '*':
                    xyz_local += [ np.array([terms[2], terms[4], terms[6]]).astype(np.float64) ]
                else:
                    xyz_local += [ np.array(terms[2], terms[3], terms[5]).astype(np.float64) ]
                sym += [ sym_dict[terms[1]] ]
                k += 1
            xyz[:] = np.array(xyz_local).flatten()
            mass = []
            mass = np.zeros(3*len(sym))
            mdict = {'C':12.01, 'H':1.008, 'N':14.00, 'O': 16.00}
            for x in range(len(mass)):
                mass[x] = mdict[sym[x//3]] * 1850
            hess['sym'] = sym
            hess['mass'] = mass
            i = k + 1
        elif "State  1,  Mult. 1,  E-E(1)=  0.00000" in lines[i] and Emin == 0:
            Emin = float(lines[i].split()[-2])
            i += 1
        elif "GRADIENT NORM =" in lines[i]:
            norm = float(lines[i].split()[3])
            i += 1
            if norm > 10.0:
                # print('Hessian is not used under equilibrium!')
                pass
        elif 'EIGENVECTORS OF THE MASS-WEIGHTED' in lines[i]:
            k = i+3
            colidx = np.array(lines[k].split()).astype(np.int32) - 1
            k += 2
            freq[colidx] = np.array(lines[k].split()).astype(np.float64)
            k += 2
            while 'CARTESIAN DISPLACEMENT' not in lines[k]: 
                if lines[k].strip() == '': 
                    k += 1
                    continue
                terms = lines[k].split()
                if len(terms) == len(colidx) + 1:
                    rowidx = int(terms[0]) - 1
                    Tmod[rowidx, colidx] = np.array(terms[1:]).astype(np.float64)
                    k += 1
                if len(terms) <= len(colidx):
                    colidx = np.array(terms).astype(np.int32) - 1
                    k += 2
                    terms = lines[k].split()
                    freq[colidx] = np.array(terms).astype(np.float64)
                    k += 1
            find = True
            hess['vpes'] = Emin
            hess['x0'] = xyz
            hess['w'] = freq
            hess['Tmod'] = Tmod
            hess['hess'] = np.einsum('ik,k,kj->ij', Tmod, freq, Tmod.T)
            i = k + 1
            break 
        else:
            i += 1
    return find, hess

def qm_job(qm_data, args):
    qm_config = qm_data["qm_config"]
    mndo_config = qm_config["QM"]["MNDO"]
    natom = qm_data["natom"]
    znumber = qm_data["znumber"]
    xyz = qm_data["geom_xyz"]
    
    job_str = ""
    if args.task == 0 and "toplines" in mndo_config:
        job_str += mndo_config["toplines"]  + '\n'
    elif args.task ==1 and "toplines_1" in mndo_config:
        job_str += mndo_config["toplines_1"]  + '\n'
    elif args.task ==2 and "toplines_2" in mndo_config:
        job_str += mndo_config["toplines_2"]  + '\n'
    elif args.task ==3:
        f = open(directory + '/stat.dat', 'w')
        f.write('1\n')
        f.write('FINALLY FAILED in LEVEL 3\n')
        f.flush()
        f.close()
        return
    elif "first_keywords" in mndo_config:
        for k in mndo_config["first_keywords"]:
            if mndo_config["first_keywords"][k] == True:
                job_str += mndo_config["first_keywords"][k] + ' '
            else:
                job_str += k + '=' + mndo_config["first_keywords"][k] + ' '
        job_str += '+\n'

        if "second_keywords" not in mndo_config:
            print('error')
            exit(0)

        for k in mndo_config["second_keywords"]:
            if mndo_config["second_keywords"][k] == True:
                job_str += mndo_config["second_keywords"][k] + ' '
            else:
                job_str += k + '=' + mndo_config["second_keywords"][k] + ' '
        job_str += '\n'
    elif "keywords" in mndo_config:
        for k in mndo_config["keywords"]:
            if mndo_config["keywords"][k] == True:
                job_str += mndo_config["keywords"][k] + ' '
            else:
                job_str += k + '=' + mndo_config["keywords"][k] + ' '
        job_str += '\n'

    job_str += "Automatically generated by MNDO.py \n[%s]\n" % datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    for i in range(natom):
        job_str += "{:<3d} {: 12.8e} {: 3d} {: 12.8e} {: 3d} {: 12.8e} {: 3d}\n".format(
                znumber[i], xyz[i][1], 0, xyz[i][2], 0, xyz[i][3], 0)
    job_str += "{:<3d} {: 12.8e} {: 3d} {: 12.8e} {: 3d} {: 12.8e} {: 3d}\n".format(
        0, 0, 0, 0, 0, 0, 0)

    if "bottomlines" in mndo_config:
        job_str += mndo_config["bottomlines"]  + '\n'

    qm_config["QM"]["env"] = {
        "input_is_ready":True,
        "generated": "QM.run.MNDO",
        "directory": args.directory,
        "output": args.output,
    }

    directory = qm_config['QM']['env']['directory']
    if not os.path.exists(directory): os.makedirs(directory)

    f = open(directory + '/' + qm_config['QM']['env']['generated'], 'w')
    f.write(job_str)
    f.flush()
    f.close()

    exe_str = '%s < %s > %s'%(
        qm_config['QM']['MNDO']['path'],
        directory + '/' + qm_config['QM']['env']['generated'],
        directory + '/' + qm_config['QM']['env']['output']
    )
    os.system(exe_str)

    parse_result(
        qm_data,
        directory + '/' + qm_config['QM']['env']['output']
    )

def parse_result(qm_data, log_file):
    try:
        natom = qm_data['natom']
        qm_config = qm_data['qm_config']
        F = int(qm_config['QM']['MNDO']['F'])
        N = int(qm_config['QM']['MNDO']['N'])
        iroot = F # int(qm_config['QM']['MNDO']['keywords']['iroot'])

        f = open(log_file, 'r')
        lines = f.readlines()

        finde, ERROR_MSG = _get_errormsg_from_lines(lines)
        findu, UNABLE_MSG = _get_unable_from_lines(lines)
        find0, energy = _get_energy_from_lines(lines)
        find1, gradient = _get_gradient_from_lines(lines)
        findc, nacvector = _get_nac_from_lines(lines)
        findf, fstrength = _get_oscstrength_from_lines(lines)
    except (KeyError, IOError):
        print(format_exc())
        
    stat = 0
    if finde or findu:
        stat = 1

    eig = energy / QMutils.au_2_kcal_1mea
    dE  = np.zeros((N,F))
    nac = np.zeros((N,F,F))
    for i in range(F):
        dE[:,i] = gradient['%d'%i].flatten() / QMutils.au_2_kcal_1mea_per_ang
    for i in range(F):
        for k in range(i+1,F):
            nac[:,i,k] = nacvector['%d-%d'%(i,k)].flatten() / (1.0e0 / QMutils.au_2_ang)
            nac[:,k,i] = nacvector['%d-%d'%(k,i)].flatten() / (1.0e0 / QMutils.au_2_ang)

    qmout = QMutils.QMout(natom=natom,
        energy=eig,
        gradient=dE,
        nac=nac
    )

    f = open(qm_config['QM']['env']['directory'] + '/interface.ds', 'w')
    f.write('interface.stat\n')
    f.write(f'kids_int {1}\n')
    f.write(f'{stat}\n\n')

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

    f.write('interface.strength\n')
    f.write(f'kids_real {F}\n')
    for i in range(F):
        f.write('{: 12.8e} '.format(fstrength[i]))
    f.write('\n')
    f.close()

    #pprint(qmout)
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
            fo.write('{: 12.8e} \n\n'.format(hess['vpes'] / QMutils.au_2_ev))
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
            val = hess['w'].flatten() / QMutils.au_2_wn
            fo.write(f'kids_real {len(val)}\n');
            for i in range(len(val)):
                fo.write('{: 12.8e} '.format(val[i]))
            fo.write('\n\n')
        if 'x0' in hess:
            fo.write('model.x0\n');
            val = hess['x0'].flatten() / QMutils.au_2_ang
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
            #print(hess['x0'])

            w = np.array(hess['w']).flatten()
            x0 = np.array(hess['x0']).flatten() # * QMutils.au_2_ang
            Tmod = np.array(hess['Tmod'])
            mass = np.array(hess['mass']).flatten()
            sym = hess['sym']
            #print(sym, mass)
            Natom = len(w)//3
            ft = open('nma-%d.xyz'%i, 'w')
            beta = 3.15e5 / 300 
            sigma = 1 /(np.sqrt(beta) * w[i] / QMutils.au_2_wn)
            print(sigma)
            for t in np.linspace(0, np.pi*10/w[i], 100):
                ft.write('%d\n\n'%Natom)
                dx = x0*0 
                dx[i] = 1*np.sin(w[i]*t)
                xt = x0 +  sigma / np.sqrt(mass) * np.dot(Tmod, dx)

                for ia in range(Natom):
                    #print(ia, ia//3,  sym[ia//3])
                    ft.write('%s %12.8f %12.8f %12.8f\n'%(sym[ia], xt[3*ia], xt[3*ia+1], xt[3*ia+2]))
            ft.close()
        # pprint(hess)

    # debug
    # pprint(_get_unable_from_lines(lines))
    # pprint(_get_errormsg_from_lines(lines))
    # pprint(_get_oscstrength_from_lines(lines))
    # pprint(_get_energy_from_lines(lines))
    # pprint(_get_gradient_from_lines(lines))
    # pprint(_get_nac_from_lines(lines))
    # pprint(_get_hessian_from_lines(lines))

    main()
