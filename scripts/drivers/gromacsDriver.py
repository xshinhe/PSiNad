#!/usr/bin/env python3
#   Coding=utf-8

#   KIDS SCRIPTS
#   Author: xshinhe
#   
#   Copyright (c) 2024 PeKing Univ. - GNUv3 License

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


## importing external functions and modules
import numpy as np
import os
import shutil
import shlex
import subprocess
import math
import time

import kids_env
from kids_log import Log

def prepare(geometry, ks_config):
    """ prepare the input files for MM gromacs calculations"""
    gro_real, gro_modelH = ks_config.args.grofile.split(',')
    top_real, top_modelH = ks_config.args.topology.split(',')
    
    dict1 = {}
    dict2 = {}

    f0new = top_modelH.replace('.top', '-noc.top')
    lines = open(top_modelH, 'r', encoding='utf8').readlines()
    fout = open(f0new, 'w')
    for line in lines:
        if 'qtot' in line:
            terms = line.split()

            key = terms[3]+'-'+terms[4]
            value = terms[6]
            dict1[key] = value

            terms[6] = "0.0000000"
            terms[10] = "0.0000000"
            nl = '  '.join(terms)
            fout.write(nl+'\n')
        else:
            fout.write(line)
    fout.close()

    f0new = top_real.replace('.top', '-modelnoc.top')
    lines = open(top_real, 'r', encoding='utf8').readlines()
    fout = open(f0new, 'w')
    for line in lines:
        if 'qtot' in line:
            terms = line.split()
            key = terms[3]+'-'+terms[4]
            value = terms[6]
            dict2[key] = value

            if key in dict1:
                terms[6] = "0.0000000"
                terms[10] = "0.0000000"
                nl = '  '.join(terms)
                fout.write(nl+'\n')
            else:
                fout.write(line)                
        else:
            fout.write(line)
    fout.close()

    CRG_real = []
    CRG_model_H = []

    lines = open(gro_real, 'r', encoding='utf8').readlines()
    natom = int(lines[1])
    for i in range(2, natom+2):
        terms = lines[i].split()
        key = ''.join([i for i in terms[0] if not i.isdigit()])
        key += '-' + terms[1]
        if key in dict1:
            CRG_model_H += [float(dict1[key])]

        if key in dict2:
            CRG_real += [float(dict2[key])]
        else:
            print(key)
            print('error')
            exit(0)
    return CRG_real, CRG_model_H


def MMcalculations(geometry, ks_config):

    gro_real, gro_modelH = ks_config.args.grofile.split(',')
    top_real, top_modelH = ks_config.args.topology.split(',')

    dict_xyz = {gro_real: geometry.cartesian, gro_modelH: geometry.modelH}

    for grofile in [gro_modelH, gro_real]:
        f1 = open(grofile, 'r')
        lines1 = f1.readlines()
        natom = int(lines1[1])
        fo = open(grofile.replace('.gro', '-gen.gro'), 'w')
        for i in range(2):
            fo.write(lines1[i])
        for i in range(2, 2+natom):
            idx = i-2
            terms = lines1[i].split()
            terms[3] = (dict_xyz[grofile][0][idx] / 10)
            terms[4] = (dict_xyz[grofile][1][idx] / 10)
            terms[5] = (dict_xyz[grofile][2][idx] / 10)
            # print(terms)
            nl = "%5d%-5s%5s%5d%8.4f%8.4f%8.4f" % (
                idx+1,
                terms[0][-3:],
                terms[1],
                int(terms[2]),
                float(terms[3]),
                float(terms[4]),
                float(terms[5]),
            )
            fo.write(nl+'\n')
        for i in range(2+natom, len(lines1)):
            fo.write(lines1[i])
        fo.close()

    mdp = ks_config.get_nested('MM.GROMACS.level0', '')
    if mdp != '':
        f = open('gromacs.mdp', 'w')
        f.write(mdp)
        f.close()
    
    #####################
    # gmx full MM single point on real
    E_real, Fxyz_real = gmx_real(geometry, ks_config)
    #####################
    # gmx modelnoc MM single point on real w ZERO charges on model-H
    E_real_modelnoc, Fxyz_real_modelnoc = gmx_real_modelnoc(ks_config)

    #####################
    # gmx modelnoc MM single point on totally zero charges 
    #E_totalnoc, Fxyz_totalnoc = sandertotalnoc(ks_config)
    #####################
    # run model-H MM single point with ZERO charges
    E_modelH, Fxyz_modelH = gmx_modelH(geometry,ks_config)

    return E_real, Fxyz_real, E_real_modelnoc, Fxyz_real_modelnoc, E_modelH, Fxyz_modelH


def gmx_real(geometry, ks_config):
    if geometry.calculationType == 'M' or geometry.calculationType == 'ML':
        Log.writeLog("Perform a single point calculation and collect energy and gradient ...")
    
    comm1 = shlex.split("gmx grompp -f gromacs.mdp -c real-gen.gro -p real.top -o real.tpr -backup no")
    subprocess.call(comm1)

    comm2 = shlex.split("gmx mdrun -s real.tpr -o real.trr -e real.edr -g real.log -backup no")
    subprocess.call(comm2)

    Fxyz_real = read_forcesA('real.trr')
    E_real = read_energyMM('real.log')

    if geometry.calculationType == 'M' or geometry.calculationType == 'ML':
        Log.writeLog(' done!\n')
        Log.writeLog('MM energy is {0} Hartree\n'.format(E_sander_real))

    return [E_real, Fxyz_real]

def gmx_real_modelnoc(ks_config):
    Log.writeLog("\nPerform a single point calc. of the entire system with ZERO charge on the \nHigh layer, and collect energy and gradient ...")

    comm1 = shlex.split("gmx grompp -f gromacs.mdp -c real-gen.gro -p real-modelnoc.top -o real-modelnoc.tpr -backup no")
    subprocess.call(comm1)

    comm2 = shlex.split("gmx mdrun -s real-modelnoc.tpr -o real-modelnoc.trr -e real-modelnoc.edr -g real-modelnoc.log -backup no")
    subprocess.call(comm2)
    
    Fxyz_modelnoc = read_forcesA('real-modelnoc.trr')
    E_modelnoc = read_energyMM('real-modelnoc.log')

    Log.writeLog(' done!\n')
    Log.writeLog('MM energy is {0} Hartree\n'.format(E_modelnoc))
#    Log.writeLog(kids_log.matrix_prettystring(np.array(Fxyz_modelnoc), ".6f"), 2)

    return [E_modelnoc,Fxyz_modelnoc]

def gmx_modelH(geometry,ks_config):
    Log.writeLog("\nPerform a single point calc. of the High layer, with H-saturated bonds and \nZERO charge, and collect energy and gradient ...")

    comm1 = shlex.split("gmx grompp -f gromacs.mdp -c model-H-gen.gro -p model-H-noc.top -o model-H-noc.tpr -maxwarn 1 -backup no")
    subprocess.call(comm1)

    comm2 = shlex.split("gmx mdrun -s model-H-noc.tpr -o model-H-noc.trr -e model-H-noc.edr -g model-H-noc.log -backup no")
    subprocess.call(comm2)

    Fxyz_modelH = read_forcesA('model-H-noc.trr')
    E_modelH = read_energyMM('model-H-noc.log')

    Log.writeLog(' done!\n')
    Log.writeLog('MM energy is {0} Hartree\n'.format(E_modelH))
#    Log.writeLog(kids_log.matrix_prettystring(np.array(Fxyz_modelH), ".6f"), 2)

    return [E_modelH,Fxyz_modelH]


def read_forcesA(filename):  # BUG TODO
    c1=2625.49963 # 1 au = 2625.49963 kJ/mol.
    c2=0.05291772083 # 1 au = 0.0529 nm
    c3=c2/c1
    c4=c1/c2
    c5=1/(c1*c2)

    o1 = subprocess.check_output('gmx dump -f %s | head -n2'%filename, shell=True)
    length =  int(o1.split(b'\n')[1].strip().split()[1])
    # print(length)

#    print length
    filein = subprocess.check_output('gmx dump -f %s | grep -A%s "f ("'%(filename,str(length)), shell=True)
    filein=filein.split(b'\n')
    # print(filein)
    fx=[]
    fy=[]
    fz=[]
    for i in range(1,length+1):
        row=filein[i].split(b'{')[1]
        # print(row)
        row = row.split(b'}')[0]
        row=row.split(b',')
        fx.append(float(row[0])*c3) # convert to au
        fy.append(float(row[1])*c3)
        fz.append(float(row[2])*c3)
    Fx=np.array(fx)
    Fy=np.array(fy)
    Fz=np.array(fz)
    Fxyz=[Fx,Fy,Fz]
    return Fxyz

def read_energyMM(filename):
    c1=2625.49963 
    filein=subprocess.check_output('grep "Potential Energy" ' + filename, shell=True)
    filein=filein.split()
    ene1=float(filein[3])/c1
    return ene1

def clean():
    """ clean up the run directory from all the files that have been used to run Amber and that
        are no longer needed at the end of the calculation """

    # list of files/directory to remove
    toRemove = ['real-noc.xml', "model-H-noc.xml", "model-H.rst"]
    if False:
        for f in toRemove:
            # if f is an existing file, remove it
            if os.path.isfile(f):
                os.remove(f)
            # if f is an existing directory, clean tree
            elif os.path.isdir(f):
                shutil.rmtree(f)

    # done
    return
