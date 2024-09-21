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

## importing external functions and modules
import numpy as np
import os
import shutil
import shlex
import subprocess
import CBF
import logwrt
import qmmmenv
import math

def prepare(geometry, cobcom, command):
    """ prepare the input files for MM amber calculations"""

    # check that the environment is properly set
    envDefined, errorMsg = qmmmenv.checkAmberEnv()
    if not envDefined: logwrt.fatalerror(errorMsg)

    # set the input files
    fileinp_4="model-H.top"
    fileinp_5="real.top"

    # check their existence
    if not os.path.isfile(fileinp_5): logwrt.fatalerror( 'mandatory input file {0} not found!'.format(fileinp_5))

    ## get real sander input from cobram.command
    low_sander_input=CBF.ReadCobramCommand(cobcom,'sander','real-sander.inp')
    if len(low_sander_input) == 0:
        logwrt.fatalerror( "missing !sander ... ?sander section in the command file" )

    # check whether the option "ntxo = 1" is present, and in case add it to the amber input options
    if not "ntxo" in [ line.split("=")[0].strip() for line in low_sander_input ]:
        low_sander_input.insert(2, "ntxo = 1,") # ntxo=1, ASCII format other than NetCDF

    if geometry.calculationType == 'M' or geometry.calculationType == 'ML':
        #logwrt.writelog('No topology check needed: model-H.top is not needed')
        pass
    else:
        ## check the consistency of the atom type in real.top e model-H.top
        if not os.path.isfile(fileinp_4): logwrt.fatalerror( 'mandatory input file {0} not found!'.format(fileinp_4))
        checktoptop(geometry,fileinp_5,fileinp_4)

    ## Read the charge of real from real.top
    CRG_real=read_crg_from_top(fileinp_5)

    if geometry.calculationType == 'HML' or geometry.calculationType == 'HL' or geometry.calculationType == 'ML':
        ## get free mm atoms from cobram.command
        freeMM=CBF.ReadCobramCommand(cobcom,'free-mm','free-mm.dat')

        ## get bellymask from cobram.command
        bellymask=CBF.ReadCobramCommand(cobcom,'bellymask','bellymask.dat')

        ## complete the real-sander.inp and freeze HIGH and MEDIUM
        new_sander_input=free_MM(freeMM,geometry,low_sander_input,bellymask)

    if (geometry.calculationType in ['HML' ,'HL' ,'HM']):
        ## set to 0 the charges of the model part
        chargezeromodel("real.top","real-modelnoc.top",geometry)
        chargezero("model-H.top","model-H-noc.top")

    ## create the real second input
    cretesandersecond(low_sander_input,geometry,command)
    return CRG_real

def MMcalculations(step, geometry, command, cobcom):

    # amber optimization of low layer
    ifskip = skip_MMfirst(step, command, cobcom)
    if "L" in geometry.calculationType and ifskip == 0:  # low optimization is not skipped
        sanderfirst(command) # optimization of low
    else:  # low optimization is skipped
        shutil.copy('real.crd','real.rst')
    # update the real geometry
    geometry.updatereal('real.rst')

    #####################
    # amber full MM single point on real
    E_second, Fxyz_second = sandersecond(geometry, command)
    #####################
    # amber modelnoc MM single point on real w ZERO charges on model-H
    E_modelnoc, Fxyz_modelnoc = sandermodelnoc(command)
    #####################
    # run model-H MM single point with ZERO charges
    E_modelH, Fxyz_modelH = sandermodelH(geometry,command)

    return Fxyz_second, E_second, E_modelnoc, Fxyz_modelnoc, Fxyz_modelH, E_modelH


def greptop(starting,ending,filename):
    filein=open(filename)
    top=[]
    x=[]
    ## grep a part from a topology file
    for n, line in enumerate(filein, 1):
        top.append(line)
        try:
            point=line.split()
            if point[1].strip() == starting:
                start=n
            if point[1].strip() == ending:
                end=n
        except:
            pass
    filein.close()
    ## tranform in a list
    for i in range(start+1,end-1):
        x=x+(top[i].strip().split())
    return x

def checktoptop(geometry,realtop1,modelhtop1):
    ## check the consistency of the atom type in real.top e model-H.top
    realat=greptop('AMBER_ATOM_TYPE','TREE_CHAIN_CLASSIFICATION','real.top')
    modelhat=greptop('AMBER_ATOM_TYPE','TREE_CHAIN_CLASSIFICATION','model-H.top')
    ## check the number of QM atoms
    if (geometry.NatomQM+geometry.NsubH)!=(len(modelhat)):
        logwrt.fatalerror('the number of QM atoms in model-H.top ('+str(len(modelhat))+
                       ') and in real_layers.xyz ('+str(len(modelhat)+geometry.NsubH)+') are different\n')
    ## check the correspondance between atom types
    for i in range(geometry.NatomQM):
        if (modelhat[i]!=realat[geometry.list_HIGH[i]-1]):
            logwrt.fatalerror('atom types do not correspond in real.top \t'+realat[geometry.list_HIGH[i]-1]+'\n'+
                           '\tand model-H.top \t\t\t\t'+modelhat[i]+' atom '+str(i+1)+'\n'+
                           '\tcheck your input files\n')
    #logwrt.writelog('atom type check done!')
    return


def free_MM(frozenMM,geometry,low_sander_input,bellymask):
    ## create the list of frozen atoms in real sander and produce the
    ## real-sander.inp
    ## run the ambmask
    if (geometry.calculationType == 'HM'):
        return
    else:
        if (frozenMM!=[] and bellymask!=[]):
            logwrt.fatalerror('It is not possible to define freeMM and bellymask atoms in a single run\n')
        if (frozenMM==[]):
            frozen_MM=geometry.list_LOW
            startfrozen=[]
            endfrozen=[]
            frozen_MM2=[]
            ## condesate the atoms in the amber format
            for i in range(len(frozen_MM)-1):
                if (frozen_MM[i] != frozen_MM[i+1]-1):
                    startfrozen.append(frozen_MM[i+1])
                    endfrozen.append(frozen_MM[i])
            startfrozen.insert(0,frozen_MM[0])
            endfrozen.insert(len(endfrozen),frozen_MM[len(frozen_MM)-1])
            for i in range(len(startfrozen)):
                if (startfrozen[i] == endfrozen[i]):
                    frozen_MM2.append(str(startfrozen[i]))
                else:
                    frozen_MM2.append(str(startfrozen[i])+'-'+str(endfrozen[i]))
            ## play a little bit to write correctly the sander input
            for i in range(len(frozen_MM2)-1):
                frozen_MM2[i]=str(frozen_MM2[i])+',\n '
            low_sander_input.pop(len(low_sander_input)-1)
            low_sander_input[0]=low_sander_input[0].strip()+'\n'
            for i in range(1,len(low_sander_input)):
                low_sander_input[i]='   '+low_sander_input[i].strip()+'\n'
            low_sander_input.insert(len(low_sander_input)," bellymask='@ ")
            new_input=low_sander_input+list(frozen_MM2)
            try:
                for i in range(len(bellymask)):
                    if ( i == 0):
                        new_input.insert(len(new_input),'\n & ')
                    new_input.insert(len(new_input),str(bellymask[i].strip())+'\n')
            except:
                pass
            new_input.insert(len(new_input)," '\n")
            new_input.insert(len(new_input)," /")
            sander=open('real-sander-first.inp','w')
            ## write the final real-sander.inp
            for i in range(len(new_input)):
                sander.write(new_input[i])
            sander.write('\n')
            sander.close()
        else:
            atomnumber=[]
            frozenPOP=[]
            frozen_MM=[]
            subprocess.call( 'ambmask -p real.top -c real.crd -find '+"'"+frozenMM[0].strip()+"'"+
               ' -out pdb '+'> frozenmmatoms.pdb', shell=True)
            filein=open('frozenmmatoms.pdb')
            ## grep the atom numer from pdb
            for line in filein:
                row=line.split()
                if (row[0]=="ATOM"):
                    atomnumber.append(int(row[1]))
            total=list(atomnumber)+list(geometry.list_LOW)
            ## get the movable atoms
            for j in range(len(atomnumber)):
                n=total.count(total[j])
                if (str(n) == '2'):
                    frozenPOP.append(int(total[j]))
            frozen_MM=np.array(frozenPOP)
            startfrozen=[]
            endfrozen=[]
            frozen_MM2=[]
            ## condesate the atoms in the amber format
            for i in range(len(frozen_MM)-1):
                if (frozen_MM[i] != frozen_MM[i+1]-1):
                    startfrozen.append(frozen_MM[i+1])
                    endfrozen.append(frozen_MM[i])
            startfrozen.insert(0,frozen_MM[0])
            endfrozen.insert(len(endfrozen),frozen_MM[len(frozen_MM)-1])
            for i in range(len(startfrozen)):
                if (startfrozen[i] == endfrozen[i]):
                    frozen_MM2.append(str(startfrozen[i]))
                else:
                    frozen_MM2.append(str(startfrozen[i])+' '+str(endfrozen[i]))
            ## play a little bit to write correctly the sander input
            for i in range(len(frozen_MM2)):
                frozen_MM2[i]='ATOM '+str(frozen_MM2[i])+'\n'
            low_sander_input.pop(len(low_sander_input)-1)
            low_sander_input[0]=low_sander_input[0].strip()+'\n'
            for i in range(1,len(low_sander_input)):
                low_sander_input[i]=' '+low_sander_input[i].strip()+'\n'
            slash=["/\n"]
            title=["TITLE\n"]
            new_input=low_sander_input+slash+title+list(frozen_MM2)
            try:
                for i in range(len(bellymask)):
                    if ( i == 0):
                        new_input.insert(len(new_input),'\n & ')
                    new_input.insert(len(new_input),str(bellymask[i].strip())+'\n')
            except:
                pass
            new_input.insert(len(new_input),"END\nEND\n")
            sander=open('real-sander-first.inp','w')
            ## write the final real-sander.inp
            for i in range(len(new_input)):
                sander.write(new_input[i])
            sander.write('\n')
            sander.close()
        #logwrt.writelog('MM Frozen atoms check done!')
    return new_input[i]

def cretesandersecond(low_sander_input,geometry,command):
    ## search and replace in the real-sander.inp
    if len(low_sander_input)==0:
# (OW) Do the pmemd trick, compute one MD step in dielektrikum to obtain velocities for GPU-code
        if command[20]=='1':
            low_sander_input= ["minimizzazione","&cntrl","imin = 1,","ncyc = 1,","maxcyc = 0,","ntb = 0,","igb = 1,","extdiel = 1","ntr = 0,","DRMS = 0.0005,","ntxo=1,","cut = "+command[41],"/"]
        else:
            low_sander_input= ["minimizzazione","&cntrl","imin = 1,","ncyc = 1,","maxcyc = 0,","ntb = 0,","igb = 0,","ntr = 0,","DRMS = 0.0005,","ntxo=1,","cut = "+command[41],"/"]



    ## commands and write a correct second input
    p='belly'
    low_sander_input=list(filter(lambda low_sander_input:low_sander_input.find(p) == -1,low_sander_input))
#    p='iwforc'
#    low_sander_input=filter(lambda low_sander_input:low_sander_input.find(p) == -1,low_sander_input)
    p='maxcyc'
    low_sander_input=list(filter(lambda low_sander_input:low_sander_input.find(p) == -1,low_sander_input))
    if command[20]!='0':
        p='imin'
        low_sander_input=list(filter(lambda low_sander_input:low_sander_input.find(p) == -1,low_sander_input))
    p='/'
    low_sander_input=list(filter(lambda low_sander_input:low_sander_input.find(p) == -1,low_sander_input))
    p='cut'
    low_sander_input=list(filter(lambda low_sander_input:low_sander_input.find(p) == -1,low_sander_input))
    low_sander_input2=[]
    for i in range(len(low_sander_input)):
        low_sander_input2.append(low_sander_input[i].strip())
    a=low_sander_input2.index('&cntrl')
    low_sander_input.insert(a+1,' maxcyc = 0,')

## write the final real-sander-second.inp
    sander=open('real-sander-second.inp','w')
    if geometry.calculationType in ['HML','HL','HM','M','ML'] :
        for i in range(len(low_sander_input)):
            sander.write(' '+low_sander_input[i]+'\n')
# do the pmemd-trick to get forces from velocities (OW)
        sander.write('  cut = '+command[41]+',\n')
        if command[20]!='0':
           sander.write('  nstlim = 1,\n')
           sander.write('  dt     = 0.01,\n')
           sander.write('  ntwr   = 1,\n')
           sander.write('  ntt=0,\n')
           sander.write('  ntxo=1,\n')
# (OW)
        else:
           sander.write('  ntxo=1,\n')
           sander.write(' /\n')
           sander.write(' DUMP FORCES\n')
           sander.write(' &debugf\n')
           sander.write(' do_debugf=1,\n')
           sander.write(' dumpfrc=1,\n')
        sander.write(' /\n')
        sander.close()
    if geometry.calculationType in ['HML','HL','HM'] :
        sander=open('model-H-sander.inp','w')
        ## write the final model-H-sander.inp
        for i in range(len(low_sander_input)):
            sander.write(' '+low_sander_input[i]+'\n')
        sander.write('  cut = '+command[41]+',\n')

# do the pmemd-trick (OW)
        if command[20]!='0':
           sander.write('  nstlim = 1,\n')
           sander.write('  dt     = 0.01,\n')
           sander.write('  ntwr   = 1,\n')
           sander.write('  ntt=0,\n')
           sander.write('  ntxo=1,\n')
# (OW)
        else:
           sander.write('  ntxo=1,\n')
           sander.write(' /\n')
           sander.write(' DUMP FORCES\n')
           sander.write(' &debugf\n')
           sander.write(' do_debugf=1,\n')
           sander.write(' dumpfrc=1,\n')
        sander.write(' /\n')
        sander.close()
        return

def read_crg_from_top(filename):
    Camber=math.sqrt(332.0)
    crg_real_tmp=greptop('CHARGE','MASS',filename)
    crg_real=[]

    try:
       for i in range(len(crg_real_tmp)):
           c=float(crg_real_tmp[i])/Camber
           crg_real.append(c)
    except:
           crg_real=[]
           logwrt.writelog('Amber version > 11 detected!\n')
           crg_real_tmp=greptop('CHARGE','ATOMIC_NUMBER',filename)
           for i in range(len(crg_real_tmp)):
              c=float(crg_real_tmp[i])/Camber
              crg_real.append(c)

    CRG_real=np.array(crg_real)
    return CRG_real


def checksanderdone(filename, sandertype):
    """ Check if a sander calculation eneded correctly, if not exit from program """

    # initialize error flag
    SanderError = False
    # process the content of the output file, if the "Error" string is there set SanderError to true
    with open(filename) as filein:
        for line in filein:
            if str(line.strip()) == 'Error': SanderError = True
    # if SanderError is true, terminate COBRAMM with fatalerror
    if SanderError: logwrt.fatalerror('error termination in '+sandertype+' calculation!')
    return


def sanderfirst(command):
    """ Run first AMBER calculation, with the full model, for the optional optimization of the system with with HIGH and MEDIUM fixed"""

    # get the command for sander parallel execution
    para_exe = os.getenv('PARA_EXE')
    if para_exe is None:
       para_exe='mpirun -np'

    logwrt.writelog("Perform a full calculation of the entire system ...")

    # run sander
    if command[20]=='1':
        _command = shlex.split( 'pmemd.cuda -O -i real-sander-first.inp -o real-sander-first.out -p real.top -c real.crd -r real.rst' )
    elif command[20]=='2':
        _command = shlex.split( para_exe+' '+str(command[7])+' sander.MPI -O -i real-sander-first.inp -o real-sander-first.out -p real.top -c real.crd -r real.rst' )
    else:
        _command = shlex.split( 'sander -O -i real-sander-first.inp -o real-sander-first.out -p real.top -c real.crd -r real.rst' )
    subprocess.call( _command )

    # check the correct end of sander
    checksanderdone("real-sander-first.out",'real-sander-first')
    logwrt.writelog(' done!\n')

    return


def sandersecond(geometry,command):
    """ Run second AMBER calculation, with the SP calculation and the full model"""

    # get the command for sander parallel execution
    para_exe = os.getenv('PARA_EXE')
    if para_exe is None:
       para_exe='mpirun -np'

    if geometry.calculationType == 'M' or geometry.calculationType == 'ML':
        logwrt.writelog("Perform a single point calculation and collect energy and gradient ...")

    # run the second real sander without constrains
    if command[20]=='1':
        _command = shlex.split( 'pmemd.cuda -O -i real-sander-second.inp -o real-sander-second.out -p real.top -c real.rst -r real-second.rst')
    elif command[20]=='2':
        _command = shlex.split(para_exe+' '+str(command[7])+' sander.MPI -O -i real-sander-second.inp -o real-sander-second.out -p real.top -c real.rst -r real-second.rst')
    else:
        _command = shlex.split('sander -O -i real-sander-second.inp -o real-sander-second.out -p real.top -c real.rst -r real-second.rst')
    subprocess.call( _command )

    # check the correct end of sander
    checksanderdone("real-sander-second.out",'real-sander-second')

    # get energy and gradient from real second sander
    if command[20]!='0' :
         Fxyz_second=read_forcesA_pmemd('real-second')
         E_sandersecond=read_energyMM_out('real-sander-second.out')
    else:
         Fxyz_second=read_forcesA('forcedump.dat')
         E_sandersecond=read_energyMM()

    if geometry.calculationType == 'M' or geometry.calculationType == 'ML':
        logwrt.writelog(' done!\n')
        logwrt.writelog('MM energy is {0} Hartree\n'.format(E_sandersecond))

    return [E_sandersecond,Fxyz_second]


def sandermodelnoc(command):
    """ Run AMBER calculation, with the full model, for the optional optimization of the system with with HIGH and MEDIUM fixed"""

    # get the command for sander parallel execution
    para_exe = os.getenv('PARA_EXE')
    if para_exe == None:
       para_exe='mpirun -np'

    logwrt.writelog("\nPerform a single point calc. of the entire system with ZERO charge on the \nHigh layer, and collect energy and gradient ...")

    # run sander
    if command[20]=='1':
        _command = shlex.split('pmemd.cuda -O -i real-sander-second.inp -o real-sander-modelnoc.out -p real-modelnoc.top -c real.rst -r real-modelnoc.rst')
    elif command[20]=='2':
        _command = shlex.split(para_exe+' '+str(command[7])+' sander.MPI -O -i real-sander-second.inp -o real-sander-modelnoc.out -p real-modelnoc.top -c real.rst -r real-modelnoc.rst')
    else:
        _command = shlex.split('sander -O -i real-sander-second.inp -o real-sander-modelnoc.out -p real-modelnoc.top -c real.rst -r real-modelnoc.rst')
    subprocess.call( _command )

    # get energy and gradient from real modelnoc sander
    checksanderdone("real-sander-modelnoc.out",'real-sander-modelnoc')

    if command[20]!='0':
         Fxyz_modelnoc=read_forcesA_pmemd('real-modelnoc')
         E_modelnoc=read_energyMM_out('real-sander-modelnoc.out')
    else:
         Fxyz_modelnoc=read_forcesA('forcedump.dat')
         E_modelnoc=read_energyMM()

    logwrt.writelog(' done!\n')
    logwrt.writelog('MM energy is {0} Hartree\n'.format(E_modelnoc))
#    logwrt.writelog(logwrt.matrix_prettystring(np.array(Fxyz_modelnoc), ".6f"), 2)

    return [E_modelnoc,Fxyz_modelnoc]

def sandermodelH(geometry,command):

    logwrt.writelog("\nPerform a single point calc. of the High layer, with H-saturated bonds and \nZERO charge, and collect energy and gradient ...")

    # make model-H crd file
    # makemodelHcrd('model-H.crd',geometry)
    geometry.makemodelHcrd('model-H.crd')

    # run sander
    if (command[20])=='1':
      _command = shlex.split('pmemd.cuda -O -i model-H-sander.inp -o model-H.out -p model-H-noc.top -c model-H.crd -r model-H.rst')
    elif command[20]=='2': ## model-H cannot be run in parallel, there is only 1 residue!
      _command = shlex.split('sander -O -i model-H-sander.inp -o model-H.out -p model-H-noc.top -c model-H.crd -r model-H.rst')
    else:
      _command = shlex.split('sander -O -i model-H-sander.inp -o model-H.out -p model-H-noc.top -c model-H.crd -r model-H.rst')
    subprocess.call( _command )

    ## get energy and gradient from model-H sander
    checksanderdone("model-H.out",'model-h_sander')

    if command[20]!='0':
        Fxyz_modelH=read_forcesA_pmemd('model-H')
        E_modelH=read_energyMM_out('model-H.out')
    else:
        Fxyz_modelH=read_forcesA('forcedump.dat')
        E_modelH=read_energyMM()

    logwrt.writelog(' done!\n')
    logwrt.writelog('MM energy is {0} Hartree\n'.format(E_modelH))
#    logwrt.writelog(logwrt.matrix_prettystring(np.array(Fxyz_modelH), ".6f"), 2)

    return [E_modelH,Fxyz_modelH]


def chargezeromodel(filename,filename2,geometry):
    ## put to zero the charges of the high part of the sander topology
    filein=open(filename)
    top=[]
    ## grep a part from a topology file
    start=0
    end=0
    amb12=0
    for n, line in enumerate(filein, 1):
        top.append(line)
        try:
            point=line.split()
            if point[1].strip() == 'CHARGE':
                start=n
            if point[1].strip() == 'ATOMIC_NUMBER':
#                print "This is Amber12"
                amb12=1
                end=n
            if point[1].strip() == 'MASS' and amb12==0:
                end=n
        except:
            pass
    filein.close()
    ## transform in a list
    charge=[]
    for i in range(start+1,end-1):
        charge=charge+top[i].split()
    for i in range(geometry.NatomQM):
        charge[int(geometry.list_HIGH[i])-1]='0.00000000E+00'
    ## write the modelnoc topology
    modelnoc=open(filename2,'w')
    for i in range(start+1):
        modelnoc.write(top[i])
    for i in range(0,len(charge),5):
        try:
            modelnoc.write('%16.8e' % float(charge[i]))
            modelnoc.write('%16.8e' % float(charge[i+1]))
            modelnoc.write('%16.8e' % float(charge[i+2]))
            modelnoc.write('%16.8e' % float(charge[i+3]))
            modelnoc.write('%16.8e' % float(charge[i+4])+'\n')
        except:
            modelnoc.write('\n')
    for i in range(end-1,len(top)):
        modelnoc.write(top[i])
    modelnoc.close()
    return

def chargezero(filename,filename2):
    filein=open(filename)
    top=[]
    ## grep a part from a topology file
    amb12=0
    start=0
    end=0
    for n, line in enumerate(filein, 1):
        top.append(line)
        try:
            point=line.split()
            if point[1].strip() == 'CHARGE':
                start=n
            if point[1].strip() == 'ATOMIC_NUMBER':
#                print "This is Amber12!"
                amb12=1
                end=n
            if point[1].strip() == 'MASS' and amb12==0:
                end=n
        except:
            pass
    filein.close()
    ## tranform in a list
    for i in range(start+1,end-1):
        top[i]=top[i].replace("1","0")
        top[i]=top[i].replace("2","0")
        top[i]=top[i].replace("3","0")
        top[i]=top[i].replace("4","0")
        top[i]=top[i].replace("5","0")
        top[i]=top[i].replace("6","0")
        top[i]=top[i].replace("7","0")
        top[i]=top[i].replace("8","0")
        top[i]=top[i].replace("9","0")
    modelnoc=open(filename2,'w')
    ## write the final real-modelnoc.top
    for i in range(len(top)):
        modelnoc.write(top[i])
    modelnoc.close()
    return

def makemodelHcrd(filename,geometry):
    ## check the presence of the real.crd file in the working directory
    ## and make it if not found
    AXYZ_real=geometry.modelH
    ## check the existance of the file
    crd=open(filename,'w')
    crd.write('\n'+str(len(AXYZ_real[0]))+'\n')
    j=0
    for i in range(int(len(AXYZ_real[0])/2+1)):
        ## write the real.crd
        try:
            crd.write('%12.7f' % AXYZ_real[0][j]+'%12.7f' % AXYZ_real[1][j]+'%12.7f' % AXYZ_real[2][j])
            try:
                crd.write('%12.7f' % AXYZ_real[0][j+1]+'%12.7f' % AXYZ_real[1][j+1]+'%12.7f' % AXYZ_real[2][j+1]+'\n')
                j=j+2
            except:
                pass
        except:
            pass
    crd.close()
    #logwrt.writelog(filename+' created!')
    return


def read_crgA():
    ## take thecharges from model-H.top
    Camber=math.sqrt(332.0)
    ch=greptop('CHARGE','MASS','model-H.top')
    charges=[]
    try:
       for i in range(len(ch)):
          charges.append(float(ch[i])/Camber)
    except:
          logwrt.writelog('Amber version > 11 detected!\n')
          charges=[]
          ch=greptop('CHARGE','ATOMIC_NUMBER','model-H.top')
          for i in range(len(ch)):
             charges.append(float(ch[i])/Camber)
    CRG=np.array(charges)
    return CRG


def save_MM_amber_step(step,command):
    """ save gaussian output and chk file for last single point QM calculation
        in a common directory where QM data is stored """

    # for a frequency calculation there is no amber output to save
    if command[1] == 'freqxg': return

    # name of the gaussian calculation (common to chk and log file)
    logName = 'real-sander-first.out'
    # directory and files where to store QM results
    dirName = 'MM_data'
    allName = os.path.join(dirName,'real-sander-first.all')

    # at first step, remove directory from previous calculations and do nothing else
    if step == 0:
        if dirName in os.listdir('.'): shutil.rmtree(dirName)
        return

    # at second step, also create a new directory
    # (there is no need at first step because SP calculation do not run real-sander-first calculation)
    if step == 1: os.mkdir(dirName)

    # check if the file exists, otherwise print a warning to screen
    if not os.path.exists(logName):
        logwrt.writewarning( "Amber file {0} cannot be found: the {1} file will not be updated".format(logName, allName) )

    else:
        # write the content of the QM single point calc in the ALL file, decorated with
        # comments that highligh the step number of the QM single point
        with open(allName,'a') as ambtot:

            ambtot.write('='*80 + '\n')
            ambtot.write("Start full AMBER calculation of STEP : "+str(step)+"\n")
            ambtot.write(' || '*20 + '\n')
            ambtot.write(' \/ '*20 + '\n')
            with open(logName,"r") as amberOut:
                ambtot.write( amberOut.read() )
            ambtot.write(' /\ '*20 + '\n')
            ambtot.write(' || '*20 + '\n')
            ambtot.write("End full AMBER calculation of STEP : "+str(step)+"\n")
            ambtot.write('='*80 + '\n')
            ambtot.write('\n')


def read_forcesA(filename):
    ##read forces from fort.77 file
    c1=627.5095 # 1 au = 627.5095 kcal/mol.
    c2=0.5291772083 # 1 au = 0.529 A
    c3=c2/c1
    c4=c1/c2
    c5=1/(c1*c2)


    filename='forcedump.dat'
    length=int(subprocess.check_output('head -n1 '+filename, shell=True))
#    print length
    filein=subprocess.check_output('grep -A'+str(length)+' Total '+filename, shell=True)
    filein=filein.split(b'\n')
#    print filein
    fx=[]
    fy=[]
    fz=[]
    for i in range(1,length+1):
        row=filein[i].split()
        fx.append(float(row[0])*c3) # convert to au
        fy.append(float(row[1])*c3)
        fz.append(float(row[2])*c3)
    Fx=np.array(fx)
    Fy=np.array(fy)
    Fz=np.array(fz)

#    print "THE MM FORCES, 1st and last"
#    print fx[0]
#    print fz[length-1]
#    close=(filename)
#    os.remove(filename)
    Fxyz=[Fx,Fy,Fz]
    return Fxyz


### New routine to get forces when using pmend GPU code (OW,08/2015)
def read_forcesA_pmemd(filename):
    c1=627.5095
    c2=0.5291772083
    c3=c2/c1
    c4=c1/c2
    c5=1/(c1*c2)

    top=[]
    atom=[]

    if filename=='real-second':
         filein=input('real.top')
    else:
         filein=input(filename+'.top')
    for line in filein:
       top.append(line)
       try:
           point=line.split()
           if (point[1].strip() == 'ATOM_NAME'):
               start=filein.lineno()
           if (point[1].strip() == 'CHARGE'):
               end=filein.lineno()
       except:
           pass
    filein.close()
    for i in range(start+1,end-1):
        tmp1=top[i].strip()
        for j in range(0,len(tmp1),4):
            if tmp1[j:j+2]=='Na' or tmp1[j:j+2]=='Cl' or tmp1[j:j+2]=='Mg':
                atom.append(tmp1[j:j+2])
            else:
                atom.append(tmp1[j][:1])
    print( atom)
    print( "We found: ",len(atom)," Atoms")
    top=[]
    velox,veloy,veloz=[],[],[]
    filein=input(filename+'.rst')
    for line in filein:
       top.append(line)
    if len(atom) % 2 == 0 :
       pos=int((len(atom))/2)+2
    else:
       pos=int(round((len(atom))/2)+3)
    print( "POS=",pos,len(atom))
    for i in range(pos,pos*2-2):
       tmp1=top[i].split()
       print( tmp1,i)
       for j in range(0,len(tmp1),3):
          velox.append(float(tmp1[j]))
          veloy.append(float(tmp1[j+1]))
          veloz.append(float(tmp1[j+2]))
    print( "We found", len(atom)," Atoms and have ", len(velox), "velocity vectors")
    force=[]
    count=0
    # now compute the Forces from velocities
    forcex,forcey,forcez=[],[],[]
    factor=9.77756   # this factor for a timestep of 0.01 ps in amber!
    for i in range(len(atom)):
       if atom[i]=='C':
           forcex.append(velox[i]*factor*12.01*c3)
           forcey.append(veloy[i]*factor*12.01*c3)
           forcez.append(veloz[i]*factor*12.01*c3)
       elif atom[i]=='N':
           forcex.append(velox[i]*factor*14.01*c3)
           forcey.append(veloy[i]*factor*14.01*c3)
           forcez.append(veloz[i]*factor*14.01*c3)
       elif atom[i]=='H':
           forcex.append(velox[i]*factor*1.008*c3)
           forcey.append(veloy[i]*factor*1.008*c3)
           forcez.append(veloz[i]*factor*1.008*c3)
       elif atom[i]=='O':
           forcex.append(velox[i]*factor*16.0*c3)
           forcey.append(veloy[i]*factor*16.0*c3)
           forcez.append(veloz[i]*factor*16.0*c3)
       elif atom[i]=='S':
           forcex.append(velox[i]*factor*32.06*c3)
           forcey.append(veloy[i]*factor*32.06*c3)
           forcez.append(veloz[i]*factor*32.06*c3)
       elif atom[i]=='P':
           forcex.append(velox[i]*factor*30.97*c3)
           forcey.append(veloy[i]*factor*30.97*c3)
           forcez.append(veloz[i]*factor*30.97*c3)
       elif atom[i]=='Na':
           forcex.append(velox[i]*factor*22.99*c3)
           forcey.append(veloy[i]*factor*22.99*c3)
           forcez.append(veloz[i]*factor*22.99*c3)
       elif atom[i]=='Cl':
           forcex.append(velox[i]*factor*35.45*c3)
           forcey.append(veloy[i]*factor*35.45*c3)
           forcez.append(veloz[i]*factor*35.45*c3)
       elif atom[i]=='Mg':
           forcex.append(velox[i]*factor*24.305*c3)
           forcey.append(veloy[i]*factor*24.305*c3)
           forcez.append(veloz[i]*factor*24.305*c3)
       elif atom[i]=='F':
           forcex.append(velox[i]*factor*19.00*c3)
           forcey.append(veloy[i]*factor*19.00*c3)
           forcez.append(veloz[i]*factor*19.00*c3)
       else:
           logwrt.fatalerror("Atom "+atom[i]+" not implemented in forces_pmemd routine, please change the source code!")
    for i in range (len(atom)):
       print( forcex[i],forcey[i],forcez[i])

    Fx=np.array(forcex)
    Fy=np.array(forcey)
    Fz=np.array(forcez)

    Fxyz=[Fx,Fy,Fz]

    #now restore the orinigal rst file, delete the velocity entries
    if filename[:4]=='real':
        shutil.copyfile('real.crd', filename+'.rst')
    else:
        shutil.copyfile('model-H.crd', filename+'.rst')

    return Fxyz


def read_energyMM():
    ## get energies from forcedump file
    c1=627.5095
    filename='forcedump.dat'
    filein=subprocess.check_output('grep -A1 Energies '+filename, shell=True)
    filein=filein.split(b'\n')

    ene1=float(filein[1])/c1
#    filein.close()
    os.remove('forcedump.dat')
    return ene1


def read_energyMM_out(filename):
    ## (OW) get energies from the mdout file
    ## to be used with the pmemd code, where no forcedump.dat is created
    c1=627.5095
    filein=subprocess.check_output('grep -A1 Etot '+filename+'|head -n1', shell=True)
    print( "trying to energy from out-file ",filename,", this is the result:")
    print( filein,filein.split())
    filein=filein.split()
    ene1=float(filein[2])/c1
    return ene1


def skip_MMfirst(step, command, cobcom):
    # command=lists[6]
    ifskip=0
    if command[1] == 'freqxg':
        ifskip=1  # low optimization is skipped
        logwrt.writelog('This is a frequency run: the initial full MM calculation is skipped\n')
    elif step==0:
        ifskip=1  # low optimization is skipped
        logwrt.writelog('This is the 0th step: the initial full MM calculation is skipped\n')
    else:
        #logwrt.writelog('Cobram is looking for some coordinates to be numerically differentiated\n')
        Dnumber=0
        frozenQM=CBF.ReadCobramCommand(cobcom,'redundant','redundant.dat')
        for i in range(len(frozenQM)):
            splitted=frozenQM[i].split()
            if splitted[-1] in ['d','D']:
                Dnumber=Dnumber+1
        if Dnumber!=0:
            logwrt.writelog("Numerical 2nd derivative is requested for "+str(Dnumber)+" coordinates\n")
            logwrt.writelog('so for the first '+str(Dnumber)+' steps the low layer is not optimized\n')
            if step<Dnumber+1:
                logwrt.writelog('step= '+str(step)+': is in the range [1 , '+str(Dnumber)+']\n')
                logwrt.writelog('low layer is NOT optimized\n')
                ifskip=1
            else:
                logwrt.writelog('step= '+str(step)+': is out of the range [1 , '+str(Dnumber)+']\n')
                logwrt.writelog('low layer is optimized\n')
                ifskip=0
    return ifskip


def clean_MM_amber():
    """ clean up the run directory from all the files that have been used to run Amber and that
        are no longer needed at the end of the calculation """

    # list of files/directory to remove
    toRemove = ['real-sander-first.inp', "real-sander-first.out", "real-sander-second.inp", "real-sander-second.out",
                "real.rst", "real-second.rst", 'mdinfo', 'restrt',
                "real-sander-modelnoc.out", "real-modelnoc.top", "real-modelnoc.rst",
                "model-H-sander.inp", "model-H.out", "model-H-noc.top", "model-H.crd", "model-H.rst"]

    for f in toRemove:
        # if f is an existing file, remove it
        if os.path.isfile(f):
            os.remove(f)
        # if f is an existing directory, clean tree
        elif os.path.isdir(f):
            shutil.rmtree(f)

    # done
    return


def chargerebuildsander(filename,CRG_real):
    ## build the new real.top with the QM charges
    filein=open(filename)
    top=[]
    ## grep a part from a topology file
    for line in filein:
        top.append(line)
        try:
            point=line.split()
            if (point[1].strip() == 'CHARGE'):
                start=filein.lineno()
            if (point[1].strip() == 'MASS'):
                end=filein.lineno()
        except:
            pass
    filein.close()
    ## tranform in a list
    Camber=math.sqrt(332.0)
    charge=list(CRG_real*Camber)
    ## write the modelnoc topology
    modelnoc=open(filename,'w')
    for i in range(start+1):
        modelnoc.write(top[i])
    for i in range(0,len(charge),5):
        try:
            modelnoc.write('%16.8e' % float(charge[i]))
            modelnoc.write('%16.8e' % float(charge[i+1]))
            modelnoc.write('%16.8e' % float(charge[i+2]))
            modelnoc.write('%16.8e' % float(charge[i+3]))
            modelnoc.write('%16.8e' % float(charge[i+4])+'\n')
        except:
            modelnoc.write('\n')
    for i in range(end-1,len(top)):
        modelnoc.write(top[i])
    modelnoc.close()
    return


def modelHcharges(QM_Results):## put QM charges in model-H.top
    filein=input('model-H.top')
    top=[]
    ## grep a part from a topology file
    for line in filein:
        top.append(line)
        try:
            point=line.split()
            if (point[1].strip() == 'CHARGE'):
                start=filein.lineno()
            if (point[1].strip() == 'MASS'):
                end=filein.lineno()
        except:
            pass
    filein.close()
    ## tranform in a list
    Camber=math.sqrt(332.0)
    charge=QM_Results[2]*Camber
    ## write the model-H topology
    model_H=open('model-H.top','w')
    for i in range(start+1):
        model_H.write(top[i])
    for i in range(0,len(charge),5):
        try:
            model_H.write('%16.8e' % float(charge[i]))
            model_H.write('%16.8e' % float(charge[i+1]))
            model_H.write('%16.8e' % float(charge[i+2]))
            model_H.write('%16.8e' % float(charge[i+3]))
            model_H.write('%16.8e' % float(charge[i+4])+'\n')
        except:
            model_H.write('\n')
    for i in range(end-1,len(top)):
        model_H.write(top[i])
    model_H.close()
    return
