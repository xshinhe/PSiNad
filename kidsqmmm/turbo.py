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


import numpy as np
from math import *
from fileinput import *
from time import *
import os
import shutil
import gzip
from os import system
from sys import exit
import CBF
# import commands
import turbo
import logwrt

def prepare(lists,step):
    thisdir=os.listdir('.')
    k=0
    if step==0:
        if 'turboinp' in thisdir:
            logwrt.writelog('The turboinp directory was found ')
            turboindir=os.listdir('./turboinp')
            if len(turboindir)!=0:
                system('cp turboinp/* .')
                system('cp control control_init')
                system('mv turboinp turboinp_0')
                logwrt.writelog('The files in turboinp have been copied in the main directory')
            else:
                logwrt.fatalerror('The turboinp directory is empty: please provide!')
        else:
            logwrt.fatalerror('The turboinp directory is not present in the main directory: please provide!')
        if 'turbo_orb' in thisdir:
            system('rm -rf  turbo_orb')
        system('mkdir turbo_orb')
    else:
        system('cp control_init control ')

## read gen or cards before charges from cobram.command
    gen=[]
    keys=[]
    prepared=[keys,gen]
    return prepared

def makeinp(dftb_input,splitted,CRG_emb,gen,lists,step):
    covcost=0.52917720830008323378
    command=lists[6]
    guess=command[6]
    Nproc=command[7]
    chMED=[]
    AXYZ=splitted[1]
    natom=len(splitted[1][0])

    if command[1] not in ['optxgp','cip','tsp','freqxgp','mdvp','ircp']:
        logwrt.writelog('\n'+' '*12+'*************************\n'+
                          ' '*12+'ENTER INPUT GENERATOR FOR\n'+
                          ' '*12+'        TURBOMOLE        \n'+
                          ' '*12+'*************************\n')

    logwrt.writelog('Constructing of the coord file')
    coord_file=open('coord','w')
    coord_file.write('$coord\n')
    for i in range(natom):
        a=AXYZ[0][i].strip().lower()
        tmp=(' %16.8f' % (AXYZ[1][i]/covcost)+ ' %16.8f' %(AXYZ[2][i]/covcost)+' %16.8f' % (AXYZ[3][i]/covcost)+' %5s' %a)
        coord_file.write(tmp+'\n')
    coord_file.write('$end')
    coord_file.close()

    if len(CRG_emb)!=0 and lists[8] not in ['H','M','ML']:
        logwrt.writelog('Appending the atomic point charges to the control file')
        system('sed -i \'$ d\' control')
        control_file=open('control','a')
        control_file.write('$point_charges selfenergy \n')
        for i in range(len(CRG_emb)):
            xcoord=float(splitted[2][1][i])/covcost
            ycoord=float(splitted[2][2][i])/covcost
            zcoord=float(splitted[2][3][i])/covcost
            control_file.write('%14.8f' % xcoord +'%14.8f' % ycoord+'%14.8f' % zcoord+'%16.8f\n' % CRG_emb[i])
#            control_file.write('%12.7f' % splitted[2][1][i]+'%12.7f' % splitted[2][2][i]+'%12.7f' % splitted[2][3][i]+'%16.8f\n' % CRG_emb[i])

        # add movable point charges to control file (OW 3/2019)

        if len(CRG_emb)!=0 and command[60] != '1':
              chMED=[]
              control_file.write('$pointval fld fmt=xyz geo=point \n')
              for i in range(len(CRG_emb)):
                    if lists[4][i] in lists[1]:
                         xcoord=float(splitted[2][1][i])/covcost
                         ycoord=float(splitted[2][2][i])/covcost
                         zcoord=float(splitted[2][3][i])/covcost
                         control_file.write('%14.8f' % xcoord+'%14.8f' % ycoord+'%14.8f\n' % zcoord)
                         chMED.append(CRG_emb[i])
        control_file.write('$end')
        control_file.close()
        control_file.close()

    if command[1] not in ['optxgp','cip','tsp','freqxgp','mdvp','ircp']:
        logwrt.writelog('\n'+' '*14+'********************\n'+
                             ' '*14+'EXIT INPUT GENERATOR\n'+
                             ' '*14+'********************\n\n')

    return chMED

def launch(lists,step):
    command=lists[6]
    guess=command[6]
    Nproc=command[7]
    # OW avoid stupid input error and prevent from computing rubbish
    if lists[6][190]=='1' and lists[6][191]!='0':
       logwrt.fatalerror('The use of egrad in combination with rdgrad or ricc2 does not make sense, please choose the correct program to compute the gradient!')

    if lists[6][191]=='2':
       system('dscf >turbo.out')
       system('ricc2 >>turbo.out')
#       system('rdgrad >>turbo.out')
    if lists[6][191]=='1':
       system('ridft >turbo.out')
       system('rdgrad >>turbo.out')
    if lists[6][191]=='0':
       print( "SOME VARIABLES")
       system('echo $PARNODES')
       system('echo $PARA_ARCH')
       system('which egrad')
       system('echo $PATH')
       system('dscf >turbo.out')
       if lists[6][190]=='1':
         system('egrad >>turbo.out')
       else:
         system('grad >>turbo.out')


#    system('cp test.out turbo.out')  # Olli Extreme Debugging

def turboEneGradCrg(lists,splitted,CRG_emb_scaled,step,chMED):

    logwrt.writelog('\n'+' '*6+'***************************************\n'+
                         ' '*6+'ENTER COLLECT ENERGIES, GRADIENTS, etc.\n'+
                         ' '*6+'             TURBOMOLE                 \n'+
                         ' '*6+'***************************************\n\n')

    natom=len(splitted[1][0])

    if (step == 0):
        system('cat turbo.out >turboALL.out')
    else:
        system('cat turbo.out >>turboALL.out')

## take energy
##takegradient
    enefile=input('gradient')
    tmp=[]
    gradch=[[],[],[]]
    for line in  enefile:
        tmp.append(line.strip())
    enefile.close()
    system('rm gradient energy')
    tmp1=tmp[1].split()
#    if lists[6][190]=='0' and lists[6][191]=='0':   #OW get DSCF ground state energy

    try:
           turboenergy=[]
           grp=commands.getoutput('grep "Total energy" turbo.out')
           grp=grp.split('\n')
           for i in range(len(grp)):
                grp2=grp[i].split()
                turboenergy.append(float(grp2[2]))
           logwrt.writelog('Found excited state energies from EGRAD, we will continue with those')
    except:
           logwrt.writelog('Taking energy from the Turbomole energy-file')
           turboenergy=[]
           turboenergy.append(float(tmp1[6]))

    x=[]
    y=[]
    z=[]
    for i in range(2+natom,2+2*natom):
        gradxyz=tmp[i].strip().split()
        x.append(-float(gradxyz[0].replace('D','E')  ))
        y.append(-float(gradxyz[1].replace('D','E')))
        z.append(-float(gradxyz[2].replace('D','E')))
    gradient=[np.array(x),np.array(y),np.array(z)]
##takecharge  (put to 0.1 to avoid Division by zero in CBF charge-rally (OW 2013)
    charges=[]
    for i in range(natom):
        charges.append(0.1)
    ch=np.array(charges)

#get the charge-gradients of movable region from tf.xyz (OW 3/2019)

    if len(CRG_emb_scaled)!=0 and lists[6][60] != '1':
       tffile=input('tf.xyz')
       tf=[]
       for line in  tffile:
           tf.append(line.strip())
       tffile.close()
       for i in range(len(tf)):
           if tf[i].strip().startswith('# cartesian coordinates')==1:
               NatomQM=len(lists[3])
               NatomMedium=len(splitted[4][1])
               a=0
               for j in range(NatomMedium):
                   element=tf[i+j+1].split()
#                   print "ELEMENT=",element
                   gradch[0].append(float(element[3].replace('D', 'E'))*(-chMED[a]))
                   gradch[1].append(float(element[4].replace('D', 'E'))*(-chMED[a]))
                   gradch[2].append(float(element[5].replace('D', 'E'))*(-chMED[a]))
                   a=a+1


## take dipole in Debye
    tout=input('turbo.out')
    tmpo=[]
    for line in  tout:
        tmpo.append(line.strip())
    enefile.close()

    if lists[6][191]=="2":

      for i in range(len(tmpo)):
            if tmpo[i].strip().startswith('Analysis of relaxed')==1:
                count=i

      dipx=float(tmpo[count+7].strip().split()[1]   )
      dipy=float(tmpo[count+8].strip().split()[1]   )
      dipz=float(tmpo[count+9].strip().split()[1]   )
      norm=float(tmpo[count+11].strip().split()[5]   )

    else:

      for i in range(len(tmpo)):
            if tmpo[i].strip().startswith('| dipole moment | =')==1:
                count=i

      dipx=float(tmpo[count-4].strip().split()[3]   )
      dipy=float(tmpo[count-3].strip().split()[3]   )
      dipz=float(tmpo[count-2].strip().split()[3]   )
      norm=float(tmpo[count].strip().split()[5]   )

    dipole=[dipx,dipy,dipz,norm]

    Selfenergy=0.0

    turboResults=[turboenergy,gradient,ch,Selfenergy,dipole,gradch]
#    print "THE TURBOMOLE RESULTS:"
#    print turboResults
    logwrt.writelog('\n'+' '*6+'**************************************\n'+
                      ' '*6+'EXIT COLLECT ENERGIES, GRADIENTS, etc.\n'+
                      ' '*6+'**************************************\n')

    return turboResults

def clean(step):
    system('mkdir turboinp')
    list_save=['mos','alpha','beta','coord','basis','auxbasis']
    list_remove=['control','turbo.out','tf.xyz']
    thisdir=os.listdir('.')
    for i in range(len(list_save)):
        if list_save[i] in thisdir:
            system('cp '+list_save[i]+' turboinp')
    system('cp control_init turboinp/control')
    system('cp mos coord basis auxbasis turboinp')
    system('cp control turboinp/control_full')      # easier for later DFT/MRCI computations
    system('tar -czvf turboinp.tar.gz turboinp')
    system('mv turboinp.tar.gz turbo_orb/turboinp_'+str(step)+'.tar.gz')
#    system('rm -rf turboinp')
    for i in range(len(list_remove)):
        if list_remove[i] in thisdir:
            system('rm '+list_remove[i])
    system('rm -rf control.old*')




