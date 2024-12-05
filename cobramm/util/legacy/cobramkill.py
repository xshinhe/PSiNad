#! /usr/bin/env python
from fileinput import *
from math import *
from Numeric import * 
from time import *
from os import system
from os import listdir
import os
import shutil
import gzip
from sys import *
from sys import exit

print ""
print "======================================="
print "The cobram job running here will be "
print "          !! KILLED !!"



def GZIP(file):
    r_file = open(file, 'r')
    w_file = gzip.GzipFile(file + '.gz', 'w', 9)
    w_file.write(r_file.read())
    w_file.flush()
    w_file.close()
    r_file.close()
    os.unlink(file) #We don't need the file now
    
def G03_getpid(namefile):
    PID=''
    try:
        file=open(namefile,'r')
        file.close()
    except:
        print'The file '+namefile+' is not present in this directory'
        pass
    os.system('head -10 "'+namefile+'" |grep "PID=" > kill.tmp')
    file=open('kill.tmp','r')
    for line in file:
        if line.find('PID=')!=-1:
            tmp=line.split()
            for i in range(len(tmp)):
                if tmp[i]=='PID=':
                    k=i
            PID=int(tmp[k+1].strip('.'))
    file.close()
    os.remove('kill.tmp')
    return PID
########################################################
tobekilled=[]
PID=G03_getpid('cobram.out')
print "\nThe PID of the current process is "+str(PID)+"\n"
tobekilled.append(str(PID))

if 'geometry.log' in os.listdir('.'):
    PID_geom=G03_getpid('geometry.log')
    tobekilled.append(str(PID_geom))
if 'gaussian-QM.log' in os.listdir('.'):
    PID_gaussian_QM=G03_getpid('gaussian-QM.log')
    tobekilled.append(str(PID_gaussian_QM))
cont=0
while cont==0:
    answ1=raw_input("Do you want to proceed? (y/n) [y]")
    if answ1=='':
        cont=1
    elif answ1=='y':
        cont=1
    elif answ1=='n':
        cont=2
    else:
        cont=0
        
if cont==2:
    print "\nThe process will NOT be terminated"
    
##    exit("=======================================")
    
elif cont==1:
    for i in range(len(tobekilled)):
        system('kill -9 '+tobekilled[i])

cont=0
while cont==0:
    answ1=raw_input("\nDo you want to remove temporary files? (y/n) [y]")
    if answ1=='':
        cont=1
    elif answ1=='y':
        cont=1
    elif answ1=='n':
        cont=2
    else:
        cont=0
        
if cont==1:

    toremove=['cobram.command','real_layers.xyz','real.top','model-H.top', 'model-H-sander.inp','model-H.top','model-H-noc.top','real.crd','real-sander-second.inp','real-sander-first.inp','real-sander-first.out','gaussian-QM.rwf','gaussian-QM.d2e','gaussian-QM.int','gaussian-QM.chk','gaussian-QM.log','gaussian-QM.com','mdinfo','restrt','real.crd','real.rst','real-second.rst','real-sander-second.out','real-sander-modelnoc.out','model-H.out','model-H.crd','control.dat','geometry.com','geometry.d2e','geometry.int','geometry.rwf','gradient.dat','fort.44','fort.77','fort.99','cobram.Err','frozenmmatoms.pdb','real-modelnoc.top','optxg_control.dat','optxg_control2.dat','test.inp','test.out']
    for i in range(len(toremove)):
        if toremove[i]  in os.listdir('.'):
            os.remove(toremove[i])
    tmpdir=os.listdir('.')
    for i in range(len(tmpdir)):
        if  tmpdir[i].startswith('Gau-') and tmpdir[i].endswith('.scr'):
            os.remove(tmpdir[i])
    if 'tmpgarbage' in os.listdir('.'):
        shutil.rmtree('tmpgarbage')
    if 'geometry.chk' in os.listdir('.'):
        GZIP('geometry.chk')
print "\nThe process was terminated" 
print "and all files have been removed"
print "======================================="
        
        
