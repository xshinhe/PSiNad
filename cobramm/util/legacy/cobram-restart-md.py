#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" Utility to restart COBRAM MDV 
"""

__author__ = "Piero Altoè <piero.altoe@gmail.com>"
__date__ = "18th July 2008"
__version__ = ': 0.1 $'
__credits__ = """Guido van Rossum, for an excellent programming language."""


##------------------------------------------------------------------------------
## Functions
##------------------------------------------------------------------------------

from fileinput import *
from os import system
from os import path


##------------------------------------------------------------------------------
## Files names
##------------------------------------------------------------------------------

file1='out/dynamic.out'
file2='tmpgarbage/real.crd'
if path.isfile(file2):
    print "Found actual real.crd in tmpgarbage"
else:
    file2='real.crd'
    print "Using real.crd from current directory"
file3='real_layers.xyz'
file4='restart/velocity.dat'
file5='restart/real.crd'
file6='restart/AtomLink.xyz'

##------------------------------------------------------------------------------
## Variables
##------------------------------------------------------------------------------

a=0
Gstart,Gstop,Vstart,Vstop,all=[],[],[],[],[]   
x,y,z=[],[],[]
vx,vy,vz=[],[],[]
num=[]
X,Y,Z=[],[],[]

##------------------------------------------------------------------------------
## Main
##------------------------------------------------------------------------------

GV=input(file1)
for line in GV :   
        con=line.split()
        all.append(con)
        try:
                if con[0] == '!Geometry' :
                        Gstart.append(a)
                if con[0] == '?Geometry' :
                        Gstop.append(a)
                if con[0] == '!Velocity' :
                        Vstart.append(a)
                if con[0] == '?Velocity' :
                        Vstop.append(a)
        except:
                pass
        a=a+1
#print "Gstart", Gstart
#print "Gstop",Gstop
#print "Vstart", Vstart
#print "Vstop", Vstop
GV.close()

natom=Gstop[len(Gstop)-2]-Gstart[len(Gstart)-2]
print "lenVstart", len(Vstart)
print "lenGstart", len(Gstart)
if len(Vstart)!=len(Gstart):
    print "ERROR - inconsistency in dynamics.out: the number of geometries does not match"
    print "        the number of veloctiy entries, please check!"
    exit()
for i in range(len(all)):
        if i == Gstart[len(Gstart)-1] :
                for j in range(1,natom):
                        x.append(float(all[i+j][1]))
                        y.append(float(all[i+j][2]))
                        z.append(float(all[i+j][3]))
        if i == Vstart[len(Vstart)-1] :
                for j in range(1,natom):
                        vx.append(all[i+j][0])
                        vy.append(all[i+j][1])
                        vz.append(all[i+j][2])

a=0

LAYER=input(file3)
for line in LAYER :
        con=line.split()
        if (con[5] == 'M') or (con[5] == 'H'):
                num.append(a)
        a=a+1
LAYER.close()

a=0
CRD=input(file2)
for line in CRD:
        if a != 1 and a != 0:
                con=line.split()
                X.append(float(con[0]))
                Y.append(float(con[1]))
                Z.append(float(con[2]))
                try:
                        X.append(float(con[3]))
                        Y.append(float(con[4]))
                        Z.append(float(con[5]))
                except:
                        pass
        a=a+1
CRD.close()

system('mkdir restart')
system('cp -r real_layers.xyz cobram.command WORK.wfu work.wfu WFU/WORK.wfu WFU/work.wfu wchk real.top model-H.top cobram-sef imomap.dat fort.11 restart >/dev/null')
system('echo "1" >restart/NSTEP')
system('mkdir restart/out;awk -vRS="STEP" \'NR>1{print s RT} {s=$0}\' ORS="" out/dynamic.out >restart/out/dynamic.out')
system('echo " RESTART ----------" >>restart/out/dynamic.out')
VEL=open(file4,'w')
for i in range(len(vx)):
        VEL.write('%16.10f' %float(vx[i])+'%16.10f' %float(vy[i])+'%16.10f' %float(vz[i])+'\n')
VEL.close() 

AL=open(file6,'w')
for i in range(len(num),natom-1):
        AL.write(' H   %12.7f' %x[i]+'%12.7f' %y[i]+'%12.7f\n' %z[i])
AL.close()

a=0
NCRD=open(file5,'w')
NCRD.write('\n'+str(len(X))+'\n')
for i in range(0,len(X),2):
        if i in num:
                X[i]=x[a]
                Y[i]=y[a]
                Z[i]=z[a]
                a=a+1
        if i+1 in num:
                try:
                        X[i+1]=x[a]
                        Y[i+1]=y[a]
                        Z[i+1]=z[a]
                        a=a+1
                except:
                        pass
        NCRD.write('%12.7f' %X[i]+'%12.7f' %Y[i]+'%12.7f' %Z[i])
        try:
                NCRD.write('%12.7f' %X[i+1]+'%12.7f' %Y[i+1]+'%12.7f\n' %Z[i+1])
        except:
                NCRD.write('\n')
NCRD.close()

print 'DONE'

