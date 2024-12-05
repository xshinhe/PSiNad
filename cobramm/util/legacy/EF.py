#! /usr/bin/env python

## importing external functions and modules
from Numeric import *
from fileinput import *
from math import *
from time import *
from os import system
import copy

k=0
try:
    attempt=open('cobram.out','r')
    attempt.close()
    k=1
except:
    pass
    k=0
    print '\nThe file named cobram.out is not present in the present directory\n'
    
if k==1: 
    system('grep "E12 "  -A2 cobram.out >  cobram.out.tmp')
    system('grep "I1 "   -A2 cobram.out >> cobram.out.tmp')
    system('grep "I2 "   -A2 cobram.out >> cobram.out.tmp')    
    file_cobram_out=input('cobram.out.tmp')
    tmp=[]
    E12=[]
    I1=[]
    I2=[]
    for line in file_cobram_out:
        tmp.append(line)
        try:
            splitted_line=line.split() 
            if (splitted_line[0].strip() == 'E12'):
                E12.append(file_cobram_out.lineno())
            if (splitted_line[0].strip() == 'I1'):
                I1.append(file_cobram_out.lineno())
            if (splitted_line[0].strip() == 'I2'):
                I2.append(file_cobram_out.lineno())
        except:
            pass
    file_cobram_out.close()  
    
    OUT=open('ef.dat','w')
          
    print ('%4s %10s  %10s %10s' %('step','Etot','Fmax','Frms'))
    OUT.write('%4s %10s  %10s %10s\n' %('step','Etot','Fmax','Frms'))

    for i in range(len(E12)):
        tmpE=tmp[E12[i]].strip()
        tmpFmax=tmp[I1[i]].strip() 
        tmpFrms=tmp[I2[i]].strip() 

        try:
            toprint=('%4u   %10f %10s %10s' %(i,float(tmpE),tmpFmax,tmpFrms))
        except:
            toprint=('%4u   %10s %10s %10s' %(i,tmpE,tmpFmax,tmpFrms))
        print toprint
        OUT.write(toprint+'\n')
    OUT.close()
    system('rm cobram.out.tmp') 
