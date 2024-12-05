#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" Utility to get results from COBRAM MDV 
"""

## OW April 2018: added functionality for traj. occupaptions and populations

__author__ = "Piero Altoè <piero.altoe@gmail.com>"
__date__ = "5th February 2008"
__version__ = ': 0.2 $'
__credits__ = """Guido van Rossum, for an excellent programming language."""


##------------------------------------------------------------------------------
## Functions
##------------------------------------------------------------------------------

from os import remove

##------------------------------------------------------------------------------
## Functions
##------------------------------------------------------------------------------

def writeout(write,file):
    out=open(file,'a')
    out.write(write)
    out.close()

##------------------------------------------------------------------------------
## Files names
##------------------------------------------------------------------------------

file1='dynamic.out'

##------------------------------------------------------------------------------
## Required output
##------------------------------------------------------------------------------

req=raw_input('Required data : g=geometry v=velocity e=energy (tot,pot,kin)\n'+
              '                gr=gradient eAll=QM pot energy of all states\n'+
              '                eq=QMMM Energy of all states p=electronic \n'+
              '                populations o=trajectory state occupations \n'  )

##------------------------------------------------------------------------------
## variable names
##------------------------------------------------------------------------------

if req.upper() == 'G':
    choice='Geometry' 
elif req.upper() == 'V':
    choice='Velocity'
elif req.upper() == 'E':
    choice='ETotEPotEKin'
elif req.upper() == 'G+Q':
    choice='Charges'
elif req.upper() == 'GR':
    choice='Gradient'
elif req.upper() == 'EALL':
    choice='ALLstateEnergy'
elif req.upper() == 'EQ':
    choice='ALLstateEQMMM'
elif req.upper() == 'P':
    choice='Population'
elif req.upper() == 'O':
    choice='Occupation'
else:
    exit('No resonable choice\n')


start,stop,time=[],[],[]
t=0.0
a=0

##------------------------------------------------------------------------------
## Main program
##------------------------------------------------------------------------------

## open input file
try:
    filein=open(file1)
except:
    exit('ERROR: File not found !!!\n \n')
input=filein.read()
filein.close()
## split in lines
input=input.split('\n')
## grep the positions of the informations
for i in range(len(input)-1):
    try:
       line=input[i].strip().split()
       if line[0] == '!'+choice :
          start.append(i+1)
       if line[0] == '?'+choice :
          stop.append(i)
       if line[0] == 'timestep:' :
          if a == 0 :
              pass
          else:
              t=t+float(line[1])/float('41.341373337')
          time.append(t)
          a=a+1
    except:
       pass
## remove old file if present
try:
    remove(req.upper()+'.dat')
except:
    pass
## write the output file
for i in range(len(start)):
    if req.upper() == 'G' or req.upper() == 'GR' or req.upper() == 'V':
        writeout(str(stop[0]-start[0])+'\n\n',req.upper()+'.dat')
    if req.upper() == 'EALL' or req.upper() == 'EQ' :
        if i == 0 :
            writeout(' %8.3f ' %time[i]+'   ',req.upper()+'.dat')
            for j in range(start[i],stop[i]):
                line=input[j].strip().split()
                writeout(line[1]+'     ',req.upper()+'.dat')
            writeout('\n',req.upper()+'.dat') 
        else:
            if input[start[i]].strip().split() == input[start[i-1]].strip().split():
                pass
            else:
                try:
                    writeout(' %8.3f ' %time[i]+'   ',req.upper()+'.dat')
                except:
                    writeout(' %8.3f ' %time[i-1]+'   ',req.upper()+'.dat')
                for j in range(start[i],stop[i]):
                    line=input[j].strip().split()
                    writeout(line[1]+'     ',req.upper()+'.dat')
                writeout('\n',req.upper()+'.dat') 

    elif req.upper() == 'O' or req.upper() == 'P' :
        if i == 0 :
            writeout(' %8.3f ' %time[i]+'   ',req.upper()+'.dat')
            for j in range(start[i],stop[i]):
                line=input[j].strip().split()
                writeout(line[1]+'     ',req.upper()+'.dat')
            writeout('\n',req.upper()+'.dat')
        else:
                try:
                    writeout(' %8.3f ' %time[i]+'   ',req.upper()+'.dat')
                except:
                    writeout(' %8.3f ' %time[i-1]+'   ',req.upper()+'.dat')
                for j in range(start[i],stop[i]):
                    line=input[j].strip().split()
                    writeout(line[1]+'     ',req.upper()+'.dat')
                writeout('\n',req.upper()+'.dat')



    else:
        if req.upper() == 'E':
            for j in range(start[i],stop[i]):
                writeout(' %8.3f ' %time[i]+'   '+input[j]+'\n',req.upper()+'.dat')
        else:
            for j in range(start[i],stop[i]):
                writeout(input[j]+'\n',req.upper()+'.dat')


