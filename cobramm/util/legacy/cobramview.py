#! /usr/bin/env python
## utility to create molden crd and pdb files from cbram.out
## Version 1.0 31/1/2006, Piero Altoe'

##create molden & vmd output
from fileinput import *
import Numeric
from time import *
from os import system
from sys import *

## starting time
start=time()
print ""
print "====================================================="
print "COBRAM molden & vmd creator starts"
print ""
## input filename
##filename=raw_input('Insert pdb filename (default cobram.molden)\n')

## default filename
##if ( filename == ''):
filename='cobram.molden'
filename2='geometry.log'
##print "_______________________________________"
##print ""
print "The file cobram.molden will be generated"
print ""

output=[]
listH=[]
listM=[]
listL=[]
axyz=[]
fxyz=[]
ene=[]
ch=[]
Fmax=[]
Frms=[]
Dmax=[]
Drms=[]
OT=''
## read all cobram.out
filein=open('cobram.out')
output=filein.read()
filein.close()
## split in lines
output=output.split('\n')
i=0
for i in range(len(output)):
    element=output[i].split(':')
    if (output[i]. startswith("CalculationType")) :
        CT=output[i].split()[1]
    if (output[i]. startswith("OptimizationType")) :
        OT=output[i].split()[1]
    if (element[0] == 'R5a '):
        topstart=i+2
    if (element[0] == 'list_HIGH '):
        topend=i-4
    ## list definition
    if (element[0] == 'list_HIGH '):
        j=0
        for j in range(int(element[1])):
            listH.append(int(output[i+j+1].strip()))
        i=i+j
    elif (element[0] == 'list_MEDIUM '):
        j=0
        for j in range(int(element[1])):
            listM.append(int(output[i+j+1].strip()))
        i=i+j
    elif (element[0] == 'list_LOW '):
        j=0
        for j in range(int(element[1])):
            listL.append(int(output[i+j+1].strip()))
        i=i+j
    ## real geometry
    elif ( element[0] == 'G1 '):
        stepn=element[1]
        stepi=int(stepn.split()[1].strip())
        if stepi!=0:
            atoms=int(element[2])
            for j in range(int(element[2])):
                axyz.append(output[i+j+1])
            i=i+j
    ## forces medium high
    elif (( element[0] == 'F10 ') or ( element[0] == 'F7 ') )and stepi!=0:
        j=0
        for j in range(int(element[2])):
            fxyz.append(output[i+j+1])
        i=i+j
    ## total energy, max and rms forces and displacements
    elif ( element[0] == 'E12 ')and stepi!=0:
        ene.append(output[i+1])
    elif ( element[0] == 'I1 ')and stepi!=0:
        Fmax.append(output[i+1])
    elif ( element[0] == 'I2 ')and stepi!=0:
        Frms.append(output[i+1])
    elif ( element[0] == 'I3 ')and stepi!=0:
        Dmax.append(output[i+1])
    elif ( element[0] == 'I4 ')and stepi!=0:
        Drms.append(output[i+1])
    ## charges real
    elif ( element[0] == 'C1 ')and stepi!=0:
        j=0
        for j in range(int(element[2])):
            ch.append(output[i+j+1])
        i=i+j

## number of cicles, the 0th data are discarder
step=int(stepn.split()[1].strip())

## create temporary topology
top=open('tmp.top','w')
for i in range(topstart,topend):
    top.write(output[i]+'\n')
top.close()

output=[0]

x=[]
y=[]
z=[]
j=0
## create temporary crd file
for i in range(atoms):
    coord=axyz[i].split()
    x.append(float(coord[1]))
    y.append(float(coord[2]))
    z.append(float(coord[3]))
AXYZ_real=[x,y,z]
crd=open('tmp.crd','w')
crd.write('****\n')
crd.write(str(atoms)+'\n')
for i in range(atoms/2+1):
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
coord=[]
vmd=[]
k=0.0
for i in range(len(axyz)):
    coord=axyz[i].split()
    vmd.append(float(coord[1]))
    vmd.append(float(coord[2]))
    vmd.append(float(coord[3]))
crd=open('cobram.crd','w')
crd.write('\n')
for j in range(len(vmd)/(atoms*3)):
    for i in range(atoms*3):
        crd.write('%8.3f' % vmd[i+(atoms*3*j)])
        k=(i+1)%10
        if ( k == 0 ):
            crd.write('\n')
    crd.write('\n')

## create pdb
if CT!='H':
    print "creating cobram.pdb cobram.crd"
    system('ambpdb -p tmp.top < tmp.crd > cobram.pdb') 
system('rm tmp.top tmp.crd')


## write geometry of medium high in molden file
if OT == 'irc':
    print "writing only converged points along the IRC"
    atom_names=[]
    tmp=[]
    lists=Numeric.sort(Numeric.array((listH+listM)))-1
    natom=len(lists)
    blocks=int((3*natom+2)/5)
    if (3*natom+2) % 5 != 0:
        blocks=blocks+1
    gaussian=open(filename2,'r')
    for line in gaussian:
        tmp.append(line)
    gaussian.close()
    for i in range(len(tmp)):
        if tmp[i].find('Z-MATRIX') != -1:
            for j in range(natom):
                atom_names.append(tmp[i+j+3].split()[2])
            break
    for i in range(len(tmp)):
        if tmp[i].find('# OF POINTS ALONG THE PATH') != -1:
            irc_points=int(tmp[i].split()[7])+1
        if tmp[i].find('Summary of reaction path following') != -1 :
            molden=open(filename,'w')
            for j in range(irc_points):
                atom_counter=0
                energy=float(tmp[i+4+j].split()[1])
                molden.write(str(natom)+'\n'+str(energy)+'\n')
                for k in range(blocks):
                    tmp1=tmp[i+4+j+k*(irc_points+1)].split()
                    if k == 0:
                        molden.write(atom_names[atom_counter]+' %14.8f' %float(tmp1[3])+' %14.8f' %float(tmp1[4])+' %14.8f' %float(tmp1[5])+'\n')
                        atom_counter=atom_counter+1
                    elif (k-1) % 3 == 0 and k != 0:
                        molden.write(atom_names[atom_counter]+' %14.8f' %float(tmp1[1])+' %14.8f' %float(tmp1[2])+' %14.8f' %float(tmp1[3])+'\n')
                        atom_counter=atom_counter+1
                        if atom_counter == natom:
                           break
                        molden.write(atom_names[atom_counter]+' %14.8f' %float(tmp1[4])+' %14.8f' %float(tmp1[5]))
                        atom_counter=atom_counter+1
                        if atom_counter == natom:
                           molden.write(' %14.8f' %float(tmp1[1])+'\n')
                           break
                    elif (k-1) % 3 == 1:
                        molden.write(' %14.8f' %float(tmp1[1])+'\n'+atom_names[atom_counter]+' %14.8f' %float(tmp1[2])+' %14.8f' %float(tmp1[3])+' %14.8f' %float(tmp1[4])+'\n')
                        atom_counter=atom_counter+1
                        if atom_counter == natom:
                           break
                        molden.write(atom_names[atom_counter]+' %14.8f' %float(tmp1[5]))
                        atom_counter=atom_counter+1
                        if atom_counter == natom:
                           molden.write(' %14.8f' %float(tmp1[1])+' %14.8f' %float(tmp1[2])+'\n')
                           break
                    elif (k-1) % 3 == 2:
                        molden.write(' %14.8f' %float(tmp1[1])+' %14.8f' %float(tmp1[2])+'\n'+atom_names[atom_counter]+' %14.8f' %float(tmp1[3])+' %14.8f' %float(tmp1[4])+' %14.8f' %float(tmp1[5])+'\n')
                        atom_counter=atom_counter+1
			if atom_counter == natom:
                           break
            molden.close()
else:
    molden=open(filename,'w')
    molden.write(' [MOLDEN FORMAT]\n'+
                ' [GEOMETRIES] (XYZ)\n')
    passo=len(axyz)/step
    lists=Numeric.sort(Numeric.array((listH+listM)))-1
    natom=len(lists)
    #print passo
    #print len(axyz),step
    
    for i in range(0,len(axyz),passo):
        molden.write(str(natom)+'\n\n')
    
        for j in range(natom):
            if len(ch)!=0:
                molden.write(str(axyz[lists[j]+i]).strip()+'   '+str(ch[lists[j]+i])+'\n')
            else:
                molden.write(str(axyz[lists[j]+i]).strip()+'   '+'0.0000\n')
                
    ## write energy max and rms forces and displacement in molden file
    molden.write(' [GEOCONV]\n'+
                 ' energy\n')
    for i in range(len(ene)):
        molden.write(ene[i]+'\n')
    molden.write(' max-force\n')
    for i in range(len(Fmax)):
        molden.write(Fmax[i]+'\n')
    molden.write(' rms-force\n')
    for i in range(len(Frms)):
        molden.write(Frms[i]+'\n')
    molden.write(' max-step\n')
    for i in range(len(Dmax)):
        molden.write(Dmax[i]+'\n')
    molden.write(' rms-step\n')
    for i in range(len(Drms)):
        molden.write(Drms[i]+'\n')
    
    ## write forces in molden file
    molden.write(' [FORCES]\n')
    passo=len(fxyz)/step
    k=1
    for i in range(0,len(fxyz),passo):
        molden.write('point '+str(k)+'\n')
        molden.write(str(natom)+'\n')
        k=k+1
        for j in range(natom):
            molden.write(str(fxyz[j+i]).strip()+'\n')
    
    molden.close()

end=time()
total=end-start
print ""
print "Total time: "+'%4i' % total+" seconds"
print "COBRAM molden & vmd creator ends"
print "====================================================="
print ""
system('molden -l -A cobram.molden &')
