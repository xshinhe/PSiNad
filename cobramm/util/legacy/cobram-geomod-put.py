#! /usr/bin/env python
##puts modified HIGH-MEDIUM geometry in real one and write Molden and Amber formats 
from fileinput import *
from math import *
from Numeric import * 
from time import *
from os import system
from sys import *

print ""
print "======================================="
print "COBRAM geometry-modifier (put) starts"
print ""


fileinp1=input('list.dat')
fileinp2=input('unmodified_real.xyz')
fileinp3=input('HM.com')

fileout1=open('HM.xyz','w')
fileout2=open('new_real.xyz','w')
fileout3=open('new_real.crd','w')


tmp1=[]
for line in fileinp1:
    tmp1.append(line)
    try:
        splitted_line=line.split() 
        if (str(splitted_line[0]) == "list_HIGH_MEDIUM") :
            start_list_HM=int(fileinp1.lineno())
            num_list_HM=int(splitted_line[2])  
    except:
        pass

tmp2=[] 
for line in fileinp2:
    tmp2.append(line)
        

tmp3=[]
for line in fileinp3:
    tmp3.append(line)
    try:
        splitted_line=line.split() 
        if (str(splitted_line[0]) == "HIGH-MEDIUM") :
            start_new_HM=int(fileinp3.lineno())
    except:
        pass

list_HM=[]
for i in range(start_list_HM,start_list_HM+num_list_HM):
    line_splitted=tmp1[i].split()
    list_HM.append(int(line_splitted[0]) )
    
    
AT=[]
x=[]
y=[]
z=[]
for i in range(2,len(tmp2)):
    row=tmp2[i].split()
    AT.append(row[0])
    x.append(float(row[1]))
    y.append(float(row[2]))
    z.append(float(row[3]))
X=array(x)
Y=array(y)
Z=array(z)
AXYZ_new=[AT,X,Y,Z]
AXYZ_old=(AT,X,Y,Z)


AT=[]
x=[]
y=[]
z=[]
for i in range(start_new_HM+2,start_new_HM+2+len(list_HM)):
    row=tmp3[i].split()
    AT.append(row[0])
    x.append(float(row[2]))
    y.append(float(row[3]))
    z.append(float(row[4]))
X=array(x)
Y=array(y)
Z=array(z)
AXYZ_HM=[AT,X,Y,Z]


fileout1.write(str(len(list_HM))+'\n')
fileout1.write('\n')
for i in range(len(list_HM)):
    fileout1.write( str(AXYZ_HM[0][i])+' '+str(AXYZ_HM[1][i])+' '+str(AXYZ_HM[2][i])+' '+str(AXYZ_HM[3][i])+'\n')


for i in range(len(list_HM)):
    j=list_HM[i]-1
    AXYZ_new[0][j]=AXYZ_HM[0][i]
    AXYZ_new[1][j]=AXYZ_HM[1][i]
    AXYZ_new[2][j]=AXYZ_HM[2][i]
    AXYZ_new[3][j]=AXYZ_HM[3][i]
    
fileout2.write(str(len(AXYZ_new[0]))+'\n')
fileout2.write('\n')
for i in range(len(AXYZ_new[0])):
    fileout2.write( str(AXYZ_new[0][i])+' '+str(AXYZ_new[1][i])+' '+str(AXYZ_new[2][i])+' '+str(AXYZ_new[3][i])+'\n')
   
    
    
    
crd=fileout3
crd.write('\n'+str(len(AXYZ_new[1]))+'\n')
j=0
for i in range(len(AXYZ_new[1])/2+1):
    ## write the real.crd
    try:
        crd.write('%12.7f' % AXYZ_new[1][j]+'%12.7f' % AXYZ_new[2][j]+'%12.7f' % AXYZ_new[3][j])
        try:
            crd.write('%12.7f' % AXYZ_new[1][j+1]+'%12.7f' % AXYZ_new[2][j+1]+'%12.7f' % AXYZ_new[3][j+1]+'\n')
            j=j+2
        except:
            pass
    except:
        pass
crd.close() 



fileinp1.close()
fileinp2.close()
fileinp3.close()
fileout1.close()
fileout2.close()
fileout3.close()




print ""
print "COBRAM geometry-modifier (put) ended"
print "'new_real.crd' was successfully created"
print "======================================="

print ""
