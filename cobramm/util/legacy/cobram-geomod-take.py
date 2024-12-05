#! /usr/bin/env python
##Takes HIGH-MEDIUM geometry for you to modify it
from fileinput import *
from math import *
from Numeric import * 
from time import *
from os import system
from sys import *


print ""
print "======================================="
print "COBRAM 3.0 geometry-modifier (taker) starts"
print ""
control=0
while control==0:
    print ""
    kontrol=0
    while kontrol==0:
        in_wanted_step=raw_input('What step do you need (-1 for last one)? [-1]')
        if in_wanted_step=="":
            in_wanted_step="-1"
        try:
            kontrol=1
            wanted_step=int(in_wanted_step)
            if wanted_step !=-1:
                str_wanted_step=str(wanted_step)
                namedir="geomod_inputs_from_step_"+str_wanted_step
            else:
                str_wanted_step="last_step"
                namedir="geomod_inputs_from_"+str_wanted_step
            system("mkdir "+namedir)
        except:
            print "Please an integer number!"
            print" "
            kontrol=0
            
    fileout1=open(namedir+'/cobram.command','w')
##    fileout2=open(namedir+'/real-sander-first.inp','w')
    fileout3=open(namedir+'/real_layers.xyz','w')
    fileout6=open(namedir+'/unmodified_real.crd','w')
    fileout7=open(namedir+'/unmodified_real.xyz','w')
    fileout8=open(namedir+'/unmodified_HM.xyz','w')
    fileout9=open(namedir+'/unmodified_HM.com','w')
    fileout10=open(namedir+'/list.dat','w')

    
    fileinp=input('cobram.out')
    tmp=[]
    STEP=[]
    top_c_start=[]
    top_c_end=[]
    for line in fileinp:
        tmp.append(line)
        try:
            splitted_line=line.split() 
            if (splitted_line[0].strip() == "CalculationType") :
                CT=splitted_line[1].strip()
            if (splitted_line[0].strip() == "R1a") :
                start1=(int(fileinp.lineno()))
            if (splitted_line[0].strip() == "R1b" ) :
                end1=(int(fileinp.lineno()))
            if (splitted_line[0].strip() == "R2a") :
                start2=(int(fileinp.lineno()))
            if (splitted_line[0].strip() == "R2b" ) :
                end2=(int(fileinp.lineno()))
            if (splitted_line[0].strip() == "R3a") :
                start3=(int(fileinp.lineno()))
            if (splitted_line[0].strip() == "R3b" ) :
                end3=(int(fileinp.lineno()))
            if (splitted_line[0].strip() == "R4a") :
                start4=(int(fileinp.lineno()))
            if (splitted_line[0].strip() == "R4b" ) :
                end4=(int(fileinp.lineno()))
            if (splitted_line[0].strip() == "R5a") :
                start5=(int(fileinp.lineno()))
            if (splitted_line[0].strip() == "R5b" ) :
                end5=(int(fileinp.lineno()))
            if (splitted_line[0].strip() == "%FLAG") and (splitted_line[1].strip() == "CHARGE") :
                top_c_start.append(int(fileinp.lineno()))
            if (splitted_line[0].strip() == "%FLAG") and (splitted_line[1].strip() == "MASS") :
                top_c_end.append(int(fileinp.lineno()))                
            if (splitted_line[0].strip() == "list_HIGH" ) :
                start_list_H=int(fileinp.lineno())
                num_list_H=int(splitted_line[2])          
            if (splitted_line[0].strip() == "list_MEDIUM" ) :
                start_list_M=int(fileinp.lineno())
                num_list_M=int(splitted_line[2])  
            if (splitted_line[0].strip() == "STEP" ) :
                STEP.append(int(splitted_line[1].strip()))
            if (splitted_line[0].strip() == 'G1' ) and (splitted_line[3].strip() == str(STEP[wanted_step]))  :
                start=int(fileinp.lineno())
                num=int(splitted_line[5])
            if (splitted_line[0].strip() == 'G8' ) and (splitted_line[3].strip() == str(STEP[wanted_step]))  :
                start_HM=int(fileinp.lineno())
                num_HM=int(splitted_line[5])
            if (splitted_line[0].strip() == 'C1' ) and (splitted_line[3].strip() == str(STEP[wanted_step]))  :
                start_c_real=int(fileinp.lineno())
                num_c_real=int(splitted_line[5])  
            if (splitted_line[0].strip() == 'C3' ) and (splitted_line[3].strip() == str(STEP[wanted_step]))  :
                start_c_model_H=int(fileinp.lineno())
                num_c_model_H=int(splitted_line[5])
        except:
            pass
    fileinp.close()
    
    if wanted_step > STEP[-1]:
        print "The step you want does not exist: "
        print " the number you entered before inot in range[0-"+str(STEP[-1])+"]"
        system("rm -rf "+namedir)
        control=0
    else:
        print "The step you want is in range[0-"+str(STEP[-1])+"]"
        print " OK: I will write datas from step "+str(STEP[wanted_step])
        control=1

totSTEPs=len(STEP)-1
for i in range(start1+1,end1-2):
    fileout1.write(tmp[i])
fileout1.close()
##for i in range(start2+1,end2-2):
##    fileout2.write(tmp[i])
##fileout2.close()
for i in range(start3+1,end3-2):
    fileout3.write(tmp[i])
fileout3.close()

if CT=='HM' or CT=='HML' or CT=='HL':
    fileout4=open(namedir+'/model-H.top','w')
    ## topology of model-H.top
    crg_model_H=[]
    for i in range(start_c_model_H,start_c_model_H+num_c_model_H):
        row=tmp[i].split()
        crg_model_H.append(float(row[0])*math.sqrt(332.0))
    for i in range(start4+1,top_c_start[0]+1):
        fileout4.write(tmp[i])
    count=divmod(len(crg_model_H),5)
    for i in range(count[0]):
        j=5*i
        ## write the model_H charges
        fileout4.write('%16.8E' % crg_model_H[j]+'%16.8E' %crg_model_H[j+1]\
            +'%16.8E' % crg_model_H[j+2]+'%16.8E' % crg_model_H[j+3]\
            +'%16.8E' % crg_model_H[j+4]+'\n')
    j=count[0]*5
    if count[1]==4:
        fileout4.write('%16.8E' % crg_model_H[j]+'%16.8E' %crg_model_H[j+1]\
            +'%16.8E' % crg_model_H[j+2]+'%16.8E' % crg_model_H[j+3]+'\n')
    elif count[1]==3:
        fileout4.write('%16.8E' % crg_model_H[j]+'%16.8E' %crg_model_H[j+1]\
            +'%16.8E' % crg_model_H[j+2]+'\n')
    elif count[1]==2:
        fileout4.write('%16.8E' % crg_model_H[j]+'%16.8E' %crg_model_H[j+1]+'\n')
    elif count[1]==1:
        fileout4.write('%16.8E' % crg_model_H[j]+'\n')
        
    for i in range(top_c_end[0]-1,end4-2):
        fileout4.write(tmp[i])
    fileout4.close()

if CT=='HM' or CT=='HML' or CT=='HL' or CT=='ML' or CT=='M':
    fileout5=open(namedir+'/real.top','w')
    ## topology of real.top
    crg_real=[]
    for i in range(start_c_real,start_c_real+num_c_real):
        row=tmp[i].split()
        crg_real.append(float(row[0])*math.sqrt(332.0))
    for i in range(start5+1,top_c_start[-1]+1):
        fileout5.write(tmp[i])
    count=divmod(len(crg_real),5)
    for i in range(count[0]):
        j=5*i
        ## write the real charges
        fileout5.write('%16.8E' % crg_real[j]+'%16.8E' %crg_real[j+1]\
            +'%16.8E' % crg_real[j+2]+'%16.8E' % crg_real[j+3]\
            +'%16.8E' % crg_real[j+4]+'\n')
    j=count[0]*5
    if count[1]==4:
        fileout5.write('%16.8E' % crg_real[j]+'%16.8E' %crg_real[j+1]\
            +'%16.8E' % crg_real[j+2]+'%16.8E' % crg_real[j+3]+'\n')
    elif count[1]==3:
        fileout5.write('%16.8E' % crg_real[j]+'%16.8E' %crg_real[j+1]\
            +'%16.8E' % crg_real[j+2]+'\n')
    elif count[1]==2:
        fileout5.write('%16.8E' % crg_real[j]+'%16.8E' %crg_real[j+1]+'\n')
    elif count[1]==1:
        fileout5.write('%16.8E' % crg_real[j]+'\n')
    for i in range(top_c_end[-1]-1,end5-2):
        fileout5.write(tmp[i])
    fileout5.close()


## finds in tmp the wanted real geometry and converts it in crd format
AT=[]
x=[]
y=[]
z=[]
for i in range(start,start+num):
    row=tmp[i].split()
    AT.append(row[0])
    x.append(float(row[1]))
    y.append(float(row[2]))
    z.append(float(row[3]))
X=array(x)
Y=array(y)
Z=array(z)
AXYZ=[AT,X,Y,Z]
crd=fileout6
crd.write('\n'+str(len(AXYZ[1]))+'\n')
j=0
for i in range(len(AXYZ[1])/2+1):
    ## write the real.crd
    try:
        crd.write('%12.7f' % AXYZ[1][j]+'%12.7f' % AXYZ[2][j]+'%12.7f' % AXYZ[3][j])
        try:
            crd.write('%12.7f' % AXYZ[1][j+1]+'%12.7f' % AXYZ[2][j+1]+'%12.7f' % AXYZ[3][j+1]+'\n')
            j=j+2
        except:
            pass
    except:
        pass
crd.close()   


## write real geometry in Molden readable format
fileout7.write(str(num)+'\n')
fileout7.write(' \n')
for i in range(start,start+num):
    fileout7.write(tmp[i])
fileout7.close()


## write a list file 
list_HM=[]
list_H=[]
list_M=[]
k=0
for i in range(start_list_H-1,start_list_H+num_list_H):
    fileout10.write(tmp[i])
    line_splitted=tmp[i].split()

    if k>0:
        list_HM.append(int(line_splitted[0]) )
        list_H.append(int(line_splitted[0]) )
    k=k+1
    
fileout10.write(' \n')

k=0
for i in range(start_list_M-1,start_list_M+num_list_M):
    fileout10.write(tmp[i]) 
    line_splitted=tmp[i].split()
    if k>0:
        list_HM.append(int (line_splitted[0]) )
        list_M.append(int(line_splitted[0]) )
    k=k+1
        
fileout10.write(' \n')   

list_HM.sort()

fileout10.write('list_HIGH_MEDIUM : '+str(len(list_HM))+'\n') 
for i in range (len(list_HM)):
    fileout10.write(str(list_HM[i])+' \n')
    
fileout10.write(' \n')     
fileout10.close()   


## write wanted HIGH-MEDIUM geometry in Molden readable format
fileout8.write(str(num_HM)+'\n')
fileout8.write(' \n')
for i in range(start_HM,start_HM+num_HM):
    fileout8.write(tmp[i])
fileout8.close()
    
## write wanted HIGH-MEDIUM geometry in Gaussview readable format
fileout9.write('%chk=HM.chk \n')
fileout9.write('%mem=100MB \n')
fileout9.write('#p oniom(hf:am1) \n')
fileout9.write(' \n')
fileout9.write('HIGH-MEDIUM geometry \n')
fileout9.write(' \n')
fileout9.write('0 1 0 1 0 1\n')
k=0
for i in range(start_HM,start_HM+num_HM):
    line_splitted=tmp[i].split()
    
    if list_HM[k] in list_H:
        label="H"
    else:
        label="M"
        
    fileout9.write(line_splitted[0]+'   0   '+line_splitted[1]+'   '+line_splitted[2]+'   '+line_splitted[3]+'   '+label+'\n')
    k=k+1
fileout9.write(' ')  
fileout9.close()


system('cp '+namedir+'/unmodified_HM.com '+namedir+'/HM.com')



print ""
print "COBRAM geometry-modifier (taker) ended"
print ""
print "After you modify 'HM.com' use ''cobram-geomod-put.py'' to create 'new_real.crd'"
print "======================================="
print ""
