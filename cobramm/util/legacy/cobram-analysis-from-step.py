#! /usr/bin/env python
from fileinput import *
from math import *
from Numeric import * 
from time import *
from os import system
from os import listdir
from sys import *


def pdb_read(namepdb):
    pdbinp=input(namepdb)
    tmp_pdb=[]
    RESname_vec_all=['zero']
    RESnumber_vec_all=['zero']
    ATOMnumber_vec_all=['zero']
    RESname_vec=['zero']
    RESnumber_vec=['zero']
    ATOMnumber_vec=['zero']
    for line in pdbinp:
            tmp_pdb.append(line)
            splitted_pdb=line.split()
            ATOMnumber_vec_all.append(splitted_pdb[1]) 
            RESname_vec_all.append(splitted_pdb[3]) 
            RESnumber_vec_all.append(splitted_pdb[4]) 
    
    len_vec_all=len(ATOMnumber_vec_all)
    for i in range(1,len_vec_all):
        if RESnumber_vec_all[i]!=RESnumber_vec_all[i-1]:
            ATOMnumber_vec.append(ATOMnumber_vec_all[i])
            RESname_vec.append(RESname_vec_all[i])
            RESnumber_vec.append(RESnumber_vec_all[i])
            
    pdbinp.close()
    PDB_parsed=[ATOMnumber_vec_all,RESname_vec_all,RESnumber_vec_all,ATOMnumber_vec,RESname_vec,RESnumber_vec]
    return PDB_parsed


print ""
print "======================================="
print "COBRAM analysis_from_step starts"
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
                namedir="analysis_from_step_"+str_wanted_step
            else:
                str_wanted_step="last_step"
                namedir="analysis_from_"+str_wanted_step
            system("mkdir "+namedir)
        except:
            print "Please an integer number!"
            print" "
            kontrol=0
            
    fileout1=open(namedir+'/cobram.command','w')
    fileout3=open(namedir+'/real_layers.xyz','w')
    fileout4=open(namedir+'/model-H.top','w')
    fileout5=open(namedir+'/real.top','w')
    fileout6=open(namedir+'/real.crd','w')
    
    fileinp=input('cobram.out')
    tmp=[]
    STEP=[]
    top_c_start=[]
    top_c_end=[]
    gen_pres1='n'
    gen_pres2='n'
    
    for line in fileinp:
        tmp.append(line)
        try:
            splitted_line=line.split() 
            
            if (splitted_line[0].strip() == "CalculationType") :
                CT=splitted_line[1].strip()
                print CT
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
            if (splitted_line[0].strip() == "STEP" ) :
                STEP.append(int(splitted_line[1].strip()))
            if (splitted_line[0].strip() == 'G1' ) and (splitted_line[3].strip() == str(STEP[wanted_step]))  :
                start=int(fileinp.lineno())
                num=int(splitted_line[5])
            if (splitted_line[0].strip() == 'G3' ) and (splitted_line[3].strip() == str(STEP[wanted_step]))  :
                start_g_model_H=int(fileinp.lineno())
                num_g_model_H=int(splitted_line[5])
            if (splitted_line[0].strip() == 'C1' ) and (splitted_line[3].strip() == str(STEP[wanted_step]))  :
                start_c_real=int(fileinp.lineno())
                num_c_real=int(splitted_line[5])  
            if (splitted_line[0].strip() == 'C3' ) and (splitted_line[3].strip() == str(STEP[wanted_step]))  :
                start_c_model_H=int(fileinp.lineno())
                num_c_model_H=int(splitted_line[5])
            if (splitted_line[0].strip() == 'C8' ) and (splitted_line[3].strip() == str(STEP[wanted_step]))  :
                start_c_real_fromembedder=int(fileinp.lineno())
                num_c_real_fromembedder=int(splitted_line[5])
                
            if (splitted_line[0].strip() == 'list_HIGH' )   :
                start_list_HIGH=int(fileinp.lineno())
                num_list_HIGH=int(splitted_line[2])
                
            if (splitted_line[0].strip() == '!gaussian' ):  
                start_gaussian=int(fileinp.lineno())
            if (splitted_line[0].strip() == '?gaussian' ):  
                end_gaussian=int(fileinp.lineno())
            if (splitted_line[0].strip() == '!gen' ):  
                start_gen_keyword=int(fileinp.lineno())
                gen_pres1='y'
            if (splitted_line[0].strip() == '?gen' ):  
                end_gen_keyword=int(fileinp.lineno())
                gen_pres2='y'
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

list_HIGH=[]
for i in range(start_list_HIGH,start_list_HIGH+num_list_HIGH):
    row=tmp[i].split()
    list_HIGH.append( int(row[0]))


gaussian_tmp=[]
for i in range(start_gaussian+1,end_gaussian-1):
    tmp[i]=tmp[i].lower()
    splitted_gaussian=tmp[i].split()
    try: 
        if splitted_gaussian[0].find("%") ==-1 :
            gaussian_tmp.append(tmp[i].replace("force",""))
    except:
        gaussian_tmp.append('###')
        
gaussian=[]
for i in range(len(gaussian_tmp)):
    if gaussian_tmp[i]!='###':
        gaussian.append(gaussian_tmp[i])
    else:
        break
for i in range(len(gaussian_tmp)):
    if gaussian_tmp[i]!='###':
        cm=gaussian_tmp[i]
  
genkey=[]
if gen_pres1=='y' and gen_pres2=='y':
    for i in range(start_gen_keyword,end_gen_keyword-1):
        genkey.append(tmp[i]) 

## finds in tmp charges of 'real' and converts them in amber format
crg_real=[]
crg_real_readable=[]
for i in range(start_c_real,start_c_real+num_c_real):
    row=tmp[i].split()
    crg_real.append(float(row[0])*math.sqrt(332.0))
    crg_real_readable.append(float(row[0]))

CRG_real_fromembedder=[]
for i in range(start_c_real_fromembedder,start_c_real_fromembedder+num_c_real_fromembedder):
    row=tmp[i].split()
    CRG_real_fromembedder.append(float(row[0]))    
    
## finds in tmp charges of 'model-H' and converts them in amber format
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

for i in range(start5+1,top_c_start[1]+1):
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

for i in range(top_c_end[1]-1,end5-2):
    fileout5.write(tmp[i])
fileout5.close()


    
## finds in tmp the wanted geometry and converts it in crd format
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

AT=[]
x=[]
y=[]
z=[]
for i in range(start_g_model_H,start_g_model_H+num_g_model_H):
    row=tmp[i].split()
    AT.append(row[0])
    x.append(float(row[1]))
    y.append(float(row[2]))
    z.append(float(row[3]))
X=array(x)
Y=array(y)
Z=array(z)
AXYZ_model_H=[AT,X,Y,Z]


system('ambpdb -p '+namedir+'/real.top < '+namedir+'/real.crd>'+namedir+'/real.pdb')
system("grep -v 'TER' "+namedir+"/real.pdb|grep -v 'REMARK' |grep -v 'END' > "+namedir+"/temporary.pdb")
pdbinp=(namedir+"/temporary.pdb")
PDB_parsed=pdb_read(pdbinp)
ATOMnumber_vec_all=PDB_parsed[0]
RESname_vec_all=PDB_parsed[1]
RESnumber_vec_all=PDB_parsed[2]
ATOMnumber_vec=PDB_parsed[3]
RESname_vec=PDB_parsed[4]
RESnumber_vec=PDB_parsed[5]
  
print ' '
req_mem=raw_input('How many memory do you need? [500MB]')
if req_mem=='':
    req_mem='500MB'
req_proc=raw_input('How many processors do you need? [1]')
if req_proc=='':
    req_proc='1'
print ' '
  
system('mkdir '+namedir+'/analys_electro_DFP')
system('mkdir '+namedir+'/analys_electro_RFP')
fileout_script_DFP=open(namedir+'/analys_electro_DFP/electro_DFP_script4run.sh','w')
fileout_script_RFP=open(namedir+'/analys_electro_RFP/electro_RFP_script4run.sh','w')
fileout_script_DFP.write('## number_of_residues: '+str(len(RESnumber_vec))+'\n')
fileout_script_RFP.write('## number_of_residues: '+str(len(RESnumber_vec))+'\n')
for j in range(len(RESnumber_vec)):
    fileout7=open(namedir+'/electrostatic_DFP_tmp.com','w')
    fileout9=open(namedir+'/electrostatic_RFP_tmp.com','w')
    if j==0:
        commento_DFP='single point with no charges (DFP)\n'
        commento_RFP='single point with all charges (RFP)\n'
    else:
        RESname=RESname_vec[j]
        RESnumber=RESnumber_vec[j]
        commento_DFP='single point for residue '+RESname+' ' +str(RESnumber)+' (DFP)\n'
        commento_RFP='single point for residue '+RESname+' ' +str(RESnumber)+' (RFP)\n'
    fileout7.write('%chk=./electrostatic_DFP\n')
    fileout9.write('%chk=./electrostatic_RFP\n')
    fileout7.write('%rwf=./electrostatic_DFP\n')
    fileout9.write('%rwf=./electrostatic_RFP\n')
    fileout7.write('%int=./electrostatic_DFP\n')
    fileout9.write('%int=./electrostatic_RFP\n')
    fileout7.write('%d2e=./electrostatic_DFP\n')
    fileout9.write('%d2e=./electrostatic_RFP\n')
    fileout7.write('%mem='+req_mem+'\n')
    fileout9.write('%mem='+req_mem+'\n')
    fileout7.write('%nproc='+req_proc+'\n')
    fileout9.write('%nproc='+req_proc+'\n')
    for i in range(len(gaussian)):
        fileout7.write(gaussian[i])
        fileout9.write(gaussian[i])
    fileout7.write('charge\n')
    fileout9.write('charge\n')
    if j!=0:
        fileout7.write('guess=read\n')
        fileout9.write('guess=read\n')
    fileout7.write('\n')
    fileout9.write('\n')
    fileout7.write(commento_DFP)
    fileout9.write(commento_RFP)
    fileout7.write('\n')
    fileout9.write('\n')
    fileout7.write(cm)
    fileout9.write(cm)
    for i in range(len(AXYZ_model_H[0])):
        fileout7.write('%3s' % AXYZ_model_H[0][i]+' '+'%12.7f' % AXYZ_model_H[1][i]+' '+'%12.7f' % AXYZ_model_H[2][i]+' '+'%12.7f' % AXYZ_model_H[3][i]+' '+ '\n') 
        fileout9.write('%3s' % AXYZ_model_H[0][i]+' '+'%12.7f' % AXYZ_model_H[1][i]+' '+'%12.7f' % AXYZ_model_H[2][i]+' '+'%12.7f' % AXYZ_model_H[3][i]+' '+ '\n')
    fileout7.write('\n')
    fileout9.write('\n')
    if len(genkey) !=0:
        for i in range(len(genkey)):
            fileout7.write(genkey[i])
            fileout9.write(genkey[i])
        fileout7.write('\n')
        fileout9.write('\n')

    
    if j!=0:
        start=int(ATOMnumber_vec[j])-1
        init=int(ATOMnumber_vec[j])
        if RESnumber==RESnumber_vec[-1]:
            end=int(ATOMnumber_vec[-1])
            link_end='\n'
        else:
            end=int(ATOMnumber_vec[j+1])-1
            

    if j!=0 :
        for k in range(start,end):
            if k+1 not in list_HIGH:
                string7='%12.7f' % AXYZ[1][k]+' '+'%12.7f' % AXYZ[2][k]+' '+'%12.7f' % AXYZ[3][k]+' '+'%16.8f\n' % CRG_real_fromembedder[k]
                fileout7.write(string7) 
        fileout7.write('\n')
        for k in range(len(CRG_real_fromembedder)):
            string9='%12.7f' % AXYZ[1][k]+' '+'%12.7f' % AXYZ[2][k]+' '+'%12.7f' % AXYZ[3][k]+' '+'%16.8f\n' % CRG_real_fromembedder[k]
            if  k+1 not in list_HIGH and (k+1<init or k+1>end):
                fileout9.write(string9) 
        fileout9.write('\n')
    else:
        string7='%12.7f' % 0.0+' '+'%12.7f' % 0.0+' '+'%12.7f' % 0.0+' '+'%16.8f\n' % 0.0
        fileout7.write(string7) 
        fileout7.write('\n') 
        for k in range(len(CRG_real_fromembedder)):
            if k+1 not in list_HIGH :
                string9='%12.7f' % AXYZ[1][k]+' '+'%12.7f' % AXYZ[2][k]+' '+'%12.7f' % AXYZ[3][k]+' '+'%16.8f\n' % CRG_real_fromembedder[k]
                fileout9.write(string9) 
        fileout9.write('\n')
    fileout7.close()
    fileout9.close()
    system('mv '+namedir+'/electrostatic_DFP_tmp.com '+namedir+'/analys_electro_DFP/electrostatic_DFP_'+str(j)+'.com')
    system('mv '+namedir+'/electrostatic_RFP_tmp.com '+namedir+'/analys_electro_RFP/electrostatic_RFP_'+str(j)+'.com')
    fileout_script_DFP.write('g03 <electrostatic_DFP_'+str(j)+'.com> electrostatic_DFP_'+str(j)+'.log\n')
    fileout_script_RFP.write('g03 <electrostatic_RFP_'+str(j)+'.com> electrostatic_RFP_'+str(j)+'.log\n')

fileout_script_DFP.write('rm -rf electrostatic_DFP.* \n')
fileout_script_RFP.write('rm -rf electrostatic_RFP.* \n')
fileout_script_DFP.close()
fileout_script_RFP.close()
system('chmod a+rwx  '+namedir+'/analys_electro_DFP/electro_DFP_script4run.sh')
system('chmod a+rwx  '+namedir+'/analys_electro_RFP/electro_RFP_script4run.sh')

##vdw anal
system('mkdir '+namedir+'/analys_vdw')
fileout_script=open(namedir+'/analys_vdw/vdw_script4run.sh','w')
fileout_script.write('## number_of_residues: '+str(len(RESnumber_vec))+'\n')

for i in range(1,len(RESnumber_vec)):
    fileout_script.write('anal -O -i vdw_'+str(i)+'.inp -o vdw_'+str(i)+'.out -p real.top  -c real.crd -ref real.crd \n')
    fileout_vdw_tmp=open('vdw_tmp.inp','w')
    fileout_vdw_tmp.write('Analysis_of_the_system\n')
    fileout_vdw_tmp.write('   1 0 0 0 600 1\n')
    fileout_vdw_tmp.write('   0 0.0 0.0 0.0 0.0\n')
    fileout_vdw_tmp.write('   1 0 3 1 100 1\n')
    fileout_vdw_tmp.write('   888.0 2.0 1.2 1\n')
    fileout_vdw_tmp.write('   1 99.0 99.0 99.0 1.0 99.0 1.0 99.0 1.0 99.0\n')
    fileout_vdw_tmp.write('ENERGY\n')
    fileout_vdw_tmp.write('COBRAM VDW Interaction Energy between HIGH and '+RESnumber_vec[i]+' ' + RESname_vec[i]+'\n')
    fileout_vdw_tmp.write('RES '+ str(i)+' '+ str(i) +' \n')
    fileout_vdw_tmp.write('END\n')
    fileout_vdw_tmp.write('HIGH region (no H link)\n')
    for j in range(len(list_HIGH)):
        fileout_vdw_tmp.write('ATOM '+str(list_HIGH[j])+' \n') 
    fileout_vdw_tmp.write('END\n')
    fileout_vdw_tmp.write('END\n')
    fileout_vdw_tmp.write('STOP\n')
    fileout_vdw_tmp.close()
    system('mv vdw_tmp.inp '+namedir+'/analys_vdw/vdw_'+str(i)+'.inp')

system('cp '+namedir+'/real.top '+namedir+'/real.crd '+namedir+'/analys_vdw')
fileout_script.close()
system('chmod a+rwx  '+namedir+'/analys_vdw/vdw_script4run.sh')
print ""

maskcontrol1=0
while maskcontrol1==0:
    asked=raw_input('Do you want analyze a particular selection? n/y [n]')
    if asked=='':
        asked='n'
    if asked=='y':
        maskcontrol2=0
        maskcounter=1
        while maskcontrol2==0:
            print "Insert the mask you want in the ambmask format"
            mask=raw_input('mask: ')
            mask=mask.strip("'")
            system("ambmask -p "+namedir+"/real.top -c "+namedir+"/real.crd -out pdb -find '"+mask+"' > "+namedir+"/masked_"+str(maskcounter)+".tmp")
            
            masktmpfile=input(namedir+"/masked_"+str(maskcounter)+".tmp")
            err=0
            rem=0
            msktmp=[]
            msklab=[]
            for line in masktmpfile:
##                print line
                msktmp.append(line)
                if rem==0:
                    msklab.append(line)
                try:
                    if line.split()[0]=='Error' or (line.split()[0]=='Warning:' and line.split()[1]=='no' and line.split()[2]=='atoms' and line.split()[3]=='selected'):
                        err=1
                except:
                    pass
                try:
                    if line.split()[0]=='REMARK':
                        rem=1
                except:
                    pass
            masktmpfile.close()
            
            if err==1:
                print'Error:'
                for i in range(len(msktmp)):
                    print msktmp[i].strip()
                print ' '
                maskcontrol2=0
            else:
                system('grep -A '+str(len(ATOMnumber_vec_all))+' "REMARK" '+namedir+"/masked_"+str(maskcounter)+".tmp |grep -v 'REMARK' |grep -v 'END' > "+namedir+"/masked_"+str(maskcounter)+".pdb")
                mkspdbinp=(namedir+"/masked_"+str(maskcounter)+".pdb")
                print mkspdbinp
                mksPDB_parsed=pdb_read(mkspdbinp)
                mksATOMnumber_vec_all=mksPDB_parsed[0]
                mksRESname_vec_all=mksPDB_parsed[1]
                mksRESnumber_vec_all=mksPDB_parsed[2]
                mksATOMnumber_vec=mksPDB_parsed[3]
                mksRESname_vec=mksPDB_parsed[4]
                mksRESnumber_vec=mksPDB_parsed[5]
                mdkdat=open(namedir+'/masked_'+str(maskcounter)+'.dat','w')
                for i in range(len(mksRESnumber_vec)):
                    if mksRESnumber_vec[i]!='zero':
                        mdkdat.write(mksRESnumber_vec[i])
                mdkdat.close()
                
                fileout_script_mskDFP=open(namedir+'/analys_electro_DFP/mask_'+str(maskcounter)+'_electro_DFP_script4run.sh','w')
                fileout_script_mskRFP=open(namedir+'/analys_electro_RFP/mask_'+str(maskcounter)+'_electro_RFP_script4run.sh','w')
                for i in range(len(msklab)):
                    fileout_script_mskDFP.write('## '+msklab[i]+'\n')
                    fileout_script_mskRFP.write('## '+msklab[i]+'\n')
                
                fileout_script_mskDFP.write('g03 <electrostatic_DFP_0.com> electrostatic_DFP_0.log\n')
                fileout_script_mskRFP.write('g03 <electrostatic_RFP_0.com> electrostatic_RFP_0.log\n')
                for i in range(len(RESnumber_vec)):
                    if RESnumber_vec[i] in mksRESnumber_vec and RESnumber_vec[i]!='zero':
                        j=RESnumber_vec[i]
                        fileout_script_mskDFP.write('g03 <electrostatic_DFP_'+str(j)+'.com> electrostatic_DFP_'+str(j)+'.log\n')
                        fileout_script_mskRFP.write('g03 <electrostatic_RFP_'+str(j)+'.com> electrostatic_RFP_'+str(j)+'.log\n')
                
                fileout_script_mskDFP.close()
                fileout_script_mskRFP.close()
                print 'The scripts for the chosen mask has been created'
                maskcounter=maskcounter+1
                maskcontrol3=0
                while maskcontrol3==0:
                    askredo=raw_input('Do you want analyze another selection? n/y [n]')
                    if askredo=='':
                        askredo='n'
                    if askredo=='n':
                        maskcontrol1=1
                        maskcontrol2=1
                        maskcontrol3=1
                    elif askredo=='y':
                        maskcontrol2=0
                        maskcontrol3=1
                    else:
                        maskcontrol3=0
    elif asked=='n':
        maskcontrol1=1
    else:
        maskcontrol1=0
    
    
    
    
    
    
    
    
    
    
    
    







print "COBRAM analysis_from_step ended"
print "======================================="
print ""



    
    
