#! /usr/bin/env python
from fileinput import *
from math import *
from Numeric import * 
from time import *
from os import system
from os import listdir
from sys import *
from sys import exit
import fnmatch

def log_check(filename):
    fileinp=input(filename)
    ck=0
    tmp=[]
    for line in fileinp:
        tmp.append(line)
        try:
            splitted_line=line.split() 
            if (splitted_line[0].strip() == "Normal") and  (splitted_line[1].strip() == "termination"):
                ck=1
        except:
            pass
    fileinp.close()
    return ck

def log_read(filename):
    fileinp=input(filename)
    RESname='ref'
    RESnumber='0'
    tmp=[]
    for line in fileinp:
        tmp.append(line)
        try:
            splitted_line=line.split() 
            if (splitted_line[0].strip() == "single") and  (splitted_line[1].strip() == "point") and (splitted_line[2].strip() != "with") :
                RESname=(splitted_line[4])
                RESnumber=(int(splitted_line[5]))
            if (splitted_line[0].strip() == "SCF") and (splitted_line[1].strip() == "Done:"):
                SCF_ene=(float(splitted_line[4]))
            if (splitted_line[0].strip() == "Self") and (splitted_line[1].strip("energy") == "") and (splitted_line[2].strip() == "of") and (splitted_line[3].strip() == "the"):
                SELF_En_of_Ch=(float(splitted_line[6]))
        except:
            pass
    results=[RESname,RESnumber,SCF_ene,SELF_En_of_Ch]
    fileinp.close()
    return results

def out_read(filename):
    fileinp=input(filename)
    tmp=[]
    for line in fileinp:
        tmp.append(line)
        try:
            splitted_line=line.split() 
            if (splitted_line[0].strip() == "COBRAM") and  (splitted_line[1].strip() == "VDW") :
                RESname=splitted_line[8].strip()
                RESnumber=int(splitted_line[7].strip())
            if (splitted_line[0].strip() == "VDW") and  (splitted_line[1].strip() == "(N-B") and (splitted_line[2].strip() == "+") \
                    and (splitted_line[3].strip() == "1-4)") and (splitted_line[4].strip() == "INTERACTION")  \
                    and (splitted_line[5].strip() == "ENERGY") and (splitted_line[6].strip() == "MATRIX"):
                init=int(fileinp.lineno())
        except:
            pass
    fileinp.close()
    
    temp_value=[]
    for j in range(init,init+30):
        splitted_tmp=tmp[j].split() 
        try:
            if splitted_tmp[0].strip()=='2':
                temp_value.append(splitted_tmp[1])
        except:
            pass
    value=float(temp_value[0])
    results=[RESname,RESnumber,value]
    return results
    
def script_read(scriptname):
    fileinp=input(scriptname)
    tmp=[]
    for line in fileinp:
        tmp.append(line)
        try:
            splitted_line=line.split() 
            if (splitted_line[0].strip() == "##") and  (splitted_line[1].strip() == "number_of_residues:"):
                Rn=(splitted_line[2])
        except:
            pass
    fileinp.close()
    return int(Rn)

print ""
print "======================================="
print "COBRAM analysis-results starts "
print ""



dir_elec_DFP='analys_electro_DFP'
dir_elec_RFP='analys_electro_RFP'
dir_vdw='analys_vdw'

accepted=['0','1','2','3']
controller=0
while controller==0:
    print '//////////////////////////////////////////////////////////////////'
    print '// [0] Analysis of all contributes                              //'
    print '// [1] Analysis of electrostatic Direct   Finger Print  (DFP)   //'
    print '// [2] Analysis of electrostatic Rreverse Finger Print  (RFP)   //'
    print '// [3] Analysis of van der Waals          Finger Print   (FP)   //'
    print '//////////////////////////////////////////////////////////////////'
    choise=raw_input('What do you want to do? [0]')
    if choise=='':
        choise='0'
        controller=1
    else:
        if choise in accepted:
            controller=1
        else:
            controller=0
            
            
    Rn_DFP=int(script_read(dir_elec_DFP+'/electro_DFP_script4run.sh'))
    Rn_RFP=int(script_read(dir_elec_RFP+'/electro_RFP_script4run.sh'))    
    Rn_VDW=int(script_read(dir_vdw+'/vdw_script4run.sh'))
    if Rn_DFP==Rn_RFP and Rn_DFP==Rn_VDW:
        tot_Rn=Rn_DFP
        print "The whole system is composed by "+str(tot_Rn)+" residues"
    else:
        print "DFP residue number is ",Rn_DFP
        print "RFP residue number is ",Rn_RFP
        print "VDW residue number is ",Rn_VDW
        print "They must have the same value:"
        print " calculation aborted"
        exit()
            
    dir=listdir('.')
    masks=[]
    for i in range(len(dir)):
        if fnmatch.fnmatch(dir[i],'masked_*.dat')==1:
            masks.append(dir[i])

    print "You can restrict the analysis over "+str(len(masks))+ " masks"
    prnt=  '%s%3s%s  %s'  % ('(', 0, ')', 'All the system')
    print prnt
    if len(masks)!=0:
        maskdict={}
        for i in range(len(masks)):
            maskdict[i+1]=masks[i]
        
        for i in range(len(masks)):
            prnt=  '%s%3s%s  %s'  % ('(', i+1, ')', maskdict[i+1])      
            print prnt
        print ""
    
    mskchoisecontrol=0
    while mskchoisecontrol==0:
        ask=raw_input('What is your choice?  [0]')
        if ask=='':
            ask=0
            mskchoisecontrol=1
        else:
            if len(masks)==0:
                print 'No mask is avalialble: the whole system will be considered'
                ask=0
                mskchoisecontrol=1
            else:
                try:
                    if maskdict.has_key(int(ask))==1:
                        print 'The requested mask will be processed'
                        ask=int(ask)
                        mskchoisecontrol=1
                    else:
                        print 'The requested mask is not available'
                        mskchoisecontrol=0
                except:
                    print 'The requested mask is not available'
                    mskchoisecontrol=0
                    pass
            
if ask==0:
    list_res=range(tot_Rn+1)
    Rn=tot_Rn
else:
    maskinp=input(maskdict[ask])
    list_res=[0]
    for line in maskinp:
        list_res.append(int(line.strip()))
    maskinp.close()
    Rn=len(list_res)

system("mkdir re-run")##to repeat aborted calculations, if present

replace_script=open('re-run/replace-script.sh','w')
replace_script.write('rm *.rwf *.int *.d2e *.chk \n')
replace_script.write('cp *DFP* ../analys_electro_DFP \n')
replace_script.write('cp *RFP* ../analys_electro_RFP \n')
replace_script.close()
rerun_script=open('re-run/rerun-script.sh','w')
rerun_script.write('##lauch me please\n')
err1=0
err2=0

if choise=='0' or choise=='1':
## analysis of electrostatic DFP
    DFP_fileout=open('electrostatic_DFP.dat','w')

    DFP_RESname=[]
    DFP_RESnumber=[]
    DFP_SCF_ene=[]
    DFP_SELF_En_of_Ch=[]
    DFP_ENE_cleaned=[]
    DFP_ENE_rel_h=[]
    DFP_ENE_rel_kcal=[]
    j=0
    for i in range(tot_Rn):
        if i in list_res:
            filename1=dir_elec_DFP+'/electrostatic_DFP_'+str(i)+'.log'
            ck=log_check(filename1)
            if ck==1:
                log_results=log_read(filename1)
                DFP_RESname.append(log_results[0])
                DFP_RESnumber.append(log_results[1])
                DFP_SCF_ene.append(log_results[2])
                DFP_SELF_En_of_Ch.append(log_results[3])
                etmp1=DFP_SCF_ene[j]-DFP_SELF_En_of_Ch[j]
                DFP_ENE_cleaned.append(etmp1)
                etmp2=etmp1-DFP_ENE_cleaned[0]
                DFP_ENE_rel_h.append(etmp2)
                etmp3=etmp2*627.5
                DFP_ENE_rel_kcal.append(etmp3)
                DFP_fileout.write('%4s' % DFP_RESnumber[j]+' '+'%3s' % DFP_RESname[j]+' '+'%12.7f' % DFP_SCF_ene[j]+' ' +'%12.7f' % DFP_SELF_En_of_Ch[j] +' '+'%12.7f' % DFP_ENE_cleaned[j]+' '+'%12.7f' % DFP_ENE_rel_h[j]+' '+'%12.7f' % DFP_ENE_rel_kcal[j]   +'\n' )
            else:
                DFP_RESname.append('****')
                DFP_RESnumber.append('****')
                DFP_SCF_ene.append('****')
                DFP_SELF_En_of_Ch.append('****')
                DFP_ENE_cleaned.append('****')
                DFP_ENE_rel_h.append('****')
                DFP_ENE_rel_kcal.append('****')
                DFP_fileout.write('**************************************************************************\n' )
                print 'Please check the file electrostatic_DFP_'+str(i)+'.log'
                system('cp analys_electro_DFP/electrostatic_DFP_'+str(i)+'.com re-run')
                rerun_script.write('g03<electrostatic_DFP_'+str(i)+'.com>electrostatic_DFP_'+str(i)+'.log\n')
                err1=err1+1
            j=j+1
    DFP_fileout.close() 


if choise=='0' or choise=='2':
    ## analysis of electrostatic RFP
    RFP_fileout=open('electrostatic_RFP.dat','w')

    RFP_RESname=[]
    RFP_RESnumber=[]
    RFP_SCF_ene=[]
    RFP_SELF_En_of_Ch=[]
    RFP_ENE_cleaned=[]
    RFP_ENE_rel_h=[]
    RFP_ENE_rel_kcal=[]
    j=0
    for i in range(tot_Rn):
        if i in list_res:
            filename1=dir_elec_RFP+'/electrostatic_RFP_'+str(i)+'.log'
            ck=log_check(filename1)
            if ck==1:
                log_results=log_read(filename1)
                RFP_RESname.append(log_results[0])
                RFP_RESnumber.append(log_results[1])
                RFP_SCF_ene.append(log_results[2])
                RFP_SELF_En_of_Ch.append(log_results[3])
                etmp1=RFP_SCF_ene[j]-RFP_SELF_En_of_Ch[j]
                RFP_ENE_cleaned.append(etmp1)
                etmp2=etmp1-RFP_ENE_cleaned[0]
                RFP_ENE_rel_h.append(etmp2)
                etmp3=etmp2*627.5
                RFP_ENE_rel_kcal.append(etmp3)
                RFP_fileout.write('%4s' % RFP_RESnumber[j]+' '+'%3s' % RFP_RESname[j]+' '+'%12.7f' % RFP_SCF_ene[j]+' ' +'%12.7f' % RFP_SELF_En_of_Ch[j] +' '+'%12.7f' % RFP_ENE_cleaned[j]+' '+'%12.7f' % RFP_ENE_rel_h[j]+' '+'%12.7f' % RFP_ENE_rel_kcal[j]   +'\n' )
            else:
                RFP_RESname.append('****')
                RFP_RESnumber.append('****')
                RFP_SCF_ene.append('****')
                RFP_SELF_En_of_Ch.append('****')
                RFP_ENE_cleaned.append('****')
                RFP_ENE_rel_h.append('****')
                RFP_ENE_rel_kcal.append('****')
                RFP_fileout.write('**************************************************************************\n' )
                print 'Please check the file electrostatic_RFP_'+str(i)+'.log'
                system('cp analys_electro_RFP/electrostatic_RFP_'+str(i)+'.com re-run')
                rerun_script.write('g03<electrostatic_RFP_'+str(i)+'.com>electrostatic_RFP_'+str(i)+'.log\n')
                err2=err2+1
            j=j+1
    RFP_fileout.close()    
        
rerun_script.close()
if err1+err2==0:
    system('rm -rf re-run' )
else:
    system('chmod a+rwx re-run/*.sh' )

if choise=='0' or choise=='3':
    ## analysis of electrostatic VDW
    VDW_fileout=open('VDW_FP.dat','w')
    
    VDW_RESname=['ref']
    VDW_RESnumber=['0']
    VDW_SCF_ene=[0.0]
    VDW_SELF_En_of_Ch=[0.0]
    VDW_ENE_cleaned=[0.0]
    VDW_ENE_rel_h=[0.0]
    VDW_ENE_rel_kcal=[0.0]
    i=0
##    VDW_fileout.write('%4s' % VDW_RESnumber[i]+' '+'%3s' % VDW_RESname[i]+' '\
##        +'%12.7f' % VDW_SCF_ene[i]+' ' +'%12.7f' % VDW_SELF_En_of_Ch[i] +' '\
##        +'%12.7f' % VDW_ENE_cleaned[i]+' '+'%12.7f' % VDW_ENE_rel_h[i]+' '\
##        +'%12.7f' % VDW_ENE_rel_kcal[i]   +'\n' )
    j=0
    for i in range(1,tot_Rn):
        if i in list_res:
            filename1=dir_vdw+'/vdw_'+str(i)+'.out'
            out_results=out_read(filename1)
            VDW_RESname.append(out_results[0])
            VDW_RESnumber.append(out_results[1])
            VDW_SCF_ene.append(0.0)
            VDW_SELF_En_of_Ch.append(0.0)
            etmp1=0.0
            VDW_ENE_cleaned.append(etmp1)
            etmp2=0.0
            VDW_ENE_rel_h.append(etmp2)
            etmp3=0.0
            VDW_ENE_rel_kcal.append(out_results[2])
        
            VDW_fileout.write('%4s' % VDW_RESnumber[j]+' '+'%3s' % VDW_RESname[j]+' '\
                +'%12.7f' % VDW_SCF_ene[j]+' ' +'%12.7f' % VDW_SELF_En_of_Ch[j] +' '\
                +'%12.7f' % VDW_ENE_cleaned[j]+' '+'%12.7f' % VDW_ENE_rel_h[j]+' '\
                +'%12.7f' % VDW_ENE_rel_kcal[j]   +'\n' )
            j=j+1
    VDW_fileout.close() 

print ""
print "COBRAM analysis-results ended "
print "======================================="
