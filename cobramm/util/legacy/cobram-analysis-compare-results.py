#! /usr/bin/env python
from fileinput import *
from math import *
from Numeric import * 
from time import *
from os import system
from os import listdir
from sys import *

def read_file(nomefile):
    fileinp=input(nomefile)
    tmp=[]
    row0=[]
    row1=[]
    row2=[]
    row3=[]
    row4=[]
    row5=[]
    row6=[]
    for line in fileinp:
        tmp.append(line)
        try:
            splitted_line=line.split() 
            row0.append(int(splitted_line[0]))
            row1.append(splitted_line[1])
            row2.append(float(splitted_line[2]))
            row3.append(float(splitted_line[3]))
            row4.append(float(splitted_line[4]))
            row5.append(float(splitted_line[5]))
            row6.append(float(splitted_line[6]))
        except:
            pass
    datas=[row0,row1,row2,row3,row4,row4,row5,row6]
    return datas
    
def comparator(datas_A,datas_B,column,shifter):
    new_row0=[]
    new_row=[]
    shifted=[]
    for i in range(len(datas_A[0])):
        delta=datas_B[column][i]-datas_A[column][i]
        new_row.append(delta)
        shifted.append(datas_A[0][i]+shifter)
        new_row0.append(shifted)
    shifted[0]=0
    results=[shifted,datas_A[1],new_row]
    return results
    
def writer(nomefile_C,results):
    fileout=open(nomefile_C,'w')
    for i in range(len(results[0])):
        fileout.write('%4s' % results[0][i]+' '+'%3s' % results[1][i]+' '+'%12.7f' % results[2][i]+'\n')
    fileout.close()

print ""
print "======================================="
print "COBRAM analysis-compare-results starts "
print ""
print "You will be asked to type the filename of two analysis"
print "Comparison will be of type B-A"
print ""
nomefile_A=raw_input('Please enter the name of file A ')
nomefile_B=raw_input('Please enter the name of file B ')
nomefile_C=raw_input('Please enter the name of output file ')
offset       =raw_input('Please enter offset value [0] ')
if offset=='':
    offset='0'
offset=int(offset)
    
datas_A=read_file(nomefile_A)
datas_B=read_file(nomefile_B)   

column=7
results=comparator(datas_A,datas_B,column,offset)
writer(nomefile_C,results)
    
    
print ""
print "======================================="
print "COBRAM analysis-compare-results ended "
print ""
    
