#! /usr/bin/env python
# -*- coding: utf-8 -*-
##########################################################################
#                                                                        #
#  mod_real.py                                                           #
#                                                                        #
#  MODIFY real.crd  --> insert XYZ                                       #
#                                                                        #
#  this program is part of the COBRAM-beta_5.7 suite                     #
#                                                                        #
#  (C) O.Weingart February 2015                                          #
#  HHU Duesseldorf                                                       #  
#                                                                        #
##########################################################################
#
#
# insert given XYZ coordinates into the real.crd
# needs real.crd and realmod.in file specifying which range of 
# (consecutive) atoms to replace
# each line holds a range of atoms
# e.g. 
# 1207 1219
# 1221 1240
# etc.

try:
   filein=open('real.crd')
except:
   print "ERROR, no real.crd file!"
   stop
text=filein.read()
mod=text.split('\n')
atnum=int(mod[1])
x=[]
y=[]
z=[]
for i in range(2,len(mod)):
   tmp=mod[i].split()
   if tmp!=[]:
      x.append(tmp[0])
      y.append(tmp[1]) 
      z.append(tmp[2]) 
      try: 
        x.append(tmp[3])
        y.append(tmp[4]) 
        z.append(tmp[5]) 
      except:
        pass
try:
   filein=open('realmod.in')
except:
   print "ERROR, no realmod.in file!"
   stop
text=filein.read()
inp=text.split('\n')
print ""
print "Read",len(inp), "lines of input"
begin=[]
ending=[]
line=0
for i in range(len(inp)):
   tmp=inp[i].split()
   print "tmp=", tmp
   if tmp!=[]:
     begin.append(tmp[0])
     line=line+1
     try:
       ending.append(tmp[1])
     except:
       print "single atom in this line!"
       pass
print ""
print "Will modify real.crd atom numbers:"
print ""
total=0
for i in range(line):
      print begin[i]
      try:
        print "to"
        print ending[i]
        total=total+int(ending[i])-int(begin[i])+1
        print ""
      except:
        print "Exception"
        pass

xyzname=raw_input('XYZ-input-file to insert? ')
print xyzname   
filein=open(xyzname)
text=filein.read()
xyzdata=text.split('\n')
print xyzdata[0]
atxyz=int(xyzdata[0])
print "We have ",atxyz," atoms in the XYZ file"
print "and ", total , " Atoms to modify"

# now check if the user did not do any stupid input

modx=[]
mody=[]
modz=[]

if (int(total)!=atxyz):
   print "ERROR - The number of atoms in XYZ and the number of atoms to modify does not match!"
if (atxyz != len(xyzdata)-3):
   print "ERROR - Wrong number of atoms in XYZ"
   print len(xyzdata)
for i in range(2,len(xyzdata)):
   try:
      tmp=xyzdata[i].split()
      modx.append(tmp[1])
      mody.append(tmp[2])
      modz.append(tmp[3])
   except:
      pass

# now replace:
count=0
for i in range(len(begin)):
   for j in range(int(begin[i]),int(ending[i])+1):
       print "Replacing number ",j, " with ", count+1
       x[j-1]=modx[count]          
       y[j-1]=mody[count]          
       z[j-1]=modz[count]                        
       count=count+1
print ""
print "Writing new .crd file"

coord_file=open(xyzname+'.crd','w')
coord_file.write('\n')
coord_file.write(str(atnum)+'\n')
i=0
while i<atnum:
   tmp=('%12.7f' % float(x[i])+'%12.7f' % float(y[i])+'%12.7f' % float(z[i]))
   coord_file.write(tmp)
   i=i+1
   try:
       tmp=('%12.7f' % float(x[i])+'%12.7f' % float(y[i])+'%12.7f' % float(z[i]))
       coord_file.write(tmp+'\n')
       i=i+1
   except:
       pass
print ""
print "********  All done! ***********"
print ""
print "Please carefully check with "
print ""
print 'ambpdb -p real.top <'+xyzname+'.crd>'+xyzname+'.pdb'
print ""
print "if the conversion worked properly!"
print ""
coord_file.close()

                                            
