#!/usr/bin/env python3
# coding=utf-8

#    COBRAMM
#    Copyright (c) 2019 ALMA MATER STUDIORUM - Università di Bologna

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import os
import CBF
import math
import shelve
import logwrt
import constants


def product(vec,mat):
    pro=[]

    for i in range(len(vec)):
        tmp=0
        for j in range(len(vec)):
           tmp=tmp+vec[j]*mat[i][j]
        pro.append(tmp)

    res=np.array(pro)
    return res

#def decoherence(At,tstep,DEarray,state):
#    ekin=float(os.popen('cat EKIN').read())
#    ASum=0
#    if ekin == 0.0 :
#       ekin=0.8
#    for i in range(len(At)):
#       if state != i:
#          elev=-tstep/((1.0/(abs(DEarray[state][i]/627.51)))*(1.0+0.1/ekin))
#          At[i]=At[i]*math.exp(elev)
#          ASum += abs(At[i])**2
#    At[state] = At[state]*math.sqrt((1.0-ASum)/(abs(At[state])**2))
#    sef=shelve.open("cobram-sef")
#    sef['AM']=At
#    SItime=sef['SItime']
#    sef.close()
#    out=open('Amplitudes.dat','w')
#    for i in range(len(At)):
#       out.write('%16.12f' %At[i].real+' %16.12f ' %+At[i].imag)
#    out.close()
#    #append amplitudes from this step to list of amplitdues
#    out=open('AmplitudesALL.dat','a')
#    out.write('%10.4f' %SItime+'   '),
#    for i in range(len(At)):
#        out.write('%16.12f' %At[i].real+' %16.12f ' %At[i].imag+'   '),
#    out.write('\n')
#    out.close()
#    return  At 

def angle(x1,y1,z1,x2,y2,z2):
    Nfact1=1/math.sqrt(np.add.reduce(x1**2+y1**2+z1**2))
    Nfact2=1/math.sqrt(np.add.reduce(x2**2+y2**2+z2**2))
    ## Normalize velocities and der coupling
    Nx2=Nfact2*x2
    Ny2=Nfact2*y2
    Nz2=Nfact2*z2
    Nx1=Nfact1*x1
    Ny1=Nfact1*y1
    Nz1=Nfact1*z1
    N2=math.sqrt(np.add.reduce(Nx2**2+Ny2**2+Nz2**2))
    N1=math.sqrt(np.add.reduce(Nx1**2+Ny1**2+Nz1**2))
    try:
       angle=math.acos(np.add.reduce(Nx2*Nx1+Ny2*Ny1+Nz2*Nz1))*57.2958
    except:
       print( "possible round-off error, angles are either exactly orthogonal or parallel, trying to fix this:")
       toround=round(np.add.reduce(Nx2*Nx1+Ny2*Ny1+Nz2*Nz1),7)
       angle=math.acos(toround)*57.2958
    res=[angle,Nfact1,Nfact2]
    return res


#def velscale(geometry,command,xdc,ydc,zdc,vx,vy,vz,DE,tstep,cobcom):
#    amu=1822.8880
#    if command[14] == '1':
#        at = [geometry.atomLabel[i - 1] for i in geometry.list_MEDIUM_HIGH]
#        model_MEDIUM_HIGH = geometry.MEDIUM_HIGH
#        x = np.array(model_MEDIUM_HIGH[0])
#        y = np.array(model_MEDIUM_HIGH[1])
#        z = np.array(model_MEDIUM_HIGH[2])
#    else:
#        modelH = geometry.modelH
#        labelModelH = [geometry.atomLabel[i - 1] for i in geometry.list_QM]
#        for i in range(geometry.NsubH):
#            labelModelH.append('H')
#        at = labelModelH
#        x = np.array(modelH[0])
#        y = np.array(modelH[1])
#        z = np.array(modelH[2])
#    if int(command[2]) > 1:
#        print( "Length of vectors", len(x))
#
#    masses = []
#    for i in range(len(at)):
#        atommass = constants.atommass(at[i])
#        if atommass is None:
#            logwrt.fatalerror('no mass for ' + at[i] + ', please modify the code')
#        else:
#            masses.append(atommass)
#
#    ISO=CBF.ReadCobramCommand(cobcom,'isotopes','isotopes')
#    if (len(ISO))==0 and int(command[2]) > 0:
#        print( "NO ISOTOPE section, continuing")
#    else:
#        try:
#          oldmass=[]
#          for i in range(len(ISO)):
#             isotopes=ISO[i].split()
#             atnum=int(isotopes[0])
#             print( "modifying the mass of atom number:", atnum)
#             print( "current mass   :",float(masses[atnum]/amu))
#             oldmass.append(float(masses[atnum]/amu))
#             print( "new mass       :",float(isotopes[1]))
#             masses[atnum]=float(float(isotopes[1])*amu)
#             print( "in atomic units:",masses[atnum])
#        except:
#          logwrt.fatalerror('severe PROBLEM with isotope section (empty line ??)')
#
#    ang=angle(xdc,ydc,zdc,vx,vy,vz)[0]
#    print( "Angle between DC(GD) and velocity:", ang)
#    if ang > 90 :
#        print( "DC(GD) will be inverted!")
#        xdc=-xdc
#        ydc=-ydc
#        zdc=-zdc
#    b=0.0
#    a=0.0
#    c=-abs(DE/627.51)  ## use negative value to check!
#
#    for i in range(len(x)):
#        dcvec=[xdc[i],ydc[i],zdc[i]]
#        vvec=[vx[i],vy[i],vz[i]]
#        b=b+np.add.reduce(np.array(dcvec)*np.array(vvec))
#        a=a+np.add.reduce((np.array(dcvec)**2)/masses[i])
#    print( "Checking if HOP is valid")
#    print( "b = vel x DC(GD)", b)
#    print( "a = DC(GD)**2/mass", a)
#    print( "c = Energy diff (negative)", c)
#    print( "Following value must be positive to permit hopping sqrt(b^2 + 2*c*a)", b**2+2*c*a)
#    try:
#        scale=(b+math.sqrt(b**2+2*c*a))/(0.5*a*tstep)
#        if int(command[2]) > 1:
#            print( 'The positive scaling factor is:  ', scale*627.51)
#        condition='true'
#    except:
#        logwrt.writelog('-----------------------------')
#        logwrt.writelog('frustrated hop!!!')
#        logwrt.writelog('-----------------------------')
#        condition='false'
#    return condition               

def velscalev(masses,geometry,command,E_pot,E_kin,etot,tstep,tstep_OLD,vx,vy,vz,cobcom):
    vx0,vy0,vz0=[],[],[]
    for i in range(len(vx)):
        vx0.append(vx[i])
        vy0.append(vy[i])
        vz0.append(vz[i])
    vx0,vy0,vz0=np.array(vx0),np.array(vy0),np.array(vz0)
    c=-(E_pot+E_kin-etot)
    ekin=float(os.popen('cat EKIN').read())
    if int(command[2]) > 0:
        print( "###########\nEnergetics after hopping\n###########")
        print( "Kinetic energy at this step:",E_kin)
        print( "Potential energy at this step:",E_pot)
        print( "Total energy at this step:",E_kin+E_pot)
        print( "Kinetic energy at previous step:",ekin)
        print( "Potential energy at previous step:",etot-ekin)
        print( "Total energy at previous step:",etot)
        print( 'Total energy (this step) - Total energy (previous step):',c)
    if c < 0:
        print( "Total energy increases after hop!")

    amu=1822.8880
    # F: old scheme = in case of dynamics with TDC, rescale velocity of QM+MM parts (model MH), in case of dynamics with spatil NAC, rescale only modelH components
    #if command[14] == '1':
    #    at = [geometry.atomLabel[i - 1] for i in geometry.list_MEDIUM_HIGH]
    #    model_MEDIUM_HIGH = geometry.MEDIUM_HIGH
    #    x = np.array(model_MEDIUM_HIGH[0])
    #    y = np.array(model_MEDIUM_HIGH[1])
    #    z = np.array(model_MEDIUM_HIGH[2])
    #else:
    #    modelH = geometry.modelH
    #    labelModelH = [geometry.atomLabel[i - 1] for i in geometry.list_QM]
    #    for i in range(geometry.NsubH):
    #        labelModelH.append('H')
    #    at = labelModelH
    #    x = np.array(modelH[0])
    #    y = np.array(modelH[1])
    #    z = np.array(modelH[2])
    # F: new scheme = always rescale only modelH components (as in previous "else")
    modelH = geometry.modelH
    labelModelH = [geometry.atomLabel[i - 1] for i in geometry.list_QM]
    for i in range(geometry.NsubH):
        labelModelH.append('H')
    at = labelModelH
    x = np.array(modelH[0])
    y = np.array(modelH[1])
    z = np.array(modelH[2])
    if int(command[2]) > 1:
        print( "Length of vectors", len(x))

    masses = []
    for i in range(len(at)):
        atommass = constants.atommass(at[i])
        if atommass is None:
            logwrt.fatalerror('no mass for ' + at[i] + ', please modify the code')
        else:
            masses.append(atommass)

    ISO=CBF.ReadCobramCommand(cobcom,'isotopes','isotopes')
    if (len(ISO))==0 and int(command[2]) > 0:
        print( "NO ISOTOPE section, continuing")
    else:
        try:
          oldmass=[]
          for i in range(len(ISO)):
             isotopes=ISO[i].split()
             atnum=int(isotopes[0])
             print( "modifying the mass of atom number:", atnum)
             print( "current mass   :",float(masses[atnum]/amu))
             oldmass.append(float(masses[atnum]/amu))
             print( "new mass       :",float(isotopes[1]))
             masses[atnum]=float(float(isotopes[1])*amu)
             print( "in atomic units:",masses[atnum])
        except:
          logwrt.fatalerror('severe PROBLEM with isotope section (empty line ??)')

    print( "Trying to rescale the velocity vector!")
    try:
        #c = Ei_pot(t-dt) + E_kin(t-dt) - Ej_pot(t) - E_kin(t) = E_total(t-dt) - E_total(t)
        #scale = sqrt( (E_total(t-dt) - Ej_pot(t))/E_kin(t) ) = sqrt( (E_total(t-dt) - Ej_pot(t) - E_kin(t))/E_kin(t) + E_kin(t)/E_kin(t) ) = sqrt( (E_total(t-dt) - E_total(t))/E_kin(t) + 1 ) = sqrt(c//E_kin(t) + 1)
        scale=math.sqrt(c/E_kin + 1.0)
        print( "Scale", scale)
        for i in range(len(at)):
            vx[i]*=scale
            vy[i]*=scale
            vz[i]*=scale
        ang=angle(vx0,vy0,vz0,vx,vy,vz)[0]
        if int(command[2]) > 0:
            print( "Scaled kinetic energy at this step:",E_kin*scale*scale)
            print( "Potential energy at this step:",E_pot)
            print( "Total energy at previous step:",etot)
            print( 'Total energy with rescaled kinetic contribution (this step) - Total energy (previous step):',E_kin*scale*scale + E_pot - etot)
            print( "Angle between original and scaled velocities", round(ang,3))
        veloc=[vx,vy,vz]
    except:
        logwrt.writelog('-----------------------------')
        logwrt.writelog('Frustrated hop!!!')
        logwrt.writelog('-----------------------------')
        veloc='false'
    return veloc

def velrescale(masses,geometry,command,xdc,ydc,zdc,E_pot,E_kin,etot,tstep,tstep_OLD,vx,vy,vz,cobcom):
    vx0,vy0,vz0=[],[],[]
    vx1,vy1,vz1=[],[],[]
    for i in range(len(vx)):
        vx0.append(vx[i])
        vy0.append(vy[i])
        vz0.append(vz[i]) 
        vx1.append(vx[i])
        vy1.append(vy[i])
        vz1.append(vz[i])
    vx0,vy0,vz0=np.array(vx0),np.array(vy0),np.array(vz0)
    vx1,vy1,vz1=np.array(vx1),np.array(vy1),np.array(vz1)
    c=-(E_pot+E_kin-etot)  
    ekin=float(os.popen('cat EKIN').read())
  
    if int(command[2]) > 0:
        print( "\n###########\nEnergetics after hopping\n###########" )
        print( "Kinetic energy at this step:",E_kin)
        print( "Potential energy at this step:",E_pot)
        print( "Total energy at this step:",E_kin+E_pot)
        print( "Kinetic energy at previous step:",ekin)
        print( "Potential energy at previous step:",etot-ekin)
        print( "Total energy at previous step:",etot)
        print( 'Total energy (this step) - Total energy (previous step):',c)
    if c < 0:
        print( "Total energy increases after hop!")

    ####################################

    if command[206] == '2':
        print( "\n###########\nScaling velocity along GD vector\n###########")
    else:
        print( "\n###########\nScaling velocity along NAC vector\n###########")

    b=0.0
    a=0.0
    ccc=0
    for i in range(geometry.NatomHM):
        #DC GD are available only for HIGH layer
        if geometry.list_MEDIUM_HIGH[i] in geometry.list_HIGH:
            dcvec=[xdc[ccc],ydc[ccc],zdc[ccc]]
            vvec=[vx0[i],vy0[i],vz0[i]]
            b=b+np.add.reduce(np.array(dcvec)*np.array(vvec))
            a=a+np.add.reduce((np.array(dcvec)**2)/masses[i])
            ccc=ccc+1

    try:
        scale=-(b+math.sqrt(b**2+2*c*a))/(0.5*a*tstep_OLD)
    except:
        print( "HOP rejected!")
        veloc='false'
        return veloc

    if int(command[2]) > 1:
        print( 'The scaling is: ',scale)

    E_kin=0.0
    ccc=0
    for i in range(geometry.NatomHM):
        # DC and GD are available only for HIGH layer
        if geometry.list_MEDIUM_HIGH[i] in geometry.list_HIGH:
            vscx=(scale*0.5*xdc[ccc]*tstep/masses[i])
            vscy=(scale*0.5*ydc[ccc]*tstep/masses[i])
            vscz=(scale*0.5*zdc[ccc]*tstep/masses[i])
            vx[i]=vx[i]+vscx
            vy[i]=vy[i]+vscy
            vz[i]=vz[i]+vscz
            ccc=ccc+1
        E_kin=E_kin+0.5*(vx[i]**2+vy[i]**2+vz[i]**2)*masses[i]

    logwrt.writelog('Error in the energy redistribution = %16.10f\n' % abs(E_pot+E_kin-etot))
    diff=-(E_pot+E_kin-etot)
    ang1=angle(vx0,vy0,vz0,vx,vy,vz)
    if int(command[2]) > 0:
        print("Scaled kinetic energy at this step:", E_kin)
        print("Potential energy at this step:", E_pot)
        print("Total energy at previous step:", etot)
        print('Total energy with rescaled kinetic contribution (this step) - Total energy (previous step):', diff)
        print("Angle between original and scaled velocities", ang1[0])
    

    #################################

    if command[206] == '2':
        print( "\n###########\nScaling velocity along -GD vector\n###########")
    else:
        print( "\n###########\nScaling velocity along -NAC vector\n###########")
    xdc,ydc,zdc=-xdc,-ydc,-zdc
    b=0.0
    a=0.0
    ccc=0
    for i in range(geometry.NatomHM):
        # DC and GD are available only for HIGH layer
        if geometry.list_MEDIUM_HIGH[i] in geometry.list_HIGH:
            dcvec=[xdc[ccc],ydc[ccc],zdc[ccc]]
            vvec=[vx0[i],vy0[i],vz0[i]]
            b=b+np.add.reduce(np.array(dcvec)*np.array(vvec))
            a=a+np.add.reduce((np.array(dcvec)**2)/masses[i])
            ccc=ccc+1

    scale=-(b+math.sqrt(b**2+2*c*a))/(0.5*a*tstep_OLD)
    if int(command[2]) > 1:
        print( 'The scaling is: ',scale)

    E_kin=0.0
    ccc=0
    for i in range(geometry.NatomHM):
        # DC and GD are available only for HIGH layer
        if geometry.list_MEDIUM_HIGH[i] in geometry.list_HIGH:
            vscx=(scale*0.5*xdc[ccc]*tstep/masses[i])
            vscy=(scale*0.5*ydc[ccc]*tstep/masses[i])
            vscz=(scale*0.5*zdc[ccc]*tstep/masses[i])
            vx1[i]=vx1[i]+vscx
            vy1[i]=vy1[i]+vscy
            vz1[i]=vz1[i]+vscz
            ccc=ccc+1
        E_kin=E_kin+0.5*(vx1[i]**2+vy1[i]**2+vz1[i]**2)*masses[i]

    logwrt.writelog('Error in the energy redistribution = %16.10f\n' % abs(E_pot+E_kin-etot))
    ang2=angle(vx0,vy0,vz0,vx1,vy1,vz1)
    diff=-(E_pot+E_kin-etot)
    if int(command[2]) > 0:
        print( "Scaled kinetic energy at this step:",E_kin)
        print( "Potential energy at this step:",E_pot)
        print( "Total energy at previous step:",etot)
        print( 'Total energy with rescaled kinetic contribution (this step) - Total energy (previous step):',diff)
        print( "Angle between original and scaled velocities", ang2[0])

    if ang2 < ang1 :
        if command[206] == '2':
            print( "\nDecision taken: Scaling velocity along -GD vector")
        else:
            print( "\nDecision taken: Scaling velocity along -NAC vector")
        vx,vy,vz=vx1,vy1,vz1
    else:
        if command[206] == '2':
            print( "\nDecision taken: Scaling velocity along GD vector")
        else:
            print( "\nDecision taken: Scaling velocity along NAC vector")
    veloc=[vx,vy,vz]   
    return veloc
