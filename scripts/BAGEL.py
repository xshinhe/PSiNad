import argparse
import datetime
import math
import os
import re
import shutil
import subprocess
import sys
import time
import toml
from copy import deepcopy
from multiprocessing import Pool
from pprint import pprint
from traceback import format_exc
from socket import gethostname
import numpy as np

parser = argparse.ArgumentParser(description='Execute MNDO Calculation')
parser.add_argument('integers', metavar='N', type=int, nargs='+', 
    help='an integer for the accumulator')
parser.add_argument('--sum', dest='accumulate', action='store_const', const=sum, default=max,
    help='sum the integers (default: find the max)')

class PeriodicTable:
    """Hold a list of element names and a dictionary of atomic numbers."""

    def __init__(self):
        self.symbols = ["-", "H", "He",
                        "Li", "Be", "B", "C", "N", "O", "F", "Ne",
                        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
                        "K", "Ca",
                        "Sc", "Ti", "V", "Cr", "Mn",
                        "Fe", "Co", "Ni", "Cu", "Zn",
                        "Ga", "Ge", "As", "Se", "Br", "Kr",
                        "Rb", "Sr",
                        "Y", "Zr", "Nb", "Mo", "Tc",
                        "Ru", "Rh", "Pd", "Ag", "Cd",
                        "In", "Sn", "Sb", "Te", "I", "Xe",
                        "Cs", "Ba",
                        "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
                        "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
                        "Hf", "Ta", "W", "Re",
                        "Os", "Ir", "Pt", "Au", "Hg",
                        "Tl", "Pb", "Bi", "Po", "At", "Rn",
                        "Fr", "Ra",
                        "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
                        "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
                        "Rf", "Db", "Sg", "Bh",
                        "Hs", "Mt", "Ds", "Rg"]
        self.mass = ["-", 1.00794, 4.002602, 6.941, 9.0121831, 10.811, 12.0107, 14.0067, 15.9994, 18.998403163, 20.1797]
        self.atomicNumbers = {}
        for i in range(1, len(self.symbols)):
            self.atomicNumbers[self.symbols[i]] = i

    def getSymbol(self, atomicNumber):
        """Return an atomic symbol given an atomic number"""
        if atomicNumber < 1 or atomicNumber >= len(self.symbols):
            return None
        else:
            return self.symbols[atomicNumber]

    def getZ(self, symbol):
        """Return the atomic number corresponding to given atomic symbol"""
        if symbol in self.atomicNumbers:
            return self.atomicNumbers[symbol]
        else:
            raise KeyError('not valid atom number')

class RawFile():
    def __init__(self, fileName):
        self.filename = fileName
        try:
            f = open(fileName, 'r')
            self.fileContents = f.readlines()
            f.close()
        except IOError:
            sys.stderr.write("Error: could not read in %s\n" % fileName)
            sys.exit(1)  

    def searchForString(self, searchString, findFinal, startLine=-1):
        """Search for a string within the file.

        findFinal - if true, return the last occurence, else find the first
        startLine - start the search from this line (numbered from zero)

        Returns the list index of the matching line (i.e. the line number
        starting from zero),
        Returns -1 if the string wasn't found.
        """
        foundLine = -1
        lineNumber = -1
        for line in self.fileContents:
            lineNumber += 1
            if line.find(searchString) != -1 and lineNumber >= startLine:
                foundLine = lineNumber
                if not findFinal:
                    break
        return foundLine    

class OutputFile(RawFile):    
    """Hold an MNDO output file and analyse it."""
    
    def __init__(self, fileName):
        RawFile.__init__(self, fileName)
        
    def getPotential(self, num):
        '''
        Return energy of specific state in eV unit.
        '''
        frontline = self.searchForString(' SUMMARY OF MULTIPLE CI CALCULATIONS', False)      
        stateline = self.searchForString(f'    {num}    {num}', False, frontline)
        if stateline == -1:
            raise Exception( 'No such state', num)
        
        linelist = self.fileContents[stateline].split()
        energy = linelist[2]
        return float(energy)
    
    def getGrd(self, num):
        '''
        Return energy gradient of specific state in KCAL/(MOL*ANGSTROM) unit.
        '''
        grd = []
        
        stateline = self.searchForString(f'CI CALCULATION FOR STATE:  {num}', True)
        if stateline==-1:
            raise Exception( 'No such state', num)

        grdline = self.searchForString('GRADIENTS (KCAL/(MOL*ANGSTROM))', False, stateline)
        
        passline = 4
        for i in self.fileContents[grdline+passline: ]:
            splitline = i.split()
            
            if splitline == []:
                break
            grd.append(float(splitline[-3]))  
            grd.append(float(splitline[-2]))  
            grd.append(float(splitline[-1]))  
        return np.array(grd)
    
    def getCp(self, state1, state2):
        '''
        Return coupling terms between state1 and state2,
        E1 > E2, otherwise raise error
        '''
        if state1 < state2:
            state1, state2 = state2, state1
        coupling = []
        stateline = self.searchForString( 
                f' CI CALCULATION FOR INTERSTATE COUPLING OF STATES:  {state1}   {state2}' 
                ,True)
        
        if stateline==-1:
            raise Exception( 'No couple found between', state1, state2)
    
        cpline = self.searchForString('GRADIENTS (KCAL/(MOL*ANGSTROM))', False, stateline)
 
        passline = 4    
        for i in self.fileContents[cpline+passline: ]:
            splitline = i.split()
            
            if splitline == []:
                break
            coupling.append(float(splitline[-3]))  
            coupling.append(float(splitline[-2]))  
            coupling.append(float(splitline[-1]))  
        
        return np.array(coupling)
    
    def get_allPotential(self, F):
        E = []
        for i in range(1,F+1):
            E.append(self.getPotential(i))          
        return np.array(E)
    
    def get_allGrd(self, F):
        dE = []
        for i in range(1,F+1):
            dE.append(self.getGrd(i))          
        return np.array(dE)

class Loadinit_File(RawFile):
    '''
    New version!!!
    Loading initial positions and masses
    '''
    def __init__(self, filename, input_type='mndo'):
        super(Loadinit_File, self).__init__(filename)
        self.pt = PeriodicTable()
        if input_type == 'mndo':
            self._mndoinput()
        elif input_type == 'gjf':
            self._gjfinput()
        elif input_type == 'xyz':
            self._xyzinput()
            
    def _xyzinput(self):
        self._elements = []
        self._mass = []
        self._position = []

        for line in self.fileContents[2:]:
            i = line.split()
            try:
                Z = self.pt.getZ(i[0])
            except:
                Z=int(i[0])
            self._elements.append(Z)
            self._mass.append([self.pt.mass[Z]]*3)
            self._position.extend(i[1:4])
        
        
    def _mndoinput(self):
        self._elements = []
        self._mass = []
        self._position = []
        stateline = self.searchForString('+', True) 
        if stateline==-1:
            raise Exception( 'Unable to analyze input file')

        for line in self.fileContents[4+stateline:]:
            i = line.split()
            if i[0] == '0':
                break
            # Z = self.pt.getZ(i[0])
            # print(i[0])
            Z = int(i[0])
            # print(i)
            self._elements.append(Z)
            self._mass.append([self.pt.mass[Z]]*3)
            self._position.extend(i[1::2])
    
    def _gjfinput(self):
        self._elements = []
        self._mass = []
        self._position = []
        stateline = self.searchForString('Title', False)
        if stateline==-1:
            raise Exception( 'Unable to analyze input file')      
        passline = 3
        for line in self.fileContents[passline+stateline:]:
            i = line.split()
            if i == []:
                break
            Z = self.pt.getZ(i[0])
            self._elements.append(Z)
            self._mass.append([self.pt.mass[Z]]*3)
            self._position.extend(i[1:4])
        
            
    def get_mass(self):
        return np.array(self._mass).reshape(-1)
    
    def get_position(self):
        return np.array(self._position, dtype=np.float64)
    
    def get_elements(self):
        return self._elements

     

class setInputFile():
    '''
    Turn a gaussview gjf file into a mndo99 input file
    
    ****
    Old version for loading initial position, do not use it.
    Use Loadinit_File instead !
    ****
    
    '''
    def __init__(self, fileName, optimize=False, setting='default'):
        
        self.pt = PeriodicTable()
        self.fileName = fileName
        if optimize:
            self.opt = 1
        else:
            self.opt=0
        
        if setting=='default':
            with open(fileName, 'w') as f:
                f.write(
    'JOP=-2 IOP=-6 IGEOM=1 IFORM=1 ICUTS=-1 ICUTG=-1  +\n\
    ISCF=9 IPLSCF=9 DPREC=1D-8 DSTEP=1D-5 IPRINT=1 +\n\
    NCIGRD=2 IEF=1 IPREC=100 +\n\
    IMULT=1 IMOMAP=1 +\n\
    KCI=5 IOUTCI=1 MPRINT=1 ICROSS=7 +\n\
    MOVO=0 ICI1=10 ICI2=8 NCIREF=3 MCIREF=0 LEVEXC=2 IROOT=3 KITSCF=2000 \n\n'
                        )
                f.write('OM2\n')
        else:
            with open(fileName, 'w') as f:
                f.write(setting)
                f.write('OM2\n')

        
    def load_mndoinput(self, mndoinput):
        self._elements = []
        self._mass = []
        self._position = []
        
        file = RawFile(mndoinput)
        
        stateline = file.searchForString('+', True) 
        if stateline==-1:
            raise Exception( 'Unable to analyze input file')
        
        with open(self.fileName, 'a') as f:
            for line in file.fileContents[4+stateline:]:
                i = line.split()
                if i[0] == '0':
                    break
                # Z = self.pt.getZ(i[0])
                # print(i[0])
                Z = int(i[0])
                # print(i)
                self._elements.append(Z)
                self._mass.append([self.pt.mass[Z]]*3)
                self._position.extend(i[1::2])
                f.write(
    f'      {Z}  {i[1]}  {self.opt}  {i[2]}  {self.opt}  {i[3]}  {self.opt}  \n   '         
                        )
                
            f.write(
    '   0     0.0000000000 0     0.0000000000 0     0.0000000000 0\n'
                    )
            
            f.write('1 2 3')
            
    
    def loadgeom(self, gjffile):
        '''
        Tranform a gjf file
        '''
        self._elements = []
        self._mass = []
        self._position = []
        
        with open(gjffile, 'r') as f:
            contents = f.readlines()
        
        lineNumber = -1
        stateline = -1
        for line in contents:
            lineNumber += 1
            if line.find('Title') != -1:
                stateline = lineNumber
                break
        
        if stateline==-1:
            raise Exception( 'Unable to analyze input file')
        
        passline = 3
        #print(stateline)
        with open(self.fileName, 'a') as f:
            for line in contents[passline+stateline:]:
                i = line.split()
                if i == []:
                    break
                Z = self.pt.getZ(i[0])
                self._elements.append(Z)
                self._mass.append([self.pt.mass[Z]]*3)
                self._position.extend(i[1:4])
                f.write(
    f'      {Z}  {i[1]}  0  {i[2]}  0  {i[3]}  0  \n   '         
                        )
                
            f.write(
    '   0     0.0000000000 0     0.0000000000 0     0.0000000000 0\n'
                    )
            
            f.write('1 2 3')
    def load_xyz(self, xyzfile):
        self._elements = []
        self._mass = []
        self._position = []

        with open(xyzfile, 'r') as f:
            contents = f.readlines()
        
        with open(self.fileName, 'a') as f:
            for line in contents[2:]:
                i = line.split()
                try:
                    Z = self.pt.getZ(i[0])
                except:
                    Z=int(i[0])
                self._elements.append(Z)
                self._mass.append([Z]*3)
                self._position.extend(i[1:4])
                f.write(
        f'      {Z}  {i[1]}  {self.opt}  {i[2]}  {self.opt}  {i[3]}  {self.opt}  \n   '         
                            )
                    
            f.write(
    '   0     0.0000000000 0     0.0000000000 0     0.0000000000 0\n'
                    )
                
            f.write('1 2 3') 
        
    def get_mass(self):
        return np.array(self._mass).reshape(-1)
    
    def get_position(self):
        return np.array(self._position, dtype=np.float64)
    
    def get_elements(self):
        return self._elements
            
class My_InputFile(setInputFile):
    def __init__(self, filename):
        super(My_InputFile, self).__init__(filename)
    
    def writein(self, position, elements, couple_state):
        with open(self.fileName, 'a') as f:
            for i in range(len(position)//3):
                f.write(
    f'      {elements[i]}  {position[i*3]}  0  {position[i*3+1]}  0  {position[i*3+2]}  0  \n   '             
                        )
            f.write(
    '   0     0.0000000000 0     0.0000000000 0     0.0000000000 0\n'                    
                    )
            if couple_state == 3:
                f.write('1 2 3')
            elif couple_state == 2:
                f.write('1 2')
    
    def update(self, out_filename):
        p = os.system(f'mndo99 < {self.fileName} > {out_filename}')
        if p == 1:
            raise ValueError('mndo99 command unsuccessful, your computer may not install mndo99 program')
    
class Read_xyz(RawFile):
    def __init__(self, filename):
        super(Read_xyz, self).__init__(filename)
        self.pt = PeriodicTable()
        self.filename = filename
        self._process()
        
    def _process(self):
        n = int(self.fileContents[0].split()[0])
        self.n = n
        var = []
        time = []
        self._elements = []
        self.mass = []
        for i in range(0,len(self.fileContents),n+2):
            try: time.append(float(self.fileContents[i+1].split()[2]))
            except: pass
            for j in range(n):
                line = self.fileContents[i+j+2].split()
                var.extend(line[1:4])
                if i == 0:
                    try: Z=self.pt.getZ(line[0])
                    except: Z = int(line[0])
                    self._elements.append(Z)
                    self.mass += ([self.pt.mass[Z]]*3)
            # if i%(n+2)==0:
            #     p = []
            #     var.append(p)
            # elif i%(n+2)==1: 
            #     # print(self.fileContents[i].split())
            #     try:
            #         time.append(float(self.fileContents[i].split()[2]))
            #     except: continue
            # else:
            #     line = self.fileContents[i].split()
            #     p.extend(line[1:4])
            #     if i <= n+2 :
            #         try:
            #             Z = self.pt.getZ(line[0])
            #         except:
            #             Z=int(line[0])
            #         self._elements.append(Z)    
                    
        self._var = np.array(var, dtype=np.float64).reshape(-1, n*3)
        self._time = np.array(time, dtype=np.float64)
        self._step = self._var.shape[0] 
        self.mass = np.array(self.mass)
    
    def rewrite(self, to_where):
        if to_where == -1:
            return 0
        
        with open(self.filename, 'w') as f:
            f.writelines(self.fileContents[:(self.n+2)*(to_where+1)])
    
    def xyz_to_zmatrix(self, frame=0):
        ...
        
    def get_dihedral(self, *args):
        a, b, c, d = args
        
        var = self._var.reshape(-1, self.n, 3)
        # print(var)
        
        vec1 = var[:,b] - var[:,a]
        vec2 = var[:,c] - var[:,b]
        vec3 = var[:,d] - var[:,c]
        
        n1 = np.cross(vec1, vec2)
        n2 = np.cross(vec2, vec3)
        
        sgn = np.sign( np.sum(np.cross(n1, n2) * vec2, axis=1))
        
        theta = np.sum(n1*n2, axis=1) / (np.linalg.norm(n1, axis=1) * np.linalg.norm(n2, axis=1))
        return sgn * (np.arccos(theta) * 180 / np.pi)
    
    def get_distance(self, *args):
        a, b = args
        var = self._var.reshape(-1, self.n, 3)
        
        dist = np.linalg.norm(var[:,a]-var[:,b], axis=1)
        return dist
    
    def get_angle(self, *args):
        a, b, c = args
        var = self._var.reshape(-1, self.n, 3)
        
        vec1 = var[:, a] - var[:, b]
        vec2 = var[:, c] - var[:, b]
        theta = np.sum(vec1*vec2, axis=1) / (np.linalg.norm(vec1, axis=1) * np.linalg.norm(vec2, axis=1))
        return np.arccos(theta) * 180 / np.pi
    
    def get_variable(self):

        return self._var
        

        # self._var = []
        
        # for line in self.fileContents[2:]:
        #     splitline = line.split()
        #     self._var.extend(splitline[1:4])
            
        # return np.array(self._var, dtype=np.float64)
        
    def to_mndo99(self, name):
        Input = My_InputFile(f'{name}.input')
        Input.writein(self._var, self._elements, 2)
        self._output_name = f'{name}.output'
        Input.update(self._output_name)


def calc_dist(R, *args):
    if len(R.shape) != 2:
        raise ValueError('dimension must be 2')
    a, b = args
    dist = np.sqrt( np.sum( (R[a]-R[b])**2 ) )
    return dist

def calc_angle(R, *args):
    if len(R.shape) != 2:
        raise ValueError('dimension must be 2')
    a, b, c = args
    vec1 = R[a] - R[b]
    vec2 = R[c] - R[b]
    theta = np.sum(vec1*vec2) / np.sqrt(np.sum(vec1**2) * np.sum(vec2**2))
    return math.acos(theta) * 180 / math.pi

def calc_dihedral(R, *args):
    if len(R.shape) != 2:
        raise ValueError('dimension must be 2')    
    a, b, c, d = args
    vec1 = R[b] - R[a]
    vec2 = R[c] - R[b]
    vec3 = R[d] - R[c]    
    n1 = np.cross(vec1, vec2)
    n2 = np.cross(vec2, vec3)  
    sgn = np.sign(np.sum(np.cross(n1, n2) * vec2))
    theta = np.sum(n1*n2, axis=1) / np.sqrt(np.sum(n1**2) * np.sum(n2**2))
    return sgn * (math.acos(theta) * 180 / math.pi)      

def cart2zmat(coord):
    '''
    -------------------
    Copy from pyscf
    -------------------
    >>> c = numpy.array((
    (0.000000000000,  1.889726124565,  0.000000000000),
    (0.000000000000,  0.000000000000, -1.889726124565),
    (1.889726124565, -1.889726124565,  0.000000000000),
    (1.889726124565,  0.000000000000,  1.133835674739)))
    >>> print(cart2zmat(c))
    1
    1 2.67247631453057
    1 4.22555607338457 2 50.7684795164077
    1 2.90305235726773 2 79.3904651036893 3 6.20854462618583
    '''
    zstr = []
    zstr.append('1')
    if len(coord) > 1:
        r1 = coord[1] - coord[0]
        nr1 = np.linalg.norm(r1)
        zstr.append('1 %.15g' % nr1)
    if len(coord) > 2:
        r2 = coord[2] - coord[0]
        nr2 = np.linalg.norm(r2)
        a = np.arccos(np.dot(r1,r2)/(nr1*nr2))
        zstr.append('1 %.15g 2 %.15g' % (nr2, a*180/np.pi))
    if len(coord) > 3:
        o0, o1, o2 = coord[:3]
        p0, p1, p2 = 1, 2, 3
        for k, c in enumerate(coord[3:]):
            r0 = c - o0
            nr0 = np.linalg.norm(r0)
            r1 = o1 - o0
            nr1 = np.linalg.norm(r1)
            a1 = np.arccos(np.dot(r0,r1)/(nr0*nr1))
            b0 = np.cross(r0, r1)
            nb0 = np.linalg.norm(b0)

            if abs(nb0) < 1e-7: # o0, o1, c in line
                a2 = 0
                zstr.append('%d %.15g %d %.15g %d %.15g' %
                            (p0, nr0, p1, a1*180/np.pi, p2, a2))
            else:
                b1 = np.cross(o2-o0, r1)
                nb1 = np.linalg.norm(b1)

                if abs(nb1) < 1e-7:  # o0 o1 o2 in line
                    a2 = 0
                    zstr.append('%d %.15g %d %.15g %d %.15g' %
                                (p0, nr0, p1, a1*180/np.pi, p2, a2))
                    o2 = c
                    p2 = 4 + k
                else:
                    if np.dot(np.cross(b1, b0), r1) < 0:
                        a2 = np.arccos(np.dot(b1, b0) / (nb0*nb1))
                    else:
                        a2 =-np.arccos(np.dot(b1, b0) / (nb0*nb1))
                    zstr.append('%d %.15g %d %.15g %d %.15g' %
                                (p0, nr0, p1, a1*180/np.pi, p2, a2*180/np.pi))

    return '\n'.join(zstr)


def xyz_to_zmat(R, elements):
    if len(R.shape) != 2:
        raise ValueError('dimension must be 2')
    recap = [1 if i==1 else 0 for i in elements] #把H标注为1， 其余标为0，方便排序
    zmat = []
    links = []
    for i in range(len(elements)):
        if i==0:
            zmat.extend([0,0,0])
            links.extend([0,0,0])
        elif i==1:
            dist = calc_dist(R, 0, 1) 
            zmat.extend([dist, 0, 0])
            links.extend(1,0,0)
        elif i==2:
            dist = [calc_dist(R, 2, i) for i in [0,1]]
            sortlist = np.argsort(dist)
            links.extend(sortlist[0]+1, sortlist[1]+1, 0)
            angle = calc_angle(R, 2, sortlist[0], sortlist[1])
            zmat.extend([dist[sortlist[0]], angle, 0])
        else:
            pre_atom = list(range(i))
            dist = [calc_dist(R, i, j) for j in pre_atom] #计算前n-1个与第n个的距离
            sortlist = np.lexsort(dist, recap[:i])
            a = sortlist[0]
            pre_atom.pop(a)
            dist2 = [calc_dist(R, a, j) for j in pre_atom]
            sortlist2 = np.lexsort(dist2, recap[:i])


def delete_line(filename, to_where):
    with open(filename, 'r+') as f:
        lines = f.readlines()
        f.seek(0)
        f.truncate()
        f.writelines(lines[:to_where+1])

        

def read_xyz(xyzfile):
    with open(xyzfile, 'r') as f:
        contents = f.readlines()
    n = int(contents[0].split()[0])
    var = []
    for i in range(len(contents)):
        if i%(n+2)==0:
            p = []
            var.append(p)
        elif i%(n+2)==1: continue
        else:
            p.extend(contents[i].split()[1:4])
    return np.array(var, dtype=np.float64)        
   
class Read_molden(RawFile):
    def __init__(self, filename):
        super(Read_molden, self).__init__(filename)
        self.pt = PeriodicTable()
    
    def findfreq(self):
        startline = self.searchForString('[FREQ]', False)
        if startline == -1:
            raise KeyError('No FREQ attribute on molden file!')
        
        num = int(self.fileContents[startline].split()[1]) - 6
        self.freq = []
        for i in range(num):
            self.freq.append( float( self.fileContents[startline+1+i].split()[0] ) )
        
        self.freq = np.array(self.freq)
        return self.freq
    
    def findNormalMode(self):
        startline = self.searchForString('[FR-NORM-COORD]', False)
        if startline == -1:
            raise KeyError('No FR-NORM-COORD attribute on molden file!')
        
        atom_num = int(self.fileContents[startline].split()[1])//3
        mode_num = atom_num*3 - 6
        self.normalMode = []
        for i in range(mode_num):
            vib = []
            self.normalMode.append(vib)
            for j in range(atom_num):
                vib.extend(self.fileContents[startline+i*(atom_num+1)+j+2].split())
        
        self.normalMode = np.array(self.normalMode, dtype=float)
        return self.normalMode
    
    def findoptGeom(self, frame=-1):
        '''
        frame specify which frame as optimized geometry
        default -1 (the final one)
        '''
        startline = self.searchForString('[GEOMETRIES]', False)
        if startline == -1:
            raise KeyError('No GEOMETRIES attribute on molden file!')
        
        num = int(self.fileContents[startline+1].split()[0])
        # frameline = startline+1
        if frame == -1:
            self.opt_xyz = self.fileContents[-num-2:]
            # return self.opt_xyz
            return ''.join(self.opt_xyz)
    
    def get_mass(self):
       return self.mass
        
    
    def get_opt_cord(self):
        try: self.opt_xyz
        except: self.findoptGeom()

        self.elements = []
        self.mass = []
        self.opt_position = []

        for line in self.opt_xyz[2:]:
            i = line.split()
            try:
                Z = self.pt.getZ(i[0])
            except:
                Z=int(i[0])
            self.elements.append(Z)
            self.mass.append([self.pt.mass[Z]]*3)
            self.opt_position.extend(i[1:4])
        
        self.mass = np.array(self.mass).reshape(-1)
        # n = int(self.opt_xyz[0].split()[0])
        # for i in range(len(self.opt_xyz)):
        #     if i%(n+2)==0:
        #         p = []
        #     elif i%(n+2)==1: continue
        #     else:
        #         p.extend(self.opt_xyz[i].split()[1:4])
        return np.array(self.opt_position, dtype=np.float64)        
    
    def optGeom2xyz(self, xyzName):
        with open(xyzName, 'w') as f:
            try: self.opt_xyz
            except: self.findoptGeom()
            f.writelines(self.opt_xyz)
        

def array_to_xyz(filename, position, elements):
    n = len(elements)
    position = position.reshape(-1,3)
    with open(filename, 'w') as f:
        f.write(f'  {n}\n \n')
        for i in range(n):
            f.write(
    f' {elements[i]}    {position[i][0]}    {position[i][1]}    {position[i][2]}\n'
                )
        

class FormchkInterface:
    '''
    copy from pyxDH
    https://github.com/ajz34/Py_xDH/blob/master/pyxdh/Utilities/formchk_interface.py
    '''
    def __init__(self, file_path):
        self.file_path = file_path
        self.natm = NotImplemented
        self.nao = NotImplemented
        self.nmo = NotImplemented
        self.initialization()

    def initialization(self):
        self.natm = int(self.key_to_value("Number of atoms"))
        self.nao = int(self.key_to_value("Number of basis functions"))
        self.nmo = int(self.key_to_value("Number of independent functions"))

    def key_to_value(self, key, file_path=None):
        if file_path is None:
            file_path = self.file_path
        flag_read = False
        expect_size = -1
        vec = []
        with open(file_path, "r") as file:
            for l in file:
                if l[:len(key)] == key:
                    try:
                        expect_size = int(l[len(key):].split()[2])
                        flag_read = True
                        continue
                    except IndexError:
                        try:
                            return float(l[len(key):].split()[1])
                        except IndexError:
                            continue
                if flag_read:
                    try:
                        vec += [float(i) for i in l.split()]
                    except ValueError:
                        break
        if len(vec) != expect_size:
            raise ValueError("Number of expected size is not consistent with read-in size!")
        return np.array(vec)

    def total_energy(self, file_path=None):
        if file_path is None:
            file_path = self.file_path
        return self.key_to_value("Total Energy", file_path)

    def grad(self, file_path=None):
        if file_path is None:
            file_path = self.file_path
        return self.key_to_value("Cartesian Gradient", file_path).reshape((self.natm, 3))

    def dipole(self, file_path=None):
        if file_path is None:
            file_path = self.file_path
        return self.key_to_value("Dipole Moment", file_path)

    @staticmethod
    def tril_to_symm(tril: np.ndarray):
        dim = int(np.floor(np.sqrt(tril.size * 2)))
        if dim * (dim + 1) / 2 != tril.size:
            raise ValueError("Size " + str(tril.size) + " is probably not a valid lower-triangle matrix.")
        indices_tuple = np.tril_indices(dim)
        iterator = zip(*indices_tuple)
        symm = np.empty((dim, dim))
        for it, (row, col) in enumerate(iterator):
            symm[row, col] = tril[it]
            symm[col, row] = tril[it]
        return symm

    def hessian(self, file_path=None):
        if file_path is None:
            file_path = self.file_path
        return self.tril_to_symm(self.key_to_value("Cartesian Force Constants", file_path))

    def polarizability(self, file_path=None):
        if file_path is None:
            file_path = self.file_path
        # two space after `Polarizability' is to avoid `Polarizability Derivative'
        return self.tril_to_symm(self.key_to_value("Polarizability  ", file_path))

    def dipolederiv(self, file_path=None):
        if file_path is None:
            file_path = self.file_path
        return self.key_to_value("Dipole Derivatives", file_path).reshape(-1, 3)
#file = OutputFile('cis_azobene2_trajectory_-8/_step0.output')
# a = setInputFile('testinput')
# a.load_mndoinput('tpe_3/_step0.input')
#a.loadgeom('testmodel.gjf')
#
#b = a.get_elements()
#c = a.get_position()
#d = inputFile('testinput1')
#d.writein(c,b)


def main():
    print("test")


if __name__ == '__main__':
    main()