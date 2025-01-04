#!/bin/env python3
#!filename: om.py

import numpy as np
import numpy.linalg as la
import sys, os, copy, string
np.set_printoptions(linewidth=400)

PeriodicTable = '''      -----                                                               -----
1 | H |                                                               |He |
  |---+----                                       --------------------+---|
2 |Li |Be |                                       | B | C | N | O | F |Ne |
  |---+---|                                       |---+---+---+---+---+---|
3 |Na |Mg |3B  4B  5B  6B  7B |    8B     |1B  2B |Al |Si | P | S |Cl |Ar |
  |---+---+---------------------------------------+---+---+---+---+---+---|
4 | K |Ca |Sc |Ti | V |Cr |Mn |Fe |Co |Ni |Cu |Zn |Ga |Ge |As |Se |Br |Kr |
  |---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---|
5 |Rb |Sr | Y |Zr |Nb |Mo |Tc |Ru |Rh |Pd |Ag |Cd |In |Sn |Sb |Te | I |Xe |
  |---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---|
6 |Cs |Ba |LAN|Hf |Ta | W |Re |Os |Ir |Pt |Au |Hg |Tl |Pb |Bi |Po |At |Rn |
  |---+---+---+------------------------------------------------------------
7 |Fr |Ra |ACT|
  -------------
              -------------------------------------------------------------
   Lanthanide |La |Ce |Pr |Nd |Pm |Sm |Eu |Gd |Tb |Dy |Ho |Er |Tm |Yb |Lu |
              |---+---+---+---+---+---+---+---+---+---+---+---+---+---+---|
   Actinide   |Ac |Th |Pa | U |Np |Pu |Am |Cm |Bk |Cf |Es |Fm |Md |No |Lw |
              -------------------------------------------------------------'''


at_list = ['X', 'H', 'He',
'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']

at_dict = {'X':0, 'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6, 'N':7, 'O':8, 'F':9, 'Ne':10,
'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P':15, 'S':16, 'Cl':17, 'Ar':18, 'K':19, 'Ca':20,
'Sc':21, 'Ti':22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30,
'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36}

# relative atomic masses of elements (in atomic mass units [amu]) from
# "CRC Handbook" 84th ed, ed Lide, pgs 1-12 - 1-14
at_mass = {    'H' : 1.00794, 'C' : 12.0107, 'O' : 15.9994, 'N' : 14.0067,
  'F' : 18.9984, 'P' : 30.9738, 'S' : 32.0650, 'Cl': 35.4530, 'Br': 79.9040,
  'I' : 126.904, 'He': 4.00260, 'Ne': 20.1797, 'Ar': 39.9480, 'Li': 6.94100,
  'Be': 9.01218, 'B' : 10.8110, 'Na': 22.9898, 'Mg': 24.3050, 'Al': 26.9815,
  'Si': 28.0855, 'K' : 39.0983, 'Ca': 40.0780, 'Sc': 44.9559, 'Ti': 47.8670,
  'V' : 50.9415, 'Cr': 51.9961, 'Mn': 54.9380, 'Fe': 55.8450, 'Co': 58.9332,
  'Ni': 58.6934, 'Cu': 63.5460, 'Zn': 65.4090, 'Ga': 69.7230, 'Ge': 72.6400,
  'As': 74.9216, 'Se': 78.9600, 'Kr': 83.7980, 'X' :  0.0000}

## CONSTANTS ##
bondcut = 2.0 # Anstrong
rad2deg = 180.0 / np.pi
deg2rad = np.pi / 180.0

# read file data into a 2-d array
def get_file_string_array(file):
    try:
        file = open(file, "r")
    except IOError:
        print('Error: file (%s) not found!\n' % (file))
        sys.exit()
    lines = file.readlines()
    file.close()
    array = []
    for line in lines:
        array.append(line.split())
    return array


################### VEC3D MATH FUNCTIONS##################
'''
    using vectorization programs via numpy
'''
# outer product
def vec3d_outer(V1, V2):
    """ input numpy array
    V1 ~ (3) or (:,3)
    V2 ~ (3) or (:,3)
    """
    if(len(np.shape(V1))== 1):
        A = V1.reshape((len(V1)//3, 3))
        B = V2.reshape((len(V2)//3, 3))
    elif(len(np.shape(V1))== 2):
        A = V1; B = V2
    C = np.zeros(np.shape(A))
    C[:,0] = A[:,1]*B[:,2] - A[:,2]*B[:,1]
    C[:,1] = A[:,2]*B[:,0] - A[:,0]*B[:,2]
    C[:,2] = A[:,0]*B[:,1] - A[:,1]*B[:,0]
    return C.reshape(np.shape(V1))

# inner product
def vec3d_inner(V1, V2):
    """ input numpy array
    V1 ~ (3) or (:,3)
    V2 ~ (3) or (:,3)
    """
    if(len(np.shape(V1))== 1):
        A = V1.reshape((len(V1)//3, 3))
        B = V2.reshape((len(V2)//3, 3))
        C = np.zeros((len(V1)//3))
    elif(len(np.shape(V1))== 2):
        A = V1; B = V2
        C = np.zeros((len(V1[:,0])))
    C[:] = A[:,0]*B[:,0] + A[:,1]*B[:,1] + A[:,2]*B[:,2]
    if(len(np.shape(V1))== 1):
        return C[0]
    elif(len(np.shape(V1))== 2):
        return C

# bond distance
def vec3d_bond(V1, V2):
    """ input numpy array
    V1 ~ (3) or (:,3)
    V2 ~ (3) or (:,3)
    """
    return np.sqrt( vec3d_inner(V1-V2, V1-V2) )


# unit vectors
def vec3d_uvec(V):
    """ input numpy array
    V ~ (3) or (:,3)
    """
    if(len(np.shape(V))== 1):
        return V / vec3d_bond(V, np.zeros(np.shape(V)))
    elif(len(np.shape(V))== 2):
        return V / np.outer( vec3d_bond(V, np.zeros(np.shape(V))), np.ones((3)) )

# cos of angles for point sets V1, V2, V3
def vec3d_cosangle(V1, V2, V3):
    """ input numpy array
    V1 ~ (3) or (:,3)
    V2 ~ (3) or (:,3)
    V3 ~ (3) or (:,3)
    """
    b1 = vec3d_bond(V1, V2)
    b2 = vec3d_bond(V3, V2)
    return vec3d_inner(V1-V2, V3-V2) / (b1*b2)

# angles for point sets V1, V2, V3
def vec3d_angle(V1, V2, V3):
    """ input numpy array
    V1 ~ (3) or (:,3)
    V2 ~ (3) or (:,3)
    V3 ~ (3) or (:,3)
    """
    return np.arccos(vec3d_cosangle(V1, V2, V3))

# cos of dihedrals for point sets V1, V2, V3, V4
def vec3d_cosdihedral(V1, V2, V3, V4):
    """ input numpy array
    V1 ~ (3) or (:,3)
    V2 ~ (3) or (:,3)
    V3 ~ (3) or (:,3)
    V4 ~ (3) or (:,3)
    """
    X1 = V2 - V1; X2 = V3 - V2; X3 = V4 - V3
    Y1 = vec3d_outer(X1, X2)
    Y2 = vec3d_outer(X2, X3)
    Z = np.zeros(np.shape(Y1))
    return vec3d_cosangle(Y1+Z, Z, Y2+Z)

# dihedrals for point sets V1, V2, V3, V4
def vec3d_dihedral(V1, V2, V3, V4):
    """ input numpy array
    V1 ~ (3) or (:,3)
    V2 ~ (3) or (:,3)
    V3 ~ (3) or (:,3)
    V4 ~ (3) or (:,3)
    """
    X1 = V2 - V1; X2 = V3 - V2; X3 = V4 - V3
    Y1 = vec3d_outer(X1, X2)
    Y2 = vec3d_outer(X2, X3)
    s = np.sign( vec3d_inner(vec3d_outer(Y1, Y2), X2) )
    Z = np.zeros(np.shape(Y1))
    return np.arccos( vec3d_cosangle(Y1+Z, Z, Y2+Z) ) * s

# local axis system from V2, V3, V4 (to decide V1's zmatrice)
def vec3d_axis(V2, V3, V4):
    """ input numpy array @TODO
    V2 ~ (3) or (:,3), bond: V1-V2
    V3 ~ (3) or (:,3), angle: V1-V2-V3
    V4 ~ (3) or (:,3), dihedral: V1-V2-V3-V4
    """
    # left hand local axis
    Zuvec = vec3d_uvec( V3-V2 )                      # z axis
    Yuvec = vec3d_uvec( - vec3d_outer(V3-V2,V4-V3) ) # y axis
    Xuvec = vec3d_uvec( vec3d_outer(Zuvec, Yuvec)  ) # x axis
    return [Xuvec, Yuvec, Zuvec]

# displace with d
def vec3d_displace(V, d):
    """ input numpy array
    V ~ (3) or (:,3)
    d ~ (3)
    """
    if(len(np.shape(V))== 1):
        A = V.reshape((len(V)//3, 3))
    elif(len(np.shape(V1))== 2):
        A = V
    dc = d.reshape((3))
    return A + np.outer( np.ones((len(V)//3)), dc )


# rotate with theta aound n
def vec3d_rotate(V, n, theta):
    """ input numpy array
    V ~ (3) or (:,3)
    n ~ (3)
    theta: scalar
    """
    if(len(np.shape(V))== 1):
        A = V.reshape((len(V)//3, 3))
    elif(len(np.shape(V1))== 2):
        A = V
    nc = vec3d_uvec(n).reshape((3)) # normalization of direction cosine
    C = np.zeros((3,3)); nouter = np.zeros((3,3))
    nouter[0,1] = -nc[2]; nouter[0,2] = nc[1]; nouter[1,0] = nc[2];
    nouter[1,2] = -nc[0]; nouter[2,0] = -nc[1]; nouter[2,1] = nc[0];
    C = np.outer(nc,nc)*(1-np.cos(theta)) + np.eye(3)*np.cos(theta) + nouter*np.sin(theta)
    return (np.dot(C, A.T)).T


# molecule class for molecular data
class mol_class:
    _before = ''
    _after  = ''
    _has_zmat = False
    _has_cart = False

    # constructor
    def __init__(self, file, type='zmat', fmt='mndo'):
        self.file = file
        self.type = type
        self.fmt  = fmt
        self.name = self.file.split('/')[-1].split('.')[0]
        try:
            if self.type == 'zmat':
                self.get_zmat(self.file)
                self._has_zmat = True
            elif self.type == 'cart':
                self.get_cart(self.file)
                self._has_cart = True
            else:
                raise ValueError('type: zmat or cart')
        except ValueError as e:
            print(repr(e))

    def get_cart(self, file):
        if self.fmt == 'mndo':
            array = get_file_string_array(file)
            isbefore = True; isafter = False
            self.data = []
            for i in range(len(array)):
                if array[i] == []:
                    if isbefore:
                        self._before += '\n'
                    if isafter:
                        self._after += '\n'
                    continue
                if array[i][0] == 'OM2':
                    self._before += ' '.join(array[i])
                    isbefore = False
                    continue
                if(array[i][0] == '0'):
                    self._after += ' '.join(array[i]) + '\n'
                    isafter = True
                    continue
                if isbefore:
                    self._before += ' '.join(array[i]) + '\n'
                elif isafter:
                    self._after += ' '.join(array[i]) + '\n'
                else:
                    self.data.append(array[i])
            tmp = np.array(self.data)
            self.at_idx = tmp[:,0].astype(np.int)
            self.cart = tmp[:,1:6:2].astype(np.float)
            self.cart_fix = tmp[:,2::2].astype(np.int)
            self.natom = len(self.at_idx)
        elif self.fmt == 'gauss':
            array = get_file_string_array(file)
            isbefore = True; isafter = False
            self.data = []
            for i in range(len(array)):
                if array[i] == []:
                    if isbefore:
                        self._before += '\n'
                    elif isafter:
                        self._after += '\n'
                    else:
                        isafter = True
                        self._after += '\n'
                    continue
                if array[i][0] == '0' and array[i][1] == '1':
                    self._before += ' '.join(array[i])
                    isbefore = False
                    continue
                if isbefore:
                    self._before += ' '.join(array[i]) + '\n'
                elif isafter:
                    self._after += ' '.join(array[i]) + '\n'
                else:
                    self.data.append(array[i])
            tmp = np.array(self.data)
            self.at_flag = tmp[:,0].astype(np.str); self.at_idx = np.zeros(np.shape(self.at_flag), dtype=np.int)
            self.cart = tmp[:,1:].astype(np.float)
            self.cart_fix = np.ones(np.shape(self.cart))
            self.natom = len(self.at_flag)
            for i in range(self.natom):
                self.at_idx[i] = at_dict.get(self.at_flag[i].strip(string.digits), 0)
        elif self.fmt == 'xyz':
            array = get_file_string_array(file)
            isbefore = True; isafter = False
            self.data = []
            self.natom = int(array[0][0])
            for i in range(2,2+self.natom):
                self.data.append(array[i])
            tmp = np.array(self.data)
            self.at_flag = tmp[:,0].astype(np.str); self.at_idx = np.zeros(np.shape(self.at_flag), dtype=np.int)
            self.cart = tmp[:,1:].astype(np.float)
            self.cart_fix = np.ones(np.shape(self.cart))
            self.natom = len(self.at_flag)
            for i in range(self.natom):
                self.at_idx[i] = at_dict.get(self.at_flag[i].strip(string.digits), 0)

    # read in z-matrix from zmat file
    def get_zmat(self, file):
        if self.fmt == 'mndo':
            array = get_file_string_array(file)
            isbefore = True; isafter = False
            self.data = []
            for i in range(len(array)):
                if array[i] == []:
                    if isbefore:
                        self._before += '\n'
                    if isafter:
                        self._after += '\n'
                    continue
                if array[i][0] == 'OM2':
                    self._before += ' '.join(array[i])
                    isbefore = False
                    continue
                if(array[i][0] == '0'):
                    self._after += ' '.join(array[i]) + '\n'
                    isafter = True
                    continue
                if isbefore:
                    self._before += ' '.join(array[i]) + '\n'
                elif isafter:
                    self._after += ' '.join(array[i]) + '\n'
                else:
                    self.data.append(array[i])
            tmp = np.array(self.data)
            self.at_idx = tmp[:,0].astype(np.int)
            self.zmat = tmp[:,1:6:2].astype(np.float)
            self.zmat[:,1:] *= deg2rad # deg2rad
            self.zmat_fix = tmp[:,2:7:2].astype(np.int)
            self.zmat_ref = tmp[:,7:].astype(np.int) - 1
            self.natom = len(self.at_idx)
        elif self.fmt=='gauss':
            array = get_file_string_array(file)
            isbefore = True; isafter = False
            self.data = []
            for i in range(len(array)):
                if array[i] == []:
                    if isbefore:
                        self._before += '\n'
                    elif isafter:
                        self._after += '\n'
                    else:
                        isafter = True
                        self._after += '\n'
                    continue
                if array[i][0] == '0' and array[i][1] == '1':
                    self._before += ' '.join(array[i])
                    isbefore = False
                    continue
                if isbefore:
                    self._before += ' '.join(array[i]) + '\n'
                elif isafter:
                    self._after += ' '.join(array[i]) + '\n'
                else:
                    self.data.append(array[i])
            self.data[0] += ['0']*7
            self.data[1] += ['0']*5
            self.data[2] += ['0']*3
            tmp = np.array(self.data)
            self.at_flag = tmp[:,0].astype(np.str); self.at_idx = np.zeros(np.shape(self.at_flag), dtype=np.int)
            self.zmat = tmp[:,2::2].astype(np.float); self.zmat[:,1:] *= deg2rad # deg2rad
            self.zmat_ref = tmp[:,1:-1:2].astype(np.int) - 1
            self.zmat_fix = np.ones(np.shape(self.zmat))
            self.natom = len(self.at_flag)
            for i in range(self.natom):
                self.at_idx[i] = int( at_dict.get(self.at_flag[i], 0) )

    def topology(self):
        self.dist = np.zeros((self.natom, self.natom))
        for i in range(self.natom):
            for j in range(i):
                self.dist[i,j] = vec3d_bond(self.cart[i,:], self.cart[j,:])
                self.dist[j,i] = self.dist[i,j]
        self.blist = []
        for i in range(self.natom):
            self.blist.append( np.where( (self.dist[i,:]-1e-8)*(self.dist[i,:]-bondcut) < 0 )[0] )
        # for i in range(self.natom):
        #     print(self.blist[i], self.dist[i,self.blist[i]])
        self.topmat = np.zeros((self.natom, self.natom))
        for i in range(self.natom):
            for j in self.blist[i]:
                if j>=i:
                    continue
                self.topmat[i,j] = 1
                self.topmat[j,i] = 1

    # obtain cartesian xyz-cartinates from z-matrix values
    def zmat2cart(self):
        self.cart = np.zeros((self.natom, 3))
        if self.natom >= 1:
            self.cart[0,:] = 0
        if self.natom >= 2:
            self.cart[1,0] = self.zmat[1,0]
        if self.natom >= 3:
            if self.zmat_ref[2,0] == 1:
                self.cart[2,0] = self.zmat[1,0] - self.zmat[2,0] * np.cos(self.zmat[2,1])
                self.cart[2,1] = self.zmat[2,0] * np.sin(self.zmat[2,1])
            elif self.zmat_ref[2,0] == 0:
                self.cart[2,0] = self.zmat[2,0] * np.cos(self.zmat[2,1])
                self.cart[2,1] = self.zmat[2,0] * np.sin(self.zmat[2,1])
        for i in range(3, self.natom):
            carts2 = self.cart[self.zmat_ref[i,0],:]
            carts3 = self.cart[self.zmat_ref[i,1],:]
            carts4 = self.cart[self.zmat_ref[i,2],:]
            [X, Y, Z] = vec3d_axis(carts2, carts3, carts4)
            disp = ( X * self.zmat[i,0] * np.sin(self.zmat[i,1]) * np.cos(self.zmat[i,2])
                   + Y * self.zmat[i,0] * np.sin(self.zmat[i,1]) * np.sin(self.zmat[i,2])
                   + Z * self.zmat[i,0] * np.cos(self.zmat[i,1]) )
            carts1 = carts2 + disp
            self.cart[i,:] = carts1
        self.cart_fix = np.ones(np.shape(self.cart)).astype(np.int)
        self._has_cart = True

    # obtain cartesian xyz-cartinates from z-matrix values
    def cart2zmat(self):
        self.topology()
        self.zmat = np.zeros((self.natom, 3))
        self.zmat_ref = np.zeros((self.natom, 3), dtype=np.int)
        if self.natom >= 1:
            self.zmat[0,:] = 0
        if self.natom >= 2:
            self.zmat_ref[1,0] = 0
            self.zmat[1,0] = self.dist[1,0]
        if self.natom >= 3:
            if self.dist[2,0] < self.dist[2,1]:
                self.zmat_ref[2,0] = 0; self.zmat_ref[2,1] = 1;
                self.zmat[2,0] = self.dist[2,0]
                self.zmat[2,1] = vec3d_angle(self.cart[2,:], self.cart[0,:], self.cart[1,:])
            else:
                self.zmat_ref[2,0] = 1; self.zmat_ref[2,1] = 0;
                self.zmat[2,0] = self.dist[2,1]
                self.zmat[2,1] = vec3d_angle(self.cart[2,:], self.cart[1,:], self.cart[0,:])
        for i in range(3, self.natom):
            bone = [0,0,0,0]; rcdtmp = 999
            for s2 in self.blist[i]:
                if s2>=i:
                    continue
                for s3 in self.blist[s2]:
                    if s3>=i or s3==s2:
                        continue
                    for s4 in self.blist[s3]:
                        if s4>=i or s4==s2 or s4==s3:
                            continue
                        sumbond = self.dist[i,s2]+self.dist[s2,s3]+self.dist[s3,s4]
                        if(sumbond < rcdtmp):
                            bone = [i,s2,s3,s4]
                            rcdtmp = sumbond
            if bone[0] == i:
                self.zmat_ref[i,:] = bone[1:]
                carts1 = self.cart[i,:]
                carts2 = self.cart[self.zmat_ref[i,0],:]
                carts3 = self.cart[self.zmat_ref[i,1],:]
                carts4 = self.cart[self.zmat_ref[i,2],:]
                self.zmat[i,0] = vec3d_bond(carts1, carts2)
                self.zmat[i,1] = vec3d_angle(carts1, carts2, carts3)
                self.zmat[i,2] = vec3d_dihedral(carts1, carts2, carts3, carts4)

        self.zmat[:,1:] = self.zmat[:,1:]%(2*np.pi)
        self.zmat_fix = np.ones(np.shape(self.zmat)).astype(np.int)
        self.zmat_fix[0,:] = 0; self.zmat_fix[1,1:] = 0; self.zmat_fix[2,2:] = 0
        self._has_zmat = True

    def _string(self, fmt='xyz', type='cart'):
        istr = ''
        if type=='zmat' and self._has_zmat:
            if fmt=='mndo':
                if(self.fmt=='mndo'):
                    istr += self._before+'\n'
                else:
                    istr += 'IOP=-6 JOP=0\n\nOM2' + '\n'
                for i in range(self.natom):
                    if i==0:
                        istr += (' %d\t%12.10f\t%d\t%12.10f\t%d\t%12.10f\t%d\t%d\t%d\t%d\n'%
                            (self.at_idx[i], self.zmat[i,0], 0,
                            self.zmat[i,1]*rad2deg, 0,
                            self.zmat[i,2]*rad2deg, 0,
                            0, 0, 0) )
                    elif i==1:
                        istr += (' %d\t%12.10f\t%d\t%12.10f\t%d\t%12.10f\t%d\t%d\t%d\t%d\n'%
                            (self.at_idx[i], self.zmat[i,0], self.zmat_fix[i,0],
                            self.zmat[i,1]*rad2deg, 0,
                            self.zmat[i,2]*rad2deg, 0,
                            self.zmat_ref[i,0]+1, 0, 0) )
                    elif i==2:
                        istr += (' %d\t%12.10f\t%d\t%12.10f\t%d\t%12.10f\t%d\t%d\t%d\t%d\n'%
                            (self.at_idx[i], self.zmat[i,0], self.zmat_fix[i,0],
                            self.zmat[i,1]*rad2deg, self.zmat_fix[i,1],
                            self.zmat[i,2]*rad2deg, 0,
                            self.zmat_ref[i,0]+1, self.zmat_ref[i,1]+1, 0) )
                    else:
                        istr += (' %d\t%12.10f\t%d\t%12.10f\t%d\t%12.10f\t%d\t%d\t%d\t%d\n'%
                            (self.at_idx[i], self.zmat[i,0], self.zmat_fix[i,0],
                            self.zmat[i,1]*rad2deg, self.zmat_fix[i,1],
                            self.zmat[i,2]*rad2deg, self.zmat_fix[i,2],
                            self.zmat_ref[i,0]+1, self.zmat_ref[i,1]+1, self.zmat_ref[i,2]+1) )
                if(self.fmt=='mndo'):
                    istr += self._after + '\n'
                else:
                    istr += ' %d\t%12.10f\t%d\t%12.10f\t%d\t%12.10f\t%d\t%d\t%d\t%d\n'%(0,0,0,0,0,0,0,0,0,0)
            if fmt=='gauss':
                if(self.fmt=='gauss'):
                    istr += self._before + '\n'
                else:
                    istr += '# sp hf/STO-3G\n\nTitle Card Required\n\n0 1\n'
                if self.natom >= 1:
                    istr += ' %s\n'%at_list[self.at_idx[0]]
                if self.natom >= 2:
                    istr += (' %s\t%d\t%12.10f\n'%
                        (at_list[self.at_idx[1]], self.zmat_ref[1,0]+1, self.zmat[1,0]) )
                if self.natom >= 3:
                    istr += (' %s\t%d\t%12.10f\t%d\t%12.10f\n'%
                        (at_list[self.at_idx[2]], self.zmat_ref[2,0]+1, self.zmat[2,0],
                        self.zmat_ref[2,1]+1, self.zmat[2,1]*rad2deg ) )
                for i in range(3,self.natom):
                    istr += (' %s\t%d\t%12.10f\t%d\t%12.10f\t%d\t%12.10f\t%d\n'%
                        (at_list[self.at_idx[i]], self.zmat_ref[i,0]+1, self.zmat[i,0],
                        self.zmat_ref[i,1]+1, self.zmat[i,1]*rad2deg,
                        self.zmat_ref[i,2]+1, self.zmat[i,2]*rad2deg, 0) )
                if(self.fmt=='gauss'):
                    istr += self._after + '\n'
                else:
                    istr += '\n\n\t'
        elif type=='cart' and self._has_cart:
            if fmt == 'mndo':
                if(self.fmt=='mndo'):
                    istr += self._before + '\n'
                else:
                    istr += 'IOP=-6 JOP=0\n\nOM2\n'
                for i in range(self.natom):
                    istr += (' %d\t%12.10f\t%d\t%12.10f\t%d\t%12.10f\t%d\n'%
                        (self.at_idx[i], self.cart[i,0], self.cart_fix[i,0],
                        self.cart[i,1], self.cart_fix[i,1],
                        self.cart[i,2], self.cart_fix[i,2] ) )
                if(self.fmt=='mndo'):
                    istr += self._after + '\n'
                else:
                    istr += ' %d\t%12.10f\t%d\t%12.10f\t%d\t%12.10f\t%d\n'%(0,0,0,0,0,0,0)
            elif fmt == 'gauss':
                if(self.fmt=='gauss'):
                    istr += self._before + '\n'
                else:
                    istr += '# sp hf/STO-3G\n\nTitle Card Required\n\n0 1\n'
                for i in range(self.natom):
                    istr += (' %s\t%12.10f\t%12.10f\t%12.10f\n'%
                        (at_list[self.at_idx[i]], self.cart[i,0], self.cart[i,1],
                        self.cart[i,2] ) )
                if(self.fmt=='gauss'):
                    istr += self._after + '\n'
                else:
                    istr += '\n\n\t'
            elif fmt == 'xyz':
                istr += '%d\n'%self.natom
                istr += self.name + '\n'
                for i in range(self.natom):
                    istr += (' %s\t%12.10f\t%12.10f\t%12.10f\n'%
                        (at_list[self.at_idx[i]], self.cart[i,0], self.cart[i,1],
                        self.cart[i,2]) )
        else:
            print("\033[1;35mERROR\033[0m: don't have type of (%s)"%type)
        while istr[-1]=='\n':
            istr = istr[0:-1]
        return istr

    def print(self, fmt='xyz', type='cart'):
        print(self._string(fmt=fmt, type=type))

    def save(self, file='default.save', fmt='xyz', type='cart'):
        f=open(file, 'w')
        f.write(self._string(fmt=fmt, type=type))
        f.close()

    def mix(self, other, per):
        self.zmat[:,2] = (self.zmat[:,2]+np.pi)%(2*np.pi) - np.pi
        other.zmat[:,2] = (other.zmat[:,2]+np.pi)%(2*np.pi) - np.pi
        
        obj = copy.copy(other)
        obj.zmat = other.zmat * per + (1.0-per)*self.zmat
        #obj.zmat[13,2] = other.zmat[13,2] * per**6 + (1.0-per**6)*self.zmat[13,2]
        #print((np.abs(other.zmat[:,2] - self.zmat[:,2]) > np.pi).astype(np.int))

        #for i in
        shift = (np.abs(other.zmat[:,2] - self.zmat[:,2]) > np.pi).astype(np.int)
        shift *= (2*(other.zmat[:,2] < 0).astype(np.int)-1)
        obj.zmat[:,2] += per*shift * 2*np.pi
        obj.zmat2cart()
        #print(obj.zmat[:,1:]*rad2deg)
        return obj
        
    def ref_zmat(self, other):
        self.zmat2cart()
        #print(self.zmat)
        self.zmat_ref = other.zmat_ref        
        for i in range(3, self.natom):
            carts1 = self.cart[i,:]
            carts2 = self.cart[other.zmat_ref[i,0],:]
            carts3 = self.cart[other.zmat_ref[i,1],:]
            carts4 = self.cart[other.zmat_ref[i,2],:]
            self.zmat[i,0] = vec3d_bond(carts1, carts2)
            self.zmat[i,1] = vec3d_angle(carts1, carts2, carts3)
            self.zmat[i,2] = vec3d_dihedral(carts1, carts2, carts3, carts4)
            #print(i, self.zmat[i,:])
        self.zmat2cart()
        self.zmat[:,1:] = self.zmat[:,1:]%(2*np.pi)
        self.zmat_fix = np.ones(np.shape(self.zmat)).astype(np.int)
        self.zmat_fix[0,:] = 0; self.zmat_fix[1,1:] = 0; self.zmat_fix[2,2:] = 0
        self._has_zmat = True
        
        return None
    
    def mix_cart(self, other, per):
        obj = copy.copy(other)
        obj.cart = other.cart * per + (1.0-per)*self.cart
        #print((np.abs(other.zmat[:,2] - self.zmat[:,2]) > np.pi).astype(np.int))

        #for i in
        #ishift = (np.abs(other.zmat[:,2] - self.zmat[:,2]) > np.pi).astype(np.int)
        #shift *= (2*(other.zmat[:,2] < 0).astype(np.int)-1)
        #obj.zmat[:,2] += per*shift * 2*np.pi
        #obj.zmat2cart()
        #print(obj.zmat[:,1:]*rad2deg)
        return obj



if __name__ == '__main__':
    n = len(sys.argv)
    v = sys.argv
    dd = {}
    if n%2==0:
        flag_value = False
        dd['target'] = []
        for i in range(1,n):
            if v[i] == '--help':
                print("""
    [Usage]:

    python mgeom.py [--help | -arg flags] target1 target2 ...

    [-arg]  ------- [flags]
    -i  :-----|     [input-format] (see below)
    -o  :-----|     [output-format]
         mndo-zmat
         mndo-cart
         gauss-xyz
         gauss-cart
         xyz-cart
    -t  :-----|
    [                  ]:  default convert target from input-format to output-format
    [ -t mix           ]:  mix target with file1 (save files and print)
    [ -t ref -tf file1 ]:  rebuild targets' topology refering a zmat-format file1
    [ -t tpl -tf file1 ]:  use parameters in file1
    [ -t vmd           ]:  convert targets to xyz and open vmd to view it

    -------------------------------------------------------------------------------

    example 1]

        python mgeom.py -i gauss-zmat -o mndo-cart mol_zmat.gjf

        # convert gauss's zmat input to mndo's cartesian input

    example 2]

        python mgeom.py -i gauss-zmat -t vmd mol_zmat.mndoinp

        # convert mndo's zmat input to tmp.xyz and use vmd to open it

    example 3]

        python mgeom.py -i mndo-zmat -o mndo-zmat -t ref -tf refered_zmat.mndo mol_zmat.mndoinp

        # convert mndo's zmat input to new zmat input by refering topology in refered_zmat.mndo

    example 4]

        python mgeom.py -i mndo-zmat -o mndo-cart -t mix -tf mixed_zmat.mndo mol_zmat.mndoinp

        # mix mndo's zmat input with mixed_zmat.mndo, and print as mndo-cart format
    """)
                exit()
        
            if v[i][0] == '-':
                flag_value = True
                continue
            if flag_value:
                dd[v[i-1]] = v[i]
                flag_value = False
            else:
                dd['target'].append(v[i])

        fmti = dd.get('-i', 'mndo-zmat')  # supported: mndo-zmat, mndo-cart, gauss-xyz, gauss-cart, xyz-cart 
        fmto = dd.get('-o', 'xyz-cart')   # supported: mndo-zmat, mndo-cart, gauss-xyz, gauss-cart, xyz-cart 
        task = dd.get('-t', '')           # task type
        taskfn = dd.get('-tf', '')        # file for task

        tarlist = dd['target']
        # print(tarlist)

        for i in range(len(tarlist)):
            tar = tarlist[i]
            if task =='': # default
                ssi = fmti.split('-')
                sso = fmto.split('-')
                mol = mol_class(file=tar, fmt=ssi[0], type=ssi[1])
                if(ssi[1]=='zmat' and sso[1]=='cart'):
                    mol.zmat2cart()
                if(ssi[1]=='cart' and sso[1]=='zmat'):
                    mol.cart2zmat()
                mol.print(fmt=sso[0],type=sso[1])
            elif task =='vmd':
                ssi = fmti.split('-')
                mol = mol_class(file=tar, fmt=ssi[0], type=ssi[1])
                if(ssi[1]=='zmat'):
                    mol.zmat2cart()
                mol.save(file='tmp.xyz', fmt='xyz', type='cart')
                os.system('vmd tmp.xyz')
            elif task =='gv':
                ssi = fmti.split('-')
                mol = mol_class(file=tar, fmt=ssi[0], type=ssi[1])
                mol.save(file='tmp.com', fmt='gauss', type=ssi[1])
                os.system('/home/public/opt/Gaussian/gv/gview.csh tmp.com')
            elif task =='mix' and fmti=='mndo-zmat': # NOTE: don't mix cartesian inputs
                ssi = fmti.split('-')
                sso = fmto.split('-')
                mol = mol_class(file=tar, fmt=ssi[0], type=ssi[1])
                
                if i+1 < len(tarlist) :
                    tar2 = tarlist[i+1]
                else:
                    break
                mol2 = mol_class(file=tar2, fmt=ssi[0], type=ssi[1])
                for i in np.linspace(0,1,51):
                    molx = mol.mix(mol2,per=i)
                    molx.print(fmt=sso[0], type=sso[1])
                    # molx.save(file=mol2.name+'.mix%.2f'%i,fmt=sso[0], type=sso[1])
            elif task =='ref' and fmti=='mndo-zmat' and taskfn!='':
                ssi = fmti.split('-')
                sso = fmto.split('-')
                mol = mol_class(file=tar, fmt=ssi[0], type=ssi[1])
                mol2 = mol_class(file=taskfn, fmt=ssi[0], type=ssi[1])
                mol.ref_zmat(mol2)
                mol.print(fmt=sso[0], type=sso[1])
            elif task =='tpl' and taskfn!='':
                ssi = fmti.split('-')
                mol = mol_class(file=tar, fmt=ssi[0], type=ssi[1])
                molx = mol_class(file=taskfn, fmt=ssi[0], type=ssi[1])
                mol._before = molx._before
                mol._after = molx._after
                mol.print(fmt=ssi[0], type=ssi[1])
        
                    
                    
                    
                    
                    
