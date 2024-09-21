#! /usr/bin/env python3

#####################################################################################################

import sys  # System-specific parameters and functions
from math import *  # import mathematical functions

import numpy as np  # import numpy arrays and other stuff
from scipy.stats import boltzmann  # import boltzmann statistic
import constants  # physical and mathematical constants


#####################################################################################################

class HarmonicSampling:
    # define random sampling of the configurations of a set of oscillators
    # in the Harmonic approximation

    # ================================================================================

    def __init__(self, geom, equil_geom, masses, hessian, cart2int_transform, temperature, frozenmodes):
        # define the physics input parameters for the distribution:
        # equilibrium geometry, masses of the coordinates,
        # Cartesian hessian matrix at the equilibrium, temperature

        # store the equilibrium geometry, the masses, the temperature
        # coordinates in input are in Ang, inside the class bohr is needed
        self.equilGeom = np.copy(equil_geom) / constants.Bohr2Ang
        self.masses = np.copy(masses)
        # temperature in input is in K, inside this class we need the au value
        self.temperature = temperature * constants.K2AU
        self.nCoords = len(equil_geom)
        self.nAtoms = geom.NatomHM 
        self.HMlabels = geom.getAtomLabels("MEDIUM_HIGH")
        self.HMlabelIndex = []
        for i in self.HMlabels:
            self.HMlabelIndex.append(geom.PTE[i])
        # self.nCoords is defined as the size of self.equilGeom, but masses (self.nCoords) and
        # hessian (self.nCoords x self.nCoords) should be defined accordingly

        # transform the Cartesian hessian matrix in mass-scaled coordinates
        self.massScaledHessian = hessian
        for i in range(self.nCoords):
            for j in range(self.nCoords):
                self.massScaledHessian[i, j] = self.massScaledHessian[i, j] / sqrt(masses[i] * masses[j])

        # diagonalize mass scaled hessian matrix w/o removing translation & rotation DOF
        self.eigenVal, self.eigenVec = np.linalg.eig(self.massScaledHessian)
        self.eigenVec = np.transpose(self.eigenVec)
        # sort in increasing order of frequencies
        self.eigenVal, self.eigenVec = (list(i) for i in zip(*sorted(zip(self.eigenVal, self.eigenVec))))
        print("\nList of {0} frequencies w/o removing translational and rotational DOF:".format(len(self.eigenVal)))
        for i in range(len(self.eigenVal)):    
            if i > 0 and i%10 == 0:
                print()
            if self.eigenVal[i] < 0.:
                print("{0:8.2f}".format(-1.0*np.sqrt(np.abs(self.eigenVal[i]))/constants.wavnr2au), end=" ")
            else:
                print("{0:8.2f}".format(np.sqrt(self.eigenVal[i])/constants.wavnr2au), end=" ") 

        # transform Cartesian hessian matrix to internal coordinates 
        self.massScaledHessian_int = np.linalg.multi_dot([cart2int_transform,self.massScaledHessian,np.transpose(cart2int_transform)])

        # diagonalize mass scaled hessian matrix after removing translation & rotation DOF
        self.eigenVal, self.eigenVec = np.linalg.eig(self.massScaledHessian_int[6:self.nCoords+1,6:self.nCoords+1])

        # build diagonal matrix with 1 over sqrt of the mass
        massmat = np.diag(1.0/np.sqrt(masses))
        # extend the eigenvector matrix to 3N x 3N size by adding unit vetors (for translation and rotation)
        self.eigenVec = np.insert(self.eigenVec,0,np.array([[0],[0],[0],[0],[0],[0]]),axis=1)
        self.eigenVec = np.insert(self.eigenVec,0,np.zeros((6,self.nCoords)),axis=0)
        for i in range(6):
            self.eigenVec[i][i] = 1.0
        # transform eigenvectors back to Cartesian coordinates    
        self.eigenVecCart = np.linalg.multi_dot([massmat,np.transpose(cart2int_transform),self.eigenVec])

        # normalize
        for i in range(self.nCoords):
            norm = np.linalg.norm(self.eigenVecCart[i])
            self.eigenVecCart[i] = self.eigenVecCart[i]/norm

        # transpose
        self.eigenVecCart = np.transpose(self.eigenVecCart)

        # sort in increasing order of frequencies, first six vectors are translation and rotation
        self.eigenVecCart = self.eigenVecCart[6:]
        self.eigenVal, self.eigenVecCart = (list(i) for i in zip(*sorted(zip(self.eigenVal, self.eigenVecCart))))

        print("\n\nList of {0} frequencies after removing translational and rotational DOF: ".format(len(self.eigenVal)))
        for i in range(len(self.eigenVal)):
            if i > 0 and i%10 == 0:
                print()
            if self.eigenVal[i] < 0.:
                print("{0:8.2f}".format(-1.0*np.sqrt(np.abs(self.eigenVal[i]))/constants.wavnr2au), end=" ")
            else:
                print("{0:8.2f}".format(np.sqrt(self.eigenVal[i])/constants.wavnr2au),end=" ")
        print()

        # generate a fake Gaussian file which holds the Cartesian coordinates of the normal mode and can be visualized with Molden
        self.generate_fake_gaussian()

        # define list of frozen modes to be excluded from the sampling
        # self.frozen contains only modes > 0. (imaginary modes are excluded anyway)
        self.frozen = []
        if frozenmodes == -1:
            pass
        else:
            self.frozen.extend(frozenmodes)
            tmp = []
            for i in self.frozen:
                if self.eigenVal[i-1] > 0.:
                    tmp.append(i)
            self.frozen = tmp

        print("\nList of {0} frequencies used in the Wigner sampling: ".format(sum(1 for i in self.eigenVal if i > 0.)-len(self.frozen)))
        count = -1
        for i in range(len(self.eigenVal)):
            if self.eigenVal[i] > 0. and i+1 not in self.frozen:
                count += 1
                if count > 0 and count%10 == 0:
                    print()
                print("{0:8.2f}".format(np.sqrt(self.eigenVal[i])/constants.wavnr2au),end=" ")
        print()

    # ================================================================================

    def nm_to_cartesian(self, nm_coord):

        # transform from normal modes to mass weighted cartesian coordinates (referenced to the equilibrium geometry
        cart_coord = np.dot(np.transpose(nm_coord), self.eigenVecCart)
        #cart_coord = cart_coord.reshape(len(nm_coord))
      
        # scale by the square root of the masses
        cart_coord = np.array([x / sqrt(m) for x, m in zip(cart_coord[0], self.masses)])
        return cart_coord

    # ================================================================================
    
    def nm_to_cartesianV(self, nm_veloc):

        # transform from normal mode velocities to mass weighted cartesian velocities
        cart_veloc = np.dot(np.transpose(nm_veloc), self.eigenVecCart)
        #cart_veloc = cart_veloc.reshape(len(nm_veloc))
       
        # scale by the square root of the masses
        cart_veloc = np.array([x / sqrt(m) for x, m in zip(cart_veloc[0], self.masses)])
        return cart_veloc

    # ================================================================================

    def get_sample(self, distribution="wigner"):

        # cycle over normal modes and construct the phase-space sample in mass-weighted normal modes
        #nm_coord = [0,0,0,0,0,0]
        #nm_veloc = [0,0,0,0,0,0]
        nm_coord = []
        nm_veloc = []
        for w in self.eigenVal:
            #if w > 0. and w > (400. * constants.wavnr2au) ** 2:
            if w > 0. and self.eigenVal.index(w)+1 not in self.frozen:
                if distribution == "wigner":
                    p, q = self._thermal_wigner_sampling(sqrt(w), self.temperature, 1)
                elif distribution == "actionangle":
                    p, q = self._thermal_action_angle_sampling(sqrt(w), self.temperature, 1)
                else:
                    print("fatal error: sampling method ({0}) not implemented!".format(distribution))
                    sys.exit()
            else:
                q = [0.0]
                p = [0.0]
            nm_coord.append(q)
            nm_veloc.append(p)
        
        #compute potential energy for each displacement and the total potential energy
        ene_per_mode=[]
        for i in range(len(self.eigenVal)):
            if i+1 not in self.frozen and self.eigenVal[i] > 0.:
                #cart_coord=nm_coord[i]*self.eigenVecCart[:,i]
                cart_coord=nm_coord[i]*self.eigenVecCart[i][:]
                cart_coord=np.dot(cart_coord, cart_coord)
                cart_coord=0.5*cart_coord*self.eigenVal[i]
                ene_per_mode.append(cart_coord*constants.Hartree2eV)
        ePot=0.0
        for i in range(len(ene_per_mode)): 
            ePot=ePot+ene_per_mode[i]
        print("Total ePot = {0:.2f} eV".format(ePot))

        ene_per_mode=[]
        for i in range(len(self.eigenVal)):
            if i+1 not in self.frozen and self.eigenVal[i] > 0.:
                #cart_veloc=nm_veloc[i]*self.eigenVecCart[:,i]
                cart_veloc=nm_veloc[i]*self.eigenVecCart[i][:]
                cart_veloc=np.dot(cart_veloc, cart_veloc)
                cart_veloc=0.5*cart_veloc
                ene_per_mode.append(cart_veloc*constants.Hartree2eV)
        eKin=0.0
        for i in range(len(ene_per_mode)):
            eKin=eKin+ene_per_mode[i]
        print("Total eKin = {0:.2f} eV".format(eKin))

        # transform to cartesian coordinates, add the equilibrium geometry and return the results in Angstrom
        return ( (self.equilGeom + self.nm_to_cartesian(nm_coord))*constants.Bohr2Ang, self.nm_to_cartesianV(nm_veloc))

    # ================================================================================

    @staticmethod
    def _thermal_wigner_sampling(omega, temp, nsamples):
        # private function that draws coordinates and momenta according to a
        # thermal Wigner distribution

        # draw two random numbers with Gaussian distribution
        ptilde, qtilde = np.random.normal(size=(2, nsamples))

        # scale them with proper factors
        alpha = tanh(omega / 2.0 / temp)
        p = sqrt(omega / (2.0 * alpha)) * ptilde
        q = sqrt(1.0 / (2.0 * alpha * omega)) * qtilde

        return p, q

    # ================================================================================

    @staticmethod
    def _thermal_action_angle_sampling(omega, temp, nsamples):
        # private function that draws coordinates and momenta by sampling
        # the classical orbits corresponding to the harmonic oscillator eigenvalues

        # sample vibrational quantum numbers over boltzmann distribution
        n_max = - temp / omega * log(0.001)
        lamb = omega / temp
        nstates = boltzmann.rvs(lamb, n_max, size=nsamples)

        # now draw one random angle per each sample
        angles = np.random.uniform(0.0, 2.0 * pi, size=nsamples)

        # combine state and angle to compute random p and q
        ptilde, qtilde = [], []
        for n, a in zip(nstates, angles):
            ptilde.append(sin(a) * sqrt(2.0 * n + 1.0))
            qtilde.append(cos(a) * sqrt(2.0 * n + 1.0))

        # scale them with proper factors
        p = sqrt(omega) * np.array(ptilde)
        q = sqrt(1.0 / omega) * np.array(qtilde)

        return p, q

    # ================================================================================

    def generate_fake_gaussian(self):
        # generate a fake Gaussian file with Cartesian coordinates of the normal modes
        # required as Gaussian internal projection of rotations and translations gives a somewhat different result
        # and normal mode shape and ordering might be different 
        # to be visualized with Molden

        with open("gaussian.log", 'w') as gauss:
            gauss.write("""the Gaussian(R) 03 system (copyright 2003, Gaussian, Inc.),
                          Input orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------{0}""".format("\n"))
            for i in range(self.nAtoms):
                gauss.write('    {0:3d}        {1:3d}           0      {2:10.6f}  {3:10.6f}  {4:10.6f}\n'.format(i+1, self.HMlabelIndex[i], self.equilGeom[3*i]*constants.Bohr2Ang, self.equilGeom[3*i+1]*constants.Bohr2Ang, self.equilGeom[3*i+2]*constants.Bohr2Ang))
            gauss.write(""" ---------------------------------------------------------------------
 Symmetry turned off by external request.
     1 basis functions,     1 primitive gaussians,     1 cartesian basis functions
     1 alpha electrons        1 beta electrons
            Population analysis using the SCF density.
 **********************************************************************

 and normal coordinates:{0}""".format("\n"))    
            if (len(self.eigenVal) % 3 != 0):
                NMblocks = int(len(self.eigenVal) / 3) + 1
            else:
                NMblocks = int(len(self.eigenVal) / 3)
            for i in range(NMblocks):
                gauss.write('                    {0:2d}                     {1:2d}                     {2:2d}\n'.format(3*i+1, 3*i+2, 3*i+3))
                gauss.write('                     A                      A                      A\n')
                f1, f2, f3 = 1, 1, 1
                if self.eigenVal[3*i] < 0.0:
                    f1 = -1.
                if self.eigenVal[3*i+1] < 0.0:
                    f2 = -1.
                if self.eigenVal[3*i+2] < 0.0:
                    f3 = -1.
                gauss.write(' Frequencies --  {0:9.4f}              {1:9.4f}              {2:9.4f}\n'.format(f1*np.sqrt(np.abs(self.eigenVal[3*i]))/constants.wavnr2au, f2*np.sqrt(np.abs(self.eigenVal[3*i+1]))/constants.wavnr2au, f3*np.sqrt(np.abs(self.eigenVal[3*i+2]))/constants.wavnr2au))
                gauss.write(' Red. masses --   {0:8.4f}              {1:8.4f}              {2:8.4f}\n'.format(1, 1, 1))
                gauss.write(' Frc consts  --   {0:8.4f}              {1:8.4f}              {2:8.4f}\n'.format(1, 1, 1))
                gauss.write(' IR Inten    --   {0:8.4f}              {1:8.4f}              {2:8.4f}\n'.format(0, 0, 0))
                gauss.write(' IR Rel Inten--   {0:8.4f}              {1:8.4f}              {2:8.4f}\n'.format(0, 0, 0))
                gauss.write('  Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z\n')
                for j in range(self.nAtoms):
                    gauss.write('    {0:2d}   {1:2d}   {2:5.2f}  {3:5.2f}  {4:5.2f}    {5:5.2f}  {6:5.2f}  {7:5.2f}    {8:5.2f}  {9:5.2f}  {10:5.2f}\n'.format(j+1, self.HMlabelIndex[j], self.eigenVecCart[3*i][3*j], self.eigenVecCart[3*i][3*j+1], self.eigenVecCart[3*i][3*j+2], self.eigenVecCart[3*i+1][3*j], self.eigenVecCart[3*i+1][3*j+1], self.eigenVecCart[3*i+1][3*j+2], self.eigenVecCart[3*i+2][3*j], self.eigenVecCart[3*i+2][3*j+1], self.eigenVecCart[3*i+2][3*j+2] ))


