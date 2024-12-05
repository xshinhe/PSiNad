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

#####################################################################################################

# import statements of module from python standard library

import argparse    # command line arguments parser
import sys         # System-specific parameters and functions
import numpy as np
import os          # filesystem utilities
import shutil
import copy

# other modules

# hijack PYTHON path (sys.path) to find cobramm local modules
try:
    sys.path.insert(0, os.path.join(os.getenv('COBRAM_PATH'),'cobramm'))
except KeyError:
    raise RuntimeError('cannot import cobramm module.\n'
                       'Please check if $COBRAM_PATH is defined in your environment, if not run:\n'
                       'export COBRAM_PATH=/path/to/COBRAMM')

# imports of local modules

import logwrt  # write messages and output to log (std out or cobramm.log)
import cobrammenv  # management of the environment for running COBRAMM (and utils)

# imports of user-defined classes

from cobrammCalculator import CobrammCalculator  # Cobramm wrapper classes
from layers import Layers  # object to store information on the geometry and the layers definition
from harmonicSampling import HarmonicSampling  # class to handle Wigner sampling
from cobrammCalculator import CobrammOutput  # class to parse the output of a previously run COBRAMM calculation

#####################################################################################################

# SCRIPT PARAMETERS

_REALLAYERS_FNAME = "real_layers.xyz"
_REALCRD_FNAME = "real.crd"
_HMXYZ_FNAME = "geometry.xyz"
_VELOCITY_FNAME = "velocity.dat"
_REALTOP_FNAME = "real.top"
_MODELHTOP_FNAME = "model-H.top"
_COMMAND_FNAME = "cobram.command"

_SCANVALUE_FILE = "SCAN_VALUE"


######################################################################################################################

def main():

    # this new class defines a style for the parser message that join the two predefined styles
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    # Define command line arguments with argparse
    # create the top-level parser
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="COBRAMM utility script: wigner sampling with COBRAMM",
        formatter_class=CustomFormatter)

    parser.add_argument('--version', action='version',
                        version='%(prog)s is part of COBRAMM {}'.format(cobrammenv.getVersion()))

    # define subparsers depending for each of the script actions
    subparsers = parser.add_subparsers(title='valid actions',
                                       description='choose one of the steps of the scanning procedure')

    # create the parser for the "prepare" command
    parser_prepare = subparsers.add_parser('prepare', help='prepare real.crd input files for the Wigner snapshots',
                                           formatter_class=CustomFormatter)
    parser_prepare.add_argument("-n", "--snap-number", required=True, metavar="NSAMPLES", dest="nsamples", type=int,
                                help="number of required Wigner sampling snapshots")
    parser_prepare.add_argument("-fd", "--frequency-dir", required=True, metavar='FREQUENCYPATH', dest="freqdir",
                                help="path where the COBRAMM frequency calculation is stored")
    parser_prepare.add_argument("-T", "--temperature", default=300.0, metavar='TEMP', dest="temperature", type=float,
                                help="target temperature of the Wigner sampling in K")
    parser_prepare.add_argument("-FM", "--frozen-modes", default=-1, metavar='FROZENMODES', dest="frozenmod", type=int, nargs='*',
                                help="list of modes to exclude from the sampling")
    parser_prepare.set_defaults(func=prepare)

    args = parser.parse_args()

    # call the specific action requested by command argument
    try:
        args.func(args)
    except AttributeError:
        parser.print_help()


######################################################################################################################

def calcCOM(coords, masses, mass):
    '''Calculates center of mas
       coords : array containing one [x, y, z] array of coords for each atom
       masses : array of atomic masses
       mass : atomic mass'''
    return np.einsum("ij->j", np.einsum("ij,i->ij", coords, masses/mass))
    
def Gram_Schmidt(A):
    '''Gram Schmidt orthonormalization of the column vectors of matrix A'''
    Q, R = np.linalg.qr(A)
    indexes = np.where(np.any(np.abs(R) > 1e-12, axis=1))[0]
    return Q[:,indexes]

######################################################################################################################

def prepare(args):

    # source COBRAMM configuration file (if available)
    conf_file = cobrammenv.setCobrammProfile()
    # Print to log information on the environment configuration
    logwrt.writelog(cobrammenv.checkConfigFile(conf_file))

    # initialize COBRAMM wrapper classes
    CobrammCalculator()

    # ================================================================================================================

    # print header
    logwrt.writelog("###########################################\n")
    logwrt.writelog("###   COBRAMM WIGNER SAMPLING UTILITY   ###\n")
    logwrt.writelog("###########################################\n\n")

    # read the frequency calculation to define equilibrium geometry, normal modes vectors and frequencies
    reallayers = os.path.join(args.freqdir, _REALLAYERS_FNAME)
    coordfile = os.path.join(args.freqdir, _REALCRD_FNAME)

    # first define a Layers instance with the data contained in the real_layers.xyz file
    logwrt.writelog("Reading layers from COBRAMM file '{}'\n".format(reallayers))
    with open(reallayers, "r") as f:
        geometry = Layers.from_real_layers_xyz(f.read())

    # depending on the availability of the real.crd file, update coordinates
    if os.path.isfile(coordfile):
        geometry.updatereal(filename=coordfile)
        logwrt.writelog("Reading equilibrium coordinates from AMBER inpcrd file '{}'\n".format(coordfile))
    else:
        logwrt.writelog("Reading equilibrium coordinates from COBRAMM file '{}'\n".format(reallayers))

    # parse the results of the frequency calculation run using CobrammOutput
    freqresults = CobrammOutput(args.freqdir, store_files=True)
    
    # ================================================================================================================
    # prepare transformational matrices for projecting transational and rotational degrees of freedom 

    # extract equilibrium geometry and format as a simple 1D vector of 3N elements
    geomvector_1D = []
    #for x, y, z in zip(*geometry.model):
    for x, y, z in zip(*geometry.MEDIUM_HIGH):
        geomvector_1D.append(x), geomvector_1D.append(y), geomvector_1D.append(z)
    geom_xyz = np.reshape(geomvector_1D, (geometry.NatomHM, 3))

    # calculate center of mass
    massvec = np.array(freqresults.coord_masses[::3])
    mass = sum(freqresults.coord_masses)/3
    com = calcCOM(geom_xyz, massvec, mass)
    
    # shift coordinates to center of mass (3xN array)
    geomvector_com = []
    geomvector_com.append(np.array(geomvector_1D)[0::3] - com[0])
    geomvector_com.append(np.array(geomvector_1D)[1::3] - com[1])
    geomvector_com.append(np.array(geomvector_1D)[2::3] - com[2])
    geomvector_com_1D = []
    for x, y, z in zip(*geomvector_com):
        geomvector_com_1D.append(x), geomvector_com_1D.append(y), geomvector_com_1D.append(z)

    # calculate moment of inertia
    moimat = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            for k in range (geometry.NatomHM):
                pos = geomvector_com_1D[3*k:3*k+3]
                if i==j:
                    for l in range(3):
                        if l != i:
                            moimat[i,j] += massvec[k] *  (pos[l]**2)
                else:
                    moimat[i,j] -= massvec[k] * pos[i] * pos[j]

    # diagonalize moment of inertia matrix
    eigenVal, eigenVec = np.linalg.eig(moimat)
    eigenVec = np.transpose(eigenVec)

    # D matrix
    # initialize a matrix of random floats
    # columns 0,1,2 are set as molecular rigid translations in mass-weighted coords
    # columns 3,4,5 are set as molecular rigid rotations in mass-weighted coords
    D = np.random.random((3 * geometry.NatomHM, 3 * geometry.NatomHM))
    for ia in range(geometry.NatomHM):
        i = 3*ia # row index of D corresponding to atom ia
        # pos = [x, y, z] of atom ia
        pos = np.array([geomvector_com_1D[i], geomvector_com_1D[i+1], geomvector_com_1D[i+2]])
        D[i:i+3,:3] = np.identity(3) * np.sqrt(massvec[ia])
        D[i:i+3,3:6] = np.cross(eigenVec.transpose(), pos, axisa=0).transpose() * np.sqrt(massvec[ia])
    
    # obtain cart2intmat to transform cartesian normal modes to internal coordinates
    cart2intmat = Gram_Schmidt(D).transpose()
    
    # repeat in case of bad random D initialization
    while np.any(np.shape(cart2intmat) != np.shape(D)):
        cart2intmat = Gram_Schmidt(D).transpose()

    # ================================================================================================================
    print("\nList of user defined frozen modes", args.frozenmod) 
    # define harmonic sampling of the molecular oscillations
    mol_oscillator = HarmonicSampling(geometry, geomvector_1D, freqresults.coord_masses, freqresults.force_matrix, cart2intmat, args.temperature, args.frozenmod)
 
    # define directory where to store displaced geometries
    displ_dir = "sampling"
    # check the existence of the directory where to store samples, in case create the directory
    if not os.path.isdir(displ_dir):
        os.mkdir(displ_dir)

    # now randomly generate snapshots for subsequent calculations
    displ_geometries = []
    velocities = []
    for i_geom in range(args.nsamples):
        # get sample of the wigner distribution
        (newgeomvector, newvelvector) = mol_oscillator.get_sample()
        # create a new Layers object to store the new snapshot, move the H layer and store the snapshot
        new_geometry = Layers.from_real_layers_xyz(geometry.reallayertext)
        new_geometry.updateHMlayers([newgeomvector[0::3], newgeomvector[1::3], newgeomvector[2::3]])
        displ_geometries.append(new_geometry)
        velocities.append([newvelvector[0::3], newvelvector[1::3], newvelvector[2::3]])
        
    # now create the crd file for the displaced geometry in one of the sampling/sample_XXX subdirectories
    for i, new_geometry in enumerate(displ_geometries):
        # define name of the directory where to store the sample point
        path_snap = os.path.join(displ_dir, "sample_{0:04d}".format(i))
        os.mkdir(path_snap)
        # create the crd file
        new_geometry.makerealcrd(os.path.join(path_snap, _REALCRD_FNAME, ), os.path.join(path_snap,_HMXYZ_FNAME, ))
        with open(os.path.join(path_snap, _VELOCITY_FNAME), 'w') as vel:
            for j in range(len(velocities[i][0])):
                vel.write('{0:12.7f}{1:12.7f}{2:12.7f}{3}'.format(
                    velocities[i][0][j], velocities[i][1][j], velocities[i][2][j], '\n'))


######################################################################################################################

if __name__ == '__main__':
    main()
