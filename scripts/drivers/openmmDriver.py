#!/usr/bin/env python3
#   Coding=utf-8

#   KIDS SCRIPTS
#   Author: xshinhe
#   
#   Copyright (c) 2024 PeKing Univ. - GNUv3 License

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


## importing external functions and modules
import numpy as np
import os
import shutil
import shlex
import subprocess
import math
import time

from openmm.app import * 
from openmm import *
from openmm.unit import *

import kids_env
from kids_log import Log

def modified_LJ(system):
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    lorentz = CustomNonbondedForce(
        '4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)')
    lorentz.setNonbondedMethod(nonbonded_force.getNonbondedMethod())
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    system.addForce(lorentz)
    LJset = {}
    charges_set = []
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        charges_set.append(charge.value_in_unit(elementary_charge))
        LJset[index] = (sigma, epsilon)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(
            index, charge, sigma, epsilon * 0)
    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        # ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED
        # FORCE
        lorentz.addExclusion(p1, p2)
        if eps._value != 0.0:
            #print p1,p2,sig,eps
            sig14 = sqrt(LJset[p1][0] * LJset[p2][0])
            eps14 = sqrt(LJset[p1][1] * LJset[p2][1])
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
    return system


def prepare(geometry, ks_config):
    """ prepare the input files for MM openmm calculations"""
    pdb_real, pdb_modelH = ks_config.args.pdbfile.split(',')
    FORCEFIELD_PATH_LIST = ks_config.args.xmlff.split(',')

    for f0 in FORCEFIELD_PATH_LIST:
        f0new = f0.replace('.xml', '-noc.xml')
        lines = open(f0, 'r', encoding='utf8').readlines()
        fout = open(f0new, 'w')
        for line in lines:
            if 'charge' in line:
                terms = line.split()
                terms[2] = 'charge="0.000000"'
                nl = '  '.join(terms)
                fout.write(nl+'\n')
            else:
                fout.write(line)
        fout.close()

    real_topology = PDBFile(pdb_real).getTopology()
    modelH_topology = PDBFile(pdb_modelH).getTopology()
    forcefield_real = ForceField(*FORCEFIELD_PATH_LIST)
    system_real = forcefield_real.createSystem(real_topology, 
        nonbondedMethod=NoCutoff, removeCMMotion=False)
    system_modelH = forcefield_real.createSystem(modelH_topology, 
        nonbondedMethod=NoCutoff, removeCMMotion=False)

    CRG_real = []
    for force in system_real.getForces():
        if isinstance(force, openmm.NonbondedForce):
            for i in range(force.getNumParticles()):
                charge, sigma, epsilon = force.getParticleParameters(i)
                charge_value = charge.value_in_unit(openmm.unit.elementary_charge)
                if i < len(CRG_real):
                    CRG_real[i] = charge_value  
                else:
                    CRG_real.append(charge_value)

    CRG_model_H = []
    for force in system_modelH.getForces():
        if isinstance(force, openmm.NonbondedForce):
            for i in range(force.getNumParticles()):
                charge, sigma, epsilon = force.getParticleParameters(i)
                charge_value = charge.value_in_unit(openmm.unit.elementary_charge)
                if i < len(CRG_model_H):
                    CRG_model_H[i] = charge_value  
                else:
                    CRG_model_H.append(charge_value) 

    return CRG_real, CRG_model_H

def MMcalculations(geometry, ks_config):

    pdb_real, pdb_modelH = ks_config.args.pdbfile.split(',')
    FORCEFIELD_PATH_LIST = ks_config.args.xmlff.split(',')

    real_topology = PDBFile(pdb_real).getTopology()
    modelH_topology = PDBFile(pdb_modelH).getTopology()
    forcefield_real = ForceField(*FORCEFIELD_PATH_LIST)
    system_real = forcefield_real.createSystem(real_topology, 
        nonbondedMethod=NoCutoff, removeCMMotion=False)
    if ks_config.get_nested('MM.use_LJ2', False):
        system_real = modified_LJ(system_real)

    context_real = Context(system_real, CustomIntegrator(0))
    context_real.setPositions(Quantity( np.array(geometry.cartesian).T, unit=openmm.unit.angstrom))
    state  = context_real.getState(getEnergy=True, getForces=True)
    E_real = (state.getPotentialEnergy() / openmm.unit.AVOGADRO_CONSTANT_NA ).value_in_unit(openmm.unit.hartree)
    Fxyz_real = (state.getForces(asNumpy=True) / openmm.unit.AVOGADRO_CONSTANT_NA
            ).value_in_unit(openmm.unit.hartree/openmm.unit.bohr).T

    FORCEFIELD_NOC_PATH_LIST = FORCEFIELD_PATH_LIST
    for i in range(len(FORCEFIELD_NOC_PATH_LIST)):
        if FORCEFIELD_NOC_PATH_LIST[i] == 'model-H.xml':
            FORCEFIELD_NOC_PATH_LIST[i] = 'model-H-noc.xml'

    forcefield_real_modelnoc = ForceField(*FORCEFIELD_NOC_PATH_LIST)
    system_real_modelnoc = forcefield_real_modelnoc.createSystem(real_topology, 
        nonbondedMethod=NoCutoff, removeCMMotion=False)
    if ks_config.get_nested('MM.use_LJ2', False):
        system_real_modelnoc = modified_LJ(system_real_modelnoc)

    context_real_modelnoc = Context(system_real_modelnoc, CustomIntegrator(0))
    context_real_modelnoc.setPositions(Quantity( np.array(geometry.cartesian).T, unit=openmm.unit.angstrom))
    state  = context_real_modelnoc.getState(getEnergy=True, getForces=True)
    E_real_modelnoc = (state.getPotentialEnergy() / openmm.unit.AVOGADRO_CONSTANT_NA ).value_in_unit(openmm.unit.hartree)
    Fxyz_real_modelnoc = (state.getForces(asNumpy=True) / openmm.unit.AVOGADRO_CONSTANT_NA
            ).value_in_unit(openmm.unit.hartree/openmm.unit.bohr).T


    forcefield_modelH = ForceField(*FORCEFIELD_NOC_PATH_LIST)
    system_modelH = forcefield_modelH.createSystem(modelH_topology, 
        nonbondedMethod=NoCutoff, removeCMMotion=False)
    if ks_config.get_nested('MM.use_LJ2', False):
        system_modelH = modified_LJ(system_modelH)

    context_modelH = Context(system_modelH, CustomIntegrator(0))
    context_modelH.setPositions(Quantity( np.array(geometry.modelH).T, unit=openmm.unit.angstrom))
    state  = context_modelH.getState(getEnergy=True, getForces=True)
    E_modelH = (state.getPotentialEnergy() / openmm.unit.AVOGADRO_CONSTANT_NA ).value_in_unit(openmm.unit.hartree)
    Fxyz_modelH = (state.getForces(asNumpy=True) / openmm.unit.AVOGADRO_CONSTANT_NA
            ).value_in_unit(openmm.unit.hartree/openmm.unit.bohr).T

    # print(E_real)
    # print(Fxyz_real)
    # print(E_real_modelnoc)
    # print(Fxyz_real_modelnoc)
    # print(E_modelH)
    # print(Fxyz_modelH)
    # exit(0)
   
    return E_real, Fxyz_real, E_real_modelnoc, Fxyz_real_modelnoc, E_modelH, Fxyz_modelH

def clean():
    """ clean up the run directory from all the files that have been used to run Amber and that
        are no longer needed at the end of the calculation """

    # list of files/directory to remove
    toRemove = ['real-noc.xml', "model-H-noc.xml", "model-H.rst"]
    if False:
        for f in toRemove:
            # if f is an existing file, remove it
            if os.path.isfile(f):
                os.remove(f)
            # if f is an existing directory, clean tree
            elif os.path.isdir(f):
                shutil.rmtree(f)

    # done
    return
