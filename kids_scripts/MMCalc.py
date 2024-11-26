#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# KIDS SCRIPTS (adapted from COMBRAMM)
# Author: xshinhe
#
# Copyright (c) 2024 Peking Univ. - GNUv3 License
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

# import statements of module from python standard library
import os
import shutil
import gzip
import shelve

import numpy as np

# imports of local modules
from kids_log import Timing, Log
from Charge import Charge
from drivers import amberDriver
from drivers import openmmDriver
from drivers import gromacsDriver

# import turbo
# import molcas
# import molpro


def prepareCRG(geometry, ks_config):
    if ks_config.args.mmsolver == 'amber':
        return amberDriver.prepare(geometry, ks_config)
    elif ks_config.args.mmsolver == 'gromacs':
        pass
    elif ks_config.args.mmsolver == 'openmm':
        pass
    else:
        fatalError("unsupported MM software.")
    return None, None

def MM(geometry, ks_config):
    """ wrapper function for the MM part """

    # initialize MMresults with a series of dummy arguments
    MMresults = 0.0, [], 0.0, [], 0.0, []

    # when the scheme actually contains a MM part
    if "M" in geometry.calculationType or "L" in geometry.calculationType:
        Log.writeLog('Starting MM calculation ')

        # amber optimization of low layer
        geometry.updatereal('real.rst')

        # AMBER calculation (only option implemented now)
        if ks_config.args.mmsolver == 'amber':
            Log.writeLog('using AMBER ... \n')
            MMresults = amberDriver.MMcalculations(geometry, ks_config)
        elif ks_config.args.mmsolver == 'gromacs':
            MMresults = gromacsDriver.MMcalculations(geometry, ks_config)
        elif ks_config.args.mmsolver == 'openmm':
            MMresults = openmmDriver.MMcalculations(geometry, ks_config)
        else:
            fatalError("unsupported MM software.")
        Log.writeLog('\n')

    return MMresults

class MM2:
    def __init__(self, geometry, ks_config):
        """ wrapper function for the MM part """

        self.geometry = geometry
        self.config   = ks_config
        self.mmsolver = ks_config.args.mmsolver

        if self.mmsolver == 'amber':
            self.charge_real, self.charge_modelH = amberDriver.prepare(geometry, ks_config)
        elif self.mmsolver == 'gromacs':
            self.charge_real, self.charge_modelH = gromacsDriver.prepare(geometry, ks_config)
        elif self.mmsolver == 'openmm':
            self.charge_real, self.charge_modelH = openmmDriver.prepare(geometry, ks_config)
        else:
            fatalError("unsupported MM software.")            

        # update with QM's charges if possible
        if os.path.exists('laststep.charge'):
            Log.writeLog('update with laststep charge')
            lines = open('laststep.charge', 'r').readlines()
            self.charge_modelH = [float(l) for l in lines]
            shutil.move('laststep.charge', 'laststep.charge.old')

        # update geometry if possible
        if os.path.exists('real.rst'):
            Log.writeLog('update geometry')
            geometry.updatereal('real.rst')

        # build Charge class
        self.charges = Charge(geometry, self.charge_real, self.charge_modelH)
        self.charges.checkConsistency()
                
    def runMM(self):
        MMresults = 0.0, [], 0.0, [], 0.0, []
        if self.mmsolver == 'amber':
            MMresults = amberDriver.MMcalculations(self.geometry, self.config)
        elif self.mmsolver == 'gromacs':
            MMresults = gromacsDriver.MMcalculations(self.geometry, self.config)
        elif self.mmsolver == 'openmm':
            MMresults = openmmDriver.MMcalculations(self.geometry, self.config)
        else:
            fatalError("unsupported MM software.")

        self.energy_real = MMresults[0]
        self.grad_real   =-np.array(MMresults[1])
        self.energy_real_modelnoc = MMresults[2]
        self.grad_real_modelnoc   =-np.array(MMresults[3])
        self.energy_modelH = MMresults[4]
        self.grad_modelH   =-np.array(MMresults[5])
        Log.writeLog('\n')



