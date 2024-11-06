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

# import statements of module from python standard library
import os
import shutil
import gzip
import shelve

# imports of local modules
from kids_log import Timing, Log
import amberDriver
# import turbo
# import molcas
# import molpro



def prepareCRG(geometry, ks_config):
    CRG_real, CRG_model_H = None, None
    if ks_config.args.mmsolver == 'amber':
        CRG_real = amberDriver.prepare(geometry, ks_config)
        CRG_model_H = amberDriver.read_crgA()
    elif ks_config.args.mmsolver == 'gromacs':
        pass
    elif ks_config.args.mmsolver == 'openmm':
        pass
    else:
        fatalError("unsupported MM software.")
    Log.writeLog('\n')
    return CRG_real, CRG_model_H

def MM(geometry, ks_config):
    """ wrapper function for the MM part """

    # initialize MMresults with a series of dummy arguments
    MMresults = [], 0.0, 0.0, [], [], 0.0

    # when the scheme actually contains a MM part
    if "M" in geometry.calculationType or "L" in geometry.calculationType:
        Log.writeLog('Starting MM calculation ')

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
