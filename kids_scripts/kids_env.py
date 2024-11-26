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

import sys
import os
import subprocess
import shlex
from typing import Optional, Tuple


PROFILE_FILENAME = ".kids_profile"

KIDS_SCRIPTS_ENV = {}

def which(program: str) -> str:
    """
    @brief Locate a program file in the user's path.

    Locate a executable file and return its full path

    @param[in] program (str)
            The name of the executable file to locate.

    @return (str)
            The full path to the executable if found, otherwise None.
    """
    # Retrieve the PATH environment variable and split it into a list
    paths = os.environ.get("PATH", "").split(os.pathsep)
    
    # Iterate over each directory in the PATH
    for path in paths:
        # Construct the full path to the executable
        exe_file = os.path.join(path, program)
        
        # Check if the file exists and is a file
        if os.path.isfile(exe_file):
            # Check if the file is executable
            if os.access(exe_file, os.X_OK):
                # If the file exists and is executable, return its path
                return exe_file
    
    # If the executable is not found, return None
    return None

# ==============================================================================

def source(profilefile: str):
    """
    @brief Source the content of a file to define some environmental variables,
           it's a simple imitation of the "source" bash command.

    This function reads the content of a file to define some environment 
    variables, mimicking the behavior of the "source" bash command.

    @param[in] profilefile (str)
            The name of the file (with absolute or relative path) that contains
            the definitions of some environment variables and possibly some
            additional bash commands.

    @return None            

    @note This function does not modify the original profilefile.
    """

    # Split the bash source command into a list that can be used with subprocess
    command = shlex.split(f"bash -c 'source {profilefile} && env'")

    # Call the source command with subprocess and get the stdout as a pipe
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)

    # Process the standard output line by line
    for line in proc.stdout:
        # Redefine each environmental variable with the updated values
        (key, _, value) = line.decode().partition("=")
        try:
            if key not in ["BASH_FUNC_module%%", "BASH_FUNC_ml%%", "LS_COLORS"]:
                # Remove any trailing newline and set the environment variable
                os.environ[key] = value.rstrip('\n')
                KIDS_SCRIPTS_ENV[key] = value.rstrip('\n')
        except ValueError:
            # If there is a ValueError (e.g., no '=' in the line), do nothing
            pass

    # Wait for the process to terminate and get the remaining output
    proc.communicate()

# ==============================================================================

def findKIDSProfile() -> Optional[str]:
    """
    @brief find a profile file if it exists in the user's home directory or the 
        current working directory.

    This function checks for the existence of a profile file in the user's home 
    directory or the current working directory. If found, it return the file to 
    set up the environment. If no profile file is found, it is assumed that the 
    environment is already set up, and return None.

    @return (Optional[str])
            The path to the profile file if it was found, otherwise None.

    @note This function assumes that the profile file, if it exists, is properly 
    formatted for sourcing.
    """
    # Initialize the filenameconf variable to None
    filenameconf: Optional[str] = None

    # Define the list of directories to search for the profile file: 
    #       first the current working directory, then the home directory
    for directory in [os.getcwd(), os.path.expanduser('~')]:
        # Construct the full path to the profile file
        profile_path: str = os.path.join(directory, PROFILE_FILENAME)
        # Check if the profile file exists at the constructed path
        if os.path.isfile(profile_path):
            # If the file exists, set filenameconf to the path and break
            filenameconf = profile_path
            break

    # Return the path to the profile file if it was found, otherwise None
    return filenameconf


# ==============================================================================


def checkKIDSProfile(filenameconf: str) -> str:
    """
    @brief Checks if the environment for KIDS SCRIPTS execution is properly 
        defined and if all the necessary executables are available.
        If possible, source the file.

    @return (Tuple[bool, str])
            A tuple where the first element is a boolean indicating if the 
            environment is properly set up, and the second element is a string 
            containing an information/error message.

    @note This function also souce the profile if possible
    """
    logstring: str = ""
    if filenameconf:
        logstring += f"""The file {filenameconf} is present
It will be read to set the paths of third-party software\n"""
        source(filenameconf)
    else:
        logstring += f"""The {PROFILE_FILENAME} is NOT present in your home or working directory
It is assumed that the environment has been already set\n"""
    logstring += '\n'

    return logstring


# ==============================================================================


def checkKIDSEnv() -> Tuple[bool, str]:
    """
    @brief Checks if the environment for KIDS SCRIPTS execution is properly 
    defined and if all the necessary executables are available.

    This function verifies that the environment required for KIDS SCRIPTS is 
    correctly configured and that all required executables are accessible.

    @return (Tuple[bool, str])
            A tuple where the first element is a boolean indicating if the 
            environment is properly set up, and the second element is a string 
            containing an error message if any check fails.

    @note This function does not modify the environment or system state.
    """
    # Check if the environmental variable with the <KIDS SCRIPTS> path is defined
    message: str = ""
    try:
        kids_scripts_path = os.environ['KIDS_SCRIPTS_PATH']
        os.environ['PATH'] += os.pathsep + kids_scripts_path
        message = f'KIDS_SCRIPTS_PATH={kids_scripts_path}, add it to PATH'
    except KeyError:
        message = "Environment variable $KIDS_SCRIPTS_PATH is not defined"
        return False, message

    # Check if <KIDS SCRIPTS> executables are available
    for exe in ["kidsqm.py", "kidsqmmm.py"]:
        if not which(exe):
            message = f"{exe} executable is not available"
            return False, message

    # Return True if all the checks were OK
    return True, message

# ==============================================================================

def setGaussianProfile(gversion:str, gpath:str):
    """ Define enviromental variables for running a given Gaussian version

    The function gets as input a string labelling the version of Gaussian
    to setup, and the root directory that contains that gaussian version.
    The function then defines the environmental variables that are needed
    to run that version of Gaussian, using the correct Gaussian profile file.

    Keyword arguments:
    gversion -- name of the executable of the gaussian version to use (e.g. g17, g09 ...)
    gpath    -- path where the directory of Gaussian version is stored
                ( attention! give only /usr/local when Gaussian09 dir is /usr/local/g09 )
    """

    # create the scratch directory
    scratchDir = os.environ["GAUSS_SCRDIR"]
    if not os.path.exists(scratchDir): os.mkdir(scratchDir)

    # define the variable for the gaussian root directory, e.g. $g16root = /usr/local
    os.environ[gversion + "root"] = gpath

    # source the profile file in the gpath/gversion/bsd directory
    profilefile = os.path.join(gpath, gversion, "bsd", gversion + ".profile")
    source(profilefile)



def checkAmberEnv():
    """ the function checks if the environment for AMBER execution is properly defined and
        if all the necessary executables are available """

    # environmental variable with the AMBER path
    try:
        AMBERHOME = os.environ['AMBERHOME']
    except KeyError:
        message = "environment variable $AMBERHOME is not defined"
        return False, message
        
    command = shlex.split(os.path.join(AMBERHOME, "update_amber") + " -v")

    # call update_amber to check the current version of amber
    #command = shlex.split(os.path.join(AMBERHOME, "update_amber") + " -v")
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    # process the std out line-by-line
    amberVersion = ["8"]
    for line in proc.stdout:
        if b"AmberTools version" in line:
            amberVersion = line.split()[2].split(b".")
    proc.communicate()

    # define the list of executables depending on the version of amber
    if int(amberVersion[0]) <= 12:
        exeList = ["antechamber", "tleap", "sander", "ambpdb", "parmchk", "cpptraj"]
    else:
        exeList = ["antechamber", "tleap", "sander", "ambpdb", "parmchk2", "cpptraj"]

    # check if amber executables are available
    for exe in exeList:
        if not which(exe):
            message = exe + " executable is not available"
            return False, message

    # return True if all the checks were OK
    return True, ""

# ================================================================================

def checkGaussianQMEnv():
    """ the function checks if the environment for GAUSSIAN execution
       for QM calculations is properly defined and if all the
       necessary executables are available """

    # get the name of the executable and the path from the environmental variable
    try:
        gversion = os.environ["GAUSSIAN_EXE_QM"]
        gpath = os.environ["GAUSSIAN_DIR_QM"]
        _ = os.environ["GAUSS_SCRDIR"]
    except KeyError:
        return False, "environment variables for Gaussian optimization, $GAUSSIAN_EXE_QM, " \
                      "$GAUSSIAN_DIR_QM and $GAUSS_SCRDIR, are not defined"

    # check that the gversion that is passed corresponds to a supported version
    if gversion != "g03" and gversion != "g09" and gversion != "g16":
        return False, "{0} version of Gaussian is not currently supported by Cobramm".format(gversion)

    # source the configuration file for this version of gaussian
    try:
        setGaussianProfile(gversion, gpath)
    except:
        return False, "error encountered while sourcing the {0} profile script".format(gversion)

    # check if gaussian executable is available
    if not which(gversion):
        message = gversion + " executable is not available"
        return False, message

    # return True if all the checks were OK
    return True, ""

# ================================================================================

def checkMolcasEnv():
    """ the function checks if the environment for MOLCAS or OPENMOLCAS execution
        is properly defined and if all the necessary executables are available """

    # get the name of the executable and the path from the environmental variable
    try:
        molcasscript = os.environ['MOLCAS_SCRIPT']  # path of the MOLCAS launch script
        molcaspath = os.environ['MOLCAS']  # path of the MOLCAS code
    except KeyError:
        return False, "environment variables for Molcas/Openmolcas, MOLCAS_SCRIPT and MOLCAS, are not defined"

    # molcas script exists and has x permissions
    if not os.path.isfile(molcasscript) or not os.access(molcasscript, os.X_OK):
        return False, 'molcas launch script {} not found, please check the env variable MOLCAS_SCRIPT!'.format(
            molcasscript)
    # molcas directory exists
    if not os.path.isdir(molcaspath):
        return False, 'molcas directory {} not found, please check the env variable MOLCAS!'.format(molcaspath)

    # return True if all the checks were OK
    return True, ""


# ================================================================================

def checkMolproEnv():
    """ the function checks if the environment for Molpro execution
        is properly defined and if all the necessary executables are available """

    # get the name of the executable and the path from the environmental variable
    try:
        molproexe = os.environ['MOLPRO_EXE']  # path of the Molpro executable
    except KeyError:
        return False, "environment variable for Molpro: MOLPRO_EXE is not defined"

    # molcas script exists and has x permissions
    if not os.path.isfile(molproexe) or not os.access(molproexe, os.X_OK):
        return False, 'Molpro executable {} not found, please check the env variable MOLPRO_EXE!'.format(molproexe)

    # return True if all the checks were OK
    return True, ""


# ==============================================================================

def getVersion():
    """
    get last git commit-hash and date of the last commit
    :return: string: "{hash} {date}"
    """
    wordDir = os.getcwd()
    try:
        KS_Path = os.environ["KIDS_SCRIPTS_PATH"]
        os.chdir(KS_Path)
    except KeyError: 
        # try in current dir
        pass

    try:
        commit_hash = subprocess.check_output(["git", "rev-parse", "HEAD"]).strip().decode()
        date = subprocess.check_output(["git", "log", "-1", "--format=%cd", "--date=short"]
            ).strip().decode()
        version = '{} {}'.format(commit_hash, date)
    except:
        version = 'N.A.'
    os.chdir(wordDir)
    return version


# ================================================================================

def getHostname():
    """
    get the hostname of the machine where cobramm is running
    :return: string hostname
    """
    # We use socket instead of reading env variables for better compatibility across systems
    return __import__('socket').gethostname()

# ================================================================================

if __name__ == '__main__':
    from kids_log import Log
    from kids_arg import QMSolverType, MMSolverType

    Log.startSection(f'[TEST] {__file__}')
    Log.writeLog(f'gitversion is {getVersion()}\n')
    Log.writeLog(f'host name is {getHostname()}\n')
    Log.writeLog(f'profile name is {PROFILE_FILENAME}\n')
    Log.writeLog(f'test for <which>: {which("pwd")}\n')

    Log.writeLog(f'find profile is: {findKIDSProfile()}\n')
    Log.writeLog(f'source profile info: {checkKIDSProfile(findKIDSProfile())}\n')
    Log.writeLog(f'check environment info: {checkKIDSEnv()}\n')
    for k,v in KIDS_SCRIPTS_ENV.items():
        is_third = False
        for t in QMSolverType:
            if t.name in k:
                is_third = True
        for t in MMSolverType:
            if t.name in k:
                is_third = True
        if is_third:
            if os.access(v, os.X_OK):
                Log.writeLog(f'Is OK    [ENV INFO] {k}={v}\n')
            else:
                Log.writeLog(f'Invalid  [ENV INFO] {k}={v}\n')


