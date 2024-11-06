#!/usr/bin/env python3
#   Coding=utf-8

#   KIDS SCRIPTS (adapted from COMBRAMM)
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

import sys  # system commands
import os  # operating system utilities
import subprocess  # run external program as child process
import shlex  # simple lexical analysis for unix syntax
from typing import Optional, Tuple

# ================================================================================

PROFILE_FILENAME = ".kids_profile"  # file that defines the working environment

# ================================================================================

def which(program: str) -> str:
    """
    Locate a program file in the user's path.

    Parameters:
    program (str): The name of the executable file to locate.

    Returns:
    str: The full path to the executable if found, otherwise None.
    """
    # Retrieve the PATH environment variable and split it into a list using the path separator
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

# ================================================================================

def source(profilefile: str) -> str:
    """
    Source the content of a file to define some environmental variables,
    it is a simple imitation of the "source" bash command.

    Keyword arguments:
    profilefile -- the name (with absolute or relative path) of the file
                   that contains the definitions of some environmental variables
                   and possibly some additional bash commands
    """

    # Split the bash source command into a list that can be used with subprocess
    command = shlex.split(f"bash -c 'source {profilefile} && env'")

    # Call the source command with subprocess and get the stdout as a pipe
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)

    # Process the standard output line by line
    for line in proc.stdout:
        # Redefine each environmental variable (both old and new) with the updated values from the source command
        (key, _, value) = line.decode().partition("=")
        try:
            if key != "BASH_FUNC_module%%" and key != "BASH_FUNC_ml%%":
                # Remove any trailing newline characters and set the environment variable
                os.environ[key] = value.rstrip('\n')
        except ValueError:
            # If there is a ValueError (e.g., no '=' in the line), do nothing
            pass

    # Wait for the process to terminate and get the remaining output
    proc.communicate()


def setKIDSProfile() -> Optional[str]:
    """
    Sets the main profile by sourcing a profile file if it exists in the user's home directory or the current working directory. 
    If no profile file is found, it is assumed that the environment is already set up, and no action is taken.

    Returns:
    Optional[str]: The path to the profile file if it was found and sourced, otherwise None.
    """
    # Initialize the filenameconf variable to None
    filenameconf: Optional[str] = None

    # Define the list of directories to search for the profile file: first the current working directory, then the home directory
    for directory in [os.getcwd(), os.path.expanduser('~')]:
        # Construct the full path to the profile file
        profile_path: str = os.path.join(directory, PROFILE_FILENAME)
        # Check if the profile file exists at the constructed path
        if os.path.isfile(profile_path):
            # If the file exists, set filenameconf to the path and break out of the loop
            filenameconf = profile_path
            break

    # If a profile file was found, source it to set the environmental variables
    if filenameconf:
        # print('sourcing ', filenameconf)
        source(filenameconf)

    # Return the path to the profile file if it was found and sourced, otherwise None
    return filenameconf

def checkKIDSProfile(filenameconf: str) -> str:
    """
    Generates a log message about the environment configuration based on the presence of the profile file.

    Parameters:
    filenameconf (str): The path name of the profile file.

    Returns:
    str: A log message indicating whether the profile file was found and what it implies about the environment setup.
    """
    logstring: str = ""
    if filenameconf:
        logstring += f"""The file {filenameconf} is present
It will be read to set the paths of third-party software\n"""
    else:
        logstring += f"""The {PROFILE_FILENAME} is NOT present in your home or working directory
It is assumed that the environment has been already set\n"""
    logstring += '\n'

    return logstring

def checkKIDSEnv() -> Tuple[bool, str]:
    """
    Checks if the environment for KIDS SCRIPTS execution is properly defined and if all the necessary executables are available.

    Returns:
    Tuple[bool, str]: A tuple where the first element is a boolean indicating if the environment is properly set up,
                      and the second element is a string containing an error message if any check fails.
    """
    # Check if the environmental variable with the <KIDS SCRIPTS> path is defined
    try:
        kids_scripts_path = os.environ['KIDS_SCRIPTS_PATH']
        os.environ['PATH'] += os.pathsep + kids_scripts_path
    except KeyError:
        message = "Environment variable $KIDS_SCRIPTS_PATH is not defined"
        return False, message

    # Check if <KIDS SCRIPTS> executables are available
    for exe in ["kidsqm.py", "kidsqmmm.py"]:
        if not which(exe):
            message = f"{exe} executable is not available"
            return False, message

    # Return True if all the checks were OK
    return True, ""

# ================================================================================

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

def checkQMEnv(qmsolver:str):
    if qmsolver == 'bdf':
        pass
    elif qmsolver == 'bagel':
        return checkBagelQMEnv()
    elif qmsolver == 'mndo':
        return checkMNDOQMEnv()
    elif qmsolver == 'gaussian':
        return checkGaussianQMEnv()
    elif qmsolver == 'molcas':
        return checkMolcasEnv()
    elif qmsolver == 'molpro':
        return checkMolproEnv()

    return False, f"Unknown qmsolver {qmsolver}"


def checkGaussianOptEnv():
    """ the function checks if the environment for GAUSSIAN execution
       for the geometry optimization is properly defined and if all
       the necessary executables are available """

    # get the name of the executable and the path from the environmental variable
    try:
        gversion = os.environ["GAUSSIAN_EXE"]
        gpath = os.environ["GAUSSIAN_DIR"]
        _ = os.environ["GAUSS_SCRDIR"]
    except KeyError:
        return False, "environment variables for Gaussian optimization, $GAUSSIAN_EXE, " \
                      "$GAUSSIAN_DIR and $GAUSS_SCRDIR, are not defined"

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

def checkMNDOQMEnv():
    """ the function checks if the environment for MNDO execution
       for QM calculations is properly defined and if all the
       necessary executables are available """

    # get the name of the executable and the path from the environmental variable
    try:
        mndoversion = os.environ["MNDO_EXE_QM"]
        mndopath = os.environ["MNDO_DIR_QM"]
    except KeyError:
        return False, "environment variables for MNDO, $MNDO_EXE_QM, " \
                      "$MNDO_DIR_QM are not defined"

    # check that the gversion that is passed corresponds to a supported version
    #if mndoversion != "mndo99" and mndoversion != "mndo2020":
    #    return False, "{0} version of mndo is not currently supported by Cobramm".format(mndoversion)

    # check if mndo executable is available
    if not which(mndoversion):
        message = mndoversion + " executable is not available"
        return False, message

    # return True if all the checks were OK
    return True, ""

# ================================================================================

def checkBagelQMEnv():
    """ the function checks if the environment for BAGEL execution
       for QM calculations is properly defined and if all the
       necessary executables are available """

    # get the name of the executable and the path from the environmental variable
    try:
        bagelversion = os.environ["BAGEL_EXE_QM"]
        bagelpath = os.environ["BAGEL_DIR_QM"]
    except KeyError:
        return False, "environment variables for BAGEL, $BAGEL_EXE_QM, " \
                      "$BAGEL_DIR_QM are not defined"

    # check if bagel executable is available
    if not which(bagelversion):
        message = bagelversion + " executable is not available"
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


# ================================================================================

def getVersion():
    """
    get last git tag and date of the last commit
    :return: string: "{tag} {date}"
    """
    return "custom"
    wordDir = os.getcwd()
    cobrammPath = os.environ["COBRAM_PATH"]
    os.chdir(cobrammPath)
    try:
        tag = subprocess.check_output(["git", "describe", "--tags", "--abbrev=0"]).strip().decode()
        date = subprocess.check_output(["git", "log", "-1", "--format=%cd", "--date=short"]).strip().decode()
        cobramm_version = '{} {}'.format(tag, date)
    except:
        cobramm_version = 'N.A.'
    os.chdir(wordDir)
    return cobramm_version


# ================================================================================

def getHostname():
    """
    get the hostname of the machine where cobramm is running
    :return: string hostname
    """
    # We use socket instead of reading env variables for better compatibility across systems
    return __import__('socket').gethostname()

# ================================================================================
