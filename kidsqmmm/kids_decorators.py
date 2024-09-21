#!/usr/bin/env python3

from contextlib import ContextDecorator
from pprint import pprint
from typing import Dict, Union
import time

PROFILE_FILENAME = ".cobramm_profile"  # name of the bash file that defines the working environment for COBRAMM

class Timing(ContextDecorator):
    timers : Dict[str, Union[float, float, int]] = dict()

    def __init__(self, name : Union[str, None] = None):
        self._start_CPU_time = None
        self._start_wall_time = None
        self.name = name
        if self.name:
            self.timers.setdefault(self.name, [0.0, 0.0, 0])

    def __enter__(self):
        self._start_CPU_time = time.process_time()
        self._start_wall_time = time.perf_counter()
        return self

    def __exit__(self, *exc_info):
        # Calculate elapsed time
        CPU_elapsed_time = time.process_time() - self._start_CPU_time
        wall_elapsed_time = time.perf_counter() - self._start_wall_time
        # Accumulate elapsed time
        if self.name:
            self.timers[self.name][0] += CPU_elapsed_time
            self.timers[self.name][1] += wall_elapsed_time
            self.timers[self.name][2] += 1
        return CPU_elapsed_time, wall_elapsed_time

    @staticmethod
    def report(userprint=None):
        sorted_timers = sorted(Timing.timers.items(), key=lambda kv: kv[1][1], reverse=True)
        log = ""
        log += "=============================================================================================\n"
        log += " Section  | tot WallTime/ s | tot CPUTime/ s  | Nr calls | avg WallTime/ s | avg CPUTime/ s  \n"
        log += "---------------------------------------------------------------------------------------------\n"
        for section, times in sorted_timers:
            cputime, walltime, nrcalls = times
            if nrcalls == 0: continue
            log += "{0:9s} | {1:15.2f} | {2:15.2f} | {3:8d} | {4:15.2f} | {5:15.2f}\n".format(
                section, walltime, cputime, nrcalls, walltime/nrcalls, cputime/nrcalls)
        log += "=============================================================================================\n\n"
        if userprint is None:
            print(log)
        else:
            userprint(log)

class Inworkspace(ContextDecorator):
    def __init__(self, rundir : str, envs : Dict[str, str]):
        # store starting dir
        self.startdir = os.path.abspath(os.getcwd())
        self.rundir = os.path.abspath(rundir)
        self.envs = envs

    def __enter__(self):
        if !os.exists(rundir):
            os.makedirs(rundir)
        os.chdir(rundir)
        for k,v in self.envs: os.putenv(k, v)
        return self

    def __exit__(self, *exc_info):
        os.chdir(startdir)

    def clean():
        if os.path.relpath(rundir, startdir).startswith('..'): # be safe
            shutil.rmtree(rundir)

class Signal:
    def __init__(self, path:str = ".", waitingtime:float = 0.1):
        self.path = path
        self.waitingtime = waitingtime

    def __del__(self):
        del self.path
        del self.waitingtime

    def sendSignal(self, channel:int = 0):
        fileName = os.path.join(self.path, "SIGNAL-{0}".format(channel))
        with open(fileName, "w") as f:
           f.write("{0}".format(os.getpid()))

    def waitTillSignal(self, channel:int = 0):
        k = None
        fileName = os.path.join(self.path, "SIGNAL-{0}".format(channel))
        while True:
           # check whether the signal file is present
           signalStat = os.path.isfile(fileName)
           try:
              with open(fileName, "r") as f:
                 k = int(f.read())
           except:
              pass
           # break the loop only when the signal file is found AND its content is different from 
           # the current process PID (to prevent a process from sending a signal to itself)
           if signalStat and k != os.getpid(): break
               break
           # if still inside the loop, wait a bit before next attempt to read the control file
           time.sleep( self.waitingtime )
        os.remove(fileName)
        return



#=================================================================================================

def which(program):
    """ Check if a filename corresponds to an executable in one of the PATH directories

    Keyword arguments:
    program -- the name of the executable file
    """

    # loop over the directories of the PATH variable
    for path in os.environ["PATH"].split(os.pathsep):

        # check if path+program combination exists and has execution permission
        exe_file = os.path.join(path, program)
        if os.path.isfile(exe_file) and os.access(exe_file, os.X_OK):
            # when the combination is found, return it
            return exe_file

    # if no directory in PATH contains an executable program, return None
    return None


# ================================================================================

def source(profilefile):
    """ Source the content of a file to define some enviromental variables,
        it is a simple imitation of the "source" bash command.

    Keyword arguments:
    profilefile -- the name (with absolute or relative path) of the file
                   that contains the definitions of some environmental variables
                   and possibly some additional bash commands
    """

    # split the bash source command in a list that can be used with subprocess
    command = shlex.split("bash -c 'source " + profilefile + " && env'")

    # call the source command with subprocess and get the stdout as pipe
    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    # process the std out line-by-line
    for line in proc.stdout:
        # redefine each environmental variables (both old and new) with the updated values from the source command
        (key, _, value) = line.decode().partition("=")
        try:
            if key != "BASH_FUNC_module%%" and key != "BASH_FUNC_ml%%":
                os.environ[key] = value.replace('\n', '')
        except ValueError:
            pass
    proc.communicate()


# ================================================================================

def setCobrammProfile():
    """ Define environmental variables for COBRAMM execution

    If the profile file is available in the home directory of the
    user or in the current working directory, the file is used to set
    the environmental variables that are needed for the execution of COBRAMM.
    If no profile file is available, it is assumed that a
    working environment has been already set, and nothing is done.
    """

    # define the path of the available PROFILE_FILENAME file, first check PWD, then HOME
    filenameconf = None
    for d in [os.getenv('PWD'), os.getenv('HOME')]:
        fpath = os.path.join(d, PROFILE_FILENAME)
        if os.path.isfile(fpath):
            filenameconf = fpath
            break

    # if a PROFILE_FILENAME is available, source it
    if filenameconf: source(filenameconf)

    return filenameconf


# ================================================================================

def setGaussianProfile(gversion, gpath):
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


# ================================================================================

def checkConfigFile(filenameconf):
    """ Returns log message about the environment configuration

    Keyword arguments:
    filenameconf -- the path name of the COBRAMM profile file
    """

    logstring = ""
    if filenameconf:
        logstring += """The file """ + filenameconf + """ is present
It will be read to set the paths of third-party software\n"""
    else:
        logstring += """The """ + PROFILE_FILENAME + """ is NOT present in your home or working directory
It is assumed that the environment has been already set\n"""
    logstring += '\n'

    return logstring


# ================================================================================

def checkqmmmenv():
    """ the function checks if the environment for AMBER execution is properly
        defined and if all the necessary executables are available """

    # environmental variable with the COBRAMM path
    try:
        os.environ['COBRAM_PATH']
    except KeyError:
        message = "environment variable $COBRAM_PATH is not defined"
        return False, message

    # check if COBRAMM executables are available
    for exe in ["cobram.py", "cobramext", "freqext"]:
        if not which(exe):
            message = exe + " executable is not available"
            return False, message

    # return True if all the checks were OK
    return True, ""


# ================================================================================

def checkAmberEnv():
    """ the function checks if the environment for AMBER execution is properly defined and
        if all the necessary executables are available """

    # environmental variable with the AMBER path
    try:
        AMBERHOME = os.environ['AMBERHOME']
    except KeyError:
        message = "environment variable $AMBERHOME is not defined"
        return False, message
    
    ambersource = False
    try:
        AMBERSOURCE = os.environ['AMBERSOURCE']
        ambersource = True
    except KeyError:
        ambersource = False

    # call update_amber to check the current version of amber
    if ambersource == True:
        command = shlex.split(os.path.join(AMBERSOURCE, "update_amber") + " -v")
    else:
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


# ===============================================================================

class CheckEnvDecor(ContextDecorator):
    def __init__(self, name : str):
        self.name = name

    def __enter__(self):
        # get the name of the executable and the path from the environmental variable
        try:
            exe = os.environ[f"{self.name}_EXE"]
            exe_dir = os.environ[f"{self.name}_DIR"]
        except KeyError:
            return False, "environment variables for Gaussian optimization, $GAUSSIAN_EXE, " \
                          "$GAUSSIAN_DIR and $GAUSS_SCRDIR, are not defined"
        if not which(gversion):
        message = gversion + " executable is not available"
        return False, message

        return self

    def __exit__(self, *exc_info):
        # Calculate elapsed time
        CPU_elapsed_time = time.process_time() - self._start_CPU_time
        wall_elapsed_time = time.perf_counter() - self._start_wall_time
        # Accumulate elapsed time
        if self.name:
            self.timers[self.name][0] += CPU_elapsed_time
            self.timers[self.name][1] += wall_elapsed_time
            self.timers[self.name][2] += 1
        return CPU_elapsed_time, wall_elapsed_time


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
    if mndoversion != "mndo99" and mndoversion != "mndo2020":
        return False, "{0} version of mndo is not currently supported by Cobramm".format(mndoversion)

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
