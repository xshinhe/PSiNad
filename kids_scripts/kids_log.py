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

import sys  # system commands
import time  # provides various time-related functions
import random  # random number generator
from typing import Dict, Union, Optional, Any
from contextlib import ContextDecorator
from pprint import pprint
from traceback import format_exc

# imports of local objects
import softenv

# math libraries
import numpy as np  # numpy library for scientific computation

#####################################################################################################
# global module: Timing, Log
#####################################################################################################

# Define a ContextDecorator class for timing code sections
class Timing(ContextDecorator):
    # Class variable to store timing information
    timers: Dict[str, Union[float, float, int]] = dict()

    def __init__(self, name: Union[str, None] = None):
        """
        Initialize the Timing object with an optional name.
        
        :param name: The name of the code section to be timed.
        """
        self._start_CPU_time = None  # CPU time at the start of the timing
        self._start_wall_time = None  # Wall time at the start of the timing
        self.name = name  # Name of the code section
        if self.name:
            # Initialize the timer for the given section if it doesn't exist
            self.timers.setdefault(self.name, [0.0, 0.0, 0])

    def __enter__(self):
        """
        Start timing when entering the context.
        
        :return: The Timing object itself.
        """
        self._start_CPU_time = time.process_time()  # Start CPU time measurement
        self._start_wall_time = time.perf_counter()  # Start wall time measurement
        return self

    def __exit__(self, *exc_info):
        """
        Stop timing when exiting the context and accumulate the elapsed time.
        
        :param exc_info: Exception information (if any) that occurred within the context.
        :return: A tuple containing the CPU and wall elapsed times.
        """
        # Calculate elapsed time
        CPU_elapsed_time = time.process_time() - self._start_CPU_time
        wall_elapsed_time = time.perf_counter() - self._start_wall_time
        # Accumulate elapsed time
        if self.name:
            self.timers[self.name][0] += CPU_elapsed_time  # Accumulate CPU time
            self.timers[self.name][1] += wall_elapsed_time  # Accumulate wall time
            self.timers[self.name][2] += 1  # Increment the call count

        if exc_info[0]:
            pprint(exc_info[1])
            return False 
        return True

    @staticmethod
    def report(userprint=None):
        """
        Generate and print a report of the accumulated timing information.
        
        :param userprint: A function to use for printing the report. If None, print to stdout.
        """
        # Sort the timers by wall time in descending order
        sorted_timers = sorted(Timing.timers.items(), key=lambda kv: kv[1][1], reverse=True)
        log = ""
        log += "================================================================================\n"
        log += " Code section      | total WallTime/ s   | total CPUTime/ s    | No. calls      \n"
        log += "--------------------------------------------------------------------------------\n"
        for section, times in sorted_timers:
            cputime, walltime, nrcalls = times
            if nrcalls == 0: continue  # Skip sections that haven't been called
            log += " {0:14s}    | {1:15.2f}     | {2:15.2f}     | {3:8d}       \n".format(
                section, walltime, cputime, nrcalls)
        log += "================================================================================\n\n"
        if userprint is None:
            print(log)  # Print the report to stdout
        else:
            userprint(log)  # Use the provided function to print the report


class Log:
    # Class variables
    DEBUG_COBRAMM_RUN = False
    VERBOSITY_LEVEL: int = 0
    writeLogToStdOut: bool = True

    @staticmethod
    def setVerbosityLevel(levelnumber: Optional[int] = None) -> None:
        """
        Set the verbosity level for logging.

        :param levelnumber: Desired verbosity level. If None, level remains unchanged.
        """
        if levelnumber is not None:
            try:
                Log.VERBOSITY_LEVEL = int(levelnumber)
            except (ValueError, TypeError):
                print(f"Invalid verbosity level provided: {levelnumber}. It must be an integer.")
                return

    @staticmethod
    def writeLog(message: Any, level: int = 0) -> None:
        """
        Write a message to the log.

        :param message: The message to log.
        :param level: The verbosity level of the message.
        """
        if level <= Log.VERBOSITY_LEVEL:
            if Log.writeLogToStdOut:
                print(message, end='')
            else:
                try:
                    with open('exec.log', 'a') as log_file:
                        log_file.write(str(message))
                        log_file.flush()
                except IOError as e:
                    print(f"Failed to write to log file: {e}")

    @staticmethod
    def printGeom(geometry, level: int = 0):
        """
        Print molecular geometry information to the log.

        :param geometry: The molecular geometry object.
        """
        if Log.VERBOSITY_LEVEL <= 1 and level <= 1:
            return

        if geometry.atomNum > 500:
            Log.writeLog("Only HM atom coordinates will be printed here\n\n")
            atoms = zip(geometry.list_MEDIUM_HIGH, geometry.getAtomLabels("MEDIUM_HIGH"),
                        *geometry.getModel("MEDIUM_HIGH"))
        else:
            atoms = zip(range(1, geometry.atomNum + 1), geometry.atomLabel, *geometry.cartesian)

        Log.writeLog(
"""----------------------------------------------------------------
   Atom     Atomic              Coordinates (Angstroms)
     ID      Label             X           Y           Z
----------------------------------------------------------------
""")
        for at in atoms:
            Log.writeLog(" {0:6d} {1:>10s}     {2:12.6f}{3:12.6f}{4:12.6f}\n".format(*at))
        Log.writeLog("----------------------------------------------------------------\n\n")
        Log.writeLog("\n")
        return

    @staticmethod
    def printLayers(listH, listM, listL, listBA):
        """
        Print layer and link atom information to the log.

        :param listH: List of high layer atoms.
        :param listM: List of medium layer atoms.
        :param listL: List of low layer atoms.
        :param listBA: List of atom links.
        """
        if len(listH) > 0:
            Log.writeLog(" * {:6d} HIGH   layer atom(s): ".format(len(listH)))
            Log.writeLog(Log.prettyAtomsList(listH) + "\n")
        else:
            Log.writeLog(" *     no MEDIUM layer atom \n")

        if len(listM) > 0:
            Log.writeLog(" * {:6d} MEDIUM layer atom(s): ".format(len(listM)))
            Log.writeLog(Log.prettyAtomsList(listM) + "\n")
        else:
            Log.writeLog(" *     no MEDIUM layer atom \n")

        if len(listL) > 0:
            Log.writeLog(" * {:6d} LOW    layer atom(s): ".format(len(listL)))
            Log.writeLog(Log.prettyAtomsList(listL) + "\n")
        else:
            Log.writeLog(" *     no LOW    layer atom \n")

        if len(listBA[0]) > 0:
            Log.writeLog(" *   atom links  ( QM --> MM ): ")
            cnt = 0
            for A, B in zip(*listBA):
                pre = ''
                if cnt > 0: pre = ' '*32 
                Log.writeLog(pre + " {} --> {}\n".format(B, A))
                cnt += 1
            Log.writeLog("\n")
        Log.writeLog("\n")
        return

    @staticmethod
    def prettyAtomsList(atomList):
        """
        Return a compact string representation of a list of atom indices.

        :param atomList: List of atom indices.
        :return: String representation of the list.
        """
        string = ""
        seqStart = None
        seqPre = None

        maxlen = 50
        curlen = 0
        for atom in atomList:
            if seqPre is None:
                seqStart = atom
            else:
                if seqPre != atom - 1:
                    if seqStart != seqPre:
                        s = " {0}-{1}".format(seqStart, seqPre)
                        curlen += len(s)
                        if curlen > maxlen:
                            curlen = 0
                            s = '\n' + ' '*20 + s
                        string += s
                    else:
                        s = " {0}".format(seqStart)
                        curlen += len(s)
                        if curlen > maxlen:
                            curlen = 0
                            s += '\n' + ' '*20
                        string += s
                    seqStart = atom
            seqPre = atom

        if seqStart != atomList[-1]:
            string += " {0}-{1}".format(seqStart, atomList[-1])
        else:
            string += " {0}".format(atomList[-1])

        return string

    @staticmethod
    def printEnergies(modelH_QM, modelH_MM, real_MM, only_MM, chargeInt, full_QMMM, modelLabel):
        """
        Print energy information to the log.

        :param modelH_QM: List of QM energies.
        :param modelH_MM: List of MM energies.
        :param real_MM: Real MM energy.
        :param only_MM: Only MM energy.
        :param chargeInt: Charge interaction energy.
        :param full_QMMM: Full QM/MM energy.
        :param modelLabel: Label of the model.
        """
        if modelLabel == "M" or modelLabel == "ML":
            Log.writeLog(' ' * 20 + 'Energies (Hartrees)\n')
            Log.writeLog('   MM energy:  ' + '%18.8f\n' % only_MM)
        else:
            for i, E_QM, E_QMMM in zip(range(len(modelH_QM)), modelH_QM, full_QMMM):
                Log.writeLog('-' * 80 + "\n")
                Log.writeLog(' ' * 20 + 'STATE %3d' % (i + 1) + "\n")
                Log.writeLog('-' * 80 + "\n")
                Log.writeLog(' ' * 20 + 'Energies (Hartrees)\n')
                Log.writeLog('   Model-H QM:  {0:18.8f}\n'.format(E_QM) +
                             '   Real MM:     {0:18.8f}\n'.format(real_MM) +
                             '   Model-H MM:  {0:18.8f}\n'.format(modelH_MM) +
                             '   Emb-emb crg: {0:18.8f}\n'.format(chargeInt) +
                             '   QM Energy:   {0:18.8f}\n\n'.format(E_QM - chargeInt) + '\n' +
                             'E(tot)=E(Model-H QM)+(Real MM)-(Model-H MM)-(Emb-emb)= {0:18.8f}\n\n'.format(E_QMMM))

        return

    @staticmethod
    def printModelHCharges(geometry, charges, dipole):
        """
        Print model-H QM charges and dipole moment to the log.

        :param geometry: The molecular geometry object.
        :param charges: List of charges.
        :param dipole: List of dipole moment components.
        """
        Log.writeLog(""" ----------------------------------------------
    Atom     Atomic       Model+H Charges
      ID      Label           (a.u.)
 ----------------------------------------------
""")
        labelModelH = [geometry.atomLabel[i - 1] for i in geometry.list_QM]
        for i in range(geometry.NsubH):
            labelModelH.append('H')
        for n in range(len(labelModelH)):
            Log.writeLog("  {0:6d} {1:>10s}     {2:12.6f}\n".format(n + 1, labelModelH[n], charges[n]))
        Log.writeLog(" ----------------------------------------------\n")
        Log.writeLog("         {0:>10s}     {1:12.6f}\n".format("TOTAL", sum(charges)))
        Log.writeLog(" ----------------------------------------------\n\n")

        kconv_au_Deb = 0.39342215569939517797

        Log.writeLog(""" ----------------------------------------------
                         Dipole Moment
                     (a.u.)         (Debye)
 ----------------------------------------------
""")
        componentsLab = ["X", "Y", "Z", "magnitude"]
        for n in range(3):
            Log.writeLog("  {0:>10s}  {1:14.6f} {2:14.6f}\n".format(componentsLab[n], dipole[n], dipole[n] / kconv_au_Deb))
        Log.writeLog(" ----------------------------------------------\n")
        Log.writeLog("  {0:>10s}  {1:14.6f} {2:14.6f}\n".format(componentsLab[3], dipole[3], dipole[3] / kconv_au_Deb))
        Log.writeLog(" ----------------------------------------------\n\n")

        return

    @staticmethod
    def printMDinfo(actualTime, state, tStep, mdItems, potReference):
        Log.writeLog("{0:<20s} {1:f} fs\n\n".format("Simulation time:", actualTime / 41.341373337))
        Log.writeLog("{0:<30s} {1:d}\n".format("Current electronic state:", state + 1))
        Log.writeLog("{0:<30s} {1:f} fs\n\n".format("Velocity-Verlet time-step:", tStep / 41.341373337))
        Log.writeLog("{0:<20s} {1:f} eV\n".format("Kinetic energy:", mdItems[2] * 27.211385))
        Log.writeLog("{0:<20s} {1:f} eV\n".format("Potential energy:", (mdItems[1] - potReference) * 27.211385))
        Log.writeLog("{0:<20s} {1:f} eV\n".format("Total energy:", (mdItems[0] - potReference) * 27.211385))
        Log.writeLog("\n")

    @staticmethod
    def startSection(sectionName):
        """
        Write the title of a section to the log.

        :param sectionName: The name of the section.
        """
        titleLen = len(sectionName)
        blankSpaces = int((80 - titleLen) / 2)

        output = ('=' * 80 + '\n' +
                  ' ' * blankSpaces + sectionName + '\n' +
                  '-' * 80 + '\n\n')
        Log.writeLog(output)

    @staticmethod
    def startSubSection(subSectionName):
        """
        Write the title of a subsection to the log.

        :param subSectionName: The name of the subsection.
        """
        titleLen = len(subSectionName)
        blankSpaces = int((80 - titleLen) / 2)

        output = (' ' * blankSpaces + "*" * titleLen + "\n" +
                  ' ' * blankSpaces + subSectionName + '\n' +
                  ' ' * blankSpaces + "*" * titleLen + "\n\n")
        Log.writeLog(output)

    @staticmethod
    def writewarning(message):
        """
        Write a warning message to the log.

        :param message: The warning message to log.
        """
        Log.writeLog('\nWARNING! {0}\n'.format(message))

    @staticmethod
    def fatalError(message):
        """
        Log an error message and exit the program.

        :param message: The error message to log.
        """
        Log.startSection("CALCULATION TIMES")
        Timing.report(Log.writeLog)

        # Log the error message
        Log.writeLog('FATAL:  {0}\n\n'.format(message) +
                     '=' * 80 + '\n' +
                     'Ending time: {0} Hostname {1}\n'.format(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), softenv.getHostname()) +
                     '=' * 80 + '\n\n' +
                     '--- COBRAM calculation ABORTED!!! ---\n')

        # Exit with status 3
        sys.exit(3)

    @staticmethod
    def start():
        """ 
        Write the starting message
        """

        Log.startSection("EXECUTE KIDS SCRIPTS CALCULATION")
        # print final message to log
        Log.writeLog('''
+------------------------------------------------------------------------------+
|             .,;;                                                             |
|               ;.                                                             |
|    .XX          dXd                                                          |
|   :''.,,       lx0kll, .lxl                                                  |
|     ''.         kKk               mm   mmm    mmmmmm    mmmmm        mmmm    |
|        ..     :                   ##  ##"     ""##""    ##"""##    m#""""#   |
|        ';ddddo;           oOx     ##m##         ##      ##    ##   ##m       |
|:kxko   ,kkkkkkl;  .kxxk dd.c'     #####         ##      ##    ##    "####m   |
|'OOO;   .kkkkkkd   .dxxl.          ##  ##m       ##      ##    ##        "##  |
|          .dd:         .           ##   ##m    mm##mm    ##mmm##    #mmmmm#"  |
|             .          lX         Kernel   Integrated   Dynamics   Simulator |
|              ,        .,:.                                                   |
|            l:,:c.                          KIDS SCRIPTS PART                 |
|             ;::.                           Copyright    2024                 |
|                                                                              |
+------------------------------------------------------------------------------+
    Author:         Prof. Jian Liu
    Contributors:   Xin He, Haocheng Lu, Bingqi Li, Baihua Wu,
                    Xiangsong Cheng, Youhao Shang.
''')

    @staticmethod
    def end():
        """ 
        Write the final message and terminate execution normally,
        exit status is equal to 0 
        """

        Log.startSection("CALCULATION TIMES")
        Timing.report(Log.writeLog)

        # print final message to log
        Log.writeLog('\n' +
             '=' * 80 + '\n' +
             'Ending time: {0} Hostname {1}\n'.format(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), softenv.getHostname()) +
             '=' * 80 + '\n\n' +
             '--- calculation terminated normally ---\n')
        # exit with status 0
        sys.exit(0)

# ===============================================================================================================

def matrix_prettystring(matrix: np.ndarray, fmt: str = ".6f", atomLabels=None) -> str:
    """
    Return a string with a nice formatted representation of a matrix, stored in a ndarray numpy object.
    The optional argument fmt should contain a string that defines the numeric format of the
    data to represent (without the field size that will be determined dinamically.

    :param matrix: numpy ndarray with the matrix to print to string
    :param fmt: numerical format to use in the matrix representation
    :return: string that contains the representation of the matrix
    """
    # compute the size that is needed to represent each column of the matrix
    colmax = [max([len(("{:" + fmt + "}").format(x))+1 for x in col]) for col in matrix.T]
    # construct the format of the line
    rowformat = "".join(["{:" + str(colsize) + fmt + "}" for colsize in colmax])
    # print each row with the appropriate format and store the result in a string
    if atomLabels is None:
        string = "\n".join([rowformat.format(*row) for row in matrix]) + "\n"
    else:
        string = "\n".join([("{:4s}"+rowformat).format(atomLabels[irow],*matrix[irow]) for irow in range(len(matrix))]) + "\n" 
    # return the string with the matrix representation
    return string

#####################################################################################################
