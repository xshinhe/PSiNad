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

from contextlib import ContextDecorator
import time

#####################################################################################################


class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""

#####################################################################################################


class Timer(ContextDecorator):
    """The timer code accumulates the times spent in portion of the codes, storing
     the partial timings in the class dictionaries. The class function
     timingreport can be used to print all the timings accumulated.
     The timers dictionary stores a list of three elements: 1) the accumulated CPU time,
     2) the accumulated wall time and 3) the nr of calls of this timer """

    timers = dict()

    # =============================================================================================================

    def __init__(self, name=None):
        """Initialization: add timer to dict of timers"""

        self._start_CPU_time = None
        self._start_wall_time = None
        self.name = name

        if self.name:
            self.timers.setdefault(self.name, [0.0, 0.0, 0])

    # =============================================================================================================

    def start(self):
        """Start a new timer"""

        if self._start_CPU_time is not None or self._start_wall_time is not None:
            raise TimerError("Timer {0} already running. Please use .stop() to stop it".format(self.name))

        self._start_CPU_time = time.process_time()
        self._start_wall_time = time.perf_counter()

    # =============================================================================================================

    def stop(self):
        """Stop the timer, and report the elapsed time"""

        if self._start_CPU_time is None or self._start_wall_time is None:
            raise TimerError("Timer {0} not yet running. Please use .start() to start it".format(self.name))

        # Calculate elapsed time
        CPU_elapsed_time = time.process_time() - self._start_CPU_time
        wall_elapsed_time = time.perf_counter() - self._start_wall_time
        self._start_CPU_time = None
        self._start_wall_time = None

        # Accumulate elapsed time
        if self.name:
            self.timers[self.name][0] += CPU_elapsed_time
            self.timers[self.name][1] += wall_elapsed_time
            self.timers[self.name][2] += 1

        return CPU_elapsed_time, wall_elapsed_time

    # =============================================================================================================

    def __enter__(self):
        """Start a new timer as a context manager"""
        self.start()
        return self

    # =============================================================================================================

    def __exit__(self, *exc_info):
        """Stop the context manager timer"""
        self.stop()

    # =============================================================================================================

    @staticmethod
    def timingreport(userprint=None):

        # convert the dictionary to an ordered list of tuples
        sorted_timers = sorted(Timer.timers.items(), key=lambda kv: kv[1][1], reverse=True)

        log = ""
        log += "==================================================================================================\n"
        log += " Code section  | tot WallTime/ s | tot CPUTime/ s  | Nr calls | avg WallTime/ s | avg CPUTime/ s  \n"
        log += "--------------------------------------------------------------------------------------------------\n"
        for section, times in sorted_timers:
            cputime, walltime, nrcalls = times
            if nrcalls == 0: continue
            log += "{0:14s} | {1:15.2f} | {2:15.2f} | {3:8d} | {4:15.2f} | {5:15.2f}\n".format(
                section, walltime, cputime, nrcalls, walltime/nrcalls, cputime/nrcalls)
        log += "==================================================================================================\n\n"

        if userprint is None:
            print(log)
        else:
            userprint(log)
