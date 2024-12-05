#! /usr/bin/env python3
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

# import statments of module from python standard library

import os                       # operating system utilities 
import time                     # time functions 

#####################################################################################################

class extSignal:
    """ the class defines a procedure to communicate between different processes """

    ###########################################################################################

    def __init__(self, path=".", waitingtime=0.1):
        """ constructor of the extSignal class:
            simply store few input values (or their default values): 
            path        : the path where the signal file will be written or read 
            waitingtime : the number of second waited between each attemp
                          to catch the signal  """ 

        # store the input path
        self.path = path
        # store the waiting time
        self.waitingtime = waitingtime

    ###########################################################################################

    def __del__(self):
        """ destructor of the extSignal class, not much need to be done """

        # destroy stored content of self
        del self.path
        del self.waitingtime

    ###########################################################################################

    def sendSignal(self, channel=0):
        """ send a system-wide signal, by creating a file in the working dir
            that contains the pid of the current process """

        # define the path and name of the control file
        fileName = os.path.join(self.path, "COBRAMM-SIGNAL-{0}".format(channel))

        # write current process PID to the signal file
        with open(fileName, "w") as f:
           f.write("{0}".format(os.getpid()))

    ###########################################################################################

    def waitTillSignal(self, channel=0):
        """ pause the process, until the signal file is found in the working dir """

        k = None            # initialize the number that will be read from fileName
                            # k is None so that no input number can match its value

        # define the path and name of the control file
        fileName = os.path.join(self.path, "COBRAMM-SIGNAL-{0}".format(channel))

        while True:
           # check whether the signal file is present
           signalStat = os.path.isfile(fileName)

           # try to read a number from the file fileName, if an exception is raised, just go on
           try:
              with open(fileName, "r") as f:
                 k = int(f.read())
           except:
              pass

           # break the loop only when the signal file is found AND its content is different from 
           # the current process PID (to prevent a process from sending a signal to itself)
           if signalStat and k != os.getpid():
               break
           
           # if still inside the loop, wait a bit before next attempt to read the control file
           time.sleep( self.waitingtime )

        # if this point is reached, the signal has been received and can be switch off 
        os.remove(fileName)
        return

    ###########################################################################################
    

