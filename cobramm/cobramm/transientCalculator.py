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
import os  # filesystem utilities
import shutil  # filesystem utilities
import subprocess  # run external program as child process
import fileinput  # fileinput is used for inplace editing of files
import shlex  # simple lexical analysis for unix syntax
import math  # import mathematical functions
import re  # process output with regular expressions
import copy  # shallow and deep copy operations
import time  # provides various time-related functions

# imports of local modules

import logwrt  # write messages and output to log (std out or cobramm.log)
import constants  # physical and mathematical constants
import cobrammenv  # function to check if executables are available
from gaussianCalculator import GaussianOutput
from cobrammCalculator import CobrammCalculator, CobrammInput, CobrammOutput

# math libraries

import numpy as np  # numpy: arrays and math utilities

#try:
#    sys.path.append(os.path.join(os.environ['COBRAM_PATH'], 'cobramm'))
#except KeyError:
#    raise RuntimeError('cannot import cobramm module.\n'
#                       'Please check if $COBRAM_PATH is defined in your environment, if not run:\n'
#                       'export COBRAM_PATH=/path/to/COBRAMM')

#####################################################################################################

class SpectronError(Exception):
    """A custom exception used to report errors in use of Spectron"""

class MultiwfnError(Exception):
    """A custom exception used to report errors in use of Multiwfn"""

#####################################################################################################

class PumpProbeCalculator:

    def __init__(self, trjdir_list, method="TDDFT", simulation_time=500, time_step=0.5, delta_t=2, nstates=20, en_min=0.0,
                 en_max=10.0, en_width=0.15, grid=400, t_width=25):
        """
        :param trjdir_list: list of directories containing the trajectories
        :param method: quantum mechanical method
        :param simulation_time: simulation time of TSH
        :param time_step: time step employed in TSH
        :param delta_t: time delay to calculate vertical excitations along the trajectories
        :param nstates: number of states to calculate
        :param en_min: minimum energy to consider in the convolution of the TR spectrum
        :param en_max: maximum energy to consider in the convolution of the TR spectrum
        :param en_width: width of the Gaussian used to convolute the spectrum along the energies
        :param grid: number of grid points to convolute the spectrum on along the energies
        :param t_width: width of the Gaussian used to convolute the spectrum along the time

        """

        #Check and define the QM method
        method = (method.upper().replace("-", ""))
        available_methods = ["TDDFT", "CASPT2"]

        if method in available_methods:
            self.method = available_methods.index(method)
        else:
            raise KeyError("The method {} is not available!")

        #Initialise attributes
        self.traj = trjdir_list
        self.simtime = simulation_time
        self.delta_t = delta_t
        self.nsteps = int(simulation_time/time_step)
        self.delta = int(delta_t/time_step)
        self.time_step = time_step
        self.nstates = nstates
        self.min_en = en_min
        self.max_en = en_max
        self.engrid = grid
        self.enwidth = en_width
        self.twidth = t_width
        self.time, self.grid, self.spectrum = [], [], []
        self.matrix = np.zeros(self.engrid)

    ####################################################
    @staticmethod
    def check_folder(folder):
        """ Check the directories where the trajectories are stored """

        ## check if the simulation is terminated normally
        termination_files = ["FINISHED", "cobramm.xml"]

        for file in termination_files:
            if not os.path.exists("{0}/{1}".format(folder, file)):
                raise FileNotFoundError

        ##TODO: add checking functionality in case of CAS calculations

    ####################################################
    def setup_vertical_excitations(self, ncores=1, chrg=0, basis_set="6-31g*", functional="", cluster=False, submission_string="", template=""):
        """Setup of the vertical excitations"""

        init_dir = os.getcwd()
        cobcomm, cobcomm_first, cobcomm_second = "", "", ""
        sander_block = "!sander\ncomment line\n&cntrl\nimin = 1,\nmaxcyc = 0,\nntb = 0,\nigb = 0,\nntr = 0,\n" \
                                    "ibelly = 1,\ncut = 10\n/\n?sander\n"
        #Define the cobram.command for the method chosen
        if self.method == 0:
            software = "gauss"
            qm_string = "!gaussian\n#p {0} {1} nosym tda(nstates={2}) " \
                        "iop(9/40=4) gfinput\n\n" \
                           "gaussian input generated by COBRAMM\n\n" \
                           "{3} 1\n?gaussian".format(basis_set, functional, self.nstates, chrg)
            cobcomm = "!keyword\ntype=optxg\nnsteps=sp\nqm-type={0}\nnproc={1}\nricd=1\n?keyword\n\n" \
                      "{2}}\n{3}".format(software, ncores, sander_block, qm_string)

        elif self.method == 1:

            if not os.path.isfile("INPORB"):
                logwrt.fatalerror("File INPORB containing initial orbitals missing!")
            software = "molcas"

            ### read the template. Format: first molcas block, second molcas block, eventual basis set information

            with open(template, "r") as templ:
                molcas_command = templ.readlines()
            first_command, second_command, bs_block = [], [], []
            new_start = 0
            for index, line in enumerate(molcas_command):
                first_command.append(line)
                if "?molcas" in line:
                    new_start = index+1
                    break
            ## read the basis set
            append = False
            for line in molcas_command:
                if "!basisset" in line:
                    append = True
                if append:
                    bs_block.append(line)
            ## read the second blocK
            second_block = molcas_command[new_start:]
            for line in second_block:
                second_command.append(line)
            ## define the blocks
            first_qm_string, second_qm_string = "", ""
            for command in first_command:
                first_qm_string += str(command)
            for command in second_command:
                second_qm_string += str(command)
            for bs in bs_block:
                first_qm_string += str(bs)
                second_block += str(bs)
            ## save for the command
            cobcomm_first = "!keyword\ntype=optxg\nnsteps=sp\nqm-type={0}\nnproc={1}\nricd=1\nta=1\n?keyword\n\n" \
                                    "{2}\n{3}".format(software, ncores, sander_block, first_qm_string)
            cobcomm_second = "!keyword\ntype=optxg\nnsteps=sp\nqm-type={0}\nnproc={1}\nricd=1\n?keyword\n\n" \
                                     "{2}\n{3}".format(software, ncores, sander_block, second_qm_string)

        #Extract the data for each step
        with open("setup_PP.log", "w") as out:

            # get the steps
            for step in range(0, self.nsteps+1, self.delta):

                if os.path.isdir("inputs_from_step_{}".format(step)):
                    logwrt.writewarning("overwriting previous inputs_from_step_{} directory ".format(step))
                    shutil.rmtree("inputs_from_step_{}".format(step))

                command = shlex.split("cobramm-get-step.py -n {}".format(step))
                subprocess.run(command, stdout=out, stderr=subprocess.STDOUT)

            # change directory name in case last step is included
            if os.path.isdir("inputs_from_last_step"):
                shutil.move("inputs_from_last_step", "inputs_from_step_{}".format(self.nsteps))

            # run the calculations
            for step in range(0, self.nsteps+1, self.delta):

                current_dir = "inputs_from_step_{}".format(step)
                if self.method == 0:
                    os.chdir(current_dir)
                    PumpProbeCalculator.launch_cobramm(cobcomm, submission_string, cluster)
                if self.method == 1:
                    #copy in the step directory and move there
                    shutil.copy("INPORB", current_dir)
                    #create the directory for the first calculation and move there
                    reference_dir = "reference"
                    reference_path = os.path.join(current_dir, reference_dir)
                    shutil.copytree(os.path.abspath(current_dir), os.path.abspath(reference_path))
                    os.chdir(current_dir)
                    upper = os.getcwd()
                    os.chdir(os.path.abspath(reference_dir))
                    PumpProbeCalculator.launch_cobramm(cobcomm_first, submission_string)
                    os.chdir(upper)
                    shutil.copyfile("{}/vanilla.JobMix".format(os.path.abspath(reference_dir)), "vanilla.JobMix")
                    PumpProbeCalculator.launch_cobramm(cobcomm_second, submission_string, cluster)

                os.chdir(init_dir)

    ####################################################
    @staticmethod
    def launch_cobramm(input_commands,  sub_string, cluster=False):
        with open("cobram.command", "w") as inp:
            inp.write(input_commands)
        if cluster:
            with open("qsub.log", "w") as fstdout:
                command = shlex.split("{}".format(sub_string))
                subprocess.run(command, stdout=fstdout, stderr=subprocess.STDOUT)
        else:
            with open("cobramm.log", "w") as fstdout:
                subprocess.call("cobram.py", stdout=fstdout, stderr=fstdout)

    ####################################################
    @staticmethod
    def get_dipoles(method, calcDir):

        dipoles = []
        if method == 0:
            # run MultiWfn
            dipoles = MultiwfnCalculator.extract_TDM_MWFN(calcDir, "restart.chk")
        if method == 1:
            onlydipoles, nstates = PumpProbeCalculator.extract_tdm_molcas(calcDir)
            cas_energies = PumpProbeCalculator.extract_pt2_energies_molcas(calcDir)
            dipoles = PumpProbeCalculator.update_dipoles_cas(nstates, energies=cas_energies, dipoles=onlydipoles)

        return dipoles
    ####################################################
    @staticmethod
    def get_active_states(nstep: int):
        """Extract and return the active state for a step along the trajectory"""

        # collect active states
        with open("cobramm.xml", "r") as f:
            xml_text = f.readlines()

        all_active_states = []
        for index, line in enumerate(xml_text):
            if "<state>" in line:
                state = (xml_text[index+1].split()[0])
                all_active_states.append(int(state))

        active_state = all_active_states[nstep]

        return active_state
    ####################################################
    @staticmethod
    def extract_overlap_molcas(pt2states, logfile="QM_data/molcasALL.log"):

        #get the number of states
        with open("reference/{}".format(logfile), "r") as log:
            ref = log.read()
        state_string = "Number of CI roots used *([0-9]*)"
        scfstates = int(re.findall(state_string, ref)[0])

        #init_dir = os.getcwd()
        #os.chdir(os.path.abspath(calcDir))
        with open(logfile, "r") as log:
            out = log.read()
        #extract the matrix
        re_string = "OVERLAP MATRIX FOR THE ORIGINAL STATES:" + "(.*)" + "I/O STATISTICS"
        find_overlap = (re.findall(re_string, out, re.DOTALL | re.IGNORECASE)[0]).rstrip().lstrip().split()
        find_overlap.pop(-1)

        #calculate the number of elements of the matrix to extract
        scf_scf, elements, to_extract = 0, 0, 0
        for state in range(1, scfstates+1):
            scf_scf += state
        for state in range(scfstates+pt2states+1):
            elements += state
            to_extract = elements-scf_scf
        overlap = (find_overlap[-to_extract:])


        ## create a list with the overlap between S_scf and S_pt2 for each of the S_scf
        states_overlap = []
        index = 0
        for state in range(1, scfstates+1):
            state_ovl = []
            count = index
            for other_state in range(1, pt2states+1):
                state_ovl.append(abs(float(overlap[count])))
                count += other_state + scfstates
            index += 1
            states_overlap.append(state_ovl)
        #os.chdir(init_dir)

        return states_overlap
    ####################################################
    @staticmethod
    def extract_tdm_molcas(calcDir, logfile="QM_data/molcasALL.log"):

        dipoles = []
        init_dir = os.getcwd()
        os.chdir(os.path.abspath(calcDir))

        #extract number of states
        with open(logfile, "r") as log:
            ref = log.read()
        state_string = "Number of CI roots used *([0-9]*)"
        nstates = int(re.findall(state_string, ref)[0])

        #extract TDM1
        sizeblock = (nstates+6)
        nblocks = int(nstates / 4)
        if nstates % 4 != 0:
            nblocks += 1
        totlines = nblocks*3*sizeblock

        with open(logfile, "r") as log:
            out = log.read()
        re_string = "MATRIX ELEMENTS OF 1-ELECTRON OPERATORS" + "(.*)" + "I/O STATISTICS"
        find_tdm = (re.findall(re_string, out, re.DOTALL | re.IGNORECASE)[0]).rstrip().lstrip().splitlines()
        all_tdm = find_tdm[1:totlines+1]

        ### extract components/states
        init_count = 5
        tdm = []
        for comp in range(0,3):
            component = []
            for state in range(1, nstates+1):
                count = init_count+state
                state_tdm = []
                for block in range(nblocks):
                    newcount = count + block*sizeblock
                    splitted = all_tdm[newcount].split()
                    for i in range(1,5):
                        try:
                            state_tdm.append(float(splitted[i]))
                        except IndexError:
                            break
                component.append(state_tdm)
            init_count += sizeblock * nblocks
            tdm.append(component)
        ##reorganize in the same format as when reading Multiwfn
        for from_state in range(nstates):
            for to_state in range(from_state+1, nstates):
                state_to_state = []
                state_to_state.append(from_state)
                state_to_state.append(to_state)
                for comp in range(3):
                    state_to_state.append(tdm[comp][from_state][to_state])
                dipoles.append(state_to_state)

        #ƒormat of dipoles : [ [Si, Sj, x,y,z],[Si, Sj, x,y,z],[Si, Sj, x,y,z] ]

        os.chdir(init_dir)
        return dipoles, nstates
    ####################################################
    @staticmethod
    def extract_pt2_energies_molcas(calcDir, logfile="QM_data/molcasALL.log"):

        init_dir = os.getcwd()
        os.chdir(os.path.abspath(calcDir))
        energies = []

        #extract number of states
        with open(logfile, "r") as log:
            ref = log.read()
        state_string = "Number of CI roots used *([0-9]*)"
        nstates = int(re.findall(state_string, ref)[0])

        ## read the file and extract the caspt2 energies
        with open(logfile, "r") as log:
            out = log.read()
        
        method_string = "Type of calculation" + "(.*)" + "Fock operator"
        find_method = (re.findall(method_string, out, re.DOTALL | re.IGNORECASE)[0]).rstrip().lstrip().splitlines()
        if find_method[0] == "XMS-CASPT2":
            re_string = "Total XMS-CASPT2 energies:" + "(.*)" + "Eigenvectors"
            find_en = (re.findall(re_string, out, re.DOTALL | re.IGNORECASE)[0]).rstrip().lstrip().splitlines()
            final_en = find_en[0:nstates]
        else:
            re_string = "Total CASPT2 energies" + "(.*)" + " A NEW JOBIPH FILE NAMED 'JOBMIX' IS PREPARED"
            find_en = (re.findall(re_string, out, re.DOTALL | re.IGNORECASE)[0]).rstrip().lstrip().splitlines()
            final_en = find_en[-nstates:]
        ## store and return only the energies
        for state in range(nstates):
            energies.append(float((final_en[state])[-13:]))
        gs = energies.pop(0)
        for state in range(nstates-1):
            energies[state] = round((energies[state] - gs) * 27.21, 2)

        os.chdir(init_dir)

        return energies
    ####################################################
    @staticmethod
    def calculate_tdm_magnitude(x, y, z):
        mag = math.sqrt(math.pow(x, 2)+math.pow(y, 2)+math.pow(z, 2))
        return mag

    ####################################################
    @staticmethod
    def calculate_fos(mag, energy):
        ##from energy stored in eV to a.u.
        fos = round(.66*(pow(mag, 2))*(energy/27.21), 2)
        return fos

    ####################################################
    @staticmethod
    def update_dipoles_cas(nstates, energies: list, dipoles: list):

        #append energy and fos for the ground state
        for gs_to_es in range(nstates-1):
            dipoles[gs_to_es].append(round(energies[gs_to_es], 2))
            mag = PumpProbeCalculator.calculate_tdm_magnitude(dipoles[gs_to_es][2], dipoles[gs_to_es][3], dipoles[gs_to_es][4])
            fos = PumpProbeCalculator.calculate_fos(mag, dipoles[gs_to_es][-1])
            dipoles[gs_to_es].append(fos)

        # append energy and fos for the es to es transitions
        for es_to_es in range(nstates-1,len(dipoles)):
            from_state = dipoles[es_to_es][0]
            to_state = dipoles[es_to_es][1]
            energy = round((energies[to_state-1]-energies[from_state-1]), 2)
            dipoles[es_to_es].append(energy)
            mag = PumpProbeCalculator.calculate_tdm_magnitude(dipoles[es_to_es][2], dipoles[es_to_es][3], dipoles[es_to_es][4])
            fos = PumpProbeCalculator.calculate_fos(mag, dipoles[es_to_es][-1])
            dipoles[es_to_es].append(fos)

        return dipoles
        
    ####################################################
    #@staticmethod
    #def extract_en_fos(calcDir, nstates: int, active_state: int):
    #    """Extract TDM between excited states.
    #    This function is an old version of "extract_en_fos_polarization" and will be deleted"""
    #
    #    energies, intensities, dipoles = [], [], []
    #    #extract dipoles and energies for different methods
    #    if self.method == 0:
    #        #run MultiWfn
    #        dipoles = MultiwfnCalculator.extract_TDM_MWFN(calcDir, "restart.chk")
    #    if self.method == 1:
    #        onlydipoles = PumpProbeCalculator.extract_tdm_molcas()
    #        cas_energies = PumpProbeCalculator.extract_energies_molcas(nstates, calcDir)
    #        dipoles = PumpProbeCalculator.update_dipoles_cas(nstates, energies=cas_energies, dipoles=onlydipoles, fos=fos)
    #    #get only energies and fos
    #    energies.append(dipoles[active_state - 1][5])
    #    intensities.append(-float(dipoles[active_state - 1][6]))
    #    #store fos
    #    first_item = 0
    #    for i in range(0, active_state):
    #        first_item += nstates-i
    #    last_item = first_item + (nstates-active_state)
    #    for state in range(first_item, last_item):
    #        energies.append(dipoles[state][5])
    #        intensities.append(dipoles[state][6])
    #
    #    return energies, intensities

    ####################################################
    def extract_init_tdm(self):

        init_dir = os.getcwd()
        all_init_tdm = []
        if self.method == 0:

            for folder in self.traj:
                os.chdir(folder)
                active_state = PumpProbeCalculator.get_active_states(0)
                os.chdir("inputs_from_step_0")
                try:
                    dipoles = MultiwfnCalculator.extract_TDM_MWFN(".", "restart.chk")
                    init_tdm = [dipoles[active_state-1][2], dipoles[active_state-1][3], dipoles[active_state-1][4]]
                    all_init_tdm.append(init_tdm)
                except:
                    skip = ""
                    all_init_tdm.append(skip)

                os.chdir(init_dir)

        return all_init_tdm
    ####################################################
    @staticmethod
    def extract_en_fos_polarization(method, calcDir, nstates: int, active_state: int, polarization, init_tdm=[]):
        """Extract TDM between excited states"""

        energies, intensities = [], []
        #
        ##extract dipoles and energies for different methods
        #if self.method == 0:
        #    #run MultiWfn
        #    dipoles = MultiwfnCalculator.extract_TDM_MWFN(calcDir, "restart.chk")
        #if self.method == 1:
        #    onlydipoles = PumpProbeCalculator.extract_tdm_molcas(nstates, calcDir)
        #    cas_energies = PumpProbeCalculator.extract_energies_molcas(nstates, calcDir)
        #    dipoles = PumpProbeCalculator.update_dipoles_cas(nstates, energies=cas_energies, dipoles=onlydipoles)

        dipoles = PumpProbeCalculator.get_dipoles(method, calcDir=calcDir)

        ### store SE
        energies.append(dipoles[active_state - 1][5])
        if polarization == "parallel":
            current_tdm = [dipoles[active_state-1][2],dipoles[active_state-1][3],dipoles[active_state-1][4]]
            cos = np.dot(init_tdm, current_tdm)/(np.linalg.norm(init_tdm)*np.linalg.norm(current_tdm))
            pol_coeff = (1+(2*(cos*cos)))
            intensities.append(-float(dipoles[active_state-1][6])*pol_coeff)
        elif polarization == "orthogonal":
            current_tdm = [dipoles[active_state-1][2],dipoles[active_state-1][3],dipoles[active_state-1][4]]
            cos = np.dot(init_tdm, current_tdm)/(np.linalg.norm(init_tdm)*np.linalg.norm(current_tdm))
            pol_coeff = (2-(cos*cos))
            intensities.append(-float(dipoles[active_state-1][6])*pol_coeff)
        elif polarization == "none":
            intensities.append(-float(dipoles[active_state-1][6]))

        #store ESA
        first_item = 0
        for i in range(0, active_state):
            first_item += nstates-i
        last_item = first_item + (nstates-active_state)
        for state in range(first_item, last_item):
            energies.append(dipoles[state][5])
            if polarization == "parallel":
                current_tdm = [dipoles[state][2],dipoles[state][3],dipoles[state][4]]
                cos = np.dot(init_tdm, current_tdm)/(np.linalg.norm(init_tdm)*np.linalg.norm(current_tdm))
                pol_coeff = (1+2*(cos*cos))
                intensities.append(pol_coeff*dipoles[state][6])
            elif polarization == "orthogonal":
                current_tdm = [dipoles[state][2],dipoles[state][3],dipoles[state][4]]
                cos = np.dot(init_tdm,current_tdm)/(np.linalg.norm(init_tdm)*np.linalg.norm(current_tdm))
                pol_coeff = (2-(cos*cos))
                intensities.append(pol_coeff*dipoles[state][6])
            elif polarization == "none":
                intensities.append(dipoles[state][6])

        return energies, intensities

    ####################################################
    @staticmethod
    def decompose_esa(calcDir, nstates:int, active_state: int, target_state: int):
        """Extract TDM between excited states"""

        energies, intensities, dipoles = [], [], []
        if self.method == 0:
            #run MultiWfn
            dipoles = MultiwfnCalculator.extract_TDM_MWFN(calcDir, "restart.chk")
        if self.method == 1:
            onlydipoles = PumpProbeCalculator.extract_tdm_molcas(calcDir)
            cas_energies = PumpProbeCalculator.extract_pt2_energies_molcas(calcDir)
            dipoles = PumpProbeCalculator.update_dipoles_cas(nstates, energies=cas_energies, dipoles=onlydipoles)

        #store TDM
        first_item = 0
        for i in range(0, active_state):
            first_item += nstates-i
        delta = (target_state-active_state)
        energies.append(dipoles[first_item+delta-1][5])
        intensities.append(dipoles[first_item + delta - 1][6])

        return energies, intensities

    ####################################################
    def collect_values_single_time(self, nstep:int, polarization="none", all_init_tdm=[]):
        """Collect the spectrum values for a certain time"""

        init_dir = os.getcwd()
        collected_energies = []
        collected_intensities = []

        #get the values
        for folder in self.traj:
            print("Entering in folder {}".format(folder))
            os.chdir(folder)

            active_state = PumpProbeCalculator.get_active_states(nstep)
            ###check order of the state with xms

            os.chdir("inputs_from_step_{}".format(nstep))
            if polarization == "parallel" or polarization == "orthogonal":
                init_tdm = all_init_tdm[self.traj.index(folder)]
            elif polarization == "none":
                init_tdm = []
            #try:  ##### to speed up tests
            if self.method == 1:
                ovlp = PumpProbeCalculator.extract_overlap_molcas(self.nstates)
                state_ovl = ovlp[active_state]
                new_active_state = state_ovl.index(max(state_ovl))


            en, fos = PumpProbeCalculator.extract_en_fos_polarization(self.method, ".", self.nstates, new_active_state,
                                                                                    polarization=polarization,init_tdm=init_tdm)
            for state in range(len(en)):
                collected_energies.append(en[state])
                collected_intensities.append(fos[state])
            #except:
            #    pass

            os.chdir(init_dir)
        return collected_energies, collected_intensities

    ####################################################
    def collect_decomposed_values_single_time(self, nstep:int,selected_state,target_state):
        """Collect the spectrum values for a certain time"""

        init_dir = os.getcwd()
        collected_energies, collected_intensities = [], []

        #get the values
        for folder in self.traj:
            os.chdir(folder)
            active_state = PumpProbeCalculator.get_active_states(nstep)
            if active_state == selected_state:
                os.chdir("inputs_from_step_{}".format(nstep))
                try: ##### to speed up tests
                    en, fos = PumpProbeCalculator.decompose_esa(".", self.nstates, active_state, target_state)
                    for state in range(len(en)):
                        collected_energies.append(en[state])
                        collected_intensities.append(fos[state])
                except:
                    pass
            os.chdir(init_dir)

        return collected_energies, collected_intensities

    ####################################################
    def get_spectrum_single_time_all_traj(self, values, currtime=0.5):
        """Convolute the spectrum values for a certain time"""

        #define grid
        spectral_grid = np.linspace(self.min_en, self.max_en, self.engrid)
        spec_value = np.zeros(len(spectral_grid))

        #convolute spectrum
        for center, intensity in zip(*values):
            linef = np.array([1. / (self.enwidth * math.sqrt(2.0 * math.pi)) * np.exp(-(e - center) ** 2 / (2.0 * self.enwidth ** 2))
                              for e in spectral_grid])
            to_add = intensity*linef
            spec_value += to_add

        #store in matrix format
        self.matrix = np.vstack([self.matrix,spec_value])

        #store in list format
        for en, ints_e in zip(spectral_grid, spec_value):
            self.time.append(currtime)
            self.grid.append(en)
            self.spectrum.append(ints_e) #/ wavelength**2)


    ####################################################
    def get_final_spectrum(self, output_name: str):
        """Write the final file"""

        #write in list format
        comment = "#\n" \
                  "# time(fs)     energy(eV)    intensity(a.u.)\n" \
                  "#\n"

        with open(output_name, "w") as f:
            f.write(comment)
            index = 1
            for c_time, wav, ints in zip(self.time, self.grid, self.spectrum):
                if index == self.engrid:
                    f.write("  \n")
                    index = 1
                else:
                    f.write("{0:.4f} {1:16.4f} {2:16.4f}\n".format(c_time, wav, ints))
                    index += 1

        #write in a matrix format
        #with open("{}_matrix".format(output_name), "w") as m:
        #    np.savetxt(m, np.transpose(np.divide(self.matrix,100)), fmt='%16.12f')

    #####################################################
    def time_convolution(self, output="time_convolution_matrix.txt"):
        """Convolute along time the spectrum obtained"""

        #set up the grid
        index = 1
        temp_en_mat = np.zeros(self.simtime+1)
        temp_matrix = np.zeros([self.engrid, self.simtime+1])

        #convolute for each spectrum value
        for c_time, wav, ints in zip(self.time, self.grid, self.spectrum):

            tconv_time, tconv_wav, tconv_spec = [], [], []

            #add zero values to consistency of matrices' formats
            for t in range(0, int(c_time-self.twidth)):
                tconv_time.append(t)
                tconv_wav.append(wav)
                tconv_spec.append(0.0000)

            #convolute
            tgrid = np.linspace(c_time-self.twidth, c_time+self.twidth, self.twidth*2+1)
            spec_value = np.zeros(len(tgrid))

            gau = np.array([np.exp( (- (t-c_time)**2) / ( (2*(self.twidth/2.355))**2) ) for t in tgrid])

            to_add = gau*ints
            spec_value += to_add

            #add zero values to consistency of matrices' formats
            for t, ints_e in zip(tgrid, spec_value):
                if t < 0 or t > self.simtime:
                    pass
                else:
                    tconv_time.append(t)
                    tconv_wav.append(wav)
                    tconv_spec.append(ints_e)
            for t in range(int(c_time+self.twidth), self.simtime):
                tconv_time.append(t)
                tconv_wav.append(wav)
                tconv_spec.append(0.0000)

            #store in matrix format
            sp = np.array(tconv_spec)
            if index == 1:
                temp_en_mat = np.add(temp_en_mat, sp)
            else:
                temp_en_mat = np.vstack([temp_en_mat, sp])

            if index == self.engrid:
                temp_matrix = np.add(temp_matrix, temp_en_mat)
                index = 1
                temp_en_mat = np.zeros(self.simtime+1)
            else:
                index += 1

        time_zero = np.zeros([self.engrid,1])
        final_mat = np.hstack([time_zero,temp_matrix])

        #Write matrix of the spectrum contribution given by the Gaussian obtained
        with open(output, "w") as m:
            np.savetxt(m, np.divide(final_mat,self.engrid*self.simtime), fmt='%16.12f')

    #####################################################
    def write_gnuplot(self, toplot="time_convolution_spectrum.txt", output="plot_matrix.gp"):

        gp_commands = "set title 'Time Resolved PP Spectrum eV, fs'\n" \
                      "set yrange [{0}:{1}]\n" \
                      "set xrange [*:*]\n" \
                      "set ylabel 'energy (eV)'\n" \
                      "set xlabel 'time (fs)'\n" \
                      "set paletted defined (-1 '#00008B', -0.5 'blue', 0.5 'cyan', 1. 'yellow', 1.5 'orange', 2 'red', 2.5 'brown')\n" \
                      "set pm3d map\n" \
                      "sp '{2}' u 1:($2*{3}:$3) matrix title ''\n"\
                      "\n"\
                      "pause -1\n" \
                      "".format(self.min_en, self.max_en, toplot, (self.max_en-self.min_en)/self.engrid)

        with open(output, "w") as gp:
            gp.write(gp_commands)

#####################################################################################################
class SpectronCalculator:
    spectron_dir = ""
    multiwfn_dir = ""
    _SpectronCheck = False

    def __init__(self):
        """ the constructor of the SpectronCalculator class checks if the environment for SPECTRON execution
            is properly defined and if all the necessary executables are available """

        if not SpectronCalculator._SpectronCheck:

            # environmental variable with the COBRAMM path
            try:
                SpectronCalculator.spectron_dir = os.environ['SPECTRON_DIR']
            except KeyError:
                raise SpectronError("environment variable $SPECTRON_DIR is not defined")

            # check if COBRAMM executables are available
            for exe in ["spectron2", "iSPECTRONa.py", "iSPECTRONb.py"]:
                if not cobrammenv.which(exe):
                    raise SpectronError(exe + " executable is not available")
            SpectronCalculator._spectrumCheck = True

    ####################################################

    @staticmethod
    def runiSpectronA(tdm_file="dipoles.txt", energies_file="energies.txt", grad_files="gradient.txt",
                      modes="frequencies.txt", free=2, w1i=10000, w1f=60000, nw1=1000, signal="PP", tag="test",
                      diagrams="ESA GSB SE", t2=0):
        """A function to run iSPECTRONa.py to create the input file for Spectron

        :param tdm_file: input file containing dipoles
        :param energies_file: input file containing energies
        :param grad_files: input file containing gradients
        :param modes: input file containing normal modes and frequencies
        :param free: number of states
        :param w1i: Initial frequency (cm^-1)
        :param w1f: Final frequency (cm^-1)
        :param nw1: Number of samples along W1
        :param signal: Signal type
        :param tag: tag for the folder's name
        :return: output directory where to run Spectron
        """

        # defining output directory

        outputDir = signal + tag
        if os.path.isdir(outputDir):
            logwrt.writewarning("overwriting previous " + outputDir + " directory ")
            shutil.rmtree(outputDir)

        # running iSPECTRONa.py

        with open("iSpectronA.log", "w") as out:
            command = shlex.split(
                "iSPECTRONa.py {0} {1} {2} {3} -free {4} -w1i {5} -w1f {6} -nw1 {7} -sig {8} -tag {9} -diagrams {10} -t2 {11}".
                format(tdm_file, energies_file, grad_files, modes, free, w1i, w1f, nw1, signal, tag, diagrams, t2))
            subprocess.run(command, stdout=out, stderr=subprocess.STDOUT)

        return outputDir

    ####################################################

    @staticmethod
    def runSpectron(calcDir,
                    inputfile="input.com"):  # , simulation_time=0, time_step=0): #overwrite=False, nCores=1, store=False):
        """A function to run Spectron

        :param calcDir: directory where to run Spectron
        :param inputfile: input file for Spectron
        """

        start_dir = os.getcwd()

        # run spectron
        os.chdir(calcDir)
        with open("spectron.log", "w") as out:
            command = shlex.split("spectron2 -i {} ".format(inputfile))
            subprocess.run(command, stdout=out, stderr=subprocess.STDOUT)

        os.chdir(start_dir)

    ####################################################
    @staticmethod
    def run_pp(inpDir, outDir="2D_data", inputfile="input.com", simulation_time=int, time_step=int):

        if os.path.isdir(outDir):
            logwrt.writewarning("overwriting previous " + outDir + " directory ")
            shutil.rmtree(outDir)

        startDir = os.getcwd()
        workDir = os.path.abspath(outDir)
        pathfile = os.path.abspath(inpDir)

        # create a folder for each time step
        for step in range(0, simulation_time + 1, time_step):
            local_dir = shutil.copytree(pathfile, "{0}/t2_{1}".format(workDir, step))
            os.chdir(local_dir)

            # change the current time step
            with open(inputfile, "r") as f:
                input = f.read()

            input = input.replace("DEL_TIME2 0", "DEL_TIME2 {}".format(step))

            with open(inputfile, "w") as f:
                f.write(input)

            SpectronCalculator.runSpectron(calcDir=".")

            os.chdir(workDir)

        os.chdir(startDir)

        return outDir

    ####################################################
    @staticmethod
    def runiSpectronB(logdir, units="eV", sig2plot="PPheatmap"):

        outputDir = "result/{0}_{1}".format(sig2plot, logdir)
        if os.path.isdir(outputDir):
            logwrt.writewarning("overwriting previous " + outputDir + " directory ")
            shutil.rmtree(outputDir)
        # start_dir = os.getcwd()
        # os.chdir(logdir)

        with open("iSPECTRONb.log", "w") as out:
            command = shlex.split("iSPECTRONb.py {0} -units {1} -sig2plot {2}".
                                  format(logdir, units, sig2plot))
            subprocess.run(command, stdout=out, stderr=subprocess.STDOUT)

        return outputDir

    ###################################################

    @staticmethod
    def insert_monoexp_k(workDir, decay_time="slow", tau=0.0001, nstates=int, active_state=int, final_state=0):

        # if decay_time != "slow" or "medium" or "fast" or "custom":
        #    raise TypeError("Decay time must be 'slow', 'medium', 'fast' or 'custom' ")

        rate_constants = {"slow": 0.001, "medium": 0.02, "fast": 0.01, "custom": 1 / tau}

        KS = np.zeros((nstates, nstates))
        KS[final_state, active_state] = float(rate_constants[decay_time])
        KS[active_state, active_state] = -float(rate_constants[decay_time])

        np.savetxt("{}/transport_rates.txt".format(os.path.abspath(workDir)), KS, fmt='%1.3f')


#####################################################################################################

class readCobrammOut:
    """ The readCobrammOut class defines the object that store and print output informations collected from a Cobramm output file"""

    def __init__(self, energies=False, en_dir="", gs_dipoles=False, gs_dip_dir="", gradient=False, grad_dir="",
                 es_dipoles=False, tda_dir="", storeFiles=False):

        self.storeFiles = True

        if energies:
            self.energies = self.extract_energies_gau(os.path.join(en_dir, "QM_data", "qmALL.log"))
        else:
            self.energies = None

        if gs_dipoles:
            self.gs_dipoles, self.nstates = self.extract_TDM_gau(os.path.join(gs_dip_dir, "QM_data", "qmALL.log"))
        else:
            self.gs_dipoles = None

        if gradient:
            self.grad, self.natom, self.nroot, self.atom_list = self.extract_grads_gau(
                os.path.join(grad_dir, "QM_data", "qmALL.log"))
        else:
            self.grad = None

        if es_dipoles:
            self.gs_dipoles = None
            self.es_dipoles = self.extract_ES_TDM_gau(tda_dir)
        else:
            self.es_dipoles = None

    ####################################################
    def write_files(self, en_out="energies.txt", dip_out="dipoles.txt", grad_out="gradient"):

        if self.energies:
            with open(en_out, "w") as en:
                en.write("free format - energies\n0\n")
                for state in self.energies:
                    en.write("{}\n".format(int(state)))

        if self.gs_dipoles:
            with open(dip_out, "w") as tdm:
                tdm.write("free format - dipoles\n")
                for index in range(self.nstates):
                    state = self.gs_dipoles[index]
                    tdm.write("0 {0} {1} {2} {3}\n".format(index + 1, state[0], state[1], state[2]))

        if self.es_dipoles:
            with open(dip_out, "w") as tdm:
                tdm.write("free format - dipoles\n")
                for index in range(len(self.es_dipoles)):
                    dipole = self.es_dipoles[index]
                    tdm.write(
                        "{0} {1} {2} {3} {4}\n".format(int(dipole[0]), int(dipole[1]), dipole[2], dipole[3], dipole[4]))

        if self.grad:
            with open("{0}_S{1}.txt".format(grad_out, str(self.nroot)), "w") as gradfile:
                gradfile.write("free format - gradient\ngradient root {}\n".format(self.nroot + 1))
                for index in range(self.natom):
                    gradfile.write("{0} {1:.6f} {2:.6f} {3:.6f}\n"
                                   .format(constants.atomic_numbers[self.atom_list[index]], self.grad[0][index],
                                           self.grad[1][index], self.grad[2][index]))

                # return en_out, dip_out, gra_dout

    ####################################################

    @staticmethod
    def extract_energies_gau(file_name):
        """ use regexpr to extract the list of electronic state energies"""

        # open output file and read content
        with open(file_name, "r") as input_file:
            out_text = input_file.read()

        # excited state string
        re_string = "Excited State *[0-9]*: *.*? *([0-9]*\.[0-9]*) *eV *[0-9]*\.[0-9]* *nm *" \
                    "f=([0-9]*\.[0-9]*) *\<S\*\*2\>=[0-9]*\.[0-9]*"
        # extract results with regular expression
        results = re.findall(re_string, out_text)

        # define and return lists of electronic state energies (convert from eV to cm^-1)
        el_energies = [float(i[0]) * 8065.5443 for i in results]

        return el_energies

    ####################################################
    @staticmethod
    def extract_TDM_gau(file_name):
        """ extract the list of transition dipole moments"""

        # extract electric transition dipole moment strings
        tdm_string = "Ground to excited state transition electric dipole moments (Au)"

        # open output file and read content
        with open(file_name, "r") as f:
            out = f.read()
            f.seek(0)
            out_text = f.readlines()

        re_string = "nstates=*([0-9]*)"
        find_root = re.findall(re_string, out)
        nstates = int(find_root[0])

        # store dipoles
        dipoles = []
        for index, line in enumerate(out_text):
            if tdm_string in line:
                for i in range(nstates):
                    results = (out_text[index + 2 + i].split()[1:4])
                    dipoles.append([float(j) for j in results])

        return dipoles, nstates

    ####################################################
    @staticmethod
    def extract_ES_TDM_gau(calcDir):
        """ extract the list of transition dipole moments between excited states"""

        dipoles = MultiwfnCalculator.extract_TDM_MWFN(calcDir, "restart.chk")

        return dipoles

    ####################################################

    @staticmethod
    def get_bright_state(file_name):
        """ use regexpr to extract the brightest state"""

        # open output file and read content
        with open(file_name, "r") as input_file:
            out_text = input_file.read()

        # excited state string
        re_string = "Excited State *[0-9]*: *.*? *([0-9]*\.[0-9]*) *eV *[0-9]*\.[0-9]* *nm *" \
                    "f=([0-9]*\.[0-9]*) *\<S\*\*2\>=[0-9]*\.[0-9]*"
        # extract results with regular expression
        results = re.findall(re_string, out_text)

        # get the oscillator strengths and extract the index for the brightest state
        fos = [float(i[1]) for i in results]
        bright = fos.index(max(fos)) + 1

        return bright

    ####################################################
    @staticmethod
    def extract_grads_gau(file_name):

        # open output file and read content
        with open(file_name, "r") as grad:
            out = grad.read()
        output = out.splitlines()

        # extract which root
        re_string = "nstates=*[0-9]*, root=*([0-9]*)"
        find_root = re.findall(re_string, out)
        nroot = int(find_root[0])

        # extract gradient, number of atom and atom list in atomic numbers
        gradient = [[], [], []]
        atom_list = []
        natom = 0
        for i in range(len(output)):

            if output[i].strip() == 'Center     Atomic                   Forces (Hartrees/Bohr)':
                while True:
                    try:
                        element = output[i + natom + 3].split()
                        atom_list.append(int(element[1]))
                        gradient[0].append(-float(element[2]))
                        gradient[1].append(-float(element[3]))
                        gradient[2].append(-float(element[4]))
                    except (IndexError, ValueError):
                        break
                    natom += 1

        return gradient, natom, nroot, atom_list

############################################################################################################################################################
class MultiwfnCalculator:

    multiwfn_dir = ""
    _multiwfnCheck = False

    def __init__(self):
        """ the constructor of the SpectronCalculator class checks if the environment for SPECTRON execution
            is properly defined and if all the necessary executables are available """

        if not MultiwfnCalculator._multiwfnCheck:

        # check multiwfn and setup
            try:
                MultiwfnCalculator.multiwfn_dir = os.environ['Multiwfnpath']
            except KeyError:
                raise MultiwfnError("environment variable $Multiwfnpath is not defined")

            if not cobrammenv.which("Multiwfn"):
                raise MultiwfnError("Multiwfn executable is not available")
            MultiwfnCalculator._multiwfnCheck = True

    ##########################################################
    @staticmethod
    def extract_TDM_MWFN(calcDir, checkfile="restart.chk"):

        init_dir = os.getcwd()
        os.chdir(os.path.abspath(calcDir))
        logfile = "QM_data/qmALL.log"
        # find number of excited states
        with open(logfile, "r") as f:
            log = f.read()

        re_string = "nstates=*([0-9]*)"
        find_root = re.findall(re_string, log)
        nexstates = int(find_root[0])

        if os.path.isfile("transdipmom.txt") == False:
            # write input for mwfn
            input_string = "{0}\n18\n5\n{1}\n2\n0\nq".format(checkfile, logfile)
            with open("input_multiwfn", "w") as f:
                f.write(input_string)
                f.close()
            # run multiwfn and check
            os.system("Multiwfn < input_multiwfn > multiwfn.log")

        try:
            with open("transdipmom.txt", "r") as out:
                mwfnout = out.readlines()
        except IOError:
            raise MultiwfnError("\nSomething went wrong with Multiwfn...")

        groundstate = "Transition dipole moment between ground state (0) " \
                      "and excited states (a.u.)"
        excitedstates = "Transition dipole moment between excited states (a.u.):"
        dipoles = []
        for index, line in enumerate(mwfnout):
            if groundstate in line:
                for state in range(nexstates):
                    results = (mwfnout[index + 2 + state].split()[0:7])
                    dipoles.append([float(j) for j in results])
            if excitedstates in line:
                newindex = int(index) + 3
                for state in range(1, nexstates):
                    remaining_states = nexstates - state
                    for i in range(remaining_states):
                        results = (mwfnout[newindex].split()[0:7])
                        dipoles.append([float(j) for j in results])
                        newindex += 1
                    newindex += 1

        os.chdir(init_dir)

        return dipoles

    ##########################################################

    @staticmethod
    def read_TDM_MWFN(calcDir, checkfile="restart.chk"):

        init_dir = os.getcwd()
        os.chdir(os.path.abspath(calcDir))

            # define files to read
        logfile = "QM_data/qmALL.log"

            # find number of excited states

        with open(logfile, "r") as f:
            log = f.read()

        re_string = "nstates=*([0-9]*)"
        find_root = re.findall(re_string, log)
        nexstates = int(find_root[0])

        if os.path.isfile("transdipmom.txt") == False:
            # write input for mwfn
            input_string = "{0}\n18\n5\n{1}\n2\n0\nq".format(checkfile, logfile)
            with open("input_multiwfn", "w") as f:
                f.write(input_string)
                f.close()

            # run multiwfn and check
            os.system("Multiwfn < input_multiwfn > multiwfn.log")

        try:
            with open("transdipmom.txt", "r") as out:
                mwfnout = out.readlines()
        except IOError:
            raise MultiwfnError("\nSomething went wrong with Multiwfn...")

        groundstate = "Transition dipole moment between ground state (0) " \
                      "and excited states (a.u.)"
        excitedstates = "Transition dipole moment between excited states (a.u.):"
        dipoles = []
        for index, line in enumerate(mwfnout):
            if groundstate in line:
                for state in range(nexstates):
                    results = (mwfnout[index + 2 + state].split()[0:7])
                    dipoles.append([float(j) for j in results])
            if excitedstates in line:
                newindex = int(index) + 3
                for state in range(1, nexstates):
                    remaining_states = nexstates - state
                    for i in range(remaining_states):
                        results = (mwfnout[newindex].split()[0:7])
                        dipoles.append([float(j) for j in results])
                        newindex += 1
                    newindex += 1

        os.chdir(init_dir)

        return dipoles

    ##########################################################
    @staticmethod
    def read_resp_charges(calcDir, checkfile="restart.chk",state="gs"):

        init_dir = os.getcwd()
        os.chdir(os.path.abspath(calcDir))

        logfile = "QM_data/qmALL.log"

            # find number of excited states

        if state == "gs":

            if os.path.isfile("transdipmom.txt") == False:
                # write input for mwfn
                input_string = "{0}\n7\n18\n8\n1\n{1}\ny\n0\n0\nq".format(checkfile, logfile)
                with open("input_multiwfn", "w") as f:
                    f.write(input_string)
                    f.close()

                # run multiwfn and check
                os.system("Multiwfn < input_multiwfn > multiwfn.log")

            try:
                with open("restart.chg", "r") as out:
                    mwfnout = out.read().splitlines()
            except IOError:
                raise MultiwfnError("\nSomething went wrong with Multiwfn...")

            charges = []
            for q in range(len(mwfnout)):
                charges.append(mwfnout[q][40:51])

            os.chdir(init_dir)

        return charges
