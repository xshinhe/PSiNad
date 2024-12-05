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

# standard library
import sys                      # operating system utilities
import os                       # system-specific parameters and functions
import re                       # regular expressions
import subprocess               # spawn external processes
import pickle                   # persistence
import copy                     # deep copy utility
# external modules
import numpy as np              # numpy library for scientific computation
# cobramm local modules
from layers import Layers       # Layers class to manage geometries
import logwrt


class Output:
    """
    Object class to read/write cobramm.xml output file
    """

    def __init__(self, filename='cobramm.xml', parse=False):
        """
        Constructor method: init class attributes when reading an output
        Only the filename attribute is initialized when writing to a file
        :param filename: output file name
        :param parse: boolean, True: read file, False: write file (defaults to False)
        """
        self.filename = filename
        self.step = 0

        if not parse:
            # check if file already exists
            if os.path.isfile(self.filename):
                # check if this is a restart
                with open(self.filename, 'rb') as f:
                    # read last 2 lines
                    out_tail = f.readlines()[-2:]
                    if b'RESTART' in out_tail[1]:
                        self.step = int(out_tail[0].strip().split(b' ')[-1][:-1]) + 1
                    else:
                        logwrt.fatalerror('Old output file found in the working directory, '
                                          'but this is not a restart.')
            else:
                return

        # if parsing an existing output file, init empty variables
        self.raw_data = ''
        self.steps = 0
        self.version = ''
        self.calculation_type = ''
        self.optimization_type = ''
        self.qm_type = ''
        self.atom_labels = None
        self.atom_layers = None
        self.atom_links = None
        self.modelH_top = ''
        self.real_top = ''
        self.geometry = None
        self.charge = None
        self.chargeMM = None
        self.chargeQM = None
        self.normal_modes = None

        # start parsing
        self.parse_init()


    def parse_init(self):
        """
        Method needed to parse cobramm.xml file and populate class attributes
        :return: None
        """
        # read output data and total number of steps
        with open(self.filename, 'rb') as f:
            self.raw_data = f.read()
            # check for RESTART in last line
            if b'RESTART' in self.raw_data[-25:]:
                # read second-last line
                tail_response = subprocess.Popen('tail -n2 {0}| head -n1'.format(self.filename),
                                                 shell=True,
                                                 stdout=subprocess.PIPE)
            else:
                # read last line
                tail_response = subprocess.Popen(['tail', '-n1', self.filename],
                                                 stdout=subprocess.PIPE)
            tail_response.wait()
            try:
                self.steps = int( tail_response.stdout.read().split(b' ')[-1].strip()[:-1] )
            except ValueError:
                self.steps = -1
            tail_response.stdout.close()

        init_data = self.grep('INITIAL DATA', self.raw_data)
        if init_data is None:
            raise IOError('No text in file')

        # for string fields
        parse_keys = ['version', 'calculation_type', 'optimization_type', 'qm_type']
        # for pickle dumped objects
        pickle_keys = ['atom_labels', 'atom_layers', 'atom_links', 'geometry',
                       'chargeMM', 'chargeQM', 'modelH_top', 'real_top']

        # populate class attributes
        for k in parse_keys:
            exec("self.{0} = self.grep('{0}', init_data)".format(k))
        for pk in pickle_keys:
            result = self.grep('{0}'.format(pk), init_data)
            if result:
                try:
                    exec("self.{0} = pickle.loads(eval(result),encoding='bytes')".format(pk))
                except TypeError:
                    exec("self.{0} = pickle.loads(eval(result).encode(),encoding='bytes')".format(pk))
            else:
                exec("self.{0} = None".format(pk))

        for k in parse_keys + ['modelH_top', 'real_top']:
            try:
                exec("self.{0} = self.{0}.decode()".format(k))
            except (TypeError, AttributeError):
                pass
        for pk in ['atom_labels', 'atom_layers']:
            try:
                exec("self.{0} = [a.decode() for a in self.{0}]".format(pk))
            except (TypeError, AttributeError):
                pass

        # rebuild geometry as Layer object (this block will be moved in layers.py)
        X = self.geometry[0]
        Y = self.geometry[1]
        Z = self.geometry[2]
        geom_stream = ['{0} {1} {2} {3} {4} {5}'.format(self.atom_labels[i],
                                                        0,
                                                        X[i],
                                                        Y[i],
                                                        Z[i],
                                                        self.atom_layers[i]) for i in range(len(self.atom_labels))]

        # change atom_links data type from float to int
        self.atom_links = self.atom_links.astype('int64') if np.any(self.atom_links) else [[],[]]
        # append atom links (both BA and DC types) to geom_stream
        for i in range(len(self.atom_links[0])):
            geom_stream[self.atom_links[0][i] - 1] += ' {0} {1}'.format(self.atom_layers[self.atom_links[1][i] - 1],
                                                                        self.atom_links[1][i])
        self.geometry = Layers.from_real_layers_xyz(geom_stream)

        if self.optimization_type in ['freqxg', 'freqxgp']:
            normal_modes_data = self.grep('FREQ', self.raw_data)
            self.normal_modes = pickle.loads(eval(normal_modes_data))


    def freq_to_string(self):
        """
        String representation of self.normal_modes, used by cobramm-human-readable
        :return: string
        """
        NM = self.normal_modes
        freqs = sorted(NM.keys())
        text = ''
        i = 0
        while i < len(freqs):
            text += ''.join(['{0:>12.4f}'.format(f) for f in freqs[i:i+5]])
            text += '\n'
            nms = [NM[f] for f in freqs[i:i+5]]
            tnms = np.transpose(nms)
            for line in tnms:
                text += ''.join(['{0:>12.5f}'.format(nm) for nm in line])
                text += '\n'
            text += '\n'
            i += 5
        return text


    def get_step(self, step):
        """
        Method to fetch single step info
        :param step: int, step number
        :return: Step object
        """
        # check if requested step exists, if not raise IndexError
        if int(step) > self.steps:  # or self.steps == -1:
            raise IndexError

        # fetch raw data string for requested step
        step_string = self.grep('STEP {0}'.format(step), self.raw_data)
        # geometry data passed to Step constructor, needed to rebuild a Layer object
        geom_data = [self.atom_labels, self.atom_layers, self.atom_links]

        return Step(step_string, geom_data)


    @staticmethod
    def grep(grep_key, selection):
        """
        Static method to retrieve data from XML-style cobramm.xml
        :param grep_key: string, data key to fetch
        :param selection: string, data blob to be parsed
        :return: re.search result (string) or None
        """
        pattern = r'<{0}>\s+(.*?)\s+</{0}>'.format(grep_key).encode()
        result = re.search(pattern, selection, re.DOTALL)
        return result.group(1).strip() if result else None


    def grep_all(self, grep_key):
        """
        This method fetch all occurrences of "grep_key" in cobramm.xml
        and loads them as objects
        :param grep_key: string
        :return: list of objects
        """
        pattern = r'<{0}>\s+(.*?)\s+</{0}>'.format(grep_key).encode()
        results = re.findall(pattern, self.raw_data, re.DOTALL)

        # load results objects
        if grep_key in ['step', 'nIRCstep', 'state', 'optstate']:
            results = [int(res) for res in results]
        else:
            try:
                results = [pickle.loads(eval(res), encoding='bytes') for res in results]
            except TypeError:
                results = [pickle.loads(eval(res).encode(), encoding='bytes') for res in results]

        # deepcopy geometry objects and update cartesian
        if grep_key == 'geometry':
            _results = []
            for res in results:
                copy_geometry = copy.deepcopy(self.geometry)
                copy_geometry.cartesian = res
                _results.append(copy_geometry)
            results = _results

        return results


    @staticmethod
    def dump(obj):
        """
        Prepare objects to be dumped using pickle (when needed)
        :param obj: object to dump, can be anything
        :return: None for zero/empty fields
        obj for integers and short strings
        pickle.dumps(obj) everywhere else
        """
        # integers and short strings are not modified
        if isinstance(obj, int) or \
                ( isinstance(obj, str) and 0 < len(obj) < 100 ):
            data = obj

        # non-empty/non-zero lists/arrays, long strings and non-zero floats are dumped using pickle
        # np.any returns False (False) for zero-value lists/arrays
        elif (
                ( isinstance(obj, list) or isinstance(obj, np.ndarray) ) and
                    len(obj) > 0 and
                        # list/array have elements in it and contains strings or non-zero values
                        ( isinstance(obj[0], str) or np.any(obj) )
            ) or \
                ( isinstance(obj, str) and len(obj) >= 100 ) or \
                ( isinstance(obj,float) and obj ):

            data = repr(pickle.dumps(obj))

        # everything else is set to None (zero/empty values)
        else:
            data = None

        return data


    def write_init(self, version, geometry, charge, command, modelH_top, real_top):
        """
        Write initial data into cobramm.xml
        :param version: cobramm version
        :param geometry: Layers class object
        :param charge: Charge class object
        :param command: lists of commands
        :param modelH_top: topology of modelH
        :param real_top: topology of complete system
        :return: None
        """

        # if this is a restart, don't write initial data
        if self.step != 0:
            return

        qm_type = {'1': 'gaussian',
                   '6': 'molcas',
                   '7': 'molpro',
                   '11': 'sharc'}

        i_keys = {'version': version,
                  'calculation_type': geometry.calculationType,
                  'optimization_type': command[1],
                  'qm_type': qm_type[command[51]],
                  'atom_labels': geometry.atomLabel,
                  'atom_layers': geometry.atomType,
                  'atom_links': np.concatenate((geometry.atomLink_BA,
                                                geometry.atomLink_DC), axis=1),
                  'modelH_top': modelH_top,
                  'real_top': real_top,
                  'geometry': geometry.cartesian,
                  'chargeMM': charge.CRG_real,
                  'chargeQM': charge.CRG_model_H}

        # prepare text
        text = '<INITIAL DATA>\n'
        for k in i_keys:
            i_keys[k] = self.dump(i_keys[k])
            if i_keys[k] is not None:
                text += '    <{0}>\n        {1}\n    </{0}>\n'.format(k, i_keys[k])
        text += '</INITIAL DATA>\n'
        # write text to file
        with open(self.filename, 'wb') as out_file:
            out_file.write(text.encode())


    def write_step(self, step, geometry, charge, energies, forces_QMMM_HM, dyn, opt):
        """
        Write step info after every step
        :param step: step number
        :param geometry: Layers class object
        :param charge: Charge class object
        :param energies: list of energies
        :param forces_QMMM_HM: total gradient
        :param dyn: MD data
        :param opt: geometry optimization data
        :return: None
        """

        # converting gradient to forces
        if forces_QMMM_HM is not None:
            forces_QMMM_HM = [-f for f in np.array(forces_QMMM_HM)]
        else:
            forces_QMMM_HM = None

        # step base keys: step num, geometry, charges, forces
        step_keys = {'step': step,
                     'geometry': geometry.cartesian,
                     'chargeMM': charge.CRG_real,
                     'chargeQM': charge.CRG_model_H,
                     'forces_QMMM_HM': forces_QMMM_HM}

        # energies
        ene_list = ['E_MM_real', 'E_MM_nocharge', 'E_MM_modelH',
                    'E_self_emb', 'E_QM', 'E_QMMM']
        ene_keys = {}
        for i,k in enumerate(ene_list):
            ene_keys[k] = energies[i]

        # dynamic
        dyn_keys = {}
        if dyn:
            dyn_list = ['time', 'state', 'velocity',
                        'MD_Etot', 'MD_Epot', 'MD_Ekin', 'MD_temp']
            for i,k in enumerate(dyn_list):
                dyn_keys[k] = dyn[i]

        # optimization
        opt_keys = {}
        if opt:
            opt_list = ['optstate', 'nIRCstep', 'Fmax', 'Frms', 'Dmax', 'Drms']
            for i,k in enumerate(opt_list):
                opt_keys[k] = opt[i]

        # prepare text
        # keys are kept separated to write similar fields in sequence
        text = '<STEP {0}>\n'.format(self.step)
        for k in step_keys:
            step_keys[k] = self.dump(step_keys[k])
            if step_keys[k] is not None:
                text += '    <{0}>\n        {1}\n    </{0}>\n'.format(k, step_keys[k])
        for k in ene_keys:
            ene_keys[k] = self.dump(ene_keys[k])
            if ene_keys[k] is not None:
                text += '    <{0}>\n        {1}\n    </{0}>\n'.format(k, ene_keys[k])
        for k in dyn_keys:
            dyn_keys[k] = self.dump(dyn_keys[k])
            if dyn_keys[k] is not None:
                text += '    <{0}>\n        {1}\n    </{0}>\n'.format(k, dyn_keys[k])
        for k in opt_keys:
            opt_keys[k] = self.dump(opt_keys[k])
            if opt_keys[k] is not None:
                text += '    <{0}>\n        {1}\n    </{0}>\n'.format(k, opt_keys[k])
        text += '</STEP {0}>\n'.format(self.step)

        # write text to file
        with open(self.filename, 'ab') as out_file:
            out_file.write(text.encode())

        # increment self.step
        self.step += 1


    def write_normal_modes(self, NM):
        """
        Dump normal modes in cobramm.xml for frequencies calculations
        :param NM: high precision normal modes string
        """
        # FREQ field is inserted between initial data and steps data
        self.normal_modes = NM
        nmdump = repr(pickle.dumps(self.normal_modes))
        insert_string = '</INITIAL DATA>\n<FREQ>\n    {0}\n</FREQ>\n'.format(nmdump)
        with open(self.filename, 'rb') as f:
            old_string = f.read()

        old_string_split = old_string.split(b'</INITIAL DATA>\n')

        # write tmp file
        with open('{0}.tmp'.format(self.filename), 'wb') as tmp:
            tmp.write(old_string_split[0])
            tmp.write(insert_string.encode())
            tmp.write(old_string_split[1])

        os.rename(self.filename, '{0}.old'.format(self.filename))
        os.rename('{0}.tmp'.format(self.filename), self.filename)
        os.remove('{0}.old'.format(self.filename))


class Step:
    """
    Class object representing a single step, used when parsing an existing cobramm.xml
    """
    # tuple of possible attributes for each step
    step_keys = ('step', 'geometry', 'chargeMM', 'chargeQM',
                 'E_MM_real', 'E_MM_nocharge', 'E_MM_modelH', 'E_self_emb', 'E_QM', 'E_QMMM',
                 'forces_QMMM_HM',
                 'time', 'state', 'velocity',
                 'MD_Etot', 'MD_Epot', 'MD_Ekin', 'MD_temp',
                 'optstate', 'nIRCstep', 'Fmax', 'Frms', 'Dmax', 'Drms')

    def __init__(self, step_string, geom_data):
        """
        Constructor method: init step_string and geom_data attributes
        Every other attribute is computed on demand with __getaddr__
        :param step_string: string, blob data to parse for the step
        """
        self.step_string = step_string
        self.geom_data = geom_data

    def __getattr__(self, item):
        """
        Compute attributes dynamically when requested
        :param item: keyword
        :return: object
        """
        # use Output.grep to fetch wanted item
        if item in Step.step_keys:
            result = Output.grep('{0}'.format(item), self.step_string)
            if result:
                if item in ['step', 'nIRCstep', 'state', 'optstate']:
                    # step, nIRCstep and state are integers, pickling not needed
                    self.__dict__[item] = int(result)
                elif item in ['geometry', 'chargeMM', 'chargeQM', 'E_QM', 'E_QMMM', 'forces_QMMM_HM', 'velocity']:
                    # loading lists/arrays
                    try:
                        self.__dict__[item] = pickle.loads(eval(result), encoding='bytes')
                    except TypeError:
                        self.__dict__[item] = pickle.loads(eval(result).encode(), encoding='bytes')
                else:
                    # make sure floats are loaded correctly
                    try:
                        self.__dict__[item] = float(pickle.loads(eval(result), encoding='bytes'))
                    except TypeError:
                        self.__dict__[item] = float(pickle.loads(eval(result).encode(), encoding='bytes'))
            else:
                # init every attribute
                self.__dict__[item] = None
        else:
            raise AttributeError('Step object has no attribute {0}'.format(item))

        if item == 'geometry':
            # rebuild geometry as Layers object using data_stream argument
            X = self.geometry[0]
            Y = self.geometry[1]
            Z = self.geometry[2]
            atom_labels = self.geom_data[0]
            atom_layers = self.geom_data[1]
            atom_links = self.geom_data[2]
            geom_stream = ['{0} {1} {2} {3} {4} {5}'.format(atom_labels[i],
                                                            0,
                                                            X[i],
                                                            Y[i],
                                                            Z[i],
                                                            atom_layers[i]) for i in range(len(atom_labels))]
            for i in range(len(atom_links[0])):
                geom_stream[atom_links[0][i] - 1] += ' {0} {1}'.format(atom_layers[atom_links[1][i] - 1],
                                                                       atom_links[1][i])
            self.geometry = Layers.from_real_layers_xyz(geom_stream)

        return self.__dict__[item]


    def get_attribute(self, label):
        """
        Alternative method to retrieve Step attributes
        Also compute attributes not stored in cobramm.xml (gradient, occupation)
        :param label: attribute wanted
        :return: attribute or None
        """
        value = None
        try:
            # we store forces, compute gradient when needed
            if label.lower() == 'gradient':
                value = [-f for f in self.forces_QMMM_HM]
            elif label.lower() in ['population', 'occupation']:
                value = []
                for i in range(len(self.E_QM)):
                    value.append(1.0 if i+1 == self.state else 0.0)
            else:
                exec( "value = self.{0}".format(label) )
        except (AttributeError, TypeError):
            value = None
        return value
