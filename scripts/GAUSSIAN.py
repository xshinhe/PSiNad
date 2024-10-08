import argparse
import datetime
import math
import os
import re
import shutil
import subprocess
import sys
import time
import toml
from copy import deepcopy
from multiprocessing import Pool
from pprint import pprint
from traceback import format_exc
from socket import gethostname
import numpy as np
import QMutils
import typing

parser = argparse.ArgumentParser(description='Execute GAUSSIAN Calculation')
parser.add_argument('-d', '--directory', dest='directory', nargs='?', default='GAUSSIAN', type=str,
    help='work directory')
parser.add_argument('-i', '--input', dest='input', nargs='?', default='QM.in.GAUSSIAN', type=str,
    help='input file')
parser.add_argument('-o', '--output', dest='output', nargs='?', default='QM.out.GAUSSIAN', type=str,
    help='output file')
# args = parser.parse_args()

def qm_job(qm_data, args):
    qm_config = qm_data["qm_config"]
    mndo_config = qm_config["QM"]["GAUSSIAN"]
    natom = qm_data["natom"]
    znumber = qm_data["znumber"]
    xyz = qm_data["geom_xyz"]
    try:
        F = int(qm_config['QM']['GAUSSIAN']['F'])
        N = int(qm_config['QM']['GAUSSIAN']['N'])
    except KeyError:
        print(format_exc())


        def fileText(self, inpFileName, memory="500MB", nproc=1, gversion="g16"):
        """Prepare the text of gaussian input file named after the inpFileName variable (input file
        should be named <inpFileName>.com, so inpFileName should not include the extension!), using the
        input data defined in the instance of gaussianInput. The memory (var memory) and
        the number of cores (var nproc) to use in the execution and the version of Gaussian (gversion)
         are given as input because they are decided when running the Gaussian calculation and are not
         considered parameters of the QM calculation"""

    nproc = 1
    try:
        nproc = int(qm_config['QM']['nproc'])
    except KeyError:
        pass
    memory = '500MB'
    try:
        memory = int(qm_config['QM']['memory'])
    except KeyError:
        pass
    

    # initialize the text of the input file
    job_str = ""

    # add the section with the names of the auxiliary files
    job_str += '%chk={0}\n%rwf={0}\n%int={0}\n%d2e={0}\n'.format(inpFileName)
    # add the number of cores to use and the memory of the calculation
    job_str += '%Nproc={0}\n%Mem={1}\n'.format(nproc, memory)

    # add the gaussian keys using the variable self.keys up to the first empty line
    for line in self.keys:
        job_str += line.lower()
        if line.strip() == '':
            break
        else:
            job_str += "\n"

    # if forces are required in the output
    if self.otheropt["forces"]: job_str += " force"
    # if the calculation has background charges
    if self.otheropt["charges"]: job_str += " charge"
    # if computing electric field is required
    if self.otheropt["field"]: job_str += " prop=(field,read)"
    # check if a command for producing charges is there, otherwise add the option POP=(FULL,CHELPG)
    if not _optionInRoute(self.keys, ['CHELP', 'MK']): job_str += " pop=(full,chelpg)"
    # if a chk file is provided for restarting the orbitals
    if self.otheropt["restart"]: job_str += " guess=(read)"
    # options that are always needed
    if not _optionInRoute(self.keys, ['gfinput']): job_str += " gfinput"
    job_str += " scf=tight"

    # when it is requested to compute td couplings, require tamm-damkov and force printing of all excitations
    if self.otheropt["tdcouplings"]:
        if _optionInRoute(self.keys, ['td']):  # when the TD option is active, substitute it with TDA (Tamm-Damkov)
            logwrt.writewarning("Calculation of time derivative couplings is only possible with TDA...replacing TD keyword in Gaussian input\n")
            job_str = job_str.replace(" td ", " tda ")  # without subsequent options
            job_str = job_str.replace(" td(", " tda(")  # with subsequent options (syntax 1)
            job_str = job_str.replace(" td=(", " tda(")  # with subsequent options (syntax 2)
        else:
            job_str += " tda"
        job_str += " IOp(9/40=16)"  # add option to print the whole set of excitation coefficients
    job_str += "\n"

    # when it is requested to change the electronic state number
    addTDSection = None
    if self.otheropt["nstate"] is not None:
        if _optionInRoute(self.keys, ['td', 'tda']):  # when the calculation is TD or TDA
            if self.otheropt["nstate"] == 0:
                findtdgroup = re.findall("(tda?\(.*?\))", job_str)
                if findtdgroup:
                    addTDSection = findtdgroup[0]
                    job_str = job_str.replace(findtdgroup[0], "")
                else:  # in this case there is an isolated td / tda keyword
                    addTDSection = "tda"
                    job_str = job_str.replace(" td ", " ")
                    job_str = job_str.replace(" tda ", " ")
            else:
                findroot = re.findall("(root *= *[0-9]*)", job_str)  # look for an expression of type "root = N"
                if findroot:
                    if self.otheropt["nstate"] != int(findroot[0].split("=")[1]) and \
                            not self.otheropt["forcestate"]:
                        logwrt.fatalerror("Mismatch of the electronic state between command file and QM input")
                    job_str = job_str.replace(findroot[0], "root={0}".format(self.otheropt["nstate"]))
                else:
                    findtdgroup = re.findall("(tda?=?\(.*?\))", job_str)  # look for a group tda( ) or td( )
                    if findtdgroup:
                        job_str = job_str.replace(findtdgroup[0], findtdgroup[0][:-1] +
                                                      ",root={0})".format(self.otheropt["nstate"]))
                    else:  # in this case there is an isolated td / tda keyword
                        job_str = job_str.replace(" td ", " td(root={0}) ".format(self.otheropt["nstate"]))
                        job_str = job_str.replace(" tda ", " tda(root={0}) ".format(self.otheropt["nstate"]))
        else:  # this means that this is not a td/tda calculation, the "nstate" option is not implemented
            logwrt.fatalerror("cannot change the electronic state number for these Gaussian options")

    # add a blank line, the title, another blank line and then charge and multeplicity
    job_str += "\nQM single-point calculation for COBRAMM\n\n"
    job_str += self.keys[-1] + "\n"

    # write molecular structure
    for atom in zip(self.symbols, *self.coords):
        job_str += "{0} {1:16.8f} {2:16.8f} {3:16.8f}\n".format(*atom)
    job_str += "\n"

    # inserting charges
    if self.otheropt["charges"]:
        coords = self.otheropt["charges"][0]
        chrg = self.otheropt["charges"][1]
        for atom in zip(chrg, *coords):
            job_str += "{1:14.8f} {2:14.8f} {3:14.8f} {0:16.8f}\n".format(*atom)
        job_str += "\n"

    # insert basis set information only when the GEN keyword is found in the given route options
    if _optionInRoute(self.keys, ['GEN']):
        for line in self.gen:
            job_str += line + "\n"
        job_str += "\n"

    # in g16 the gaussweights section (weights of the CASSCF) should be before the EF
    if gversion == 'g16' and len(self.gaussweights) != 0:
        for line in self.gaussweights:
            job_str += line + '\n'
        job_str += '\n'

    # inserting EF for gradient computation
    if self.otheropt["field"]:
        for coords in zip(*self.otheropt["field"][0]):
            job_str += "{0:14.8f} {1:14.8f} {2:14.8f}\n".format(*coords)
        job_str += '\n'

    # in g16 the gaussweights section (weights of the CASSCF) should be after the EF
    if gversion != 'g16' and len(self.gaussadd) != 0:
        for line in self.gaussweights:
            job_str += line + '\n'
        job_str += '\n'

    # add the remaining part of the input, stored in the gaussadd variable
    for line in self.gaussadd:
        job_str += line + '\n'
    job_str += '\n'

    # in case of a TD calculation with ground state gradient, append a second input at the end
    if not self.otheropt["suppress_TDDFT"] and addTDSection is not None:
        # add the "link" separation
        job_str += "--Link1--\n"
        # add the section with the names of the auxiliary files
        job_str += '%chk={0}\n%rwf={0}\n%int={0}\n%d2e={0}\n'.format(inpFileName)
        # add the number of cores to use and the memory of the calculation
        job_str += '%Nproc={0}\n%Mem={1}\n'.format(nproc, memory)
        # add the route line
        job_str2 = "" 
        for line in self.keys:
            job_str2 += line.lower()
            if line.strip() == '':
                break
            else:
                job_str2 += "\n"
        job_str2 = job_str2.replace(" td ", " tda ")  # without subsequent options
        job_str2 = job_str2.replace(" td(", " tda(")  # with subsequent options
        job_str2 = job_str2.replace(" force ", " ")   # remove force  
        job_str2 = job_str2.replace(",eqsolv", "")   # ES not in equilibrium when dynamics in GS
        job_str2 = job_str2.replace("eqsolv,", "")   # ES not in equilibrium when dynamics in GS
        job_str2 = job_str2.replace("eqsolv", "")   # ES not in equilibrium when dynamics in GS
        # if the calculation has background charges
        if self.otheropt["charges"]: job_str2 += " charge"
        job_str2 += " gfinput"
        job_str2 += " scf=tight"
        job_str += job_str2
        # add restart options
        #job_str += " guess=(read) geom=(check) " + addTDSection + " IOp(9/40=16) \n"
        job_str += " guess=(read) geom=(check) IOp(9/40=16) \n"
        # add a blank line, the title, another blank line and then charge and multeplicity
        job_str += "\nQM single-point calculation for COBRAMM\n\n"
        job_str += self.keys[-1] + "\n\n"

        # inserting charges also for second calculation
        if self.otheropt["charges"]:
            coords = self.otheropt["charges"][0]
            chrg = self.otheropt["charges"][1]
            for atom in zip(chrg, *coords):
                job_str += "{1:14.8f} {2:14.8f} {3:14.8f} {0:16.8f}\n".format(*atom)
            job_str += "\n"

    # in case of a CI optimization (with gmean branching plane) we also need gradient of lower state
    if self.otheropt["CIopt"]:
        # set lower state and determine if it is GS (DFT only) or not (TD required)
        lower_state = self.otheropt["upper_state"] - 1
        if lower_state == 1:
            GS = True
        else:
            GS = False
        # add the "link" separation
        job_str += "--Link1--\n"
        # add the section with the names of the auxiliary files
        job_str += '%chk={0}\n%rwf={0}\n%int={0}\n%d2e={0}\n'.format(inpFileName)
        # add the number of cores to use and the memory of the calculation
        job_str += '%Nproc={0}\n%Mem={1}\n'.format(nproc, memory)
        # add the route line
        job_str2 = "" 
        for line in self.keys:
            job_str2 += line.lower()
            if line.strip() == '':
                break
            else:
                job_str2 += "\n"
        if GS:
            findtdgroup = re.findall("(tda?=?\(.*?\))", job_str2)
            if findtdgroup:
                job_str2 = job_str2.replace(findtdgroup[0], "")
            else:  # in this case there is an isolated td / tda keyword
                job_str2 = job_str2.replace(" td ", " ")
                job_str2 = job_str2.replace(" tda ", " ")  
        else:
            findroot = re.findall("(root *= *[0-9]*)", job_str2)  # look for an expression of type "root = N"
            if findroot:
                job_str2 = job_str2.replace(findroot[0], "root={0}".format(lower_state - 1))
            else:
                findtdgroup = re.findall("(tda?\(.*?\))", job_str2)  # look for a group tda( ) or td( )
                if findtdgroup:
                    job_str2 = job_str2.replace(findtdgroup[0], findtdgroup[0][:-1] +
                                                  ",root={0})".format(lower_state - 1))
                else:  # in this case there is an isolated td / tda keyword
                    job_str2 = job_str2.replace(" td ", " td(root={0}) ".format(lower_state - 1))
                    job_str2 = job_str2.replace(" tda ", " tda(root={0}) ".format(lower_state - 1))
        if self.otheropt["forces"]: job_str2 += " force"
        # if the calculation has background charges
        if self.otheropt["charges"]: job_str2 += " charge"
        # if computing electric field is required
        if self.otheropt["field"]: job_str2 += " prop=(field,read)"
        # check if a command for producing charges is there, otherwise add the option POP=(FULL,CHELPG)
        if not _optionInRoute(self.keys, ['CHELP', 'MK']): job_str2 += " pop=(full,chelpg)"
        job_str2 += " gfinput"
        job_str2 += " scf=tight"
        job_str += job_str2
        # add restart options
        job_str += " guess=(read) geom=(check)\n"
        # add a blank line, the title, another blank line and then charge and multeplicity
        job_str += "\nQM single-point calculation for COBRAMM\n\n"
        job_str += self.keys[-1] + "\n\n"

        # inserting charges also for second calculation
        if self.otheropt["charges"]:
            coords = self.otheropt["charges"][0]
            chrg = self.otheropt["charges"][1]
            for atom in zip(chrg, *coords):
                job_str += "{1:14.8f} {2:14.8f} {3:14.8f} {0:16.8f}\n".format(*atom)
            job_str += "\n"

        # insert basis set information only when the GEN keyword is found in the given route options
        if _optionInRoute(self.keys, ['GEN']):
            for line in self.gen:
                job_str += line + "\n"
            job_str += "\n"

        # in g16 the gaussweights section (weights of the CASSCF) should be before the EF
        if gversion == 'g16' and len(self.gaussweights) != 0:
            for line in self.gaussweights:
                job_str += line + '\n'
            job_str += '\n'

        # inserting EF for gradient computation
        if self.otheropt["field"]:
            for coords in zip(*self.otheropt["field"][0]):
                job_str += "{0:14.8f} {1:14.8f} {2:14.8f}\n".format(*coords)
            job_str += '\n'

        # in g16 the gaussweights section (weights of the CASSCF) should be after the EF
        if gversion != 'g16' and len(self.gaussadd) != 0:
            for line in self.gaussweights:
                job_str += line + '\n'
            job_str += '\n'

        # add the remaining part of the input, stored in the gaussadd variable
        for line in self.gaussadd:
            job_str += line + '\n'
        job_str += '\n'



    job_str = '{"bagel": [\n'

    # molecule part
    job_str += '{\n'
    job_str += '"title": "molecule",\n'
    job_str += '"basis": "%s",\n'%(qm_config['QM']['GAUSSIAN']['basis'])
    job_str += '"df_basis": "%s",\n'%(qm_config['QM']['GAUSSIAN']['df_basis'])
    job_str += '"angstrom": true,\n'
    job_str += '"geometry": [\n'
    for i in range(natom):
        end = ',' if i!=natom-1 else ''
        job_str += '{"atom": "%s", "xyz": [%12.8e, %12.8e, %12.8e]}%s\n'%(
                QMutils.element_list[znumber[i]][0], xyz[i][1], xyz[i][2], xyz[i][3], end)
    job_str += ']\n' # end geometry
    job_str += '},\n' # end molecule part

    # force part
    job_str += '{\n'
    job_str += '"title": "forces",\n'
    job_str += '"dipole": true,\n'
    job_str += '"export": true,\n'
    job_str += '"nproc": 1,\n'
    job_str += '"grads": [\n'
    for i in range(F):
        end = ',' if i!=F-1 else ''
        job_str += '{"title": "force", "target": %d},\n'%i
    for i in range(F):
        for k in range(i+1,F):
            end = '' if i==k-1 and k==F-1 else ','
            job_str += '{"title": "nacme", "target": %d, "target2": %d, "nacmtype": "full"}%s\n'%(i,k,end)
    job_str += '],\n' # end grads
    # force methods
    job_str += '"method": ['

    if qm_config['QM']['GAUSSIAN']['method'] == 'casscf' or qm_config['QM']['GAUSSIAN']['method'] == 'caspt2':
        job_str += '{\n'
        job_str += '"title": "casscf",\n'
        job_str += '"nstate": %d,\n'%(qm_config['QM']['GAUSSIAN']['nstate']) #should >= F (always ==)
        job_str += '"nact": %d,\n'%(qm_config['QM']['GAUSSIAN']['nact'])
        # job_str += '"nopen": %d,\n'%(qm_config['QM']['GAUSSIAN']['nopen'])
        job_str += '"nclosed": %d,\n'%(qm_config['QM']['GAUSSIAN']['nclosed'])
        if 'active' in qm_config['QM']['GAUSSIAN']:
            job_str += '['
            cnt = 0
            for istate in qm_config['QM']['GAUSSIAN']['active']:
                if cnt !=0: job_str += ','
                job_str += '%d'%istate
                cnt += 1
            job_str += '],\n'
        if qm_config['QM']['GAUSSIAN']['method'] == 'caspt2':
            job_str += '"simth": {\n'
            job_str += '"method": "caspt2",\n'
            if 'ms' in qm_config['QM']['GAUSSIAN']:
                job_str += '"ms": "%s",\n'%(qm_config['QM']['GAUSSIAN']['ms'])
            if 'xms' in qm_config['QM']['GAUSSIAN']:
                job_str += '"xms": "%s",\n'%(qm_config['QM']['GAUSSIAN']['xms'])
            if 'sssr' in qm_config['QM']['GAUSSIAN']:
                job_str += '"sssr": "%s",\n'%(qm_config['QM']['GAUSSIAN']['sssr'])
            job_str += '"shift": %.3f\n'%(qm_config['QM']['GAUSSIAN']['shift'])
            job_str += '},\n'

        job_str += '"charge": %d,\n'%(qm_config['QM']['GAUSSIAN']['charge'])
        job_str += '"nspin": %d\n'%(qm_config['QM']['GAUSSIAN']['nspin'])
        job_str += '}\n'
    else:
        raise ValueError("unsupport method")

    job_str += ']\n' # end method
    job_str += '}\n' # end force part

    job_str += ']}' # end of bagel parameter

    qm_config["QM"]["env"] = {
        "input_is_ready":True,
        "generated": "QM.run.GAUSSIAN.json",
        "directory": args.directory,
        "output": args.output,
    }

    directory = qm_config['QM']['env']['directory']
    if not os.path.exists(directory): os.makedirs(directory)

    f = open(directory + '/' + qm_config['QM']['env']['generated'], 'w')
    f.write(job_str)
    f.flush()
    f.close()

    exe_str = 'cd %s && %s  %s > %s && cd -'%(
        qm_config['QM']['env']['directory'],
        qm_config['QM']['GAUSSIAN']['path'],
        qm_config['QM']['env']['generated'],
        qm_config['QM']['env']['output']
        )
    print(exe_str)
    os.system(exe_str)

    parse_result(
        qm_data,
        directory
    )

def parse_result(qm_data, log_path):
    try:
        natom = qm_data['natom']
        qm_config = qm_data['qm_config']
        F = int(qm_config['QM']['GAUSSIAN']['F'])
        N = int(qm_config['QM']['GAUSSIAN']['N'])
        nciref = F # int(qm_config['QM']['GAUSSIAN']['keywords']['nciref'])
    except KeyError:
        print(format_exc())

    ERROR_MSG = ""
    eig = np.zeros((F))
    dE  = np.zeros((F, N))
    nac = np.zeros((F, F, N))
    
    with open(log_path+'/ENERGY.out', 'r') as ifs:
        for i in range(F):
            eig[i] = float(ifs.readline())

    for i in range(F):
        with open(log_path+'/FORCE_%d.out'%i, 'r') as ifs:
            ifs.readline()
            for a in range(natom):
                dE[i,3*a:3*a+3] = ifs.readline().split()[1:]

    for i in range(F):
        for k in range(i+1, F):
            with open(log_path+'/NACME_%d_%d.out'%(i,k), 'r') as ifs:
                ifs.readline()
                for a in range(natom):
                    nac[i,k,3*a:3*a+3] = ifs.readline().split()[1:]
            nac[k,i,:] = -nac[i,k,:]
    
    # convert to au unit(no need)
    
    qmout = QMutils.QMout(natom=natom, 
        energy=eig,
        gradient=dE.T,
        nac=np.einsum('ikj->jik', nac)
    )

    f = open(qm_config['QM']['env']['directory'] + '/energy.dat', 'w')
    for i in range(F):
        f.write('{: 12.8e}\n'.format(eig[i]))
    f.close()

    f = open(qm_config['QM']['env']['directory'] + '/gradient.dat', 'w')
    for j in range(N):
        for i in range(F):
            f.write('{: 12.8e} '.format(dE[i,j]))
        f.write('\n')
    f.close()

    f = open(qm_config['QM']['env']['directory'] + '/nacv.dat', 'w')
    for j in range(N):
        for i in range(F):
            for k in range(F):
                f.write('{: 12.8e} '.format(nac[i,k,j]))
        f.write('\n')
    f.close()
    
    #pprint(qmout)
    return qmout

def main():
    return

if __name__ == '__main__':
    args = parser.parse_args()
    pprint(args)
    qm_data_in = QMutils.parseQMinput(args.input)
    qm_job(qm_data_in, args)

    main()
