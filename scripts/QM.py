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
import typing
import numpy as np
import numpy.typing as npt
import QMutils
import MNDO
import BAGEL

parser = argparse.ArgumentParser(description='Execute MNDO Calculation')
parser.add_argument('-d', '--directory', dest='directory', nargs='?', default='.', type=str,
    help='work directory')
parser.add_argument('-i', '--input', dest='input', nargs='?', default='QM.in.MNDO', type=str,
    help='input file')
parser.add_argument('-o', '--output', dest='output', nargs='?', default='QM.out.MNDO', type=str,
    help='output file')
args = parser.parse_args()

def specify_qm_job(qm_data):
    # get toml config
    qm_config = qm_data["qm_config"]

    if(qm_config['exec'] == "MNDO"):
        call_MNDO_job(qm_data, )
    elif(qm_config['exec'] == "BAGEL"):
        pass
    elif(qm_config['exec'] == "ORCA"):
        pass
    elif(qm_config['exec'] == "XXX"):
        pass
    else:
        pass

if __name__ == '__main__':
    qm_data_in = QMutils.parseQMinput(args.input)
    qm_config = qm_data_in["qm_config"]["QM"]
    if(qm_config['exec'] == "MNDO"):
        qm_data_out = MNDO.qm_job(qm_data_in, args)
    elif(qm_config['exec'] == "BAGEL"):
        qm_data_out = BAGEL.qm_job(qm_data_in, args)
    elif(qm_config['exec'] == "ORCA"):
        pass
    elif(qm_config['exec'] == "XXX"):
        pass
    else:
        pass

