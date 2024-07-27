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
parser.add_argument('integers', metavar='N', type=int, nargs='+', 
    help='an integer for the accumulator')
parser.add_argument('--sum', dest='accumulate', action='store_const', const=sum, default=max,
    help='sum the integers (default: find the max)')

def specify_qm_job(qm_data):
    # get toml config
    qm_config = qm_data["qm_config"]

    if(qm_config['exec'] == "MNDO"):
        call_MNDO_job(qm_data)
    elif(qm_config['exec'] == "BAGEL"):
        pass
    elif(qm_config['exec'] == "ORCA"):
        pass
    elif(qm_config['exec'] == "XXX"):
        pass
    else:
        pass

if __name__ == '__main__':
    qm_data_in = QMutils.parseQMinput(sys.argv[1])
    qm_config = qm_data_in["qm_config"]["QM"]
    if(qm_config['exec'] == "MNDO"):
        qm_data_out = MNDO.qm_job(qm_data_in)
    elif(qm_config['exec'] == "BAGEL"):
        pass
    elif(qm_config['exec'] == "ORCA"):
        pass
    elif(qm_config['exec'] == "XXX"):
        pass
    else:
        pass