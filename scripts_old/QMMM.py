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

parser = argparse.ArgumentParser(description='Execute MNDO Calculation')
parser.add_argument('integers', metavar='N', type=int, nargs='+', 
    help='an integer for the accumulator')
parser.add_argument('--sum', dest='accumulate', action='store_const', const=sum, default=max,
    help='sum the integers (default: find the max)')

def readfile(filename):
    try:
        f = open(filename)
        out = f.readlines()
        f.close()
    except IOError:
        print('File {filename} does not exist!', format_exc())
        sys.exit(-1)
    return out

def writefile(filename, content):
    # content can be either a string or a list of strings
    try:
        f = open(filename, 'w')
        if isinstance(content, list):
            for line in content:
                f.write(line)
        elif isinstance(content, str):
            f.write(content)
        else:
            print('content {content} cannot be written to file!')
        f.close()
    except IOError:
        print('Could not write to file {filename}!', format_exc())
        sys.exit(-1)

def readQMinput(QMfilename):
    QMinlines = readfile(QMfilename)
    try:
        natom = int(QMinlines[0])
    except ValueError:
        print('First line must be the number of atoms!', format_exc())
        sys.exit(-1)

    if len(QMinlines) < natom + 4:
        print('''
            Input file must contain at least:
                natom
                comment
                geometry
                keyword "states"
                at least one task
            ''')
        sys.exit(23)
    geom = []
    velo = []
    hasveloc = True
    for i in range(2, natom + 2):
        fields = QMinlines[i].split()
        fields[0] = fields[0].title()
        for j in range(1, 4):
            fields[j] = float(fields[j])
        if len(fields) >= 7:
            for j in range(4, 7):
                fields[j] = float(fields[j])
        else:
            hasveloc = False
        symbol = fields[0]
        geom +=  [fields[0:4]]
        if hasveloc: velo +=  [fields[4:7]]
    i = natom + 2
    toml_string = '\n'.join(QMinlines[i:])
    qm_data = toml.loads(toml_string)
    pprint(geom)
    pprint(velo)
    pprint(qm_data)

if __name__ == '__main__':
    readQMinput(sys.argv[1])