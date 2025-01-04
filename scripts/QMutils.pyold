import os
import re
import sys
import toml
from copy import deepcopy
from pprint import pprint
from traceback import format_exc
import typing
import numpy as np
import numpy.typing as npt

COORD_XYZ_PH = "$COORD_XYZ"
FIELD_XYZ_PH = "$FIELD_XYZ"

def readfile(filename: str):
    try:
        f = open(filename)
        out = f.readlines()
        f.close()
    except IOError:
        pprint('File {filename} does not exist!', format_exc())
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

def parseQMinput(filename: str) -> dict:
    lines = readfile(filename)
    try:
        natom = int(lines[0])
    except ValueError:
        print('First line must be the number of atoms!', format_exc())
        sys.exit(-1)
    if len(lines) < natom + 4:
        print("too few lines")
        sys.exit(-1)

    geom_xyz = []
    velocity = [] # or momentum be better?
    znumber  = [] 
    has_velocity = True
    for i in range(2, natom + 2):
        fields = lines[i].split()
        fields[0] = fields[0].title()
        for j in range(1, 4):
            fields[j] = float(fields[j])
        if len(fields) >= 7 and has_velocity:
            for j in range(4, 7):
                fields[j] = float(fields[j])
        else:
            has_velocity = False
        try:
            znumber += [ element_list[fields[0]][1] ]
        except ValueError:
            print("Unknown element {field[0]}")
        geom_xyz +=  [fields[0:4]]
        if has_velocity: 
            velocity +=  [fields[4:7]]
    
    toml_string = ''.join(lines[natom + 2:])
    qm_config = toml.loads(toml_string)

    # pprint(znumber)
    # pprint(geom_xyz)
    # pprint(qm_config)

    return {
        "natom" : natom,
        "znumber" : znumber,
        "geom_xyz" : geom_xyz,
        "velocity" : velocity,
        "qm_config": qm_config
    }

class QMout(typing.TypedDict):
    natom: int
    energy: npt.ArrayLike
    gradient: npt.ArrayLike 
    nacv: npt.ArrayLike
