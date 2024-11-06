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

import os
import toml
from argparse import ArgumentParser
from typing import List, Dict, Any
from functools import reduce
from operator import getitem

from constants import element_list
import kids_io
from kids_log import Log
from Layers import Layers

def parseXYZ(lines: List[str]) -> Dict[str, Any]:
    """
    Parses an XYZ file content into a dictionary containing the number of atoms, atomic numbers,
    geometry coordinates, and velocities if present.

    Parameters:
    lines (List[str]): A list of strings, each representing a line from the XYZ file.

    Returns:
    Dict[str, Any]: A dictionary with keys 'natom', 'znumber', 'xyz', and 'vel'.
    """
    natom = int(lines[0].strip())  # Number of atoms
    geom_xyz: List[List[Union[str, float]]] = []
    velocity: List[List[float]] = []
    znumber: List[int] = []
    has_velocity: bool = True  # Flag to check if velocities are present

    for i in range(2, natom + 2):
        fields = lines[i].split()
        fields[0] = fields[0].title()  # Ensure element symbol is capitalized
        for j in range(1, 4):
            fields[j] = float(fields[j])  # Convert coordinates to float
        if len(fields) >= 7 and has_velocity:
            for j in range(4, 7):
                fields[j] = float(fields[j])  # Convert velocities to float
        else:
            has_velocity = False  # If no velocities are found, set the flag to False

        try:
            znumber.append(element_list[fields[0]][0])  # Assuming atomic number is the first item in the list
        except (IndexError, KeyError):
            pprint(f"Unknown element {fields[0]}")  # Print a message if the element is not found in element_list

        geom_xyz.append(fields[0:4])  # Append atom symbol and coordinates
        if has_velocity:
            velocity.append(fields[4:7])  # Append velocities if present

    return {
        "natom": natom,
        "znumber": znumber,
        "xyz": geom_xyz,
        "vel": velocity,
    }

class Config(Dict[str, Any]):
    @staticmethod
    def load(toml_file: str, args: ArgumentParser = None):
        """
        Load TOML data from a file and convert it into a Config instance.

        :param toml_string: A string containing TOML-formatted data
        :return: Config instance with the loaded data
        """
        if not os.path.exists(toml_file):
            raise RuntimeError(f'TOML-format configuration file {args.input} is required!\n')
        return Config.loads(''.join(open(toml_file, 'r').readlines()), args)

    @staticmethod
    def loads(toml_string: str, args: ArgumentParser = None):
        """
        Load TOML data from a string and convert it into a Config instance.

        :param toml_string: A string containing TOML-formatted data
        :return: Config instance with the loaded data
        """
        # Parse the TOML string into a dictionary
        parsed_toml = toml.loads(toml_string)
        # Convert the parsed dictionary into a Config
        instance = Config({k: Config(v) if isinstance(v, dict) else v for k, v in parsed_toml.items()})

        if args is not None:
            instance.update(args)
        return instance

    def update(self, args: ArgumentParser):
        
        self.args = args
        self.config_path = args.input
        self.topo1 = args.topology.split(',')[0]  # real
        self.topo2 = args.topology.split(',')[-1] # model-H
        self.layer = args.layer
        self.geom_in_toml = True

        geom: Dict[str, Any] = {}
        if os.path.exists(args.coord) and os.path.exists(args.layer):
            Log.writeLog(f'Update GEOM from [args.layer]={args.layer} & [args.coord]={args.coord}\n')    
            with open(args.layer, "r") as f:
                geometry = Layers.from_real_layers_xyz(f.read())
            geometry.updatereal(args.coord)
            geometry.makerealcrd()
            geom = geometry
            self.geom_in_toml = False
        elif 'GEOM' in self:
            Log.writeLog(f'Read GEOM from [args.input]={args.input}\n')
            lines = ''
            if 'xyz' in self['GEOM']:
                lines = self['GEOM']['xyz'].split('\n')
            elif 'read_xyz' in self['GEOM']:
                lines = kids_io.readFile(self['GEOM']['read_xyz'])
            if not lines:
                Log.writeLog('Cannot find XYZ Information\n')
                exit(0)
            geom = Layers.from_only_xyz(lines[2:])
        # Update the dictionary with the geom and config data
        self['_geom'] = geom 

    def get_nested(self, key_path, default=None):
        """
        Access nested toml data via a dot-separated string path.

        :param key_path: Dot-separated string path, e.g., "a.b.c"
        :param default: Default value to return if the path does not exist
        :return: The value corresponding to the path or the default value
        """
        keys = key_path.split('.')
        try:
            return reduce(getitem, keys, self)
        except (KeyError, TypeError):
            return default

if __name__ == '__main__':
    # Example usage
    toml_string = """
    [database]
    server = "192.168.1.1"
    port = 3306
    """

    # Load toml data and convert to Config using object_hook
    data = Config.loads(toml_string)

    # Access data
    server = data.get_nested("database.server")
    port = data.get_nested("database.port")

    # Access a non-existent key, return default value
    username = data.get_nested("database.username", "default_user")

    print(f"Server: {server}")  # Output: Server: 192.168.1.1
    print(f"Port: {port}")      # Output: Port: 3306
    print(f"Username: {username}")  # Output: Username: default_user