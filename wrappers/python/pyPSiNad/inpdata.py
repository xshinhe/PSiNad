#!/bin/env python
import json

spinboson = {
    "model": {
        "name": "SystemBath",
        "system_flag": "SB",
        "coupling_flag": "SB",
        "strength_flag": "Alpha",
        "bath_flag": "Ohmic",
        "bath_omegac": 1.0,
        "bath_strength": 0.4,
        "bath_temperature": "0.25 auK^-1",
        "bath_classic": False,
        "nbath": 1,
        "Nb": 100,
        "N": 100,
        "F": 2,
        "bias": 1.0,
        "delta": 1.0,
        "occ": 0,
        "dt": 0.01,
        "tend": 15.0
    },
    "solver": {
        "name": "NAD",
        "sampling_ele_flag": "Constraint",
        "sampling_nuc_flag": "Gaussian",
        "M": 100,
        "msize": 4,
        "sstep": 20,
        "dump": "null",
        "gamma": -1.0,
        "representation_flag": "Adiabatic",
        "inp_repr_flag": "Diabatic",
        "ele_repr_flag": "Diabatic",
        "nuc_repr_flag": "Adiabatic",
        "naforce": "EHR",
        "hopping_choose_type": 1,
        "hopping_direction_type": 2,
        "reflect": False,
        "BATH_FORCE_BILINEAR": True,
        "use_cv": True,
        "use_fssh": False,
        "offd_projected": True,
        "conserve_scale": True,
        "basis_switch": False
    },
    "record": [
        {"rule": "CPSrho<mij>:C(w, K1{@0}<m[occ][occ]>, K2<mij>)", "save": "SB2_CPSw.dat"},
    ]
}