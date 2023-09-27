import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json 
from numba import jit


# for some global parameters & variables

N = 1
F = 1
Nb = 1
nbath = 1
omegac = 1.0
alpha  = 1.0
bias = 1.0
delta = 1.0
beta = 1.0

def ForceField_parm(json_str):
    print(json_str)
    parm = json.loads(json_str)

    N = int(parm["rdim"])
    F = int(parm["fdim"])
    Nb = int(parm["Nb"])
    nbath = int(parm["nbath"])
    omegac = float(parm["omegac"])
    alpha = float(parm["strength"])
    bias = float(parm["bias"])
    delta = float(parm["delta"])
    temp = float(parm["temp"].split(' ')[0])
    beta = temp # "auK^-1"

    global wj
    global cj

    z = np.arange(Nb) + 1
    wj = - omegac * np.log(1 - z/(Nb+1))
    cj = np.sqrt(omegac * alpha / (Nb + 1)) * wj

def ForceField_init(rdim, fdim, itraj):
    M = np.ones((rdim))*1
    Qjoverb = 0.5 * wj / np.tanh(0.5 * beta * wj)
    R = np.random.normal(size=(rdim)) * np.sqrt(Qjoverb / wj**2)
    P = np.random.normal(size=(rdim)) * np.sqrt(Qjoverb)
    erho = np.eye((fdim))*(1+0.j)
    eeac = np.ones((fdim))*(1+0.j)
    eocc = 0
    return [R,P,M,erho,eeac,eocc]

@jit(nopython=True)
def ForceField_npes(V, dV, ddV, R, P, flag, rdim):
    ## note: if flag < 1,  dV is nullptr, don't access it, it will cause fatal error!!!
    ## note: if flag < 2, ddV is nullptr, don't access it, it will cause fatal error!!!

    ## do some calculation

    ## final asignment, index[0] & index[:] is necessary
    V[0] = 0
    dV[:] = wj**2*R


@jit(nopython=True)
def ForceField_epes(V, dV, ddV, R, flag, rdim, fdim):
    ## note: if flag < 1,  dV is nullptr, don't access it, it will cause fatal error!!!
    ## note: if flag < 2, ddV is nullptr, don't access it, it will cause fatal error!!!

    ## do some calculation
    v00 = bias + np.dot(cj, R)

    ## final asignment, index[0] & index[:] is necessary
    V[:] = np.array([v00, delta, delta, -v00])
    dV[:] = np.outer(cj, np.array([1,0,0,-1])).reshape((rdim*fdim*fdim))

