import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import re

N = len(sys.argv)
fn = 'corr.dat'
if N > 1:
    fn = sys.argv[1]
a = pd.read_csv(fn, sep='\s+')

flagA = 'Kocc'
if N > 2:
    flagA = sys.argv[2]

flagB = 'Kdia'
if N > 3:
    flagB = sys.argv[3]


scheme = "pop"
if N > 4:
    scheme = sys.argv[4]

pattern=r'\d+$'
rex = re.search(pattern, a.columns[-1])
F = int(rex.group()) + 1
if scheme == 'pop': F = F
if scheme == 'rho': F = int(np.sqrt(F))
if N > 5:
    F = int(sys.argv[5])

if scheme == "pop":
    for i in range(F):
        plt.plot(a['t'], a['R%s0%s%d'%(flagA, flagB, i)], label='%d'%(1+i))
    plt.legend()
    plt.show()

if scheme == "rho":
    for i in range(F):
        plt.plot(a['t'], a['R%s0%s%d'%(flagA, flagB, i*(F+1))], label='%d'%(1+i))
    plt.legend()
    plt.show()

