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

flagA = 'wK1occ'
if N > 2:
    flagA = sys.argv[2]

flagB = 'K2dia'
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

pops = np.zeros((len(a['t']), F))

if scheme == "pop":
    for i in range(F):
        pops[:,i] = a['R%s0%s%d'%(flagA, flagB, i)]

if scheme == "rho":
    for i in range(F):
        pops=a['R%s0%s%d'%(flagA, flagB, i*(F+1))]

for i in range(F):
    plt.plot(a['t'], pops[:,i]/np.sum(pops, axis=1), '.-', label='%d'%(i+1))
plt.legend()
plt.show()


