import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

a = pd.read_csv(sys.argv[1], sep='\s+')
F = len(a.values[0,:])//4
for i in range(F):
    plt.plot(a['t'], a['RwK0occ0[F0]%s'%(2*i)], label='T%d'%(1+i))
    plt.plot(a['t'], a['RwK0occ0[F0]%s'%(2*i+1)], label='R%d'%(1+i))
plt.legend()
plt.show()

