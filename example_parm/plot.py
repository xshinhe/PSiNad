import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

a = pd.read_csv(sys.argv[1], sep='\s+')
F = len(a.values[0,:])//2
for i in range(F):
    plt.plot(a['t'], a['RwK0occ0Ktdia%s'%i], label='%d'%(1+i))
plt.legend()
plt.show()

