import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

a = pd.read_csv(sys.argv[1], sep='\s+')
F = (len(a.values[0,:])//2) // 4
print(F)
F = 7
for i in range(F):
    plt.plot(a['t'], 
            a['RKA0KA%s'%(i*(1))] +
            a['RKA0KD%s'%(i*(1))] +
            a['RKD0KA%s'%(i*(1))] +
            a['RKD0KD%s'%(i*(1))]
            , label='%d'%(1+i))
plt.legend()
plt.show()

