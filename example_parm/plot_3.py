import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

a = pd.read_csv(sys.argv[1], sep='\s+')
F = (len(a.values[0,:])//2) // 4
print(F)
F = 2
for i in range(F):
    plt.plot(a['t'], 
            a['R[F0]0Kt%s'%(i*(F+1))] +
            a['R[F1]0Kt%s'%(i*(F+1))] +
            a['R[F2]0KtQ%s'%(i*(F+1))] +
            a['R[F3]0KtQ%s'%(i*(F+1))]
            , label='%d'%(1+i))
plt.legend()
plt.show()

