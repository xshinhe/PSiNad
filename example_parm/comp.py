import matplotlib.pyplot as plt
import pandas as pd

a = pd.read_csv("samp0-dia.dat", sep='\s+')
b = pd.read_csv("samp0-adia.dat", sep='\s+')

for f in ['x', 'p', 'f']:
    plt.plot(a['t'], a[f+'0'], 'r-')
    plt.plot(b['t'], b[f+'0'], 'b--')
    plt.show()
