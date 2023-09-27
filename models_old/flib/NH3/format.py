import numpy as np

a = np.loadtxt('nh3_param.inp', skiprows = 14, usecols=(0,1,2), dtype=str)
print(len(a))
[print(f'{_[0]:<9s} = {_[2]:>14s}d0') for _ in a if _[1]=='0']