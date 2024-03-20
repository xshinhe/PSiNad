import numpy as np
import pandas as pd
import sys

alist = ['C']*7 + ['O'] + ['H']*7 + ['C','H','N','H'] + ['C']*4 + ['N','H','C'] + ['H']*4
au_2_ang = 0.52918

if __name__ == '__main__':
    file = sys.argv[1]
    out  = sys.argv[2]

    df = pd.read_csv(file, sep='\s+')
    t  = df["t"]
    xyz = df.filter(regex='I0x', axis=1).values

    Nt = len(t)
    N  = len(xyz[0,:])
    A  = N//3
    
    xyz = xyz.reshape((Nt,A,3))*au_2_ang

    f = open(out, 'w')
    for i in range(Nt):
        f.write('%d\n'%A)
        f.write('t=%f fs\n'%t[i])

        for k in range(A):
            f.write('%s   %12.8e  %12.8e  %12.8e\n'%(alist[k], xyz[i,k,0], xyz[i,k,1], xyz[i,k,2]))
        #f.write('\n')
    f.close()
    #print(np.shape(xyz))
