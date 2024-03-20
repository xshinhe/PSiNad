import numpy as np
import pandas as pd
import os
import sys
import re

def getfiles(path):
    f_list_raw = os.listdir(path)
    f_list = []
    for i in f_list_raw:
        if os.path.splitext(i)[1] == '.ds':
            f_list += [i]
    return f_list

if __name__ == '__main__':
    path = sys.argv[1]
    out_path = sys.argv[2]
    from_size = int(sys.argv[3])
    samp_size = int(sys.argv[4])

    #files = getfiles(path)
    #Nf = len(files)
    Nf=from_size
    file_id = np.arange(Nf)
    np.random.shuffle(file_id)
    print(file_id)

    frp = []
    E1 = []
    E2 = []
    cnt = 0
    for i in file_id:
        print(cnt)
        cnt+=1

        frps = re.split("[ |\n]+", os.popen(
            "cat %s | grep -A2 f_rp | tail -n 1"%(path + '/samp%d.ds'%i)
            ).read().strip() )
        frp += [float(frps[1])]

        Es = re.split("[ |\n]+", os.popen(
            "cat %s | grep -A2 model.rep.E | tail -n 1"%(path + '/samp%d.ds'%i)
            ).read().strip() )
        E1 += [float(Es[0])]
        E2 += [float(Es[1])]

    frp = np.array(frp)
    E1  = np.array(E1)
    E2  = np.array(E2)
    v   = frp / (E2-E1)**2

    pd.DataFrame(np.vstack([file_id, frp, E1, E2, v]).T).to_csv(
        out_path + '_dat.csv',
        header=['id', 'frp', 'E1', 'E2', 'v'],
        index=None
    )

    maxv = np.max(v)

    samp_files = []
    while len(samp_files) < samp_size:
        print('#'*10)

        r = np.random.random(size=(Nf))
        idx = np.where(v > r*maxv)[0] 
        samp_files += list(file_id[idx])

    samp_files = samp_files[:samp_size]
    print(samp_files)

    os.popen('mkdir %s'%out_path)
    cnt = 0
    for i in samp_files:
        os.popen('cp %s %s'%(path + '/samp%d.ds'%i, out_path + '/init%d.ds'%cnt))
        cnt += 1

