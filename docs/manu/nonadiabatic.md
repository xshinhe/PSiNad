# Nonadiabatic Dynamics {#manu_nad}

This is an example for simulation of photon-deriven motor under MNDO calculation.

## Step 1: preperation of initial condition

make a new directory
```bash
mkdir samp
```

sampling initial condition

```bash
mpirun -np 8 ./kids -w -handler=sampling -d=samp -p param_sampling.json
```

Here the parameters file is like

```txt
{
    "model": "Interf_MNDO",
    "mndoinp" : "S0.inp",
    "exec_file": "mndo2020",
    "hess_log": "zp_hess.out",
    "init_nuclinp" : "#hess", # under Harmonic approximation
    "solver": "NAF-adapt",
    "F": 2,
    "N": 90,
    "dt": 0,
    "occ": 0,
    "temperature": "300.0 K",
    "classical_bath": false, # use wigner distribution
    "M": 1000,
}
```

where the S0.inp is the minimum of S0 state of the motor, while zp_hess.out are log file of Hessian calculation by this structure (making `JOP=2 KPRINT=1`). The S0.inp looks like

```txt
JOP=-2 IOP=-6 IGEOM=1 IFORM=1 ICUTS=-1 ICUTG=-1 NSAV7=4 NSAV13=2 +
ISCF=9 IPLSCF=9 DPREC=1D-8 DSTEP=1D-5 IPRINT=1 +
NCIGRD=2 IEF=1 IPREC=100 +
IMULT=1 imomap=3 mapthr=95 +
KCI=5 IOUTCI=1 MPRINT=1 IUVCD=2 ICROSS=7 +
MOVO=0 ICI1=6 ICI2=3 NCIREF=3 MCIREF=0 LEVEXC=2 IROOT=3 KITSCF=250  

This file is generated automatically 
6  0.04380338793550  0  0.35077864691677  0 -0.10882721806895  0  
6  1.27121168655916  0 -0.23979603265391  0 -0.10208501977209  0  
6  2.84226287946701  0 -1.92652319997327  0 -0.09579012490127  0  
6  2.65466322153261  0  0.34008820792743  0  0.04769625311877  0  
6 -1.66527000184424  0  2.15037595470014  0 -0.63643785262607  0  
6 -0.09037626881828  0  1.92025140395509  0 -0.29070923410677  0  
6  0.20040960553294  0  2.48113248777928  0  1.08479762064004  0  
8  2.85709734292449  0  1.62697083010599  0 -0.07327426390024  0  
1 -2.13756286310946  0  2.97755695887732  0 -0.03354811552012  0  
1 -1.79490161761576  0  2.33562930149058  0 -1.74596516994872  0  
1  0.58016298640936  0  2.54792774220182  0 -0.97734602010758  0  
1  1.17729539774138  0  2.25501127109489  0  1.43060350197030  0  
1 -0.57814464898532  0  2.05183223242283  0  1.79183498981984  0  
1  0.28460862204132  0  3.54051041950305  0  1.18549546962596  0  
1  3.09488609181207  0 -3.07973921031701  0 -0.13332221609733  0  
6  3.61236452215164  0 -0.75440449571611  0  0.04026746833597  0  
1  4.61142564115244  0 -0.52382073151543  0  0.15976722920344  0  
7  1.55843885154116  0 -1.63613870859904  0 -0.13877851677858  0  
1  0.98250807404859  0 -2.40645279632863  0 -0.52658707768494  0  
6 -3.34344788021834  0 -1.20692875542908  0 -0.02955012333173  0  
6 -3.56856116615990  0  0.22243729700008  0 -0.20050925408584  0  
6 -2.27946396293341  0  0.84912004673990  0 -0.33761888260090  0  
6 -1.27176847168875  0 -0.18466608507948  0 -0.07635502292883  0  
7 -1.91738252475071  0 -1.41686609558211  0  0.08211433345851  0  
1 -4.50249697318154  0  0.71308015177610  0 -0.60361324565497  0  
6 -1.28430653514342  0 -2.59490599037353  0  0.60434180765570  0  
1 -0.96827978921506  0 -3.20081206711561  0 -0.22748342879066  0  
1 -0.50639587335968  0 -2.34933980795280  0  1.36248208576779  0  
1 -2.06980017157749  0 -3.13442523272339  0  1.37893965343842  0  
1 -4.23683976135963  0 -1.89224882856376  0  0.09017577052527  0  
0  0.000000000000000 0  0.00000000000000  0  0.00000000000000  0
1  2
```

list `samp/` showing 
```bash
ls samp/
## fort.13  fort.7  fort.91  molden.dat  samp0.ds  samp1.ds  samp2.ds  samp3.ds  ... 
```
all `xxx.ds` is sampling by wigner distribution in Harmonic approximation. While if `classical_bath` is true, the Boltzman distribution is sampled instead.

However, the initial structures should be filtered by frequency strength, which is proceeded by a python script:
```python3
# filename select_init.py

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
    path        = sys.argv[1]
    out_path    = sys.argv[2]
    from_size   = int(sys.argv[3])
    samp_size   = int(sys.argv[4])

    #files = getfiles(path)
    #Nf = len(files)
    Nf=from_size
    file_id = np.arange(Nf)
    np.random.shuffle(file_id)
    print(file_id)

    frp = []
    E1  = []
    E2  = []
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
```

using `python3 select_init.py samp_dir filter_dir number_of_samp_ds number_of_filter_ds`, for example
```
python3 select_init.py samp init 1000 240
```
will filter 240 structures from 1000 sampled structures.

## Step 2: run nonadiabatic simulation

an example file is simulation carried out by NaF-TW:

```
// filename param_nad.json
{
    "model": "Interf_MNDO",
    "solver": "NAF-adapt",
    "cmsh_flag": "CVSH",
    "mndoinp" : "S0.inp",
    "exec_file": "mndo2020",
    "init_nuclinp" : "init/init0.ds",
    "F": 2,
    "N": 90,
    "occ": 1,
    "dt": "0.5 fs",
    "tend": "1000.0 fs",
    "time_unit": "1.0 fs",
    "msize": 1024,
    "sstep": 4,
    "M": 1,
    "representation_flag" : "Adiabatic",
    "gamma": 0.333333,
    "hopping_type1": 1,
    "hopping_type2": 2,
    "reflect": false,
    "use_cv" : true,
    "sqc_init": 0,
    "use_sqc" : true,
    "only_adjust": false,
    "conserve_scale": true,
    "result": [
        // the results you want to collect
    ]
}
```

by run
```txt
./kids -w -handler=parallel -d run0 -p param_nad.json
```

then each run will generate a `res.dat` file. And then can be treated with pandas.

## Step 3: analysis result


