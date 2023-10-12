import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--save', type=str, default="corr-std.dat")
parser.add_argument('files', nargs='+')
args = parser.parse_args()

df = pd.read_csv(args.files[0], sep='\s+')
a = df.values
d1 = len(a[:,0])
d2 = len(a[0,:])
d0 = len(args.files)

all = np.zeros((d0,d1,d2))
for i in range(d0):
    all[i,:,:] = pd.read_csv(args.files[i], sep='\s+').values

std = np.std(all, axis=0)
pd.DataFrame(std).to_csv(args.save, header=df.columns, sep='\t', float_format='%.8e')

