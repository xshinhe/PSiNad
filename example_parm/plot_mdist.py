import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

plt.rcParams.update({
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 10,
    'font.size': 12,
})

parser = argparse.ArgumentParser()
parser.add_argument('--xmin', type=float, default=-5)
parser.add_argument('--xmax', type=float, default=+5)
parser.add_argument('--Nxgrid', type=int, default=101)
parser.add_argument('--pmin', type=float, default=0)
parser.add_argument('--pmax', type=float, default=50)
parser.add_argument('--Npgrid', type=int, default=101)
parser.add_argument('--alpha', type=float, default=0.001)
parser.add_argument('--save', type=str, default="")
parser.add_argument('--flag', type=str, default="10pxgrid0")
parser.add_argument('--stripes', type=int, default=1)
parser.add_argument('files', nargs='+')
args = parser.parse_args()

if args.flag[0] == "R":
    args.stripes = 2

def get_data_from_df(df, flag0, size, stripes=1):
    i = df.columns.get_loc(flag0)
    return df.iloc[:,i:i+stripes*size:stripes].values

x = np.linspace(args.xmin, args.xmax, args.Nxgrid)
sx_list = []
for f in args.files:
    df = pd.read_csv(f, sep='\s+')
    sx = get_data_from_df(df, args.flag, args.Nxgrid, args.stripes)
    sx_list += [sx]

sx = np.vstack(sx_list)
Ns = len(sx[:,0])
p = np.linspace(args.pmin, args.pmax, args.Npgrid)
pIx = np.einsum('p,s,x->psx', p, np.ones((Ns)), x)
Isx = np.einsum('p,sx->psx', np.ones((args.Npgrid)), sx)
dist = np.einsum('psx->p', np.exp(1j*(pIx-Isx) - args.alpha*(pIx-Isx)**2))

norm = np.sum(dist.real) * (p[1]-p[0])
plt.plot(p, dist.real / norm)
plt.xlim(args.pmin, args.pmax)
plt.xlabel("Momentum (a.u.)")
plt.ylabel("Probability")

if(args.save != ""):
    plt.savefig(args.save, bbox_inches = 'tight', dpi=300)
    pd.DataFrame(np.transpose(np.vstack([p, dist.real]))).to_csv(args.save + '.dat', index=None, header=['p', 'dist'], sep='\t')
plt.show()

