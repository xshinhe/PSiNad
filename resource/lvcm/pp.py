import os, sys

fn = sys.argv[1]

with open(fn, 'r') as f:
    lines = f.readlines()
    cnt = 1
    _, F = lines[1].split()
    F = int(F)
    El = ['E_%d'%i for i in range(1,1+F)]
    Ev = ['%.4f'%float(lines[2+i]) for i in range(F)]
    for i in range(F):
        print('$%s$ & %s \\\\'%(El[i], Ev[i]))
    #print(F, ', '.join(El), ', '.join(Ev))
    _, N = lines[2+F].split()
    N = int(N)

    nc = 6
    nr = N//nc
    if N%nc!=0:
        nr += 1
    Wv = ['%.4f'%float(lines[3+F+i]) for i in range(N)]
    print("\\midrule[1pt]")
    print("\multirow{%d}{*}{$\\omega_{%d}\\sim\\omega_{%d}$}\n"%(nr,1,N), end='')
    for i in range(nr-1):
        print('& ' + ', '.join(Wv[nc*i:nc*i+nc]) + ', \\\\')
    print('& ' + ', '.join(Wv[nc*(nr-1):]) + '  \\\\')

    #Wl = ['\\omega_%d'%i for i in range(1,1+N)]
    #Wv = ['%.4f'%float(lines[3+F+i]) for i in range(N)]

    #print(', '.join(Wl), ', '.join(Wv))

