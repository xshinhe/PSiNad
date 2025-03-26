import sys
import os
import subprocess

if len(sys.argv) < 4:
    print("error")
    #exit(0)
else:
    mflag = sys.argv[1]
    sflag = sys.argv[2]
    gamma = sys.argv[3]

Mdict = {
'PYR3': (3, 2, [1], '"500.0 fs"'),
'PYR24': (24, 2, [1], '"500.0 fs"'),
'CRC2': (2, 3, [1], '"1000.0 fs"'),
'CRC5': (5, 3, [1], '"1000.0 fs"'),
'BUTA5': (5, 2, [1], '"500.0 fs"'),
'CED2': (200, 2, [1], '2200.0'),
'CED3': (200, 3, [2], '2200.0'),
'PYR2CED': (2, 4, [1], '"200.0 fs"'),
'BEN5': (5, 3, [2], '"300.0 fs"'),
}

Sdict = {
    "EHR": (
        "Focus", "0.0",
        "EHR",
        "1", "2", "false", "true", "false",
        "EHRpop<mi>:R(K0<mii>)",
    ),
    "SHA": (
        "Focus", "0.0",
        "BO",
        "2", "1", "false", "false", "true",
        "SHApop<mi>:R(KSHA<mii>)",
    ),
    "CPS0": (
        "Constraint", "0.0",
        "EHR",
        "1", "2", "false", "true", "false",
        "CPSpop<mi>:R(w, K1{@0}<m[occ][occ]>, K2<mii>)"
    ),
    "CPSw": (
        "Constraint", "-1.0",
        "EHR",
        "1", "2", "false", "true", "false",
        "CPSpop<mi>:R(w, K1{@0}<m[occ][occ]>, K2<mii>)"
    ),
    "CPSh": (
        "Constraint", "0.5",
        "EHR",
        "1", "2", "false", "true", "false",
        "CPSpop<mi>:R(w, K1{@0}<m[occ][occ]>, K2<mii>)"
    ),
    "CPS1": (
        "Constraint", "1.0",
        "EHR",
        "1", "2", "false", "true", "false",
        "CPSpop<mi>:R(w, K1{@0}<m[occ][occ]>, K2<mii>)"
    ),
    "FC": (
        "Focus", "-1.0",
        "EHR",
        "1", "2", "false", "true", "false",
        "FCpop<mi>:R(w, K1{@0}<m[occ][occ]>, K2<mii>)"
    ),
    "NAF0": (
        "Constraint", "0.0",
        "NAF",
        "1", "2", "false", "true", "false",
        "CPSpop<mi>:R(w, K1{@0}<m[occ][occ]>, K2<mii>)"
    ),
    "NAFw": (
        "Constraint", "-1.0",
        "NAF",
        "1", "2", "false", "true", "false",
        "CPSpop<mi>:R(w, K1{@0}<m[occ][occ]>, K2<mii>)"
    ),
    "NAFh": (
        "Constraint", "0.5",
        "NAF",
        "1", "2", "false", "true", "false",
        "CPSpop<mi>:R(w, K1{@0}<m[occ][occ]>, K2<mii>)"
    ),
    "TW": (
        "SQCtri", "0.333",
        "EHR",
        "1", "2", "false", "true", "false",
        "TWpop<mi>:R(w, KTWD{@0}<m[occ][occ]>, KTWD<mii>)"
    ),
    "NAFTW2": (
        "SQCtri", "0.333",
        "NAF",
        "1", "2", "false", "true", "false",
        "TWpop<mi>:R(w, KTWD{@0}<m[occ][occ]>, KTWD<mii>)"
    ),
    "NAFTW1": (
        "SQCtest01", "0.333333",
        "NAF",
        "1", "2", "false", "false", "false",
        "TWpop<mi>:R(w, KTWD{@0}<m[occ][occ]>, KTWD<mii>)"
    ),
    "MSSH": (
        "Constraint", "-2.0",
        "BO",
        "0", "0", "true", "false", "false",
        "MSSHpop<mi>:R(w, K1QD{@0}<m[occ][occ]>, K1<mii>)"
    ),
}

def get_param(mflag, S, sflag, k=0):
    return '''{
    "model": {
        "name": "LVCM",
        "lvcm_flag": "Read",
        "lvcm_file": "%s_lvcm.dat",
        "N": %d,
        "F": %d,
        "occ": %d,
        "dt": 0.01,
        "tend": %s
    },
    "solver": {
        "name": "NAD",
        "sampling_ele_flag": "%s",
        "sampling_nuc_flag": "Gaussian",
        "M": 96000,
        "msize": 4,
        "sstep": 20,
        "dump": "null",
        "gamma": %s,
        "representation_flag": "Adiabatic",
        "inp_repr_flag": "Diabatic",
        "ele_repr_flag": "Diabatic",
        "nuc_repr_flag": "Adiabatic",
        "naforce": "%s",
        "hopping_choose_type": %s,
        "hopping_direction_type": %s,
        "reflect": %s,
        "use_cv": %s,
        "use_fssh": %s,
        "offd_projected": true,
        "conserve_scale": true,
        "basis_switch": false
    },
    "record": [
        {"rule": "%s", "save": "%s"},
    ]
}
'''%(mflag, Mdict[mflag][0], Mdict[mflag][1], Mdict[mflag][2][k], Mdict[mflag][3], 
             S[0], S[1], S[2], S[3], S[4], S[5], S[6], S[7], S[8], 
             mflag+'_'+sflag+'.dat')

def subjob(idir):
    #os.system('sbatch task.sh %s'%idir)
    #subprocess.run('sbatch task.sh %s'%idir, shell=True)
    return

cnt = 0
begin=0
final=800
for s in Sdict:
    if s in ['CED2', 'CED3']:
       continue
    for m in Mdict:
        if m not in ['CED2', 'CED3']:
            continue
        for k in range(len(Mdict[m][2])):
            print(cnt, begin, final)
            if cnt < begin or cnt >= final:
                print('xxx')
                pass #print("###")
            else:
                idir='LVCMs'
                p = get_param(m, Sdict[s], s, k)
                try:
                    os.makedirs(idir)
                except IOError:
                    pass
                f = open('%s/param_%s_%s_init%d_from_read.json'%(idir,m,s,Mdict[m][2][k]), 'w')
                f.write(p)
                f.close()
                subjob(idir)
            cnt += 1

