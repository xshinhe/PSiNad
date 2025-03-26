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
'SB1': ('1.0', '0.1', '0.25'),
'SB2': ('1.0', '0.4', '0.25'),
'SB3': ('2.5', '0.1', '0.25'),
'SB4': ('2.5', '0.4', '0.25'),
'SB5': ('1.0', '0.1', '5.0'),
'SB6': ('1.0', '0.4', '5.0'),
'SB7': ('2.5', '0.1', '5.0'),
'SB8': ('2.5', '0.4', '5.0'),
}

Sdict = {
    "EHR": (
        "Focus", "0.0",
        "EHR",
        "1", "2", "false", "true", "false",
        "EHRrho<mij>:C(K0<mij>)",
    ),
    "SHA": (
        "Focus", "0.0",
        "BO",
        "2", "1", "false", "false", "true",
        "SHArho<mij>:C(KSHA<mij>)",
    ),
    "CPS0": (
        "Constraint", "0.0",
        "EHR",
        "1", "2", "false", "true", "false",
        "CPSrho<mij>:C(w, K1{@0}<m[occ][occ]>, K2<mij>)"
    ),
    "CPSw": (
        "Constraint", "-1.0",
        "EHR",
        "1", "2", "false", "true", "false",
        "CPSrho<mij>:C(w, K1{@0}<m[occ][occ]>, K2<mij>)"
    ),
    "CPSh": (
        "Constraint", "0.5",
        "EHR",
        "1", "2", "false", "true", "false",
        "CPSrho<mij>:C(w, K1{@0}<m[occ][occ]>, K2<mij>)"
    ),
    "CPS1": (
        "Constraint", "1.0",
        "EHR",
        "1", "2", "false", "true", "false",
        "CPSrho<mij>:C(w, K1{@0}<m[occ][occ]>, K2<mij>)"
    ),
    "FC": (
        "Focus", "-1.0",
        "EHR",
        "1", "2", "false", "true", "false",
        "FCrho<mij>:C(w, K1{@0}<m[occ][occ]>, K2<mij>)"
    ),
    "NAF0": (
        "Constraint", "0.0",
        "NAF",
        "1", "2", "false", "true", "false",
        "CPSrho<mij>:C(w, K1{@0}<m[occ][occ]>, K2<mij>)"
    ),
    "NAFw": (
        "Constraint", "-1.0",
        "NAF",
        "1", "2", "false", "true", "false",
        "CPSrho<mij>:C(w, K1{@0}<m[occ][occ]>, K2<mij>)"
    ),
    "NAFh": (
        "Constraint", "0.5",
        "NAF",
        "1", "2", "false", "true", "false",
        "CPSrho<mij>:C(w, K1{@0}<m[occ][occ]>, K2<mij>)"
    ),
    "TW": (
        "SQCtri", "0.333",
        "EHR",
        "1", "2", "false", "true", "false",
        "TWrho<mij>:C(w, KTWD{@0}<m[occ][occ]>, KTWD<mij>)"
    ),
    "NAFTW2": (
        "SQCtri", "0.333",
        "NAF",
        "1", "2", "false", "true", "false",
        "TWrho<mij>:C(w, KTWD{@0}<m[occ][occ]>, KTWD<mij>)"
    ),
    "NAFTW1": (
        "SQCtest01", "0.333333",
        "NAF",
        "1", "2", "false", "false", "false",
        "TWrho<mij>:C(w, KTWD{@0}<m[occ][occ]>, KTWD<mij>)"
    ),
    "MSSH": (
        "Constraint", "-2.0",
        "BO",
        "0", "0", "true", "false", "false",
        "MSSHrho<mij>:C(w, K1QD{@0}<m[occ][occ]>, K1<mij>)"
    ),
}

def get_param(mflag, S, sflag, k=0):
    return '''{
    "model": {
        "name": "SystemBath",
        "system_flag": "SB",
        "coupling_flag": "SB",
        "strength_flag": "Alpha",
        "bath_flag": "Ohmic",
        "bath_omegac": %s,
        "bath_strength": %s,
        "bath_temperature": "%s auK^-1",
        "bath_classic": false,
        "nbath": 1,
        "Nb": 100,
        "N": 100,
        "F": 2,
        "bias": 1.0,
        "delta": 1.0,
        "occ": 0,
        "dt": 0.01,
        "tend": 15.0
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
'''%(Mdict[mflag][0], Mdict[mflag][1], Mdict[mflag][2], 
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
    for m in Mdict:
        for k in range(len(Mdict[m][2])):
            print(cnt, begin, final)
            if cnt < begin or cnt >= final:
                print('xxx')
                pass #print("###")
            else:
                idir='SpinBosons'
                p = get_param(m, Sdict[s], s, k)
                try:
                    os.makedirs(idir)
                except IOError:
                    pass
                f = open('%s/param_%s_%s.json'%(idir,m,s), 'w')
                f.write(p)
                f.close()
                subjob(idir)
            cnt += 1

