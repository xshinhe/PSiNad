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
'MORSE3A': '',
'MORSE3B': '',
'MORSE3C': '',
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
        "name": "NAD1D",
        "nad1d_flag": "%s",
        "N": 1,
        "F": 3,
        "occ": 0,
        "dt": "0.01 fs",
        "tend": "80.0 fs"
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
        "basis_switch": false,
        "time_unit": "1.0 fs"
    },
    "record": [
        {"rule": "%s", "save": "%s"},
    ]
}
'''%(mflag, # Mdict[mflag][2][k],
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
        for k in range(1):
            print(cnt, begin, final)
            if cnt < begin or cnt >= final:
                print('xxx')
                pass #print("###")
            else:
                idir='NAD1D'
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

