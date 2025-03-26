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
'FMO': ('1.0', '0.1', [0,5]),
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
        "system_flag": "FMO",
        "coupling_flag": "SE",
        "strength_flag": "Lambda",
        "bath_flag": "Debye",
        "bath_omegac": "106.14 wn",
        "bath_strength": "35.0 wn",
        "bath_temperature": "300 K",
        "bath_classic": false,
        "nbath": 7,
        "Nb": 50,
        "N": 350,
        "F": 7,
        "occ": %d,
        "dt": "0.1 fs",
        "tend": "1000.0 fs"
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
'''%(Mdict[mflag][2][k],
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
                idir='FMO'
                p = get_param(m, Sdict[s], s, k)
                try:
                    os.makedirs(idir)
                except IOError:
                    pass
                f = open('%s/param_%s_%s_init%d.json'%(idir,m,s,Mdict[m][2][k]), 'w')
                f.write(p)
                f.close()
                subjob(idir)
            cnt += 1

