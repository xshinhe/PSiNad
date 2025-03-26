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
'U_CAM': (30, 8, [1,4,5]),
'U_PBE0': (30, 7, [1,4,6]),
'T_CAM': (39, 7, [1,4,5]),
'T_PBE0': (39, 7, [1,4,5]),
'C_CAM': (33, 7, [0,4]),
'C_PBE0': (33, 7, [0,3]),
'A_CAM': (39, 6, [1,2]),
'A_PBE0': (39, 7, [1,2]),
'G_CAM': (42, 9, [0,3]),
'G_PBE0': (42, 9, [1,2]),
'7HG_CAM': (42, 6, [0,4]),
'7HG_PBE0': (42, 6, [0,4]),
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
        "dt": "0.05 fs",
        "tend":"300.0 fs"
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
'''%(mflag, Mdict[mflag][0], Mdict[mflag][1], Mdict[mflag][2][k], 
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
    #if s not in ['MSSH']:
    #    continue
    for m in Mdict:
        for k in range(len(Mdict[m][2])):
            print(cnt, begin, final)
            if cnt < begin or cnt >= final:
                print('xxx')
                pass #print("###")
            else:
                idir='Nucleobases'
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

