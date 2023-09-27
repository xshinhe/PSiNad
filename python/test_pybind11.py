import libopendf as odf
import json
import copy

print(dir(odf))
print(dir(odf.models))
print(dir(odf.solvers))

m_parm = {
        "unit": ["fs", "Bohr", "kcal", "K"],
        "Nb": 50,
        "nbath": 7,
        "rdim":350,
        "fdim":7,
        "coup_flag":"#se",
        "Hsys_flag":"#fmo",
        "spec_flag":"debye",
        "lamb_flag":"lambda",
        "strength": "35.0 wn",  # cm-1
        "omegac": "106.14 wn",  # cm-1
        "temp": "77.0 K", # (300K)
        "occ": 0,
        "tend": "1 ps",
        "dt": "0.1 fs"
    }
s_parm = {
    "gamma0" : 0, # Q-rep
    "ntraj": 1,
    "sstep": 20,
   }

# m = odf.models.Nad_ForceField(json.dumps(m_parm))
m = odf.models.SystemBath_ForceField(json.dumps(m_parm))
# s = odf.solvers.Solver(json.dumps(s_parm), m)


# m = odf.models.SystemBath_ForceField(json.dumps(m_parm))
s = odf.solvers.MMD_Solver(json.dumps(s_parm), m)
s.run()

