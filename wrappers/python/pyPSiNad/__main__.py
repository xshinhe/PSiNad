import pyPSiNaD

platform = pyPSiNaD.usePlatform('CPU')
PM = pyPSiNaD.Param('input.json', pyPSiNaD.Param.fromFile)

DS = pyPSiNaD.DataSet() # build as null dataset
model = pyPSiNaD.testmodels.PyModel_SBtest('TEST 1')
# pyPSiNaD.modelfactory('systembath')

# eqv. model = pyPSiNaD.SpinbosonModel()
# system = pyPSiNaD.System(model, PM, DS)
solver1 = pyPSiNaD.onlySampling(model) # / from system
solver2 = pyPSiNaD.defaultSolverFactory('NAD', model)

# eqv. solver2 = pyPSiNaD.NADSolver(model)
# modification for Param
appl1 = pyPSiNaD.appl.PopulationTransfer(model) # can merged in Param
# appl2 = pyPSiNaD.appl.LinearAbsoptionSpectrum(model)

solver2.addApplication(appl1) # append to recorder
# solver2.addApplication(appl2)

# read Param & DataSet
solver1.setInputParam(PM) # can be omitted if solver is constructed by system
solver1.setInputDataSet(DS) # can be omitted if solver is constructed by system
solver2.setInputParam(PM)
solver2.setInputDataSet(DS)


stat = pyPSiNaD.Status()
solvers = [solver1, solver2] # sequence of solvers
context = pyPSiNaD.Context(platform, system , solvers)
context.run(stat)

print(DS)
DS.dump('show.ds')