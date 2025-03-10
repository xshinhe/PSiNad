import pyPSiNad

platform = pyPSiNad.usePlatform('CPU')
PM = pyPSiNad.Param('input.json', pyPSiNad.Param.fromFile)

DS = pyPSiNad.DataSet() # build as null dataset
model = pyPSiNad.testmodels.PyModel_SBtest('TEST 1')
# pyPSiNad.modelfactory('systembath')

# eqv. model = pyPSiNad.SpinbosonModel()
# system = pyPSiNad.System(model, PM, DS)
solver1 = pyPSiNad.onlySampling(model) # / from system
solver2 = pyPSiNad.defaultSolverFactory('NAD', model)

# eqv. solver2 = pyPSiNad.NADSolver(model)
# modification for Param
appl1 = pyPSiNad.appl.PopulationTransfer(model) # can merged in Param
# appl2 = pyPSiNad.appl.LinearAbsoptionSpectrum(model)

solver2.addApplication(appl1) # append to recorder
# solver2.addApplication(appl2)

# read Param & DataSet
solver1.setInputParam(PM) # can be omitted if solver is constructed by system
solver1.setInputDataSet(DS) # can be omitted if solver is constructed by system
solver2.setInputParam(PM)
solver2.setInputDataSet(DS)


stat = pyPSiNad.Status()
solvers = [solver1, solver2] # sequence of solvers
context = pyPSiNad.Context(platform, system , solvers)
context.run(stat)

print(DS)
DS.dump('show.ds')