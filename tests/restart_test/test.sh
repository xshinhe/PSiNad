../../build/psinad -w -p param_FMO_NAFTW1_init5.json
#./psinad -w -p param_FMO_NAFTW1_init5.json -load=default/record-dump0-100.ds  # just load and no merge! may use for reading samplings.

# load from 100 steps and continue, which will recover the total trajectory from default dir, initial time start in 0
../../build/psinad -w -p param_FMO_NAFTW1_init5.json -load=default/record-dump0-100.ds:continue -d continue

# load from 100 steps and restart with different time length, which will generate a new trajectory in restart dir with initial time start in 100 steps
../../build/psinad -w -p param_FMO_NAFTW1_init5.json -load=default/record-dump0-500.ds:restart  -d restart
../../build/psinad -w -p param_FMO_NAFTW1_init5_diff.json -load=default/record-dump0-500.ds:restart -d restart2

