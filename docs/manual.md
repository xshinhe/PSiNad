# Manual {#manu}

@subpage manu_basic
@subpage manu_models
@subpage manu_solvers
@subpage manu_md
@subpage manu_nad
@subpage manu_abinitio
@subpage manu_qmmm
@subpage manu_semi
@subpage manu_wp

[TOC]


## Inputs and Outputs

## Usage

It accepts gflags-type command flag input (see more at https://github.com/gflags/gflags). 

or by:

```bash
./opendf --help # get help (something defined with gflags)
# opendf: [-s outputer][-r restart_file][-p paramemter_file][...]
```

the basic example is given here,
```bash
./opendf -p param.json
```

Note: opendf only run one traj for each time by default (so `ntraj` keyword is not used)! And you can you can parallelize PIMD or forcefield in each trajectory. However, you may want to parallelly run a batch of trajectories, where you can open flag `-para_flag=traj_para`, it will call solvers' `run_parallel` method, which simulate `ntraj` trajectories in each task.


Here `param.json` contains the parameters describing models and solvers. Some possible tests are available in `example_parm/` directory. For example, it at least should contains following keys:

```json
{ 
    // this is a comment line
    "model" : "smallmol",       // required flag to indicate the model
    "model_parm": {             // contains parameters to describe model
        "unit": {
            "time" : "ns",      // output time in ns
        },
        "rdim": 12,
        "type": "#h2o2",
    },
    "solver": "pimd",           // required flag to indicate the solver
    "solver_parm": {            // contains parameters to describe solver
        "gamma": 0.01,          // float is treated as a.u. by defualt
        "temp": "300 K",        // float is parsed with unit
        "nbead": 64,
        "dt": "0.1 fs",
        "tend": "2 ps",
        "nthstp": 1,
        "sstep": 10,
        "not used": 0,          // no-effect with additional paramters
    },
}
```

the `-s` flag inidicates outputers (by default `-s=ENER,TRAJ,ETRAJ,CORR`). Here `TRAJ/ETRAJ` writes each trajectory into a file, for nuclear DOFs and electronic DOFs resplectively. `CORR` gives time-correlation functions.

`SAMP/ESAMP`: writes all trajectories but for each step into a file, for nuclear DOFs and electronic DOFs respectively.

For example, if you only want to calculate time-correlation function by amount of trajectories, you needn't write trajectory or save restart files actually if you are confident with you calculation (because it will cost more time in disk IO), so you can disable `TRAJ` and `ETRAJ` by just reclaration the flags as `-s=CORR`. It will save mush of time.

the `-r` flag indicates retart file as additional input. The standart format for restart file is a HDF5 file, which records most important internal variables of the simulation at a certain time (we call it a context). Each solver will generate a context called `run.h5` by default. You can re-check or restart from a specific context by using flag `-r`.

Handling of restart and context:

1) simple run a restart 

```
mpirun -np 4 ./opendf -r saved.h5 -p param.json
```


## Models

## Solvers



<div class="section_buttons">

| Previous                        |                              Next |
|:--------------------------------|----------------------------------:|
| [Installation](installation.md) | [Development](development.md)     |
</div>