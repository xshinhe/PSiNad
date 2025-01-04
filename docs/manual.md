# Manual {#manu}

@subpage manu_models
@subpage manu_solvers
@subpage manu_nad
@subpage manu_abinitio
@subpage manu_qmmm
@subpage manu_semi
@subpage manu_wp

[TOC]


## Basic Usage

This (cpp-)PSiNad project comprises three standalone files: `kids`, `libcppkids.so`, and `libpykids.so`. Detailed documentation for the backend functionality provided by `libcppkids.so` can be found in the [C++ APIs](api/cpp.md), while information regarding the frontend functionality of `libpykids.so` is available in the [Python APIs](api/python.md). This manual primarily focuses on the C++ executable `kids`.

The kids C++ frontend leverages MPI for parallel processing, a feature not available in Python. Consequently, the Python APIs are implemented natively in Python to provide similar functionality. The simulation is controlled by two primary sources of parameters. The first source is command-line arguments, parsed using gflags (more details at https://github.com/gflags/gflags), which are specific to the C++ frontend. The second source is user-defined custom JSON configuration files, commonly utilized in both the C++ and Python frontends.

Here's a quick(/basic) example demonstrating parallel simulation with MPI:
```sh
mpirun -np 32 ./kids -handler=parallel -w -d output_dir -dump final -p param.json
```

In this example, `mpirun -np 32 ...` specifies arguments passed to MPI. Depending on the queue system used, this format may vary; for instance, in Slurm, you might use `srun -N 32 ...`.

On the other hand, `-handler=parallel -w -d output_dir -dump final -p param.json` are gflags passed to kids, the C++ executable. Here, `param.json` serves as the user's configuration file for parsing.

(Note: gflags supports four types of formats for specifying flags: `-flag=value`, `-flag value`, `--flag=value`, and `--flag value`.)

The task is managed by the parallel handler (it is the default handler), ensuring efficient parallel execution. All output results, including the final dumped file specified by `-dump`, can be found under the working directory `output_dir`.


## Inputs

To explore the full list of available command-line arguments enabled by gflags, execute:

```sh
./kids --help # get help (something defined with gflags)
```

This command will display:

```
kids: Kernel Integrated Dynamics Simulator

  -backup_time : double
      Specifies the timestep for backup (/1h)
      default: -1
  -d : string
      Specifies the output directory
      default: "default"
  -dump : string
      Specifies the dataset file to dump
      default: ""
  -handler : string
      Specifies the handler type
      [parallel | single | single_mpi | sampling | help | help_param |
      help_dataset ]
      default: "parallel"
  -load : string
      Specifies the dataset file to load
      default: ""
  -p : string
      paramemter inputs
      default: "param.json"
  -timing : bool
      Enables profiling for time costs
      default: false
  -w : bool
      rewrite output
      default: false
```

These flags offer precise control over simulation parameters and behavior, giving users enhanced flexibility and customization options. Some flags correspond to keywords in the JSON configuration file, governing simulation behavior. These keywords can be manipulated using command-line arguments. Notably, when both command-line arguments and JSON file values are provided, the command-line arguments take precedence. The correspondence between flags and keywords in the JSON configuration file is as follows:

|   Keyword     |   Type        |   Default     |   Command-line Arguments  |   Meaning                                 |
|:--------------|:--------------|:--------------|:--------------------------|:------------------------------------------|
|   handler     |   <string>    |   "parallel"  |   `-handler=...`          |   Specifies the handler type              |
|   directory   |   <string>    |   "default"   |   `-d=...`                |   Specifies the output directory          |
|   load        |   <string>    |   ""          |   `-load=...`             |   Specifies the dataset file to load      |
|   dump        |   <string>    |   ""          |   `-dump=...`             |   Specifies the dataset file to dump      |
|   backup_time |   <float>     |   -1          |   `-backup_time=...`      |   Specifies the timestep for backup (/1h) |
|   timing      |   <bool>      |   false       |   `-timing`               |   Enables profiling for time costs        |

 
Among the keywords (or gflags), one of the most crucial is the `handler` for the C++ frontend, which accepts the following possible values:

|   Type        |   Meaning                                                                                             |
|:--------------|:------------------------------------------------------------------------------------------------------|
|   parallel    |   Executes parallel independent Monte Carlo simulations (supports MPI) \[default\]                    |
|   single      |   Runs a single calculation only (MPI not supported)                                                  |
|   single_mpi  |   Runs a single calculation with MPI support, often used for complex force fields                     |
|   sampling    |   Performs initial sampling code without further kernel execution (MPI supported)                     |
|   help        |   Provides guidelines for models and solvers                                                          |
|   help_param  |   Offers comprehensive guidance for full parameter inputs of models/solvers, with detailed examples   |
|   help_dataset|   Provides guidance on accessible variables defined in the dataset                                    |


However, you are also required to configure your parameter file like param.json, specifying models, solvers, and other settings. If you're unsure about kids capabilities and how to write param.json, you can use the following command:
```sh
./kids -handler help
```
This command displays all possible models and solvers, along with their combinations. From this list, you can select the configurations needed for your simulation. For instance, if you find the LVCM model and NaF dynamics solver compatible and suitable, you can prepare your initial parameter JSON file as shown below:

```txt
{
    // This is a comment line enabled by kids (not native in JSON!)
    "handler": "parallel",
    "model": "LVCM",
    "solver": "NaF"
}
```
Based on this file, you'll also need to add details to the model_param and solver_param fields to instruct kids on how to set up the parameters of your models and solvers. However, you can rely on an automatic way to generate default fields by running:

```sh
./kids -handler=help_param -p param.json
```
This command will generate a new `param_gen.json` file in your workspace directory (specified by the `-d` flag). Finally, a full parameter file may look like this:

```txt
{
    "model": "LVCM",
    "model_param" : {
        "lvcm_flag": "PYR24",
        "N": 24,
        "F": 2,
        "occ": 1,
        "dt": "0.1 fs", // can appear both model_param and solver_param fields, but former is prior
        "tend": "500.0 fs",
        "unit_dt": "1.0 fs",
        "temperature": "300 K",
        "custom defined (not used)": 0,  // no-effect with additional paramters
    },
    "solver": "CMM",
    "solver_param": {
        "use_cv" : true,
        "sstep": 10,
        "N_mc": 1000,
        "dt": "0.2 fs",
    },
    "result": [
        // ...
    ],
    // additional keywors enable by command line arguments
    "handler": "parallel",
    "dump": "final", // it will generate a lot of trash file for parallel handler
    "directory": "default",
    // ...
}
```

As demonstrated above, the `param.json` file should contain parameters describing models and solvers. You can find some typical example tests in the `example_parm/` directory. For detail paramemters of various models and solvers, please refer to [Models](man/models.md) and [Solvers](man/solvers.md), resplectively.

Note there is another important field `result`, specifies the diverse quantities as outputs, which is discussed in the following subsection.

## Outputs

The gflags `-d`, `-dump`, `-backup_time`, and the keyword result in the custom JSON configuration file determine the content of the outputs. These outputs consist of two types of file formats, both located under the working directory.

One format is the dataset file (.ds), which you can dump during or after the simulation. For more details about the (.ds) format, please refer to [About the DataSet file format](man/dataset.md). (This format can also be converted into HDF5 format). The outputs in (.ds) format may appear as follows:
```txt
...
backup.1.occ_nuc 
kids_int        1
                1

backup.1.p 
kids_real       1  
1.26779265e+02

backup.1.x 
kids_real       1  
2.22689348e+01
...
```

The other format consists of (.dat) files containing data generated by the record kernel (Kernel_Record), following the rules specified in the "result" field. The contents of the (.dat) outputs may appear as follows:
```txt
               t            stat         I0Ekin0         I0Epot0         I0Etot0        ...
  0.00000000e+00               1  0.00000000e+00  4.24824397e-01  4.24824397e-01        ...
  5.00000000e-01               1  3.47060399e-03  4.21353793e-01  4.24824397e-01        ...
  1.00000000e+00               1  1.36565516e-02  4.11167845e-01  4.24824397e-01        ...
  1.50000000e+00               1  2.99105741e-02  3.94913823e-01  4.24824397e-01        ...
  2.00000000e+00               1  5.12481759e-02  3.73576221e-01  4.24824397e-01        ...
  2.50000000e+00               1  7.64675310e-02  3.48356866e-01  4.24824397e-01        ...
  3.00000000e+00               1  1.04281374e-01  3.20543023e-01  4.24824397e-01        ...
  3.50000000e+00               1  1.33438318e-01  2.91386079e-01  4.24824397e-01        ...
  4.00000000e+00               1  1.62817428e-01  2.62006969e-01  4.24824397e-01        ...
```

For example, if you want output the coordinate and momentum of each trajectory, you can use
```txt
{
    "result": [
        "$COORDINATE",  
        "$MOMENTUM", 
        // ... other keywords by orders
    ],
}
```
The intrinsically built keywords are listed as:

|   Keywords    |   Meaning                                                                                             |
|:--------------|:------------------------------------------------------------------------------------------------------|
|   $COORDINATE |   The coordinate of the nuclear DOFs                                                                  |
|   $MOMENTUM   |   The momentum of the nuclear DOFs                                                                    |
|   $VELOCITY   |   The velocity of the nuclear DOFs                                                                    |
|   $FORCE      |   The force of the nuclear DOFs                                                                       |
|   $RHO        |   The Ehrenfest's density of the electronic DOFs                                                      |
|   $KK_CPS     |   The kernel-kernel time correlation function by CPS (i.e., CMM)                                      |
|   $KK_SQC     |   The kernel-kernel time correlation function by SQC-TW (etc.)                                        |
|   $KK_CW1     |   The kernel-kernel time correlation function by CW1                                                  |
|   $KK_MM      |   The kernel-kernel time correlation function by Meyer-Miller form (including FOCUS dynamics)         |
|   $KK_GDTWA   |   The kernel-kernel time correlation function by GDTWA                                                |

**Advanced Usage**

Keywords prefixed with "$" denote several intrinsic rules. For more advanced control over the simulation outputs, you can utilize kids's result rule syntax.
- Without the "$" prefix, you can access any variable names defined in the DataSet. The variable names follow the name syntax `NAME{field@time}<ESSHAPE>`. In this approach, if the symbol is not in the DataSet, kids will report an error. (To view available variables in the DataSet, use the `-handler=help_dataset` flag).
- You can use the `NAMEOUT (NAME1, NAME2, ..., NAMEk) = FORMULA` format (where each name obeys the name syntax) to output custom-defined variables. If `NAMEOUT` is not defined in the DataSet, kids will assist in defining it.
- This syntax also allows for more complex rules, enabling both einsum rules and formula parsing. For examples of all possible types of rules in the result, refer to the documentation.

Full examples of different rules is shown as:
```txt
{
    "result": [
        // Intrinsic keywords
        "$...",

        // WARNING: The following usages are intended for experts only.

        // Short rules for single variables
        "x",        // x variable defined in field "integrator.x"
        "x{model}", // x variable defined in field "model.x"
        "x{@0}",    // x variable defined in field "integrator.x" but only valued at t=0 time
        "x<[0]>",   // 0-th component of the variable x

        // Short rules for compound variables
        "KK<ijkl>(K{@0}<ij>, K<kl>) = _1 * _2", // Output all KK correlations
        "KK<ij>  (K{@0}<ii>, K<jj>) = _1 * _2", // Output only diagonal KK correlations

        // A full rule example
        {
            "rule": "KK<ijkl>(K{@0}<ij>, K<kl>) = _1 * _2",
            "mode": "average",      // Modes: average, each, initial, final
            "save": "CMM_KK",       // Save as CMM_KK instead of res.dat
            "format": "default",    // Supports XYZ for coordinate for real systems
        }
    ]
}
```

## Restart

In certain scenarios, users may encounter situations where it becomes necessary to restart a simulation due to corruption or to continue from a specific point in the computation. In such cases, the `-load` flag in kids proves invaluable. Whether it's to rectify errors or seamlessly resume an ongoing simulation, `-load=xxx.ds` empowers users to efficiently recalculate or restart their simulations with ease. This feature not only ensures data integrity but also provides a streamlined workflow, allowing users to pick up where they left off without unnecessary complications or delays.


<div class="section_buttons">

| Previous                        |                              Next |
|:--------------------------------|----------------------------------:|
| [Installation](installation.md) | [Development](development.md)     |
</div>