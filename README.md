# PSiNaD: Phase Space Integrated Nonadiabatic Dynamics

![](docs/img/PSiNaD.png)

![version](https://img.shields.io/badge/version-0.0.1-green)  ![language](https://img.shields.io/badge/language-C%2B%2B%2FPython-blue)  ![license](https://img.shields.io/badge/license-Peking%20University-blue)

## Overview
**PSiNaD (Phase Space Integrated Nonadiabatic Dynamics)** is an open-source computational framework specifically designed for nonadiabatic dynamics simulations, with a focus on phase space integration methodologies. It aims to provide a flexible, efficient, and extensible platform for researchers in condensed matter physics, computational chemistry, and related fields to explore complex dynamical processes involving electronic nonadiabatic effects.

Nonadiabatic dynamics, which describes transitions between electronic states, is crucial for understanding phenomena such as photoexcitation, charge transfer, and chemical reactions. PSiNaD integrates state-of-the-art theoretical approaches and computational tools to enable accurate and efficient simulations of these processes across diverse systems.


## Objectives
The primary goals of PSiNaD are:
- To provide a unified framework for implementing and testing phase space-integrated nonadiabatic dynamics methods.
- To bridge classical and quantum dynamics approaches, supporting a spectrum of methodologies from Newtonian/Hamiltonian mechanics to advanced quantum approximations.
- To facilitate easy integration with external electronic structure packages and force fields, enabling multiscale simulations.
- To offer a user-friendly interface for setting up simulations and a modular architecture for developers to extend functionality.


## Key Features
- **Comprehensive Dynamics Methods**: Supports mixed quantum-classical dynamics (e.g., CMM/NaF-TW phase-space mapping methods, Ehrenfest dynamics, surface hopping), as well as quantum dynamics approaches (quantum phase space approximations).
- **Flexible Workflow Management**: Intuitive tools for constructing and controlling simulation workflows (see [Simulation Flow Control](docs/manu/flow.md)).
- **Extensibility for Developers**: Modular design simplifies adding new algorithms, models, and solvers. Guides for extending functionality are available in [Developing Models](docs/dev/dev_models.md) and [Adding Solvers](docs/dev/dev_solvers.md).
- **Multi-Scale Integration**: Diverse interfaces to external force fields and ab initio quantum chemistry packages (e.g., BAGEL, BDF, Gaussian, ORCA, PySCF, MNDO; see full list in `scripts/psnd_arg.py`).
- **Parallel Computing Support**: MPI-enabled C++ frontend for efficient parallel simulations (not available in the Python wrapper).
- **Cross-Platform Compatibility**: Tested on Linux (Ubuntu 20.04 LTS + Intel oneAPI) and macOS (Sequoia + Open-mpi 5.0.6); Windows support is in development.
- **Python Integration**: Available as a Python package ([PyPSiNaD](docs/api/python.md)) for enhanced accessibility and scripting capabilities.


## Installation
Refer to the detailed [Installation Guide](docs/installation.md) for step-by-step instructions. Prerequisites include:
- Unix-like OS (Linux/macOS; Windows pending)
- CMake 3.17+
- C++ compiler supporting C++17 (GCC, Clang, Intel C++ Compiler)
- Optional dependencies: MPI, MKL, gflags, glog, pybind11 (for Python bindings), Catch2 (for testing)

Basic installation steps:
```bash
# Clone the repository with submodules
git clone --recurse-submodules http://path_to_repository/PSiNaD.git
cd PSiNaD

# Create and navigate to build directory
mkdir build && cd build

# Configure (customize with -D options as needed)
cmake ..

# Compile and install
make -j$(nproc)
make install

# Optional: Install Python wrapper
# pip install pybind11 Cython
make PythonInstall
```


## Basic Usage
The C++ frontend `psinad_mpi` supports MPI for parallel simulations. A typical parallel run command:
```bash
mpirun -np 32 ./psinad_mpi -w -d output_dir -dump final -p param.json
```
- `mpirun -np 32`: Specifies 32 MPI processes (adjust based on your queue system, e.g., `srun -N 32` for Slurm).
- `-w`: Overwrite existing output directory.
- `-d output_dir`: Set output directory.
- `-dump final`: Dump final results.
- `-p param.json`: Path to the parameter file.

The full description of parameters can be found in the [introduction to input file](docs/manu/input.md)

For Python usage, refer to the [PyPSiNaD API](docs/api/python.md).


## Documentation
- **User Manual**: Detailed guides on [models](docs/manu/manu_models.md), [solvers](docs/manu/manu_solvers.md), [nonadiabatic dynamics](docs/manu/manu_nad.md), and more in the [Manual](docs/manual.md).
- **Tutorials**: Get started with [tutorial.md](tutorial.md), covering model and real-system simulations.
- **APIs**: Documentation for C++ backend ([C++ APIs](docs/api/cpp.md)) and Python frontend ([Python APIs](docs/api/python.md)).
- **Development**: Guidelines for contributors in [Development Docs](docs/development.md) (code style, optimization, etc.).


## Contributors
- Xin He (<xshinhe@pku.edu.cn>, <hexin@bjzgca.edu.cn>)
- Haocheng Lu (<1800011742@pku.edu.cn>)
- Members of the Jian Liu Group, CCME, Peking University


## License
Copyright © Peking University. Maintained by the Jian Liu Group, College of Chemistry and Molecular Engineering (CCME), Peking University. For inquiries, contact [jianliupku@pku.edu.cn](mailto:jianliupku@pku.edu.cn).

The PSiNaD framework is free for use and integration, with specific algorithms potentially requiring proper attribution as per academic standards.


## Support
For questions, bug reports, or feature requests, please contact the maintainers via the group email or submit an issue on the repository.


| Read Next |
|-----------|
| [Installation](docs/installation.md) |
| [Tutorial](tutorial.md) |
| [Developer Guide](docs/development.md) |
