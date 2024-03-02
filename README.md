# KIDS

[TOC]

![version](https://img.shields.io/badge/version-0.0.1-g)  ![language](https://img.shields.io/badge/language-c++-orange) ![wrapper](https://img.shields.io/badge/wrapper-python3-red) ![platform](https://img.shields.io/badge/platform-linux-lightgrey) ![cmake](https://img.shields.io/badge/cmake->3.16-yellow) ![status](https://img.shields.io/badge/status-up-blue) ![uptime](https://img.shields.io/badge/uptime-100%25-g) ![coverage](https://img.shields.io/badge/coverage-100%25-pink) ![codecy](https://img.shields.io/badge/codecy-B-g) 

## Introduction

**KIDS (Kernel Integrated Dynamics Simulator)** offers an open-source framework tailored for simulating chemical and physical dynamics, with a primary focus on atomic and molecular scales in condensed matter. It is designed for (classical / qauntum) dynamics simulation of small system, few-body system, reduced systems, and even large many-particle (i.e. molecules & condensed matter) systems. It provides a versatile platform for the development of advanced algorithms, offering ease of use and accessibility at a minimal cost.

KIDS encompasses a diverse array of models and solvers. These include harmonic models, forcefield models, and ab initio interfaces, ensuring comprehensive coverage across various simulation needs. KIDS also caters to a wide spectrum of computational approaches, from classical dynamics rooted in Newtonian and Hamiltonian mechanics to advanced quantum dynamics featuring quantum phase space approximation, quantum trajectories, path integral techniques, influence functional methods, and semi-classical wavepackets etc. Furthermore, ongoing development efforts aim to expand the library of models and solvers, with contributions welcomed from the open community.

~~We've developed the KIDS project primarily in C++, known as **cpp-KIDS**, building upon earlier projects like OpenDF and MODA~~. 

However, we're now open to expanding its reach by integrating it as a Python package ([PyKIDS](docs/api/python.md)) and/or a Rust crate ([rust-KIDS](docs/api/rust.md)).

### License

~~The KIDS framework is free of use and integration, while some algorithms may be used with proper guarantee.~~

### Support

Author: Xin He(<xshinhe@pku.edu.cn>)

Contributors: 

The project is maintained under Jian Liu group, CCME, PeKing University.

## Features

- Utilizes high-performance and parallelized code written in modern C++ (-std=11 or later), well compiled with C-style APIs. Parallelization with MPI ([Intel OneAPI MPI](https://www.intel.cn/content/www/cn/zh/developer/tools/oneapi/toolkits.html)), OpenMP and SIMD techniques. For some functions, a GPU acceleration is also on-going and optional. Many thirdpart matrix/tensor libraries are well compiled and integrated, such as [Eigen](https://eigen.tuxfamily.org
), [Intel OneAPI MKL](https://www.intel.cn/content/www/cn/zh/developer/tools/oneapi/toolkits.html), xtensor, libtorch and python's numpy etc.
- Functions are modular and separable, all of which are intepreted as Kernel, encompassing models, solvers, and more. KIDS offers a variety of intrinsic models for testing and exploring new dynamics approaches, such as harmonic (quadratic / quartic / Morse) PES, system-bath models (e.g. spin-boson and FMO), LCVM models, light-matter models etc. KIDS also involves versatile solvers, including both Born-Oppenheimer Statistics & Dynamics and Non-adiabatic Statistics & Dynamics. See list of [what KIDS does](#autotoc_md3), while for more detail documention, refer to [intrinsic models](docs/manu/models.md) and [intrinsic solvers](docs/manu/solvers.md), respectively.
- Well-organized/seperate code is potential for custom use (most are head-only), please see [how to use KIDS code under proper LICNESE in a custom project](docs/api/libs.md). 
- A standalone dynamic library, such as `libkids.so`, is seperate with C++'s Handler and provides various useful APIs, which can be integrated into other langauges like python/rust. See [Python APIs](docs/api/python.md) and [rust APIs](docs/api/rust.md), respectively.
- It is clear and convenient for user to construct and manage the flow of the simulation. See more about the [control the flow of the simulation](docs/manu/flow.md).
- It is friendly for developer to add more advanced algorithms and organize the different versions. For integrating more models, see [how to develop models](docs/dev/dev_models.md), while for developing methods, see [how to add new solvers](docs/dev/dev_solvers.md)
- Diverse intrinsic interfaces support external force fields and ab initio calculations (as well as QM/MM). It can be interfaced to few-bodies force fields for small molecules in Fortran90 (i.e., [potlib](https://comp.chem.umn.edu/potlib)) as well as general molecular forcefield with OpenMM. Interfaces with Gaussian16 and MNDO99 for on-the-fly calculation are supported. Additionally, KIDS allows user add custom interface to any other software only based python script, without modifying the c++ code!

## What KIDS does

- MD and MD-based techniques (i.e., enhenced sampling) for classical and molecular dynamics.
- PIMD (techniques) for statistics for Boltzman particle, Boson/Fermion, and multi-electronic-state systems.
- Real-time dynamics: Full semi-classical dynamics up to the HK-IVR and LCSIVR(classical Wigner). Path-integral based RPMD, CMD, and PILD-variants.
- Quantum dynamics built on discrete basis with a dirty realization, such as DVR, multi-configuration time-dependent Hartree (MCTDH), and tensor decomposition methods.
- Quantum dynamics for system-bath problem built on influence functional, including Hierarchy equation of motion (HEOM), Stochastic Schrödinger equation (SSE), Path integral techniques (QCPI). Some approximation like Redfield equation is also supported.
- Featured: Quantum phase space approximation, classical-wavepackets and quantum hydrodynamics (TODO).
- Featured: Mixed quantum-classical dynamics based on phase-space mapping methods (CMM/NaF) and its extensions. Other approximation like Ehrenfest dynamics, multiple versions of surface hopping-based method are also available.
- Further to support ab initio calculation like Hartree-Fock/DFT realization.

<div class="section_buttons">

|                        Read Next |
|---------------------------------:|
| [Installation](docs/installation.md) |

</div>
