# Installation {#install}

[TOC]

## Basic Installation

Prerequisite:
- unix-like OS
- a C++ compiler
- Intel One API (MPI and MKL)
- cmake >= 3.16

a fast and short guide for installation

```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
# or specify compilers like: CC=g++ FC=gfortran cmake ..
make -j8
# make install
```

If you don't want compile all components, you can specify modules in file `config.json`. Note: don't revised the default `config.json`, just copy the file to `${CMAKE_BINARY_DIR}`, and it will load the later one if it exists.

> **Problems Reports:**
>
> - [x] **internal compiler error**: if it reports `internal compiler error: Segmentation fault` error, please enlarge swap file and try again as following:
>
> ```bash
> sudo dd if=/dev/zero of=/swapfile bs=1M count=1024
> sudo mkswap /swapfile
> sudo swapon /swapfile
> # try compile again
> sudo swapoff /swapfile # don't forget to close swap file
> ```
>

It needs MKL and MPI library (Here I recommend Intel's oneAPI), which should be manually configued in `CMakeLists.txt`.

```cmake
# for example
set(CUSTOMIZED_MKL_DIR "/opt/intel/oneapi/mkl/latest/")
set(CUSTOMIZED_MPI_DIR "/opt/intel/oneapi/mpi/latest/")
```
you can revised `cmake/FindMKLMod.cmake` and `cmake/FindMKLMod.cmake` to adapt to environment in your PC.

## Benchmark Test


## PyKIDS

## rust-KIDS


<div class="section_buttons">

| Previous          |                              Next |
|:------------------|----------------------------------:|
| [Home](README.md) | [Manual](manual.md)               |
</div>