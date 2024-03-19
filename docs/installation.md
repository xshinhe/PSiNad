# Installation {#install}

[TOC]

## Basic Installation

**Prerequisites**:
- Unix-like operating system
- C++ compiler
- Python intepretor (optional for `libpykids.so`)
- Intel One API (MPI and MKL)
- cmake version >= 3.16

Here's a quick guide for installation:

```bash
mkdir build
cd build
CC=gcc CXX=g++ cmake .. -DCMAKE_BUILD_TYPE=Release
make -j8
# make install
```
If you prefer not to compile all components, you can specify modules in the config.json file. Note: refrain from modifying the default config.json; simply copy the file to ${CMAKE_BINARY_DIR}, and the latter will be loaded if it exists.

> **Problem Reports:**
>
> - [x] **internal compiler error**: If you encounter the `internal compiler error: Segmentation fault` message, consider enlarging the swap file and trying again using the following commands:
> ```bash
> sudo dd if=/dev/zero of=/swapfile bs=1M count=1024
> sudo mkswap /swapfile
> sudo swapon /swapfile
> # try compile again
> sudo swapoff /swapfile # don't forget to close swap file
> ```

kids requires the MKL and MPI libraries (Intel's oneAPI is recommended), which need to be manually configured in the CMakeLists.txt file.
```cmake
# for example
set(CUSTOMIZED_MKL_DIR "/opt/intel/oneapi/mkl/latest/")
set(CUSTOMIZED_MPI_DIR "/opt/intel/oneapi/mpi/latest/")
```
You can modify `cmake/FindMKLMod.cmake` and `cmake/FindMKLMod.cmake` to suit your PC environment.

## Benchmark Test


## PyKIDS

## rust-KIDS


<div class="section_buttons">

| Previous          |                              Next |
|:------------------|----------------------------------:|
| [Home](README.md) | [Manual](manual.md)               |
</div>