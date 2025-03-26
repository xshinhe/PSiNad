# Installation {#install}

[TOC]

## Basic Installation

This document provides a comprehensive guide for compiling and installing the PSiNad library using CMake, covering essential configuration options, build types, and installation procedures.

### 1. Prerequisites

Before proceeding, ensure that the following are installed on your system:

- **Unix-like OS**. (Win platform is not supported now)
- **CMake**: Version **3.17** or higher.
- **C++ Compiler**: (>= c++11, c++17 is better)
    - **GCC** (Linux/macOS)
    - **Clang** (Linux/mOS)
    - **MSVC** (Windows, TODO)
    - **Intel C++ Compiler** (optional)
- **Python Interpreter**: Required for building Python wrappers.
- **Additional Dependencies**: 
    - **MPI** (optional, for MPI-enabled executables)
    - **MKL** (optional, beter with intel oneAPI) 
    - **gflags** and **glog** (for logging and command-line flag parsing)
    - **pybind11** (for Python bindings)
    - **Catch2** (for unit testing, if enabled)
    - **ccache** (optional, for faster compilation)
    - **compile-time-perf** (optional, for performance measurements)

### 2. Build Configuration Options

The `CMakeLists.txt` file offers several options to customize the build process. Below are the key options and their descriptions:

- **Build Type**: 
    - **Release**: Optimized build with no debugging information (`-O3 -DNDEBUG`). 
        - **Default for UNIX systems** if not specified.
    - **Debug**: Includes debugging information and disables optimizations (`-g`).
        - **Set manually** if needed.
    - **Custom Flags**: You can override the default flags by setting `CMAKE_CXX_FLAGS_DEBUG` or `CMAKE_CXX_FLAGS_RELEASE`.

- **Library Type**: 
    - **Shared Library**: Enabled by default (`PSINAD_BUILD_SHARED_LIB=ON`). 
        - Generates a shared library (`libPSiNad.so` or `PSiNad.dll`).
    - **Static Library**: Enabled by default (`PSINAD_BUILD_STATIC_LIB=ON`). 
        - Generates a static library (`libPSiNad.a` or `PSiNad.lib`).

- **Backend Builds**: 
    - **C++ Backend**: Enabled by default (`PSINAD_BUILD_CXX_BACKENDS=ON`). 
        - Builds C++ executables (`psinad`, `psidyn`, `psisamp`) and their MPI-enabled counterparts (`psinad_mpi`, `psidyn_mpi`, `psisamp_mpi`).
    - **Python Wrappers**: Enabled by default (`PSINAD_BUILD_PYTHON_WRAPPERS=ON`). 
        - Builds Python bindings using `pybind11`.

- **Additional Features**: 
    - **Verbose Output**: Disabled by default (`ENABLE_VERBOSE_MAKEFILE=OFF`). 
        - Enabling it sets `CMAKE_VERBOSE_MAKEFILE` to `ON`, providing detailed build output.
    - **ccache**: Disabled by default (`CCACHE_ENABLE=OFF`). 
        - When enabled, `ccache` is used to speed up compilation by caching object files.
    - **compile-time-perf**: Disabled by default. 
        - If enabled and found, it integrates performance measurement tools during compilation.

- **Unit Testing**: Disabled by default (`PSINAD_BUILD_TESTS=OFF`). 
    - When enabled, it builds unit tests using `Catch2`.

- **Examples**: Disabled by default (`PSINAD_BUILD_EXAMPLES=OFF`). 
    - When enabled, it builds example executables.


### 3. Compilation Instructions

1. First you should clone this project with all submodules:

```bash
git clone --recurse-submodules http://path_to_this_repository/PSiNad.git 
cd PSiNad

# if you forget the `--recurse-submodules` when you clone the repository, you can do as follows:
git clone http://path_to_this_repository/PSiNad.git
cd PSiNad
git submodule update --init
```

2. Create a build directory: 
```bash
mkdir build
cd build
```

3. Configure the build: 
- Default configuration:
```bash
cmake ..
```
- Custom Configuration: 
   Enable/disable options using `-DOPTION_NAME=VALUE`.
   Example: 
    ```bash
    cmake -DPSINAD_BUILD_SHARED_LIB=ON -DPSINAD_BUILD_STATIC_LIB=OFF -DCCACHE_ENABLE=ON ..
    ```
- **Common Options**:
    - `CMAKE_BUILD_TYPE`: `Debug` or `Release`.
    - `ENABLE_VERBOSE_MAKEFILE`: `ON` or `OFF`.
    - `CCACHE_ENABLE`: `ON` or `OFF`.
    - `PSINAD_BUILD_TESTS`: `ON` or `OFF`.
    - `PSINAD_BUILD_EXAMPLES`: `ON` or `OFF`.

4. **Build the Project**: 
```bash
make -j$(nproc) # for c++ executables
make install 
```
The `-j$(nproc)` flag utilizes all available CPU cores for faster compilation.
```bash
make PythonInstall  # install python wrapper (optional)
```
This will build python wrapper for the libPSiNad.

4. **Run Tests (Optional)**: If unit tests are enabled:
```bash
ctest
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

### 4. Installation Instructions

1. **Install the Library and Executables**:
    ```bash
    make install
    ```
    - **Default Installation Path**:
        - **Windows**: `%ProgramFiles%\PSINAD`
        - **Linux/macOS**: `/usr/local/psinad`
    - **Custom Installation Path**: 
        - Set the `CMAKE_INSTALL_PREFIX` option during configuration.
            ```bash
            cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..
            ```
        - Alternatively, set the `PSINAD_INSTALL_PREFIX` variable if `CMAKE_INSTALL_PREFIX` is not initialized.

2. **Verify Installation**:
    - The library files (`libPSiNad.so` or `PSiNad.dll` / `libPSiNad.a` or `PSiNad.lib`) will be located in the `lib` directory.
    - The header files will be located in the `include` directory.
    - Executables (`psinad`, `psidyn`, `psisamp`, etc.) will be located in the `bin` directory.

### 5. Post-Installation

- **Environment Variables**: 
    - You may need to update your `PATH` (Windows) or `LD_LIBRARY_PATH` (Linux) environment variables to include the installation directory.
    - For Python wrappers, ensure that the Python `site-packages` directory is in your `PYTHONPATH`.

- **Documentation**: 
    - Refer to the official PSiNad documentation for detailed information on using the library and its features.

### 6. Additional Notes

- **In-Source Builds**: The `CMakeLists.txt` enforces out-of-source builds by checking if the source and build directories are the same.  If they are, the build will fail with an error message.

- **Compiler Compatibility**: The build system detects the compiler type and adjusts settings accordingly.  Ensure that your compiler is compatible with the C++17 standard.

- **Dependencies**:  The build system uses `find_package` to locate dependencies.  If dependencies are not found, ensure they are installed and accessible to CMake.

- **Customization**:  Advanced users can modify the `CMakeLists.txt` file to customize the build process further, such as adding new source directories, enabling additional features, or integrating with other tools.

By following this guide, you should be able to successfully compile and install the PSiNad library with the desired configuration.

## Benchmark Test


<div class="section_buttons">

| Previous          |                              Next |
|:------------------|----------------------------------:|
| [Home](README.md) | [Manual](manual.md)               |
</div>