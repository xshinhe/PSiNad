# PSiNad Overview

PSiNad is a kernel-tree based framework for nonadiabatic molecular dynamics. Algorithms are implemented as composable Kernels, assembled by Solvers, with Models providing the underlying physics. All runtime data is stored in a central DataSet; configuration comes from Param.

## Architecture

- Kernel (`kids/Kernel.h`)
  - Non-virtual entry points: `initializeKernel(stat)`, `executeKernel(stat)`, `finalizeKernel(stat)`
  - Virtual hooks: `*_impl(stat)`; the entry points call the corresponding hook and then (if enabled) call child kernels in order
  - Supports a tree of child kernels (composition over inheritance)
- DataSet (`kids/DataSet.h`)
  - Tree-structured container of tensors; variables defined via `DS->def(DATA::...)` returning `span<T>` views
- Param (`kids/Param.h`)
  - Parameter store (JSON today), load from file or string
- Status (`kids/Status.h`)
  - Execution flags and counters (first_step, frozen, istep, icalc, etc.)
- Dimensions (`kids/vars_list.h` + `src/core/src/vars_list.cpp`)
  - Global sizes `Dimension::M/P/N/F` and derived sizes (FF, NFF, ...)
  - `Kernel_Read_Dimensions` reads sizes from Param and calls `Dimension::static_build_shapes()`

## Execution flow (cpprun/psinad.cpp)

1) Build `Param` from input; sync gflags if needed
2) Create `Model` and two `Solver`s (Sampling + main), each exposing a root `Kernel` (a kernel tree)
3) For each root kernel: `setInputParam(PM) -> setInputDataSet(DS)`
4) In Monte Carlo loop, call in order:
   - `initializeKernel(stat)` -> `executeKernel(stat)` -> `finalizeKernel(stat)`
   - Each entry calls its `*_impl` then, if enabled, recurses to children

Note: Many kernels perform a one-off compute in `initializeKernel_impl` (sometimes calling `executeKernel(stat)` once).

## Representation & transforms

- `Kernel_Representation` manages electronic representations (Diabatic/Adiabatic)
- `transform(A, T, F, from, to, SpacePolicy)` converts quantities between bases
  - `SpacePolicy::L` (Liouvillian, matrix): similar transform (e.g., `A' = T^† A T`)
  - `SpacePolicy::H` (Hilbert, vector): single-side multiply

## Running

Typical build & run:

```bash
cmake -S . -B build -DPSINAD_BUILD_SHARED_LIB=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/psinad -p params.json -d out
```

Optional CMake switches:
- `PSINAD_BUILD_SHARED_LIB` / `PSINAD_BUILD_STATIC_LIB`
- `PSINAD_BUILD_TESTS` (enable internal Catch2 tests)
- `PSINAD_BUILD_EXAMPLES`
- `PSINAD_BUILD_PYTHON_WRAPPERS`

## Extending a Kernel

Implement a class deriving from `Kernel` and override:
- `setInputParam_impl(std::shared_ptr<Param>)`
- `setInputDataSet_impl(std::shared_ptr<DataSet>)`
- `initializeKernel_impl(Status&)`
- `executeKernel_impl(Status&)`

Attach it to a solver’s kernel tree (e.g., during solver assembly via `appendChild`).

## Minimal testing approach (standalone)

For an isolated kernel test (without touching the built-in tests):
- Set `Dimension::M/P/N/F`, then call `Dimension::static_build_shapes()`
- Prepare `Param` and `DataSet`; define all required `DATA::...` variables with `DS->def`
- Create the kernel, `setInputParam(pm)`, `setInputDataSet(ds)`, and call `initializeKernel(stat)` or `executeKernel(stat)`

A small example is provided in `test_kernel_naforce/` (independent mini CMake target) that validates `Kernel_NAForce` on a 2-level system.

## Typical kernel contract (simplified)

- Inputs: Param (configuration), DataSet (tensors), Status (flow control)
- Outputs: DataSet mutated (e.g., forces, densities), Status updated
- Error modes: missing variables/dimensions raise exceptions (`kids_error`/asserts)
