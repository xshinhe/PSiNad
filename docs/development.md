# Development {#dev}

@subpage dev_basic
@subpage dev_models
@subpage dev_solvers

[TOC]

## For developers

1. For developers, it provides `add_definitions(USE_PL=ON)` to test and profile your code's numeric performance (see more at https://github.com/dfeneyrou/palanteer)  while for users but not nondevelopers, please just keep `add_definitions(USE_PL=OFF)`.

![](img/profiling.png)

2. And please use `clang-format` to standardize you code style with the config file `.clang-format`. (see more at https://clang.llvm.org/docs/ClangFormat.html)

3. Comments follow doxygen C/C++ style. (see more at https://www.doxygen.nl/manual/docblocks.html#cppblock)

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Four Basic Components

The KIDS framework is divided into the following four basic components: Paramizer, Dataizer, Algorithmizer, Taskizer.

1. **Paramizer**: The Paramizer can be implemented using JSON.

The JSON library has good implementations in C++, and we can simply encapsulate its usage. For RUST, there are also libraries available for direct encapsulation.

```c++
class Param;
```

2. **Dataizer**: We define the following classes:

For C++:

```
class Dimen;
class Shape;
class Node;
template<typename T> class Tensor;
class DataSet;
```

The Dimen class is based on dynamic dimensions of integers, such as `Dimen ndofs`. Operations implemented include bidirectional assignment, comparison operations, etc.

```c++
void operator=(const int& i, const Dimen& dim);
void operator=(const Dimen& dim, const int& i);
void operator==(const int& i, const Dimen& dim);
void operator==(const Dimen& dim, const int& i);
bool operator<(const int& i, const Dimen& dim);
bool operator<=(const Dimen& dim, const int& i);
bool operator>(const int& i, const Dimen& dim);
bool operator>=(const Dimen& dim, const int& i);
```

The Dimen class records each reference to its array shape, and updates the shape information promptly when Dimen values change.

The Shape class represents the shape of an array. It records a series of pointers to Dimens and precomputes additional information (such as the size of the primary dimension for each tensor index) for subsequent calculations. The design of the Shape class aims to facilitate the use of the Tensor class.

The Node class is an interface class designed to facilitate the hierarchical combination and management of a series of Tensors. Based on the Node class interface, memory management and tree-like directory indexing of the DataSet class (non-root nodes) and the Tensor class (root nodes) are performed using unique_ptr and polymorphism.

The DataSet provides three operations: `def()`, `at()`, and `undef()`, used to define a root Tensor based on a path in the DataSet. In addition to these, DataSet also defines functional functions such as `repr()`, `dump()`, and `load()` for structure representation, transfer, and loading.

To better pseudo-bind namespaces and DataSet structures in C++, we also introduce the following auxiliary class and macro functions:

```
class Variable<T>;
```

Here is an example of managing computational data structures:

```c++
// Pre-runtime overhead...
Dimen dim1; // Declare a dimension size (note that new operation is disabled)
Dimen dim2; // Declare another dimension size
Shape S12({&dim1, &dim2}); // Declare a Shape

namespace field3 {
Variable x("field3::x", &S12, "it is a variable under the namespace field3");
};

DataSet DS; // Declare a data set

// Runtime...

dim1 = 10; // Read value for dim1
dim2 = 20; // Read value for dim2
// S12 will be automatically updated here

int* a = DS.def<int>("field1::a", 10); // Define a block of int memory of size 10 in DS and return it
DS.def<double>("field1::field2::b", &S12);  // Define a block of double memory with shape S12 in DS
double* b = DS.def<double>("field1::field2::b", 200);  // If redefined, validate the size and return a pointer to the already defined memory directly
DS.def(field3::x); // Define memory using the Variable variable in namespace
// field3::x.data(); // Operate on memory in DS using the head pointer

std::cout << DS.repr() << "\n"; // Print DS to screen
std::cout << DS.at("field1::field2")->repr() << "\n"; // Access non-leaf node via at()
std::cout << *DS.at<double>("a::c::1") << "\n"; // Access head pointer of leaf node data memory via at()

DS.undef("field1::field2"); // Remove definition (non-leaf node or leaf node)
std::cout << DS.repr() << "\n"; // Print DS to screen

// Transfer DataSet to file
std::ofstream ofs{"test.ds"};
DS.dump(ofs);
ofs.close();

// Load file into new DataSet (note that Shape information is not retained)
std::ifstream ifs{"test.ds"};
DataSet DS2;
DS2.load(ifs);
ifs.close();
std::cout << DS2.repr() << "\n";

// We can further pseudo-bind DataSet and namespace using the following macro
#define DATASET_DEFINE_VARIABLE(type, name, shape, doc) \
    namespace name {                                    \
    VARIABLE<type> var(#name, shape, doc);              \
    };
```

The ds file format is linear, with read/write complexity of O(N*H), where N is the data volume and H is the maximum depth of the tree. It's a simple hdf5-like format.

3. **Algorithmizer**

The Algorithmizer, relative to the Dataizer and Paramizer, aims to separate algorithms from data itself.

The interface of Algorithmizer is Kernel (allowing empty Kernels to be built from scratch). A Kernel has a dynamic array of pointers to other Kernel, managed by shared_ptr for dynamic memory management, forming a tree-like calling structure. For example:
```text
[Index]      [Kernel]                                   [Time]  [Percent]
#02........: Kernel__CMM                                0.000s      0.00%
  #15......: Kernel_Load_DataSet                        0.000s      0.00%
  #16......: Kernel_Random                              0.000s      0.00%
  #17......: Kernel_Declare #01 #03 #04                 0.000s      0.00%
  #18......: Kernel_Initialize #01 #05 #11 #10          0.000s      0.00%
  #14......: Kernel_Iter                                0.000s      0.00%
    #10....: Kernel_Record                              0.000s      0.00%
    #04....: Kernel__BAOAB_Integrator                   0.000s      0.00%
      #07..: Kernel_Update_p                            0.000s      0.00%
      #08..: Kernel_Update_x                            0.000s      0.00%
      #08..: Kernel_Update_x                            0.000s      0.00%
      #01..: Model_NAD1D                                0.000s      0.00%
      #05..: Kernel_Representation                      0.000s      0.00%
      #09..: Kernel_Update_c                            0.000s      0.00%
      #11..: Kernel_Elec_CMM                            0.000s      0.00%
        #12: Kernel_Elec                                0.000s      0.00%
      #06..: Kernel_NADForce                            0.000s      0.00%
      #

07..: Kernel_Update_p                            0.000s      0.00%
      #13..: Kernel_Conserve                            0.000s      0.00%
      #03..: Kernel_Timer                               0.000s      0.00%
  #19......: Kernel_Dump_DataSet                        0.000s      0.00%
Using total time 12.0723 s
```

Kernel provides two calling functions:

```
template <typename T>
exec_kernel(T Data);

template <std::size_t level>
exec_kernel();
```

> Note: In the old approach, we divided the computation into four levels, namely `read_param(), init_data(), init_calc(), exec_kernel()` four functions, but this greatly limited the extensibility of the algorithm core. We define the execution of data transfer as the first-class calling interface, so `read_param(), init_data()` are respectively `exec_kernel<Param>()` and `exec_kernel<DataSet>()`, completing the interaction between Paramizer and Dataizer.
> The levels of calculation, such as sampling as pre-calculation, dynamics as main calculation, etc., can be controlled through the level.

Under the new interface structure, an instantiated dynamics calculation can be expressed as follows:
```c++
kernel = build_kernel(kernel_name);
kernel.exec_kernel<Param>(PM);   // Passing data parameters from top to bottom to different Algorithmizers (Algorithmizers reference these parameters)
kernel.exec_kernel<DataSet>(DS); // Passing data memory from top to bottom to different Algorithmizers (Algorithmizers reference this memory)
kernel.exec_kernel<0>(); // Execute level-0 calculation, such as initializing all algorithms (i.e., sampling)
kernel.exec_kernel<1>(); // Execute level-1 calculation, such as starting simulation
// ...
kernel.exec_kernel<1000>(); // Execute level-1000 calculation...
```

Reasons for abandoning object orientation in favor of Algorithmizer:

1) Effective separation of data and algorithms favors efficient program reuse.
2) Algorithmizer follows a functional programming approach, although its implementation in C++ also relies on polymorphism, it effectively truncates to subclass inheritance. The implementation of complex algorithms depends on the combination of Algorithmizers, rather than the hierarchical inheritance of objects. Excessive object encapsulation is not the essence of algorithm abstraction.
3) Atomicity of Algorithmizers facilitates maintenance tracking. For example, if the same algorithm has multiple implementations, `Kernel_Algo1_Impl1_1995, Kernel_Algo1_Impl2_2008`, the formal algorithm only needs to define `using Kernel_Algo1_Recommend = Kernel_Algo1_Impl2_2008`.
4) Facilitates modification, replacement, or new development. A Kernel can be replaced with a user-defined Kernel, such as

```c++
// realize customized Kernel_Algo1_Custom class
kernel = build_kernel(kernel_name); // default builder
kernel.replace(Kernel::enum::Algo1, std::shared_ptr<Kernel_Algo1_Custom>(new Kernel_Algo1_Custom()));
// ...
```

5) The construction of the algorithm tree is one-time and done before formal invocation, with minimal overhead and avoiding a large number of branch statements. It is highly efficient during formal invocation.
6) More straightforward and friendly interface with scripting languages; it can be exported to a Python module via pybind11.

#### Special Algorithmizers:
Priority container: Adjusts the priority of algorithms.
Model container: Recognizing that the computing system is based on the Model class, which is essentially a derivative of the Kernel class.
Collector container: XXX

4. **Collector and Taskizer**

As mentioned above, the collector is a special type of Kernel. It contains functions for processing data inside the DataSet and is used to create intermediate information, etc.

Taskizer, combined with the collector, can be used for general computing tasks (Applications). First, the collector determines the required intermediate information (default) based on the type of the Taskizer, and of course, additional intermediate information can be added by the user. During base calculation, information is collected.

```c++
class Collection;  // A type of Kernel, serves as the interface for the entire Kernel call tree to external Applications
class Application; // Taking hardware conditions into account, a specific computing task, such as different types of calculations, such as free energy calculation, spectrum calculation, reaction rate calculation, etc.
```

Note that due to the differences in MPI between C++ and Python, Application in C++ and Python are implemented separately.

## Plans

These components and their ideas are not limited to existing dynamic frameworks and can be extended to other fields.

### Molecular Force Field Engine Plan

Construct Kernel classes for ForceField to implement some common force fields/machine learning force fields.

### Electronic Structure Plan

Construct Kernel classes for AbInitio to build electronic structures such as DFT.

### External Interfaces

√ Gauss Interface

√ CASSCF Interface

√ MNDO Interface

DFT Interface

OpenMM Interface (only providing OpenMM interface)

### Standard for init_calc

All variable references must be defined in DataSet!!!
"

<div class="section_buttons">

| Previous                        |                              Next |
|:--------------------------------|----------------------------------:|
| [Manual](manual.md)             | [Examples](examples.md)           |
</div>