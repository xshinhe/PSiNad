# Basic {#dev_basic}

## Paramizer

The Paramizer is the container for parameter input, natively in JSON.  The JSON library has good implementations in C++, and we can simply encapsulate its usage. For Python/RUST, there are also libraries available for direct encapsulation.

```cpp
class Param;
using kids_param = Param;
```

There are two way for Param to load parameter list. One is from the string:

```cpp
auto PM = Param("{\"key\": value}", Param::fromString);
```

ans the other is loading from file:

```cpp
auto  pathtofile = "param_custom.json"
auto PM = Param(pathtofile, Param::fromFile);
```

Some keywords are used for control to simulation, which are list as follows, and all of them can be passed by command line arguments.

|   keywords    |      type      |     default    |   command line arguments  |
|:--------------|:---------------|:---------------|:--------------------------|
|   timing      |      <bool>    |     false      |   `-timing`               |
|   load        |      <string>  |     ""         |   `-load=...`             |
|   dump        |      <string>  |     ""         |   `-dump=...`             |



1) Other kind of formats

