set(CMAKE_VERBOSE_MAKEFILE ON)
enable_language(Fortran)

set(CMAKE_Fortran_FLAGS "-fPIE")

add_library(liquidne random.f90 pot.f90)
