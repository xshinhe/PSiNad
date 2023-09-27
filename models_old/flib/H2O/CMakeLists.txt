set(CMAKE_VERBOSE_MAKEFILE ON)
enable_language(Fortran)

set(CMAKE_Fortran_FLAGS "-fPIE")

add_library(h2o h2opot.f90)
