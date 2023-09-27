set(CMAKE_VERBOSE_MAKEFILE ON)
enable_language(Fortran)

set(CMAKE_Fortran_FLAGS "-fPIE")

file(GLOB_RECURSE h2o2_srcs *.c *.f *.f90)

add_library(h2o2 ${h2o2_srcs})

