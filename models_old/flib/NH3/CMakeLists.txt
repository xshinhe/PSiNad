set(CMAKE_VERBOSE_MAKEFILE ON)
enable_language(Fortran)

set(CMAKE_Fortran_FLAGS "-fPIE")

add_library(nh3 NH3-potential.f90) 
