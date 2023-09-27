set(CMAKE_VERBOSE_MAKEFILE ON)
enable_language(Fortran)

set(CMAKE_Fortran_FLAGS "-fPIE")

add_library(qspcfw error.f90 ewald_image.f90 image.f90 force.f90)