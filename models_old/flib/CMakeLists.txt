enable_language(Fortran)

add_definitions("-fPIC")

add_subdirectory(H2O)
add_subdirectory(H2O2)
add_subdirectory(NH3)
add_subdirectory(paraH2)
add_subdirectory(water)
add_subdirectory(liquidNe)

# add_custom_target(flib)
# target_link_libraries(flib h2o h2o2 nh3 parah2)
