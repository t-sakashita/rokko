set(SAMPLE frank_matrix_c2)

foreach(name ${SAMPLE})
  add_executable(${name} ${name}.c frank_matrix_wrapper_c2.cpp wrapper_c2.cpp timer_wrapper_c2.cpp)
  target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${ELPA_LIB} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES})
# add_rokko_test(${name})
endforeach(name)
