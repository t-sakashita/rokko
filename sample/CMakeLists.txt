set(SAMPLE diagonalize_frank_matrix)

foreach(name ${SAMPLE})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${ELPA_LIB} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES})
# add_rokko_test(${name})
endforeach(name)
