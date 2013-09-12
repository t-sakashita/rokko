set(SAMPLE frank)

foreach(name ${SAMPLE})
  set(target_name sample_diagonalize_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko)
# add_rokko_test(${name})
endforeach(name)

set(SAMPLE diagonalize_frank_matrix)

foreach(name ${SAMPLE})
  set(target_name sample_diagonalize_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko ${ELEMENTAL_LIBRARIES} ${ELPA_LIB} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES})
# add_rokko_test(${name})
endforeach(name)
