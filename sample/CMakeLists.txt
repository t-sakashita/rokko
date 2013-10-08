set(SAMPLE frank_dense xyz_dense)

foreach(name ${SAMPLE})
  set(target_name sample_diagonalize_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko)
# add_rokko_test(${name})
endforeach(name)

set(SAMPLE all_gather_dense_mpi gather_dense_mpi frank_dense_mpi heisenberg_dense_mpi xyz_dense_mpi)

foreach(name ${SAMPLE})
  set(target_name sample_diagonalize_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko ${ELEMENTAL_LIBRARIES} ${ELPA_LIB} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES})
# add_rokko_test(${name})
endforeach(name)
