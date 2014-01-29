set(BENCHMARKS frank_mpi heisenberg_dense_mpi)

foreach(name ${BENCHMARKS})
  set(target_name benchmark_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES})
endforeach(name)
