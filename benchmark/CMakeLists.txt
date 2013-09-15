set(BENCHMARKS frank_mpi)

foreach(name ${BENCHMARKS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES})
endforeach(name)
