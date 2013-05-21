set(BENCHMARKS frank_matrix)

foreach(name ${BENCHMARKS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko ${MPI_CXX_LIBRARIES})
# add_rokko_test(${name})
endforeach(name)
