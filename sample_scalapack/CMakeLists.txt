cmake_minimum_required(VERSION 2.8)

set(TESTS  diagonalize_scalapack diagonalize_scalapack_noeigvecs)
# frank_old diagonalize2)

foreach(name ${TESTS})
  add_executable(${name}  ${name}.cpp)
  target_link_libraries(${name}  "-lmkl_scalapack_lp64" "-lmkl_blacs_intelmpi_lp64" "-mkl=sequential")
  target_link_libraries(${name} ${MPI_CXX_LIBRARIES})
  add_rokko_test(${name})
endforeach(name)

