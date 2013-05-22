set(TESTS  diagonalize_scalapack)
#diagonalize_scalapack_noeigvecs)
# frank_old diagonalize2)

if (DEFINED SCALAPACK_LIB)
  message(STATUS "Display: SCALAPACK_LIB = ${SCALAPACK_LIB}")
  separate_arguments(SCALAPACK_LIB)

  foreach(name ${TESTS})
    add_executable(${name}  ${name}.cpp)
    target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES}
${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES})
    add_rokko_test(${name})
  endforeach(name)
endif (DEFINED SCALAPACK_LIB)

