cmake_minimum_required(VERSION 2.8)

set(TESTS_SCALAPACK  diagonalize_time_frank_scalapack)
# frank_old diagonalize2)

foreach(name ${TESTS_SCALAPACK})
  add_executable(${name}  ${name}.cpp)
  target_link_libraries(${name}  "-lmkl_scalapack_lp64" "-lmkl_blacs_intelmpi_lp64" "-mkl=sequential")
  target_link_libraries(${name} ${MPI_CXX_LIBRARIES})
  add_alps_test(${name})
endforeach(name)


set(TESTS_ELEMENTAL  diagonalize_time_frank_elemental)

if( DEFINED ELEMENTAL_LIBRARIES AND ELEMENTAL_INCLUDE_DIRS )
  include_directories(${ELEMENTAL_INCLUDE_DIRS})
  foreach(name ${TESTS_ELEMENTAL})
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} ${ELEMENTAL_LIBRARIES}  "-mkl=parallel -lifcore")
    #add_alps_test(${name})
  endforeach(name)
endif( DEFINED ELEMENTAL_LIBRARIES AND ELEMENTAL_INCLUDE_DIRS )

if( DEFINED EIGEN_S_LIB AND EIGEN_S_INC )
  message(STATUS "Display: EIGEN_S_LIB= ${EIGEN_S_LIB}")
  separate_arguments(EIGEN_S_LIB)
  set(TESTS_EIGEN_S diagonalize_time_frank_eigen_s)

  foreach(name ${TESTS_EIGEN_S})
    include_directories(${name} ${EIGEN_S_INC})
    add_executable(${name} eigen_s_wrap.f  ${name}.cpp)
    target_link_libraries(${name} ${EIGEN_S_LIB})
  endforeach(name)
endif()


if( DEFINED EIGEN_SX_LIB AND EIGEN_SX_INC )
  message(STATUS "Display: EIGEN_SX_LIB= ${EIGEN_SX_LIB}")
  separate_arguments(EIGEN_SX_LIB)
  set(TESTS_EIGEN_SX diagonalize_time_frank_eigen_sx)

  foreach(name ${TESTS_EIGEN_SX})
    include_directories(${name} ${EIGEN_SX_INC})
    add_executable(${name} ${PROJECT_SOURCE_DIR}/sample_eigenK/eigen_sx_wrap.f  ${name}.cpp)
    target_link_libraries(${name} ${EIGEN_SX_LIB})
  endforeach(name)
endif()

