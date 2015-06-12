#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

include_directories(${PROJECT_SOURCE_DIR}/rokko)
include_directories(${PROJECT_BINARY_DIR}/rokko)
include_directories(${MPI_Fortran_INCLUDE_PATH})

set(EXAMPLE frank)
foreach(name ${EXAMPLE})
  set(target_name example_fortran_${name})
  add_executable(${target_name} ${name}.f90)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko ${ELEMENTAL_LIBRARIES} ${ELPA_LIB} ${SCALAPACK_LIB} ${MPI_Fortran_LIBRARIES})
  if(ROKKO_USE_CXX_LINKER)
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
  endif(ROKKO_USE_CXX_LINKER)
endforeach(name)

if(BUILD_PARALLEL_DENSE_SOLVER OR BUILD_PARALLEL_SPARSE_SOLVER)
  set(EXAMPLE frank_mpi frank_mpi_set_global frank_mpi_set_local set_distributed_matrix)
  foreach(name ${EXAMPLE})
    set(target_name example_fortran_${name})
    add_executable(${target_name} ${name}.f90)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko ${ELEMENTAL_LIBRARIES} ${ELPA_LIB} ${SCALAPACK_LIB} ${MPI_Fortran_LIBRARIES})
    if(ROKKO_USE_CXX_LINKER)
      SET_TARGET_PROPERTIES(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
    endif(ROKKO_USE_CXX_LINKER)
  endforeach(name)

  set(EXAMPLE frank_array frank_array_all_gather frank_function frank_function_array)
  foreach(name ${EXAMPLE})
    set(target_name example_fortran_${name})
    add_executable(${target_name} ${name}.F90)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko ${ELEMENTAL_LIBRARIES} ${ELPA_LIB} ${SCALAPACK_LIB} ${MPI_Fortran_LIBRARIES})
    if(ROKKO_USE_CXX_LINKER)
      SET_TARGET_PROPERTIES(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
    endif(ROKKO_USE_CXX_LINKER)
  endforeach(name)
endif(BUILD_PARALLEL_DENSE_SOLVER OR BUILD_PARALLEL_SPARSE_SOLVER)
