#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
#                            Synge Todo <wistaria@comp-phys.org>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

include_directories(${PROJECT_SOURCE_DIR}/rokko)
include_directories(${PROJECT_BINARY_DIR}/rokko)
include_directories(${MPI_Fortran_INCLUDE_PATH})

set(SAMPLE frank frank_mpi)
foreach(name ${SAMPLE})
  set(target_name sample_fortran_${name})
  add_executable(${target_name} ${name}.f90)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko ${ELEMENTAL_LIBRARIES} ${ELPA_LIB} ${SCALAPACK_LIB} ${MPI_Fortran_LIBRARIES})
  if(NOT ROKKO_USE_Fortran_LINKER)
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES LINKER_LANGUAGE Fortran)
  endif(NOT ROKKO_USE_Fortran_LINKER)
endforeach(name)

set(SAMPLE frank_function frank_array frank_array_all_gather) # set_distributed_matrix)
foreach(name ${SAMPLE})
  set(target_name sample_fortran_${name})
  add_executable(${target_name} ${name}.F90)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko ${ELEMENTAL_LIBRARIES} ${ELPA_LIB} ${SCALAPACK_LIB} ${MPI_Fortran_LIBRARIES})
  if(NOT ROKKO_USE_Fortran_LINKER)
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES LINKER_LANGUAGE Fortran)
  endif(NOT ROKKO_USE_Fortran_LINKER)
endforeach(name)
