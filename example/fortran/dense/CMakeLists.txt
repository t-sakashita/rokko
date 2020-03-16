#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

include_directories(${PROJECT_SOURCE_DIR}/rokko)
include_directories(${PROJECT_BINARY_DIR}/rokko)

set(EXAMPLE minij minij_function frank frank_function frank_get_array_pointer helmert)
foreach(name ${EXAMPLE})
  set(target_name example_fortran_dense_${name})
  add_executable(${target_name} ${name}.f90)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko)
  if(ROKKO_USE_CXX_LINKER)
    set_target_properties(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
  endif(ROKKO_USE_CXX_LINKER)
endforeach(name)

if(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  set(EXAMPLE minij_mpi frank_array_all_gather_mpi frank_array_mpi frank_function_array_mpi frank_function_mpi frank_mpi frank_set_global_mpi frank_set_local_mpi frank_get_array_pointer_mpi frank_2groups_comm_mpi frank_group_comm_null_mpi frank_comm_split_mpi frank_cartesian_comm_mpi minij_set_local_mpi minij_set_local_f_mpi minij_set_global_mpi minij_set_global_f_mpi minij_function_mpi minij_function_f_mpi minij_get_array_pointer_mpi helmert_mpi)
  foreach(name ${EXAMPLE})
    set(target_name example_fortran_dense_${name})
    add_executable(${target_name} ${name}.f90)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko)
    if(ROKKO_USE_CXX_LINKER)
      set_target_properties(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
    endif(ROKKO_USE_CXX_LINKER)
  endforeach(name)
endif(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)

add_subdirectory(generate_matrix)
add_subdirectory(generate_vector)
