#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

include_directories(${PROJECT_SOURCE_DIR}/rokko)
include_directories(${PROJECT_BINARY_DIR}/rokko)

set(EXAMPLE frank_array matrix012)
foreach(name ${EXAMPLE})
  set(target_name example_fortran_dense_generate_matrix_${name})
  add_executable(${target_name} ${name}.f90)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko)
  if(ROKKO_USE_CXX_LINKER)
    set_target_properties(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
  endif(ROKKO_USE_CXX_LINKER)
endforeach(name)

if(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  set(EXAMPLE minij_array_mpi frank_array_mpi matrix012_mpi)
  foreach(name ${EXAMPLE})
    set(target_name example_fortran_dense_generate_matrix_${name})
    add_executable(${target_name} ${name}.f90)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko)
    if(ROKKO_USE_CXX_LINKER)
      set_target_properties(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
    endif(ROKKO_USE_CXX_LINKER)
  endforeach(name)
endif(ROKKO_HAVE_PARALLEL_DENSE_SOLVER)
