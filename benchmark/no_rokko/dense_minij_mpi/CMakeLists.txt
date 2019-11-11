#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

include_directories(${MPI_Fortran_INCLUDE_PATH})

if(ROKKO_HAVE_SCALAPACK)
set(EXAMPLE scalapack_pdsyevd scalapack_pdsyevd_fixed_block scalapack_pdsyev)
foreach(name ${EXAMPLE})
  set(target_name benchmark_dense_no_rokko_minij_mpi_${name})
  file(GLOB SOURCE "${name}.f*" "${name}.F*")
  add_executable(${target_name} ${SOURCE})
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko)
  if(ROKKO_USE_CXX_LINKER)
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
  endif(ROKKO_USE_CXX_LINKER)
endforeach(name)
endif(ROKKO_HAVE_SCALAPACK)

if(ROKKO_HAVE_ELEMENTAL)
set(EXAMPLE elemental elemental_thread_single)
foreach(name ${EXAMPLE})
  set(target_name benchmark_dense_no_rokko_minij_mpi_${name})
  file(GLOB SOURCE "${name}.cpp")
  add_executable(${target_name} ${SOURCE})
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko)
  if(ROKKO_USE_CXX_LINKER)
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
  endif(ROKKO_USE_CXX_LINKER)
endforeach(name)
endif(ROKKO_HAVE_ELEMENTAL)

if(ROKKO_HAVE_EIGENEXA)
set(EXAMPLE eigen_exa eigen_s)
foreach(name ${EXAMPLE})
  set(target_name benchmark_dense_no_rokko_minij_mpi_${name})
  file(GLOB SOURCE "${name}.f*" "${name}.F*")
  add_executable(${target_name} ${SOURCE})
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko)
  if(ROKKO_USE_CXX_LINKER)
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
  endif(ROKKO_USE_CXX_LINKER)
endforeach(name)
endif(ROKKO_HAVE_EIGENEXA)

if(ROKKO_HAVE_ELPA)
set(EXAMPLE elpa1 elpa2)
foreach(name ${EXAMPLE})
  set(target_name benchmark_dense_no_rokko_minij_mpi_${name})
  add_executable(${target_name} ${name}.F90 my_elpa_utility.F90 input_parameters.F90)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_include_directories(${target_name} PUBLIC ${ELPA_INCLUDE_DIR}/modules)  # Need .mod files in ELPA_INCLUDE_DIR
  target_link_libraries(${target_name} rokko)
  if(ROKKO_USE_CXX_LINKER)
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
  endif(ROKKO_USE_CXX_LINKER)
endforeach(name)
endif(ROKKO_HAVE_ELPA)
