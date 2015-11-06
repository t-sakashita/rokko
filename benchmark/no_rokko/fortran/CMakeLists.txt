#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

include_directories(${PROJECT_SOURCE_DIR}/rokko)
include_directories(${PROJECT_BINARY_DIR}/rokko)

set(EXAMPLE scalapack)# rokko)
foreach(name ${EXAMPLE})
  set(target_name benchmark_fortran_${name})
  file(GLOB SOURCE "${name}.f*")
  add_executable(${target_name} ${name}.f)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko)
  if(ROKKO_USE_CXX_LINKER)
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
  endif(ROKKO_USE_CXX_LINKER)
endforeach(name)

# add_executable(lapack lapack.f)
# target_link_libraries(lapack ${ROKKO_lapack_LIBRARIES})
# if(ROKKO_USE_CXX_LINKER)
#   SET_TARGET_PROPERTIES(lapack PROPERTIES LINKER_LANGUAGE CXX)
# endif(ROKKO_USE_CXX_LINKER)

# add_executable(scalapack scalapack.f)
# #target_link_libraries(scalapack ${ROKKO_LIBRARIES})
# target_link_libraries(scalapack ${ROKKO_SCALAPACKlapack_LIBRARIES})
# if(ROKKO_USE_CXX_LINKER)
#   SET_TARGET_PROPERTIES(scalapack PROPERTIES LINKER_LANGUAGE CXX)
# endif(ROKKO_USE_CXX_LINKER)

# add_executable(rokko rokko.f90)
# target_link_libraries(benchmark_rokko ${ROKKO_LIBRARIES})
# set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
# if(ROKKO_USE_CXX_LINKER)
#   SET_TARGET_PROPERTIES(benchmark_rokko PROPERTIES LINKER_LANGUAGE CXX)
# endif(ROKKO_USE_CXX_LINKER)

