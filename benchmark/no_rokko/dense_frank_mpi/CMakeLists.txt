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

# User need to specify include dir of ELPA's installation by ELPA_INCLUDE_DIR.
# It is only for no_rokko version of ELPA program.
include_directories(${ELPA_INCLUDE_DIR})

set(EXAMPLE scalapack_pdsyevd scalapack_pdsyevd_fixed_block scalapack_pdsyev eigen_exa elemental)
foreach(name ${EXAMPLE})
  set(target_name benchmark_dense_no_rokko_frank_mpi_${name})
  file(GLOB SOURCE "${name}.f*" "${name}.F*" "${name}.cpp")
  add_executable(${target_name} ${SOURCE})
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} ${ELEMENTAL_LIBRARIES} ${EIGENEXA_LIBRARIES} ${SCALAPACK_LIBRARIES} ${Boost_LIBRARIES})
  if(ROKKO_USE_CXX_LINKER)
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
  endif(ROKKO_USE_CXX_LINKER)
endforeach(name)

set(EXAMPLE elpa1_fixed_block elpa2_fixed_block elpa2)
foreach(name ${EXAMPLE})
  set(target_name benchmark_dense_no_rokko_frank_mpi_${name})
  file(GLOB SOURCE "${name}.f*" "${name}.F*" "${name}.cpp")
  add_executable(${target_name} ${SOURCE})
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} ${ELPA_LIBRARIES})
  if(ROKKO_USE_CXX_LINKER)
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
  endif(ROKKO_USE_CXX_LINKER)
endforeach(name)

