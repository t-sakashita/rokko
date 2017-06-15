#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# single process tests
set(TESTS solver)
foreach(name ${TESTS})
  set(target_name test_c_${name})
  add_executable(${target_name} ${name}.c)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  set_target_properties(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(${target_name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
  add_test(${name} ${name})
endforeach(name)
unset(TESTS)

# multi process tests
if(BUILD_PARALLEL_DENSE_SOLVER OR BUILD_PARALLEL_SPARSE_SOLVER)
  if(BUILD_PARALLEL_DENSE_SOLVER)
    set(TESTS grid)
    foreach(name ${TESTS})
      set(target_name test_c_${name})
      add_executable(${target_name} ${name}.c)
      set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
      set_target_properties(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
      target_link_libraries(${target_name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
      add_test(${name} ${name})
    endforeach(name)
    unset(TESTS)
  endif(BUILD_PARALLEL_DENSE_SOLVER)
endif(BUILD_PARALLEL_DENSE_SOLVER OR BUILD_PARALLEL_SPARSE_SOLVER)
