#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
#                            Synge Todo <wistaria@comp-phys.org>,
#               2013-2013    Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# Tests using Boost.test
# Check and build Boost.Test Unit Test Framework
find_package(Boost 1.48.0 COMPONENTS unit_test_framework)
if (Boost_FOUND)
  add_definitions(-DBOOST_TEST_DYN_LINK)
endif (Boost_FOUND)

#copy files
#set(PATTERN_FILES "heisenberg.ip")
#foreach( pattern_file ${PATTERN_FILES} )
#  message(STATUS "Display_CURRENT_BINARARY: ${CMAKE_CURRENT_BINARY_DIR}")
#  message(STATUS "Display_CURRENT_SOURCE: ${CMAKE_CURRENT_SOURCE_DIR}")
#  message(STATUS "Display_pattern: ${CMAKE_CURRENT_SOURCE_DIR}/${pattern_file}")
#  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${pattern_file} ${CMAKE_CURRENT_BINARY_DIR}/${pattern_file} COPYONLY)
#endforeach( pattern_file )

#set(TESTS frank heisenberg xyz)
set(TESTS heisenberg)

foreach(name ${TESTS})
  set(target_name generate_matrix_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES})
  add_test(${target_name} ${name})
endforeach(name)

set(TESTS heisenberg_mpi)
foreach(name ${TESTS})
  set(target_name generate_matrix_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES} ${Boost_LIBRARIES})
  add_test(${name}-2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${CMAKE_CURRENT_BINARY_DIR}/${name})
  add_test(${name}-4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_CURRENT_BINARY_DIR}/${name})
endforeach(name)
unset(TESTS)

