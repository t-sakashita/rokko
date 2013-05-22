#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
#                            Synge Todo <wistaria@comp-phys.org>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

set(TESTS frank)
# set(TESTS frank scatter gather)
foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES})
  add_rokko_test(${name})
endforeach(name)

if(BUILD_LAPACKE)
  set(TESTS dsyev)
  foreach(name ${TESTS})
    add_executable(${name} ${name}.c)
    target_link_libraries(${name} rokko ${ELEMENTAL_LIBRARIES} ${SCALAPACK_LIB} ${MPI_CXX_LIBRARIES})
    add_rokko_test(${name})
  endforeach(name)
endif(BUILD_LAPACKE)
