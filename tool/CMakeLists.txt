#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2014-2014 by Synge Todo <wistaria@comp-phys.org>,
#                            Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

set(PROG rokko_solvers)
if(BUILD_ANASAZI)
  include_directories(${ANASAZI_INCLUDE_DIR})
endif(BUILD_ANASAZI)

foreach(name ${PROG})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} rokko)
  add_rokko_test(${name})
endforeach(name)
