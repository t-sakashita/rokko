#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

set(SAMPLE parameters default_solvers solvers)
foreach(name ${SAMPLE})
  set(target_name example_c_utility_${name})
  add_executable(${target_name} ${name}.c)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  set_target_properties(${target_name} PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(${target_name} rokko)
endforeach(name)
unset(SAMPLE)

