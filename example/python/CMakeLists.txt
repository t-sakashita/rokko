#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015-2016 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

enable_testing()

if(BUILD_Python_BINDING)
  add_subdirectory(dense)
  if(BUILD_PARALLEL_SPARSE_SOLVER)
    add_subdirectory(sparse)
  endif(BUILD_PARALLEL_SPARSE_SOLVER)
endif(BUILD_Python_BINDING)
