#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2015 Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

add_subdirectory(lapack)

add_subdirectory(cxx)

if(BUILD_C_BINDING)
  add_subdirectory(c)
endif(BUILD_C_BINDING)

if(BUILD_Fortran_BINDING)
  add_subdirectory(fortran)
endif(BUILD_Fortran_BINDING)

if(BUILD_Python_BINDING)
  add_subdirectory(python)
endif(BUILD_Python_BINDING)
