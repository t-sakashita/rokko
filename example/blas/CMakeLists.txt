#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# set(PROG dgemv dgemm)
# foreach(name ${PROG})
#   set(target_name example_lapack_${name})
#   add_executable(${target_name} ${name}.c)
#   set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
#   target_link_libraries(${target_name} rokko)
#   add_test(${name} ${name})
# endforeach(name)
# unset(PROG)

set(PROG gemv gemm)
foreach(name ${PROG})
  set(target_name example_lapack_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko)
  add_test(${name} ${name})
endforeach(name)
unset(PROG)

# if(BUILD_Fortran_BINDING)
#   set(PROG dgemv_f dgemm_f)
#   foreach(name ${PROG})
#     set(target_name example_lapack_${name})
#     add_executable(${target_name} ${name}.f90)
#     set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
#     target_link_libraries(${target_name} rokko)
#     add_test(${name} ${name})
#   endforeach(name)
#   unset(PROG)
# endif(BUILD_Fortran_BINDING)
