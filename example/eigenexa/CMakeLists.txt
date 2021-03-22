#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2021 by Rokko Developers https://github.com/t-sakashita/rokko
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

set(PROG eigen_s eigen_sx)
foreach(name ${PROG})
  set(target_name example_eigenexa_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko m)
  foreach(np 1 2 3 4 5 6)
    add_test(${target_name}-${np} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${name} ${MPIEXEC_POSTFLAGS})
  endforeach(np)
endforeach(name)
unset(PROG)

set(PROG eigen_s_c eigen_sx_c)
foreach(name ${PROG})
  set(target_name example_eigenexa_${name})
  add_executable(${target_name} ${name}.c)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko m)
  foreach(np 1 2 3 4 5 6)
    add_test(${target_name}-${np} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${name} ${MPIEXEC_POSTFLAGS})
  endforeach(np)
endforeach(name)
unset(PROG)

if(BUILD_Fortran_BINDING)
  set(PROG eigen_s_f eigen_sx_f)
  foreach(name ${PROG})
    set(target_name example_eigenexa_${name})
    add_executable(${target_name} ${name}.f90)
    set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
    target_link_libraries(${target_name} rokko)
    foreach(np 1 2 3 4 5 6)
      add_test(${target_name}-${np} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${name} ${MPIEXEC_POSTFLAGS})
    endforeach(np)
  endforeach(name)
  unset(PROG)
endif(BUILD_Fortran_BINDING)
