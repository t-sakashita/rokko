#
# Rokko: Integrated Interface for libraries of eigenvalue decomposition
#
# Copyright (C) 2012-2014 by Synge Todo <wistaria@comp-phys.org>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

configure_file(xyz.dat.in ${CMAKE_CURRENT_BINARY_DIR}/xyz.dat COPYONLY)

set(SAMPLE dsyev)
foreach(name ${SAMPLE})
  set(target_name sample_diagonalize_${name})
  add_executable(${target_name} ${name}.c)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko)
  add_rokko_test(${name})
endforeach(name)

set(SAMPLE frank_dense heisenberg_dense xyz_dense all_gather_dense_mpi gather_dense_mpi frank_dense_mpi heisenberg_dense_mpi xyz_dense_mpi)
foreach(name ${SAMPLE})
  set(target_name sample_diagonalize_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} rokko)
endforeach(name)
