/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/solver.hpp>
#include <rokko/mapping_bc.hpp>
#include <rokko/rokko_dense.h>

void rokko_mapping_bc_construct(struct rokko_mapping_bc* map, int global_dim, struct rokko_grid grid, struct rokko_parallel_dense_ev solver) {
  // Fix me: All parallel dense solvers use matrix_col_major. So we omit its checking, yet...
  map->ptr = new rokko::mapping_bc<rokko::matrix_col_major>(global_dim,
							    *static_cast<rokko::grid*>(grid.ptr),
							    *static_cast<rokko::parallel_dense_ev*>(solver.ptr));
  map->major = rokko_matrix_col_major;
}

  
void rokko_mapping_bc_construct_block_size(struct rokko_mapping_bc* map, int global_dim, int block_size) {
  // Fix me: All parallel dense solvers use matrix_col_major. So we omit its checking, yet...
  map->ptr = new rokko::mapping_bc<rokko::matrix_col_major>(global_dim, block_size);
  map->major = rokko_matrix_col_major;
}

void rokko_mapping_bc_destruct(struct rokko_mapping_bc* map) {
  delete static_cast<rokko::mapping_bc<rokko::matrix_col_major>*>(map->ptr);
  map->ptr = 0;
}

