/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/parallel_sparse_ev.hpp>
#include <rokko/distributed_crs_matrix.hpp>
#include <rokko/sparse.h>

void rokko_distributed_crs_matrix_construct(struct rokko_distributed_crs_matrix* matrix, int dim1, int dim2,
					    struct rokko_parallel_sparse_ev solver) {
  matrix->ptr = new rokko::distributed_crs_matrix(dim1, dim2, *static_cast<rokko::parallel_sparse_ev*>(solver.ptr));
}

void rokko_distributed_crs_matrix_destruct(rokko_distributed_crs_matrix* matrix) {
  delete static_cast<rokko::distributed_crs_matrix*>(matrix->ptr);
  matrix->ptr = 0;
}

void rokko_distributed_crs_matrix_insert(struct rokko_distributed_crs_matrix matrix, int row, int col_size, int* cols, double* values) {
  static_cast<rokko::distributed_crs_matrix*>(matrix.ptr)->insert(row, col_size, cols, values);
}

void rokko_distributed_crs_matrix_complete(struct rokko_distributed_crs_matrix matrix) {
  static_cast<rokko::distributed_crs_matrix*>(matrix.ptr)->complete();
}

int rokko_distributed_crs_matrix_start_row(struct rokko_distributed_crs_matrix matrix) {
  return static_cast<rokko::distributed_crs_matrix*>(matrix.ptr)->start_row();
}

int rokko_distributed_crs_matrix_end_row(struct rokko_distributed_crs_matrix matrix) {
  return static_cast<rokko::distributed_crs_matrix*>(matrix.ptr)->end_row();
}

int rokko_distributed_crs_matrix_num_local_rows(struct rokko_distributed_crs_matrix matrix) {
  return static_cast<rokko::distributed_crs_matrix*>(matrix.ptr)->num_local_rows();
}

void rokko_distributed_crs_matrix_print(struct rokko_distributed_crs_matrix matrix) {
  return static_cast<rokko::distributed_crs_matrix*>(matrix.ptr)->print();
}
