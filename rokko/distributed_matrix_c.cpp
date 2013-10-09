/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>,
*                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/solver.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/rokko.h>

void rokko_distributed_matrix_construct(rokko_distributed_matrix* matrix, int dim1, int dim2,
                                        rokko_grid grid, rokko_solver solver, int matrix_major) {
  if (matrix_major == rokko_matrix_col_major)
    matrix->ptr = new rokko::distributed_matrix<rokko::matrix_col_major>(dim1, dim2,
      *static_cast<rokko::grid*>(grid.ptr), *static_cast<rokko::solver*>(solver.ptr));
  else
    matrix->ptr = new rokko::distributed_matrix<rokko::matrix_row_major>(dim1, dim2,
      *static_cast<rokko::grid*>(grid.ptr), *static_cast<rokko::solver*>(solver.ptr));
  matrix->major = matrix_major;
}
void rokko_distributed_matrix_destruct(rokko_distributed_matrix* matrix) {
  if (matrix->major == rokko_matrix_col_major)
    delete static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix->ptr);
  else
    delete static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix->ptr);
}

void rokko_distributed_matrix_generate_function(struct rokko_distributed_matrix* matrix,
  double (*func)(int i, int j)) {
  if (matrix->major == rokko_matrix_col_major)
    static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix->ptr)->generate(func);
  else
    static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix->ptr)->generate(func);
}

void rokko_distributed_matrix_print(rokko_distributed_matrix matrix) {
  if (matrix.major == rokko_matrix_col_major)
    static_cast<rokko::distributed_matrix<rokko::matrix_col_major>*>(matrix.ptr)->print();
  else
    static_cast<rokko::distributed_matrix<rokko::matrix_row_major>*>(matrix.ptr)->print();
}
