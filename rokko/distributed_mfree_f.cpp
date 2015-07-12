/*****************************************************************************
 *
 * Rokko: Integrated Interface for libraries of eigenvalue decomposition
 *
 * Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
 *                            Synge Todo <wistaria@comp-phys.org>,
 *                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
 *
 * Distributed under the Boost Software License, Version 1.0. (See accompanying
 * file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 *****************************************************************************/

#include <rokko/parallel_sparse_solver.hpp>
#include <rokko/distributed_mfree.hpp>
#include <rokko/rokko_sparse.h>

class distributed_mfree_f : public rokko::distributed_mfree {
public:
  distributed_mfree_f(void (*multiply)(int, const double*, double*), int dim, int num_local_rows)
    : multiply_(multiply), dim_(dim), num_local_rows_(num_local_rows), local_offset_(0) {
  }
  ~distributed_mfree_f() {}

  void multiply(const double* x, double* y) const {
    multiply_(num_local_rows_, x, y);
  }

  int get_dim() const { return dim_; }
  int get_local_offset() const { return local_offset_; }
  int get_num_local_rows() const { return num_local_rows_; }

private:
  void (*multiply_)(int, const double*, double*);
  int dim_;
  int num_local_rows_;
  int local_offset_;
};

void rokko_distributed_mfree_f_construct(struct rokko_distributed_mfree* matrix,
					 void (*multiply)(int, const double*, double*),
					 int dim, int num_local_rows) {
  matrix->ptr = new distributed_mfree_f(multiply, dim, num_local_rows);
}

void rokko_distributed_mfree_f_destruct(rokko_distributed_mfree* matrix) {
  delete static_cast<distributed_mfree_f*>(matrix->ptr);
}

/*
int rokko_distributed_mfree_dim(struct rokko_distributed_mfree* matrix) {
  return static_cast<rokko::distributed_mfree*>(matrix->ptr)->get_dim();
}

int rokko_distributed_mfree_num_local_rows(struct rokko_distributed_mfree* matrix) {
  return static_cast<rokko::distributed_mfree*>(matrix->ptr)->get_num_local_rows();
}

int rokko_distributed_mfree_offset(struct rokko_distributed_mfree* matrix) {
  return static_cast<rokko::distributed_mfree*>(matrix->ptr)->get_local_offset();
}
*/

