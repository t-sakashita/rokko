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

class mfree_c : public rokko::distributed_mfree {
public:
  mfree_c(void (*multiply)(const double*, double*, void*), void* vars, int dim, int num_local_rows)
    : multiply_(multiply), vars_(vars), dim_(dim), num_local_rows_(num_local_rows), local_offset_(0) {
  }
  mfree_c(void (*simple)(), int dim, int num_local_rows)
    : simple_(simple), dim_(dim), num_local_rows_(num_local_rows), local_offset_(0) {
    printf("construct2\n");
    //    double x[5], y[5];
    //    multiply2_(&x[0], &y[0]);
    //printf("efefefe\n");
    //std::cout << "fun_pointer=" << simple << std::endl;
    printf ("func_address= %p\n", simple);
  }
  ~mfree_c() {}

  void multiply(const double* x, double* y) const {
    simple_();
    printf("oooooooooo\n");
    int n = num_local_rows_;
    multiply2_(n, x, y);
  }

  int get_dim() const { return dim_; }
  int get_local_offset() const { return local_offset_; }
  int get_num_local_rows() const { return num_local_rows_; }

public:

private:
  void (*simple_)();
  void (*multiply_)(const double*, double*, void*);
  void (*multiply2_)(int, const double*, double*);
  mutable void* vars_;
  int dim_;
  int num_local_rows_;
  int local_offset_;
};

void rokko_distributed_mfree_construct(struct rokko_distributed_mfree* matrix,
				       void (*multiply)(const double*, double*, void*),
				       void* vars,
				       int dim, int num_local_rows) {
  matrix->ptr = new mfree_c(multiply, vars, dim, num_local_rows);
}

void rokko_distributed_mfree_construct2(struct rokko_distributed_mfree* matrix,
					void (*simple)(),
					int dim, int num_local_rows) {
  matrix->ptr = new mfree_c(simple, dim, num_local_rows);
}

void rokko_distributed_mfree_destruct(rokko_distributed_mfree* matrix) {
  delete static_cast<mfree_c*>(matrix->ptr);
}

int rokko_distributed_mfree_dim(struct rokko_distributed_mfree* matrix) {
  return static_cast<rokko::distributed_mfree*>(matrix->ptr)->get_dim();
}

int rokko_distributed_mfree_num_local_rows(struct rokko_distributed_mfree* matrix) {
  return static_cast<rokko::distributed_mfree*>(matrix->ptr)->get_num_local_rows();
}

int rokko_distributed_mfree_offset(struct rokko_distributed_mfree* matrix) {
  return static_cast<rokko::distributed_mfree*>(matrix->ptr)->get_local_offset();
}

