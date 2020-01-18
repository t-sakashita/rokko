/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/parallel_sparse_ev.hpp>
#include <rokko/distributed_mfree.hpp>
#include <rokko/sparse.h>

class distributed_mfree_f : public rokko::distributed_mfree_default {
public:
  distributed_mfree_f(void (*multiply)(const int*, const double *const, double *const), int dim, MPI_Comm comm = MPI_COMM_WORLD)
    : rokko::distributed_mfree_default(dim, rokko::mpi_comm{comm}), multiply_(multiply) {
  }
  ~distributed_mfree_f() {}

  void multiply(const double* x, double *const y) const {
    int num_local_rows = get_num_local_rows();
    multiply_(&num_local_rows, x, y);
  }

private:
  void (*multiply_)(const int*, const double *const, double *const);
};

void rokko_distributed_mfree_f_construct(struct rokko_distributed_mfree* matrix,
					 void (*multiply)(const int*, const double *const, double *const),
					 int dim, int comm_f) {
  MPI_Comm comm = MPI_Comm_f2c(comm_f);
  matrix->ptr = new distributed_mfree_f(multiply, dim, comm);
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

