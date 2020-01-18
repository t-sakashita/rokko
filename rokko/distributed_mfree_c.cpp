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
#include <rokko/distributed_mfree.hpp>
#include <rokko/sparse.h>

class distributed_mfree_c : public rokko::distributed_mfree_default {
public:
  distributed_mfree_c(void (*multiply)(const double *const, double *const, void*), void* vars, int dim, MPI_Comm comm)
    : rokko::distributed_mfree_default(dim, rokko::mpi_comm{comm}), multiply_(multiply), vars_(vars) {
  }
  ~distributed_mfree_c() {}

  void multiply(const double *const x, double *const y) const {
    multiply_(x, y, vars_);
  }

private:
  void (*multiply_)(const double *const, double *const, void*);
  mutable void* vars_;
};

void rokko_distributed_mfree_construct(struct rokko_distributed_mfree* matrix,
				       void (*multiply)(const double *const, double *const, void*),
				       void* vars,
				       int dim, MPI_Comm comm = MPI_COMM_WORLD) {
  matrix->ptr = new distributed_mfree_c(multiply, vars, dim, comm);
}

void rokko_distributed_mfree_destruct(rokko_distributed_mfree* matrix) {
  delete static_cast<distributed_mfree_c*>(matrix->ptr);
  matrix->ptr = nullptr;
}

int rokko_distributed_mfree_dim(struct rokko_distributed_mfree matrix) {
  return static_cast<rokko::distributed_mfree*>(matrix.ptr)->get_dim();
}

int rokko_distributed_mfree_num_local_rows(struct rokko_distributed_mfree matrix) {
  return static_cast<rokko::distributed_mfree*>(matrix.ptr)->get_num_local_rows();
}

int rokko_distributed_mfree_offset(struct rokko_distributed_mfree matrix) {
  return static_cast<rokko::distributed_mfree*>(matrix.ptr)->get_local_offset();
}

