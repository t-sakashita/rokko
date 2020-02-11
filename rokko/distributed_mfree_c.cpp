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

#include <rokko/distributed_mfree.h>
#include <rokko/distributed_mfree.hpp>

void rokko_distributed_mfree_construct(struct rokko_distributed_mfree* matrix,
				       void (*multiply)(const double *const, double *const, void*),
				       void* vars,
				       int dim, MPI_Comm comm = MPI_COMM_WORLD) {
  matrix->ptr = new rokko::distributed_mfree_holder(multiply, vars, rokko::skel::mapping_1d{dim, rokko::mpi_comm{comm}});
}

void rokko_distributed_mfree_construct_no_context(struct rokko_distributed_mfree* matrix,
				       void (*multiply)(const double *const, double *const),
				       int dim, MPI_Comm comm = MPI_COMM_WORLD) {
  matrix->ptr = new rokko::distributed_mfree_holder(multiply, rokko::skel::mapping_1d{dim, rokko::mpi_comm{comm}});
}

void rokko_distributed_mfree_destruct(struct rokko_distributed_mfree* matrix) {
  delete static_cast<rokko::distributed_mfree_holder*>(matrix->ptr);
  matrix->ptr = nullptr;
}

int rokko_distributed_mfree_dim(struct rokko_distributed_mfree matrix) {
  return static_cast<rokko::distributed_mfree_holder*>(matrix.ptr)->get_dim();
}

int rokko_distributed_mfree_num_local_rows(struct rokko_distributed_mfree matrix) {
  return static_cast<rokko::distributed_mfree_holder*>(matrix.ptr)->get_num_local_rows();
}

int rokko_distributed_mfree_start_row(struct rokko_distributed_mfree matrix) {
  return static_cast<rokko::distributed_mfree_holder*>(matrix.ptr)->start_row();
}

int rokko_distributed_mfree_end_row(struct rokko_distributed_mfree matrix) {
  return static_cast<rokko::distributed_mfree_holder*>(matrix.ptr)->end_row();
}

int rokko_distributed_mfree_offset(struct rokko_distributed_mfree matrix) {
  return static_cast<rokko::distributed_mfree_holder*>(matrix.ptr)->get_local_offset();
}
