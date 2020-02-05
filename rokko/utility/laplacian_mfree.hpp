/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_UTILITY_LAPLACIAN_MFREE_HPP
#define ROKKO_UTILITY_LAPLACIAN_MFREE_HPP

#include <rokko/distributed_mfree.hpp>

namespace rokko {

class laplacian_mfree : public rokko::distributed_mfree_default {
public:
  laplacian_mfree(int dim, MPI_Comm comm = MPI_COMM_WORLD)
    : laplacian_mfree(dim, rokko::mpi_comm{comm}) {}

  laplacian_mfree(int dim, rokko::mpi_comm const& comm)
    : distributed_mfree_default(dim, mpi_comm{comm}),
      nprocs(get_mpi_comm().get_nprocs()), myrank(get_mpi_comm().get_myrank()),
      num_local_rows_(get_num_local_rows()), end_k_(num_local_rows_ - 1),
      is_first_proc(start_row() == 0), is_last_proc(end_row() == dim) {}

  ~laplacian_mfree() = default;

  void multiply(const double *const x, double *const y) const override {
    if (num_local_rows_ == 0) return;

    double buf_m, buf_p;
    MPI_Status status_m, status_p;

    if ((!is_first_proc) && (nprocs != 1)) {
      MPI_Send(x, 1, MPI_DOUBLE, myrank-1, 0, get_comm());
      MPI_Recv(&buf_m, 1, MPI_DOUBLE, myrank-1, 0, get_comm(), &status_m);
    }

    if ((!is_last_proc) && (nprocs != 1)) {
      MPI_Recv(&buf_p, 1, MPI_DOUBLE, myrank+1, 0, get_comm(), &status_p);
      MPI_Send(&x[end_k_], 1, MPI_DOUBLE, myrank+1, 0, get_comm());
    }

    if (is_first_proc) {
      if (num_local_rows_ != 1) {
        y[0] = x[0] - x[1];
        if (nprocs != 1) y[end_k_] = - x[end_k_ - 1] + 2 * x[end_k_] - buf_p;
      } else {
        y[0] = x[0] - buf_p;
      }
    }

    if (is_last_proc) {
      if (num_local_rows_ != 1) {
        if (nprocs != 1) y[0] = - buf_m + 2 * x[0] - x[1];
        y[end_k_] = 2 * x[end_k_] - x[end_k_ - 1];
      } else {
        y[end_k_] = 2 * x[end_k_] - buf_m;
      }
    }
    if (!(is_first_proc || is_last_proc)) { // neither first or last process
      if (num_local_rows_ != 1) {
        y[0] = - buf_m + 2 * x[0] - x[1];
        y[end_k_] = - x[end_k_ - 1] + 2 * x[end_k_] - buf_p;
      } else {
        y[0] = - buf_m + 2 * x[0] - buf_p;
      }
    }
    // from 1 to end-1
    for (int k=1; k<end_k_; ++k) {
      y[k] = - x[k-1] + 2 * x[k] - x[k+1];
    }
  }

private:
  int nprocs, myrank;
  int num_local_rows_, end_k_;
  bool is_first_proc, is_last_proc;
};

} // namespace rokko

#endif // ROKKO_UTILITY_LAPLACIAN_MFREE_HPP
