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

#pragma once

#include <rokko/distributed_mfree.hpp>

namespace rokko {

class tridiagonal_toeplitz_mfree : public rokko::distributed_mfree_default {
public:
  tridiagonal_toeplitz_mfree(double a_in, double b_in, int dim, MPI_Comm comm = MPI_COMM_WORLD)
    : tridiagonal_toeplitz_mfree(a_in, b_in, dim, rokko::mpi_comm{comm}) {}

  tridiagonal_toeplitz_mfree(double a_in, double b_in, int dim, rokko::mpi_comm const& comm)
    : distributed_mfree_default(dim, mpi_comm{comm}),
      a(a_in), b(b_in),
      nprocs(get_mpi_comm().get_nprocs()), myrank(get_mpi_comm().get_myrank()),
      num_local_rows_(get_num_local_rows()), end_k_(num_local_rows_ - 1),
      need_send_recv_(get_num_local_rows() != dim),
      is_first_proc(start_row() == 0), is_last_proc(end_row() == dim),
      previous_rank_(previous_rank()), next_rank_(next_rank()) {}

  tridiagonal_toeplitz_mfree(double a_in, double b_in, int dim, int num_local_rows, rokko::mpi_comm const& comm)
    : distributed_mfree_default(dim, num_local_rows, mpi_comm{comm}),
      a(a_in), b(b_in),
      nprocs(get_mpi_comm().get_nprocs()), myrank(get_mpi_comm().get_myrank()),
      num_local_rows_(get_num_local_rows()), end_k_(num_local_rows_ - 1),
      need_send_recv_(get_num_local_rows() != dim),
      is_first_proc(start_row() == 0), is_last_proc(end_row() == dim),
      previous_rank_(previous_rank()), next_rank_(next_rank()) {}

  ~tridiagonal_toeplitz_mfree() = default;

  int next_rank() {
    const auto myrank = get_mpi_comm().get_myrank();
    const auto nprocs = get_mpi_comm().get_nprocs();

    Eigen::VectorXi target(get_mpi_comm().get_nprocs());

    MPI_Allgather(&num_local_rows_, 1, MPI_INT,
                  target.data(), 1, MPI_INT,
                  get_mpi_comm().get_comm());

    if (is_last_proc)  return -1;

    for (int rank = myrank+1; rank < nprocs; ++rank)
      if (target(rank) != 0)  return rank;  // if num_local_rows is not 0

    return -2;
  }

  int previous_rank() {
    const auto myrank = get_mpi_comm().get_myrank();
    const auto nprocs = get_mpi_comm().get_nprocs();

    Eigen::VectorXi target(get_mpi_comm().get_nprocs());

    MPI_Allgather(&num_local_rows_, 1, MPI_INT,
                  target.data(), 1, MPI_INT,
                  get_mpi_comm().get_comm());

    if (is_first_proc)  return -1;

    for (int rank = myrank-1; rank >= 0; --rank)
      if (target(rank) != 0)  return rank;  // if num_local_rows is not 0

    return -2;
  }

  void multiply(const double *const x, double *const y) const override {
    if (num_local_rows_ == 0) return;

    double buf_m, buf_p;
    MPI_Status status_m, status_p;

    if ((!is_first_proc) && need_send_recv_) {
      MPI_Send(x, 1, MPI_DOUBLE, previous_rank_, 0, get_comm());
      MPI_Recv(&buf_m, 1, MPI_DOUBLE, previous_rank_, 0, get_comm(), &status_m);
    }

    if ((!is_last_proc) && need_send_recv_) {
      MPI_Recv(&buf_p, 1, MPI_DOUBLE, next_rank_, 0, get_comm(), &status_p);
      MPI_Send(&x[end_k_], 1, MPI_DOUBLE, next_rank_, 0, get_comm());
    }

    if (is_first_proc) {
      if (num_local_rows_ > 1) {
        y[0] = a * x[0] + b * x[1];
        if (need_send_recv_) y[end_k_] = b * x[end_k_ - 1] + a * x[end_k_] + b * buf_p;
      } else {
        y[0] = a * x[0] + b * buf_p;
      }
    }

    if (is_last_proc) {
      if (num_local_rows_ > 1) {
        if (need_send_recv_) y[0] = b * buf_m + a * x[0] + b * x[1];
        y[end_k_] = a * x[end_k_] + b * x[end_k_ - 1];
      } else {
        y[end_k_] = a * x[end_k_] + b * buf_m;
      }
    }
    if (!(is_first_proc || is_last_proc)) { // neither first or last process
      if (num_local_rows_ > 1) {
        y[0] = b * buf_m + a * x[0] + b * x[1];
        y[end_k_] = b * x[end_k_ - 1] + a * x[end_k_] + b * buf_p;
      } else {
        y[0] = b * buf_m + a * x[0] + b * buf_p;
      }
    }
    // from 1 to end-1
    for (int k=1; k<end_k_; ++k) {
      y[k] = b * x[k-1] + a * x[k] + b * x[k+1];
    }
  }

  void fill_diagonal(double *const x) const override {
    if (num_local_rows_ == 0) return;

    for(int k=0; k<num_local_rows_; ++k)
      x[k] = a;
  }

private:
  double a, b;
  int nprocs, myrank;
  int num_local_rows_, end_k_;
  bool need_send_recv_;
  bool is_first_proc, is_last_proc;
  int previous_rank_, next_rank_;
};

} // namespace rokko
