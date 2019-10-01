/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_UTILITY_LAPLACIAN_MFREE_HPP
#define ROKKO_UTILITY_LAPLACIAN_MFREE_HPP

#include "mpi.h"

namespace rokko {

class laplacian_mfree : public rokko::distributed_mfree {
public:
  laplacian_mfree(int dim) : dim_(dim) {
    comm_ = MPI_COMM_WORLD;
    MPI_Comm_size(comm_, &nprocs);
    MPI_Comm_rank(comm_, &myrank);

    int tmp = dim_ / nprocs;
    int rem = dim % nprocs;
    num_local_rows_ = (dim + nprocs - myrank - 1) / nprocs;
    start_row_ = tmp * myrank + std::min(rem, myrank);
    end_row_ = start_row_ + num_local_rows_ - 1;

    is_first_proc = start_row_ == 0;
    is_last_proc = end_row_ == (dim-1);
    end_k_ = num_local_rows_ - 1;
  }
  ~laplacian_mfree() {}

  void multiply(const double* x, double* y) const {
    if (num_local_rows_ == 0) return;

    if ((!is_first_proc) && (nprocs != 1)) {
      MPI_Send((double*)&x[0], 1, MPI_DOUBLE, myrank-1, 0, comm_);
      MPI_Recv(&buf_m, 1, MPI_DOUBLE, myrank-1, 0, comm_, &status_m);
    }

    if ((!is_last_proc) && (nprocs != 1)) {
      MPI_Recv(&buf_p, 1, MPI_DOUBLE, myrank+1, 0, comm_, &status_p);
      MPI_Send((double*)&x[end_k_], 1, MPI_DOUBLE, myrank+1, 0, comm_);
    }

    if (is_first_proc) {
      if (num_local_rows_ != 1) {
        y[0] = x[0] - x[1];
        if (nprocs != 1) y[end_k_] = - x[end_k_ - 1] + 2 * x[end_k_] - buf_p;
      }
      else {
        y[0] = x[0] - buf_p;
      }
    }

    if (is_last_proc) {
      if (num_local_rows_ != 1) {
        if (nprocs != 1) y[0] = - buf_m + 2 * x[0] - x[1];
        y[end_k_] = 2 * x[end_k_] - x[end_k_ - 1];
      }
      else {
        y[end_k_] = 2 * x[end_k_] - buf_m;
      }
    }
    if (!(is_first_proc || is_last_proc)) { // neither first or last process
      if (num_local_rows_ != 1) {
        y[0] = - buf_m + 2 * x[0] - x[1];
        y[end_k_] = - x[end_k_ - 1] + 2 * x[end_k_] - buf_p;
      }
      else {
        y[0] = - buf_m + 2 * x[0] - buf_p;
      }
    }
    // from 1 to end-1
    for (int k=1; k<end_k_; ++k) {
      y[k] = - x[k-1] + 2 * x[k] - x[k+1];
    }
  }
  int get_dim() const { return dim_; }
  int get_local_offset() const { return local_offset_; }
  int get_num_local_rows() const { return num_local_rows_; }
  int get_start_row() const { return start_row_; }
  int get_end_row() const { return end_row_; }

private:
  MPI_Comm comm_;
  int nprocs, myrank;
  int dim_, local_offset_, num_local_rows_;
  int start_row_, end_row_;
  int start_k_, end_k_;
  bool is_first_proc, is_last_proc;
  mutable double buf_m, buf_p;
  mutable MPI_Status status_m, status_p;
};

} // namespace rokko

#endif // ROKKO_UTILITY_LAPLACIAN_MFREE_HPP
