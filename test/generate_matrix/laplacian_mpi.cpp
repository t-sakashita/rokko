/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <iostream>
#include <vector>

#include <rokko/solver.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_mfree.hpp>

#include <rokko/collective.hpp>

#include <rokko/rokko.hpp>
#include <rokko/utility/laplacian_matrix.hpp>

class laplacian_op : public rokko::distributed_mfree {
public:
  laplacian_op(int dim) : dim_(dim) {
    comm_ = MPI_COMM_WORLD;
    MPI_Comm_size(comm_, &nprocs);
    MPI_Comm_rank(comm_, &myrank);

    int tmp = dim_ / nprocs;
    int rem = dim % nprocs;
    num_local_rows_ = (dim + nprocs - myrank - 1) / nprocs;
    start_row_ = tmp * myrank + std::min(rem, myrank);
    end_row_ = start_row_ + num_local_rows_ - 1;

    if (start_row_ == 0)  is_first_proc = true;
    else is_first_proc = false;
    
    if (end_row_ == (dim-1))  is_last_proc = true;
    else is_last_proc = false;

    end_k_ = num_local_rows_ - 1;
    
    std::cout << "myrank=" << myrank << " start_row=" << start_row_ << " end_row=" << end_row_ << std::endl;
    std::cout << "myrank=" << myrank << " num_local_rows_=" << num_local_rows_ << std::endl;
  }
  ~laplacian_op() {}

  void multiply(const double* x, double* y) const {
    MPI_Status status;
    if (!is_first_proc) {
      //std::cout << "recv myrank=" << myrank << std::endl;
      MPI_Send(x, 1, MPI_DOUBLE, myrank-1, 0, comm_);
      MPI_Recv(&buf, 1, MPI_DOUBLE, myrank-1, 0, comm_, &status);
      y[0] = - buf + 2 * x[0] - x[1];
    }
    else { // for the first process 0
      y[0] = x[0] - x[1];
    }
    if (!is_last_proc) {
      //std::cout << "send myrank=" << myrank << std::endl;
      MPI_Send(&x[end_k_], 1, MPI_DOUBLE, myrank+1, 0, comm_);
      MPI_Recv(&buf, 1, MPI_DOUBLE, myrank+1, 0, comm_, &status);
      y[end_k_] = - x[end_k_ - 1] + 2 * x[end_k_] - buf;      
    }
    else { // for the last process
      y[end_k_] = 2 * x[end_k_] - x[end_k_ - 1];
    }
 
    for (int k=1; k<end_k_; ++k) {  //  from 1 to end-1
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
  mutable double buf;
  bool is_first_proc, is_last_proc;
};

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  rokko::parallel_dense_ev solver;
  solver.initialize(argc, argv);
  rokko::grid g(MPI_COMM_WORLD);


  int myrank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  constexpr int root = 0;

  // creating column vectors which forms a heisenberg hamiltonian.
  int N_seq = 20;
  int N = N_seq / nprocs;
  Eigen::VectorXd recv_buffer(N_seq);
  Eigen::VectorXd buffer(N);
  laplacian_op op(N_seq);
  for (int i=0; i<N_seq; ++i) {
    // sequential version
    Eigen::VectorXd v_seq(N_seq), w_seq(N_seq);
    v_seq.setZero();
    v_seq(i) = 1;
    w_seq.setZero();
    if (myrank == root) {
      rokko::laplacian_matrix::multiply(N_seq, v_seq.data(), w_seq.data());
      std::cout << "sequential version:" << std::endl;
      std::cout << w_seq.transpose() << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // MPI version
    Eigen::VectorXd v(N), w(N);
    MPI_Scatter(v_seq.data(), N, MPI_DOUBLE, v.data(), N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    w.setZero();
    op.multiply(v.data(), w.data());
    for (int proc=0; proc<nprocs; ++proc) {
      if (proc == myrank) {
        std::cout << "myrank=" << myrank << std::endl;
        std::cout << w.transpose() << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(w.data(), N, MPI_DOUBLE, recv_buffer.data(), N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0) {
      std::cout << "i=" << i << std::endl;
      std::cout << "recv=" << recv_buffer.transpose() << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Finalize();
  return 0;
}
