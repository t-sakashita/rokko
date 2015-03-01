/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014 by Synge Todo <wistaria@comp-phys.org>
*                       Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <mpi.h>
#include <iostream>
#include <vector>
#include <rokko/rokko.hpp>

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
    MPI_Status *status;
    if (!is_first_proc) {
      //std::cout << "recv myrank=" << myrank << std::endl;
      MPI_Send(&x[0], 1, MPI_DOUBLE, myrank-1, 0, comm_);
      MPI_Recv(&buf, 1, MPI_DOUBLE, myrank-1, 0, comm_, status);
      y[0] = - buf + 2 * x[0] - x[1];
    }
    else { // for the first process 0
      y[0] = x[0] - x[1];
    }
    if (!is_last_proc) {
      //std::cout << "send myrank=" << myrank << std::endl;
      MPI_Send(&x[end_k_], 1, MPI_DOUBLE, myrank+1, 0, comm_);
      MPI_Recv(&buf, 1, MPI_DOUBLE, myrank+1, 0, comm_, status);
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
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;

  rokko::grid_1d g(comm);
  int myrank = g.get_myrank();
  int root = 0;

  std::cout.precision(5);
  int nev = 5;
  int blockSize = 2;
  int maxIters = 500;
  double tol = 1.0e-2;  //6;


  rokko::parallel_sparse_solver solver("anasazi");
  laplacian_op  mat(20);

  if (myrank == root)
    std::cout << "Eigenvalue decomposition of Laplacian" << std::endl
              << "solver = LOBPCG" << std::endl
              << "dimension = " << mat.get_dim() << std::endl;

  solver.diagonalize(mat, nev, blockSize, maxIters, tol);

  if (myrank == root) {
    std::cout << "number of converged eigenpairs=" << solver.num_conv() << std::endl;
    std::cout << "smallest eigenvalues:";
    for (int i = 0; i < solver.num_conv(); ++i)
      std::cout << ' ' << solver.eigenvalue(i);
    std::cout << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  std::vector<double> eigvec;

  for (int i = 0; i < solver.num_conv(); ++i) {
    solver.eigenvector(0, eigvec);

    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank == root) {
      std::cout << "corresponding eigenvectors:";
      for (int j=0; j<eigvec.size(); ++j)
	std::cout << ' ' << eigvec[j];
      std::cout << std::endl;
    }
  }

  MPI_Finalize();
}
