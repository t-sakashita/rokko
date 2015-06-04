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
    if (num_local_rows_ == 0) return;
    
    if ((!is_first_proc) && (nprocs != 1)) {
      //std::cout << "recv myrank=" << myrank << std::endl;
      MPI_Send(&x[0], 1, MPI_DOUBLE, myrank-1, 0, comm_);
      MPI_Recv(&buf_m, 1, MPI_DOUBLE, myrank-1, 0, comm_, &status_m);
      //std::cout << "buffff=" << buf << std::endl;
    }

    if ((!is_last_proc) && (nprocs != 1)) {
      //std::cout << "send myrank=" << myrank << std::endl;
      MPI_Recv(&buf_p, 1, MPI_DOUBLE, myrank+1, 0, comm_, &status_p);
      MPI_Send(&x[end_k_], 1, MPI_DOUBLE, myrank+1, 0, comm_);
      //std::cout << "buffff=" << buf2 << std::endl;
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


int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;

  rokko::grid_1d g(comm);
  int myrank = g.get_myrank();
  int root = 0;

  std::cout.precision(5);

  rokko::parallel_sparse_solver solver("slepc");
  //  rokko::parallel_sparse_solver solver("anasazi");
  int dim = 20;
  laplacian_op  mat(dim);

  rokko::parameters params;
  params.set("Which", "LM");
  params.set("Block Size", 5);
  params.set("Maximum Iterations", 500);
  params.set("Convergence Tolerance", 1.0e-8);
  params.set("num_eigenvalues", 10);
  params.set("routine", "lanczos");
  //params.set("routine", "SimpleLOBPCG");
  //params.set("routine", "BlockDavidson");
  
  if (myrank == root)
    std::cout << "Eigenvalue decomposition of Laplacian" << std::endl
              << "solver = " << params.get_string("routine") << std::endl
              << "dimension = " << mat.get_dim() << std::endl;

  solver.diagonalize(mat, params);

  if (myrank == root) {
    std::cout << "number of converged eigenpairs=" << solver.num_conv() << std::endl;
    std::cout << "largest eigenvalues:";
    for (int i = 0; i < solver.num_conv(); ++i)
      std::cout << ' ' << solver.eigenvalue(i);
    std::cout << std::endl;
    std::cout << "theoretical eigenvalues:" << std::endl;
    for(int k=dim-1; k>=0; --k) {
      std::cout << rokko::laplacian_matrix::eigenvalue(dim, k) << " ";
    }
    std::cout << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  std::vector<double> eigvec;

  for (int i = 0; i < solver.num_conv(); ++i) {
    solver.eigenvector(i, eigvec);

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
