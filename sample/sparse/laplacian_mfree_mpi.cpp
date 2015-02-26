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
    int nprocs, myrank;
    MPI_Comm_size(comm_, &nprocs);
    MPI_Comm_rank(comm_, &myrank);
    //int n = size;

    num_local_rows_ = dim_ / nprocs;
    int rem = dim % nprocs;
    int rest_size = rem / (myrank+1);  // 0 or 1
    num_local_rows_ += rest_size;
    start_row_ = dim / nprocs * myrank + std::min(rem, myrank);
    end_row_ = start_row_ + num_local_rows_;
    //local_offset_ = num_local_rows_ * myrank;
    start_k_ = start_row_;
    if (start_k_ == 0) {
      ++start_k_;
    }
    end_k_ = end_row_;
    if (end_row_ == dim) {
      --end_k_;
    }
    std::cout << "myrank=" << myrank << " start_row=" << start_row_ << std::endl;
    //std::cout << "num_local_rows=" << num_local_rows_ << std::endl;
    //std::cout << "rest_size=" << rest_size << std::endl;
    //std::cout << "end_row=" << end_row_ << std::endl;
    int tmp = dim_ / nprocs;
    int cal = (tmp + nprocs - myrank - 1) / nprocs;
    std::cout << "myrank=" << myrank << " cal=" << cal << std::endl;
  }
  ~laplacian_op() {}

  void multiply(const double* x, double* y) const {
    y[0] += x[0] - x[1];
    for (int k=start_k_; k<end_k_; ++k) {
      y[k] += - x[k-1] + 2 * x[k] -x[k+1];
      //y[k] -= x[k-1];
    }
    //y[end_row_ - 1] += 2 * x[end_row_ - 1] - x[end_row_ - 2];
  }
  int get_dim() const { return dim_; }
  int get_local_offset() const { return local_offset_; }
  int get_num_local_rows() const { return num_local_rows_; }
  int get_start_row() const { return start_row_; }
  int get_end_row() const { return end_row_; }

private:
  MPI_Comm comm_;
  int dim_, local_offset_, num_local_rows_;
  int start_row_, end_row_;
  int start_k_, end_k_;
};

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;

  rokko::grid_1d g(comm);
  int myrank = g.get_myrank();
  int root = 0;

  std::cout.precision(5);
  int nev = 1;
  int blockSize = 1;
  int maxIters = 500;
  double tol = 1.0e-8;


  rokko::parallel_sparse_solver solver("anasazi");
  laplacian_op  mat(10);

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

  //for (int i = 0; i < solver.num_conv(); ++i) {
  solver.eigenvector(0, eigvec);

  MPI_Barrier(MPI_COMM_WORLD);

  //if (myrank == root) {
  if (myrank == 1) {
  std::cout << "corresponding eigenvectors:";
    for (int j=0; j<eigvec.size(); ++j)
      std::cout << ' ' << eigvec[j];
    std::cout << std::endl;
  }
  //}

  MPI_Finalize();
}
