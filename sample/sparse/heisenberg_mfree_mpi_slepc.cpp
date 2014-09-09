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

#include <rokko/parallel_sparse_solver.hpp>
#include <rokko/utility/heisenberg_hamiltonian_mpi.hpp>

#include <rokko/distributed_mfree.hpp>

class heisenberg_op : public rokko::distributed_mfree_slepc {
public:
  heisenberg_op(int L, const std::vector<std::pair<int, int> >& lattice) : L_(L), lattice_(lattice) {
    comm_ = MPI_COMM_WORLD;
    int nproc;
    MPI_Comm_size(comm_, &nproc);
    int n = nproc;
    int p = -1;
    do {
      n /= 2;
      ++p;
    } while (n > 0);
    local_N = 1 << (L-p);
    buffer_.assign(local_N, 0);
    dim_ = 1 << L;
  }

  ~heisenberg_op() {}

  void multiply(const double* x, double* y) const {
    rokko::heisenberg_hamiltonian::multiply(comm_, L_, lattice_, x, y, &(buffer_[0]));
  }

  void diagonal(double* x) const {
    rokko::heisenberg_hamiltonian::fill_diagonal(comm_, L_, lattice_, x);
  }

  int get_dim() const {
    return dim_;
  }

  int get_num_local_rows() const {
    return local_N;
  }

private:
  MPI_Comm comm_;
  mutable std::vector<double> buffer_;
  int L_;
  int local_N;
  std::vector<std::pair<int, int> > lattice_;
  int dim_;
};


int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;

  rokko::grid_1d g(comm);
  int myrank = g.get_myrank();
  int root = 0;

  std::cout.precision(5);
  int nev = 1; //10;
  int blockSize = 1; //5;
  int maxIters = 500;
  double tol = 1.0e-8;

  int L = 3;
  std::vector<std::pair<int, int> > lattice;
  for (int i = 0; i < L; ++i) {
    lattice.push_back(std::make_pair(i, (i+1) % L));
  }

  rokko::parallel_sparse_solver solver("slepc");
  heisenberg_op  mat(L, lattice);

  if (myrank == root)
    std::cout << "Eigenvalue decomposition of antiferromagnetic Heisenberg chain" << std::endl
              << "solver = LOBPCG" << std::endl
              << "L = " << L << std::endl
              << "dimension = " << mat.get_dim() << std::endl;

  solver.diagonalize(&mat, nev, blockSize, maxIters, tol);

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
