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


class heisenberg_op : public rokko::distributed_mfree {
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
  int nev = 1;
  int blockSize = 5;
  int maxIters = 500;
  double tol = 1.0e-8;

  int L = 8;
  int dim = 1 << L;
  std::vector<std::pair<int, int> > lattice;
  for (int i = 0; i < L; ++i) {
    lattice.push_back(std::make_pair(i, (i+1) % L));
  }

  rokko::parallel_sparse_solver solver("anasazi");

  heisenberg_op  mat(L, lattice);

  if (myrank == root)
    std::cout << "Eigenvalue decomposition of antiferromagnetic Heisenberg chain" << std::endl
              << "solver = LOBPCG" << std::endl
              << "L = " << L << std::endl
              << "dimension = " << mat.get_dim() << std::endl;

  solver.diagonalize(&mat, nev, blockSize, maxIters, tol);

  if (myrank == root) {
    std::cout << "number of converged eigenpairs=" << solver.num_conv() << std::endl;
    std::cout << "smallest eigenvalues:";
    for (int i = 0; i < solver.num_conv(); ++i)
      std::cout << ' ' << solver.eigenvalue(i);
    std::cout << std::endl;
  }

  MPI_Finalize();
}
