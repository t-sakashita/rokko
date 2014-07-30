/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <mpi.h>
#include <iostream>

#include <rokko/anasazi/core.hpp>
#include <rokko/mapping_1d.hpp>
#include <rokko/distributed_mfree.hpp>

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziSimpleLOBPCGSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_CommandLineProcessor.hpp"

#include <rokko/utility/heisenberg_hamiltonian_mpi.hpp>
#include <vector>


class heisenberg_op : public rokko::distributed_mfree {
public:
  heisenberg_op(rokko::mapping_1d const& map, int L, const std::vector<std::pair<int, int> >& lattice) : L_(L), lattice_(lattice) {
    map_ = get_mapping_1d();
    comm_ = MPI_COMM_WORLD;
    int nproc;
    MPI_Comm_size(comm_, &nproc);
    int n = nproc;
    int p = -1;
    do {
      n /= 2;
      ++p;
    } while (n > 0);
    int local_N = 1 << (L-p);
    buffer_.assign(local_N, 0);
  }

  ~heisenberg_op() {}

  void multiply(const double* x, double* y) const {
    rokko::heisenberg_hamiltonian::multiply(comm_, L_, lattice_, x, y, &(buffer_[0]));
  }
private:
  MPI_Comm comm_;
  mutable std::vector<double> buffer_;
  int L_;
  std::vector<std::pair<int, int> > lattice_;
};


int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;

  rokko::grid_1d g(comm);
  int myrank = g.get_myrank();
  int root = 0;

  std::cout.precision(5);
  int nev = 10;
  int blockSize = 5;
  int maxIters = 500;
  double tol = 1.0e-8;

  int L = 8;
  int dim = 1 << L;
  std::vector<std::pair<int, int> > lattice;
  for (int i = 0; i < L; ++i) {
    lattice.push_back(std::make_pair(i, (i+1) % L));
  }

  rokko::solver_anasazi solver;
  if (myrank == root)
    std::cout << "Eigenvalue decomposition of antiferromagnetic Heisenberg chain" << std::endl
              << "solver = LOBPCG" << std::endl
              << "L = " << L << std::endl
              << "dimension = " << dim << std::endl;

  rokko::mapping_1d map(dim, g);
  heisenberg_op  mat(map, L, lattice);

  rokko::distributed_multivector_anasazi ivec(map, blockSize);
  ivec.init_random();
  solver.diagonalize(&mat, ivec, nev, blockSize, maxIters, tol);

  if (myrank == root) {
    std::cout << "smallest eigenvalues:";
    for (int i = 0; i < solver.eigenvalues().size(); ++i)
      std::cout << ' ' << solver.eigenvalues()[i].realpart;
    std::cout << std::endl;
  }

  MPI_Finalize();
}
