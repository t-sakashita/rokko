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
#include <rokko/distributed_crs_matrix.hpp>

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

  rokko::anasazi::solver_anasazi solver;
  if (myrank == root)
    std::cout << "Eigenvalue decomposition of antiferromagnetic Heisenberg chain" << std::endl
              << "solver = LOBPCG" << std::endl
              << "L = " << L << std::endl
              << "dimension = " << dim << std::endl;

  rokko::mapping_1d map(dim, g);
  rokko::distributed_crs_matrix mat("anasazi", map);
  std::vector<double> values;
  std::vector<int> cols;
  for (int k = 0; k < map.num_rows(); ++k) {
    int row = map.rows()[k];
    cols.clear();
    values.clear();
    double diag = 0;
    for (int l = 0;  l < lattice.size(); ++l) {
      int i = lattice[l].first;
      int j = lattice[l].second;
      int m1 = 1 << i;
      int m2 = 1 << j;
      int m3 = m1 + m2;
      if (((row & m3) == m1) || ((row & m3) == m2)) {
        cols.push_back(row^m3);
        values.push_back(0.5);
        diag += -0.25;
      } else {
        diag += 0.25;
      }
    }
    cols.push_back(row);
    values.push_back(diag);
    mat.insert(row, cols, values);
  }
  mat.complete();

  rokko::distributed_multivector_anasazi ivec(map, blockSize);
  ivec.init_random();
  solver.diagonalize(mat, ivec, nev, blockSize, maxIters, tol);

  if (myrank == root) {
    std::cout << "smallest eigenvalues:";
    for (int i = 0; i < solver.eigenvalues().size(); ++i)
      std::cout << ' ' << solver.eigenvalues()[i].realpart;
    std::cout << std::endl;
  }

  MPI_Finalize();
}
