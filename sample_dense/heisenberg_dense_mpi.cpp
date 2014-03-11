/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>

#include <rokko/solver.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/collective.hpp>
#include <rokko/utility/heisenberg_hamiltonian_mpi.hpp>
#include <rokko/utility/sort_eigenpairs.hpp>

typedef rokko::matrix_col_major matrix_major;

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;

  rokko::grid g(comm, rokko::grid_col_major);
  int myrank = g.get_myrank();
  int root = 0;

  if (argc <= 1) {
    if (myrank == root)
      std::cerr << "error: " << argv[0] << " solver_name" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 34);
  }

  std::cout.precision(5);
  std::string solver_name(argv[1]);
  int L = 8;
  int dim = 1 << L;
  std::vector<std::pair<int, int> > lattice;
  for (int i = 0; i < L-1; ++i) {
    lattice.push_back(std::make_pair(i, i+1));
  }
  lattice.push_back(std::make_pair(L-1, 0));

  rokko::parallel_dense_solver solver(solver_name);
  solver.initialize(argc, argv);
  if (myrank == root)
    std::cout << "solver = " << solver_name << std::endl
              << "L = " << L << std::endl
              << "dimension = " << dim << std::endl;

  rokko::distributed_matrix<matrix_major> mat(dim, dim, g, solver);
  rokko::heisenberg_hamiltonian::generate(L, lattice, mat);
  rokko::localized_matrix<matrix_major> mat_loc(dim, dim);
  rokko::gather(mat, mat_loc, root);

  rokko::localized_vector w(dim);
  rokko::distributed_matrix<matrix_major> eigvec(dim, dim, g, solver);

  try {
    solver.diagonalize(mat, w, eigvec);
  }
  catch (const char *e) {
    if (myrank == root) std::cout << "Exception : " << e << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 22);
  }

  rokko::localized_matrix<matrix_major> eigvec_loc;
  rokko::gather(eigvec, eigvec_loc, root);
  if (myrank == root) {
    rokko::localized_vector eigval_sorted(dim);
    rokko::localized_matrix<matrix_major> eigvec_sorted(dim, dim);
    rokko::sort_eigenpairs(w, eigvec_loc, eigval_sorted, eigvec_sorted);
    std::cout << "eigenvalues:\n" << eigval_sorted.transpose() << std::endl;
    std::cout << "residual of the smallest eigenvalue/vector (A x - lambda x):" << std::endl
              << (mat_loc * eigvec_sorted.col(dim-1) - eigval_sorted(dim-1) * eigvec_sorted.col(dim-1)).transpose()
              << std::endl;
  }

  solver.finalize();
  MPI_Finalize();
  return 0;
}
