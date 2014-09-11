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
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>

#include <rokko/solver.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/collective.hpp>
#include <rokko/utility/xyz_hamiltonian_mpi.hpp>

typedef rokko::matrix_col_major matrix_major;

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;

  rokko::grid g(comm);
  int myrank = g.get_myrank();
  int root = 0;

  if (argc <= 2) {
    if (myrank == root) {
      std::cerr << "error: " << argv[0] << " solver_name lattice_file" << std::endl
                << "available solvers:";
      BOOST_FOREACH(std::string name, rokko::parallel_dense_solver::solvers())
        std::cerr << ' ' << name;
      std::cerr << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 34);
  }

  std::cout.precision(5);
  std::string solver_name(argv[1]);

  std::ifstream ifs(argv[2]);
  if (!ifs) {
    std::cout << "can't open file" << std::endl;
    exit(1);
  }
  int num_sites, num_bonds;
  std::vector<std::pair<int, int> > lattice;
  std::vector<boost::tuple<double, double, double> > coupling;
  ifs >> num_sites >> num_bonds;
  for (int i = 0; i < num_bonds; ++i) {
    int j, k;
    ifs >> j >> k;
    lattice.push_back(std::make_pair(j, k));
  }
  
  for (int i = 0; i < num_bonds; ++i) {
    double jx, jy, jz;
    ifs >> jx >> jy >> jz;
    coupling.push_back(boost::make_tuple(jx, jy, jz));
  }
  int dim = 1 << num_sites;
  if (myrank == root) {
    std::cout << "Eigenvalue decomposition of XYZ model" << std::endl
              << "solver = " << solver_name << std::endl
              << "lattice file = " << argv[2] << std::endl
              << "number of sites = " << num_sites << std::endl
              << "number of bonds = " << num_bonds << std::endl
              << "matrix dimension = " << dim << std::endl;
  }

  rokko::parallel_dense_solver solver(solver_name);
  solver.initialize(argc, argv);

  rokko::distributed_matrix<matrix_major> mat(dim, dim, g, solver);
  rokko::xyz_hamiltonian::generate(num_sites, lattice, coupling, mat);
  rokko::localized_matrix<matrix_major> mat_loc(dim, dim);
  rokko::gather(mat, mat_loc, root);

  rokko::localized_vector eigval(dim);
  rokko::distributed_matrix<matrix_major> eigvec(dim, dim, g, solver);
  try {
    solver.diagonalize(mat, eigval, eigvec);
  }
  catch (const char *e) {
    if (myrank == root) std::cout << "Exception : " << e << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 22);
  }

  rokko::localized_matrix<matrix_major> eigvec_loc(dim, dim);
  rokko::gather(eigvec, eigvec_loc, root);
  if (myrank == root) {
    std::cout << "smallest eigenvalues:";
    for (int i = 0; i < std::min(dim, 10); ++i) std::cout << ' ' << eigval(i);
    std::cout << std::endl;
    std::cout << "residual of the smallest eigenvalue/vector: |x A x - lambda| = "
              << std::abs(eigvec_loc.col(0).transpose() * mat_loc * eigvec_loc.col(0) - eigval(0))
              << std::endl;
  }

  solver.finalize();
  MPI_Finalize();
}
