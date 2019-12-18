/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
#include <rokko/collective.hpp>
#include <rokko/utility/xyz_hamiltonian_mpi.hpp>
#include <tuple>
#include <fstream>
#include <iostream>

using matrix_major = rokko::matrix_col_major;

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;
  std::string solver_name(rokko::parallel_dense_ev::default_solver());
  std::string lattice_file("xyz.dat");
  if (argc >= 2) solver_name = argv[1];
  if (argc >= 3) lattice_file = argv[2];

  rokko::grid g(comm);
  int myrank = g.get_myrank();

  std::cout.precision(5);

  std::ifstream ifs(lattice_file);
  if (!ifs) {
    throw std::runtime_error("can't open file \"" + lattice_file + "\"");
  }
  int num_sites, num_bonds;
  std::vector<std::pair<int, int>> lattice;
  std::vector<std::tuple<double, double, double>> coupling;
  ifs >> num_sites >> num_bonds;
  int dim = 1 << num_sites;
  for (int i = 0; i < num_bonds; ++i) {
    int j, k;
    ifs >> j >> k;
    lattice.push_back(std::make_pair(j, k));
  }
  for (int i = 0; i < num_bonds; ++i) {
    double jx, jy, jz;
    ifs >> jx >> jy >> jz;
    coupling.push_back(std::make_tuple(jx, jy, jz));
  }

  rokko::parallel_dense_ev solver(solver_name);
  solver.initialize(argc, argv);
  if (myrank == 0)
    std::cout << "Eigenvalue decomposition of XYZ model" << std::endl
              << "num_procs = " << g.get_nprocs() << std::endl
              #ifdef _OPENMP
              << "num_threads per process = " << omp_get_max_threads() << std::endl
              #endif
              << "solver = " << solver_name << std::endl
              << "lattice file = " << lattice_file << std::endl
              << "number of sites = " << num_sites << std::endl
              << "number of bonds = " << num_bonds << std::endl
              << "dimension = " << dim << std::endl;

  rokko::mapping_bc<matrix_major> map = solver.default_mapping(dim, g);
  rokko::distributed_matrix<double, matrix_major> mat(map);
  rokko::xyz_hamiltonian::generate(num_sites, lattice, coupling, mat);
  Eigen::MatrixXd mat_loc(dim, dim);
  rokko::gather(mat, mat_loc, 0);

  Eigen::VectorXd eigval(dim);
  rokko::distributed_matrix<double, matrix_major> eigvec(map);
  solver.diagonalize(mat, eigval, eigvec);

  Eigen::MatrixXd eigvec_loc(dim, dim);
  rokko::gather(eigvec, eigvec_loc, 0);
  if (myrank == 0) {
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
