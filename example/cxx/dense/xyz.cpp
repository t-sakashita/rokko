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
#include <rokko/utility/xyz_hamiltonian.hpp>
#include <tuple>
#include <iostream>
#include <fstream>


int main(int argc, char *argv[]) {
  std::string solver_name(rokko::serial_dense_ev::default_solver());
  std::string lattice_file("xyz.dat");
  if (argc >= 2) solver_name = argv[1];
  if (argc >= 3) lattice_file = argv[2];

  std::cout.precision(5);

  std::ifstream ifs(lattice_file);
  if (!ifs) {
    throw std::runtime_error("can't open file \"" + lattice_file + "\"");
  }
  int num_sites, num_bonds;
  std::vector<std::pair<int, int>> lattice;
  std::vector<std::tuple<double, double, double>> coupling;
  ifs >> num_sites >> num_bonds;
  for (int i = 0; i < num_bonds; ++i) {
    int j, k;
    ifs >> j >> k;
    lattice.emplace_back(std::make_pair(j, k));
  }
  for (int i = 0; i < num_bonds; ++i) {
    double jx, jy, jz;
    ifs >> jx >> jy >> jz;
    coupling.emplace_back(std::make_tuple(jx, jy, jz));
  }
  int dim = 1 << num_sites;

  rokko::serial_dense_ev solver(solver_name);
  solver.initialize(argc, argv);
  std::cout << "Eigenvalue decomposition of XYZ model" << std::endl
            << "solver = " << solver_name << std::endl
            << "lattice file = " << lattice_file << std::endl
            << "number of sites = " << num_sites << std::endl
            << "number of bonds = " << num_bonds << std::endl
            << "matrix dimension = " << dim << std::endl;

  Eigen::MatrixXd mat(dim, dim);
  rokko::xyz_hamiltonian::generate(num_sites, lattice, coupling, mat);

  Eigen::VectorXd eigval(dim);
  Eigen::MatrixXd eigvec(dim, dim);
  solver.diagonalize(mat, eigval, eigvec);
  rokko::xyz_hamiltonian::generate(num_sites, lattice, coupling, mat);

  std::cout << "smallest eigenvalues:";
  for (int i = 0; i < std::min(dim, 10); ++i) std::cout << ' ' << eigval(i);
  std::cout << std::endl;
  std::cout << "residual of the smallest eigenvalue/vector: |x A x - lambda| = "
            << std::abs(eigvec.col(0).transpose() * mat * eigvec.col(0) - eigval(0))
            << std::endl;

  solver.finalize();
}
