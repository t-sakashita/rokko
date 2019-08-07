/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <rokko/utility/xyz_lattice.hpp>
#include <rokko/utility/xyz_hamiltonian.hpp>
#include <tuple>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

typedef rokko::matrix_col_major matrix_major;

int main(int argc, char *argv[]) {
  std::string solver_name(rokko::serial_dense_ev::default_solver());
  std::string lattice_file("xyz.dat");
  if (argc >= 2) solver_name = argv[1];
  if (argc >= 3) lattice_file = argv[2];

  std::cout.precision(5);

  int num_sites;
  std::vector<std::pair<int, int> > lattice;
  std::vector<std::tuple<double, double, double> > coupling;
  rokko::read_lattice_file(lattice_file, num_sites, lattice, coupling);
  int dim = 1 << num_sites;

  rokko::serial_dense_ev solver(solver_name);
  solver.initialize(argc, argv);
  std::cout << "Eigenvalue decomposition of XYZ model" << std::endl
            << "solver = " << solver_name << std::endl
            << "lattice file = " << lattice_file << std::endl
            << "number of sites = " << num_sites << std::endl
            << "number of bonds = " << lattice.size() << std::endl
            << "matrix dimension = " << dim << std::endl;

  rokko::localized_matrix<double, matrix_major> mat(dim, dim);
  rokko::xyz_hamiltonian::generate(num_sites, lattice, coupling, mat);

  rokko::localized_vector<double> eigval(dim);
  rokko::localized_matrix<double, matrix_major> eigvec(dim, dim);
  try {
    solver.diagonalize(mat, eigval, eigvec);
  }
  catch (const char *e) {
    std::cout << "Exception : " << e << std::endl;
    exit(22);
  }
  rokko::xyz_hamiltonian::generate(num_sites, lattice, coupling, mat);

  std::cout << "smallest eigenvalues:";
  for (int i = 0; i < std::min(dim, 10); ++i) std::cout << ' ' << eigval(i);
  std::cout << std::endl;
  std::cout << "residual of the smallest eigenvalue/vector: |x A x - lambda| = "
            << std::abs(eigvec.col(0).transpose() * mat * eigvec.col(0) - eigval(0))
            << std::endl;

  solver.finalize();
}
