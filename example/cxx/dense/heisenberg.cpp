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
#include <rokko/utility/heisenberg_hamiltonian.hpp>
#include <iostream>

auto create_periodic_1dim_lattice(int L) {
  std::vector<std::pair<int, int>> lattice;
  for (auto i=0; i<L-1; ++i) {
    lattice.emplace_back(std::make_pair(i, i+1));
  }
  lattice.emplace_back(std::make_pair(L-1, 0));

  return lattice;
}

int main(int argc, char *argv[]) {
  const std::string library = (argc >= 2) ? argv[1] : rokko::serial_dense_ev::default_solver();
  const int L = (argc >= 3) ? std::stoi(argv[2]) : 8;

  std::cout.precision(5);

  const auto lattice = create_periodic_1dim_lattice(L);
  const auto dim = 1 << L;

  rokko::serial_dense_ev solver(library);
  solver.initialize(argc, argv);
  std::cout << "Eigenvalue decomposition of antiferromagnetic Heisenberg chain" << std::endl
            << "library = " << library << std::endl
            << "L = " << L << std::endl
            << "dimension = " << dim << std::endl;

  Eigen::MatrixXd mat(dim, dim);
  rokko::heisenberg_hamiltonian::generate(L, lattice, mat);

  Eigen::VectorXd eigval(dim);
  Eigen::MatrixXd eigvec(dim, dim);
  solver.diagonalize(mat, eigval, eigvec);
  rokko::heisenberg_hamiltonian::generate(L, lattice, mat);

  std::cout << "smallest eigenvalues:";
  for (int i = 0; i < std::min(dim, 10); ++i) std::cout << ' ' << eigval(i);
  std::cout << std::endl;
  std::cout << "residual of the smallest eigenvalue/vector: |x A x - lambda| = "
            << std::abs(eigvec.col(0).transpose() * mat * eigvec.col(0) - eigval(0))
            << std::endl;

  solver.finalize();
}
