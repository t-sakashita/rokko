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
#include <boost/lexical_cast.hpp>
#include <iostream>

using matrix_major = rokko::matrix_col_major;

int main(int argc, char *argv[]) {
  std::string solver_name(rokko::serial_dense_ev::default_solver());
  int L = 8;
  if (argc >= 2) solver_name = argv[1];
  if (argc >= 3) L = boost::lexical_cast<unsigned int>(argv[2]);

  std::cout.precision(5);

  std::vector<std::pair<int, int>> lattice;
  for (int i = 0; i < L-1; ++i) { lattice.push_back(std::make_pair(i, i+1)); }
  lattice.push_back(std::make_pair(L-1, 0));
  int dim = 1 << L;

  rokko::serial_dense_ev solver(solver_name);
  solver.initialize(argc, argv);
  std::cout << "Eigenvalue decomposition of antiferromagnetic Heisenberg chain" << std::endl
            << "solver = " << solver_name << std::endl
            << "L = " << L << std::endl
            << "dimension = " << dim << std::endl;

  Eigen::MatrixXd mat(dim, dim);
  rokko::heisenberg_hamiltonian::generate(L, lattice, mat);

  Eigen::VectorXd eigval(dim);
  Eigen::MatrixXd eigvec(dim, dim);
  try {
    solver.diagonalize(mat, eigval, eigvec);
  }
  catch (const char *e) {
    std::cout << "Exception : " << e << std::endl;
    exit(22);
  }
  rokko::heisenberg_hamiltonian::generate(L, lattice, mat);

  std::cout << "smallest eigenvalues:";
  for (int i = 0; i < std::min(dim, 10); ++i) std::cout << ' ' << eigval(i);
  std::cout << std::endl;
  std::cout << "residual of the smallest eigenvalue/vector: |x A x - lambda| = "
            << std::abs(eigvec.col(0).transpose() * mat * eigvec.col(0) - eigval(0))
            << std::endl;

  solver.finalize();
}
