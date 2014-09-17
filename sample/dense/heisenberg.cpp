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

#include <iostream>
#include <boost/foreach.hpp>
#include <rokko/rokko.hpp>
#include <rokko/utility/heisenberg_hamiltonian.hpp>

typedef rokko::matrix_col_major matrix_major;

int main(int argc, char *argv[]) {
  std::string solver_name(rokko::serial_dense_solver::default_solver());
  if (argc >= 2) solver_name = argv[1];

  std::cout.precision(5);
  int L = 8;
  int dim = 1 << L;
  std::vector<std::pair<int, int> > lattice;
  for (int i = 0; i < L-1; ++i) {
    lattice.push_back(std::make_pair(i, i+1));
  }
  lattice.push_back(std::make_pair(L-1, 0));

  rokko::serial_dense_solver solver(solver_name);
  solver.initialize(argc, argv);
  std::cout << "Eigenvalue decomposition of antiferromagnetic Heisenberg chain" << std::endl
            << "solver = " << solver_name << std::endl
            << "L = " << L << std::endl
            << "dimension = " << dim << std::endl;

  rokko::localized_matrix<matrix_major> mat(dim, dim);
  rokko::heisenberg_hamiltonian::generate(L, lattice, mat);

  rokko::localized_vector eigval(dim);
  rokko::localized_matrix<matrix_major> eigvec(dim, dim);
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
