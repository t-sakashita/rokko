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
#include <rokko/utility/frank_matrix.hpp>
#include <iostream>


int main(int argc, char *argv[]) {
  std::string solver_name(rokko::serial_dense_ev::default_solver());
  unsigned int dim = 10;
  if (argc >= 2) solver_name = argv[1];
  if (argc >= 3) dim = std::stoi(argv[2]);

  std::cout.precision(5);

  rokko::serial_dense_ev solver(solver_name);
  solver.initialize(argc, argv);
  std::cout << "Eigenvalue decomposition of Frank matrix" << std::endl
            << "solver = " << solver_name << std::endl
            << "dimension = " << dim << std::endl;

  Eigen::MatrixXd mat(dim, dim);
  rokko::frank_matrix::generate(mat);
  std::cout << "Frank matrix:\n" << mat << std::endl;

  Eigen::VectorXd eigval(dim);
  Eigen::MatrixXd eigvec(dim, dim);
  solver.diagonalize(mat, eigval, eigvec);
  rokko::frank_matrix::generate(mat);

  bool sorted = true;
  for (unsigned int i = 1; i < dim; ++i) sorted &= (eigval(i-1) <= eigval(i));
  if (!sorted) std::cout << "Warning: eigenvalues are not sorted in ascending order!\n";

  std::cout << "eigenvalues:\n" << eigval.transpose() << std::endl
            << "eigvectors:\n" << eigvec << std::endl;
  std::cout << "orthogonality of eigenvectors:" << std::endl
            << eigvec.transpose() * eigvec << std::endl;
  std::cout << "residual of the smallest eigenvalue/vector (A x - lambda x):" << std::endl
            << (mat * eigvec.col(0) - eigval(0) * eigvec.col(0)).transpose() << std::endl;

  solver.finalize();
}
