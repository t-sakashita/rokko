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
#include <rokko/utility/laplacian_matrix.hpp>
#include <iostream>


int main(int argc, char *argv[]) {
  const std::string library = (argc >= 2) ? argv[1] : rokko::serial_dense_ev::default_solver();
  const unsigned int dim = (argc >= 3) ? std::stoi(argv[2]) : 10;

  std::cout.precision(5);

  rokko::serial_dense_ev solver(library);
  solver.initialize(argc, argv);
  std::cout << "Eigenvalue decomposition of Laplacian matrix" << std::endl
            << "library = " << library << std::endl
            << "dimension = " << dim << std::endl;

  Eigen::MatrixXd mat(dim, dim);
  rokko::laplacian_matrix::generate(mat);
  std::cout << "Laplacian matrix:\n" << mat << std::endl;

  Eigen::VectorXd eigval(dim);
  Eigen::MatrixXd eigvec(dim, dim);
  solver.diagonalize(mat, eigval, eigvec);
  rokko::laplacian_matrix::generate(mat);

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
