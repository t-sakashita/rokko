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
#include <rokko/utility/solver_name.hpp>
#include <rokko/utility/frank_matrix.hpp>
#include <iostream>


int main(int argc, char *argv[]) {
  const std::string library_routine = (argc >= 2) ? argv[1] : rokko::serial_dense_ev::default_solver();
  const unsigned int dim = (argc >= 3) ? std::stoi(argv[2]) : 10;
  const auto [library, routine] = rokko::split_solver_name(library_routine);

  std::cout.precision(5);
  std::cout << "Eigenvalue decomposition of Frank matrix" << std::endl
            << "library:routine = " << library_routine << std::endl
	    << "library = " << library << std::endl
	    << "routine = " << routine << std::endl
	    << "dimension = " << dim << std::endl;

  rokko::serial_dense_ev solver(library);
  solver.initialize(argc, argv);

  Eigen::MatrixXd mat(dim, dim);
  rokko::frank_matrix::generate(mat);
  std::cout << "Frank matrix:\n" << mat << std::endl;

  Eigen::VectorXd eigval(dim);
  Eigen::MatrixXd eigvec(dim, dim);
  rokko::parameters params;
  params.set("routine", routine);
  params.set("upper_value", 1.2);
  params.set("lower_value", 0.1);
  //params.set("upper_index", 5);
  //params.set("lower_index", 3);
  params.set("uplow", 'L');
  //params.set("uplow", 'lower');
  params.set("verbose", true);
  solver.diagonalize(mat, eigval, eigvec, params);
  //solver.diagonalize(mat, eigval, params);
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
