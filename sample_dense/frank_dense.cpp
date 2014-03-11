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
#include <boost/lexical_cast.hpp>

#include <rokko/serial_solver.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/frank_matrix.hpp>
#include <rokko/utility/sort_eigenpairs.hpp>

typedef rokko::matrix_col_major matrix_major;

int main(int argc, char *argv[]) {
  if (argc <= 2) {
    std::cerr << "error: " << argv[0] << " solver_name matrix_size" << std::endl;
    exit(34);
  }

  std::cout.precision(5);
  std::string solver_name(argv[1]);
  unsigned int dim = boost::lexical_cast<unsigned int>(argv[2]);

  rokko::serial_solver solver(solver_name);
  solver.initialize(argc, argv);
  std::cout << "solver = " << solver_name << std::endl
            << "dimension = " << dim << std::endl;

  rokko::localized_matrix<matrix_major> mat(dim, dim);
  rokko::frank_matrix::generate(mat);
  std::cout << "Frank matrix:\n" << mat << std::endl;

  rokko::localized_vector w(dim);
  rokko::localized_matrix<matrix_major> eigvec(dim, dim);

  try {
    solver.diagonalize(mat, w, eigvec);
  }
  catch (const char *e) {
    std::cout << "Exception : " << e << std::endl;
    exit(22);
  }
  rokko::frank_matrix::generate(mat);

  rokko::localized_vector eigval_sorted(dim);
  rokko::localized_matrix<matrix_major> eigvec_sorted(dim, dim);
  rokko::sort_eigenpairs(w, eigvec, eigval_sorted, eigvec_sorted);
  std::cout << "eigenvalues:\n" << eigval_sorted.transpose() << std::endl
            << "eigvectors:\n" << eigvec_sorted << std::endl;
  std::cout << "orthogonality of eigenvectors:" << std::endl
            << eigvec_sorted.transpose() * eigvec_sorted << std::endl;
  std::cout << "residual of the largest eigenvalue/vector (A x - lambda x):" << std::endl
            << (mat * eigvec_sorted.col(0) - eigval_sorted(0) * eigvec_sorted.col(0)).transpose()
            << std::endl;

  solver.finalize();
  return 0;
}
