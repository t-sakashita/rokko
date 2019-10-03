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
#include <boost/lexical_cast.hpp>
#include <iostream>

typedef rokko::matrix_col_major matrix_major;

template<int MATRIX_MAJOR, typename VEC>
void default_diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, VEC& eigvals,
			 rokko::parameters const& params) {
  rokko::serial_dense_ev solver(rokko::serial_dense_ev::default_solver());
  int argc;
  char** xargv;
  solver.initialize(argc, xargv);
  solver.diagonalize(mat, eigvals, params);
}

template<int MATRIX_MAJOR, typename VEC>
void default_diagonalize(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, VEC& eigvals,
			 Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			 rokko::parameters const& params) {
  rokko::serial_dense_ev solver(rokko::serial_dense_ev::default_solver());
  int argc;
  char** xargv;
  solver.initialize(argc, xargv);
  solver.diagonalize(mat, eigvals, eigvecs, params);
}

int main(int argc, char *argv[]) {
  std::string library_routine(rokko::serial_dense_ev::default_solver());
  std::string library, routine;
  unsigned int dim = 10;
  if (argc >= 2) library_routine = argv[1];
  if (argc >= 3) dim = boost::lexical_cast<int>(argv[2]);
  rokko::split_solver_name(library_routine, library, routine);

  std::cout.precision(5);
  std::cout << "Eigenvalue decomposition of Frank matrix" << std::endl
            << "library:routine = " << library_routine << std::endl
	    << "library = " << library << std::endl
	    << "routine = " << routine << std::endl
	    << "dimension = " << dim << std::endl;

  rokko::serial_dense_ev solver(library);
  solver.initialize(argc, argv);

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<matrix_major>> mat(dim, dim);
  rokko::frank_matrix::generate(mat);
  std::cout << "Frank matrix:\n" << mat << std::endl;

  Eigen::VectorXd eigval(dim);
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<matrix_major>> eigvec(dim, dim);
  rokko::parameters params;
  params.set("upper_value", 1.2);
  params.set("lower_value", 0.1);
  //params.set("upper_index", 5);
  //params.set("lower_index", 3);
  params.set("uplow", 'L');
  //params.set("uplow", 'lower');
  params.set("verbose", true);
  try {
    default_diagonalize(mat, eigval, eigvec, params);
    //default_diagonalize(mat, eigval, params);
  }
  catch (const char *e) {
    std::cout << "Exception : " << e << std::endl;
    exit(22);
  }
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
