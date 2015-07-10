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
#include <boost/algorithm/string.hpp>
#include <iostream>

typedef rokko::matrix_col_major matrix_major;

int main(int argc, char *argv[]) {
  std::string library_routine(rokko::serial_dense_solver::default_solver());
  std::string library, routine;
  unsigned int dim = 4;
  if (argc >= 2) library_routine = argv[1];
  if (argc >= 3) dim = boost::lexical_cast<unsigned int>(argv[2]);
  rokko::split_solver_name(library_routine, library, routine);

  std::cout.precision(5);
  std::cout << "Eigenvalue decomposition of Frank matrix" << std::endl
            << "library:routine = " << library_routine << std::endl
	    << "library = " << library << std::endl
	    << "routine = " << routine << std::endl
	    << "dimension = " << dim << std::endl;

  rokko::serial_dense_solver solver(library);
  solver.initialize(argc, argv);

  rokko::localized_matrix<double, matrix_major> mata(dim, dim), matb(dim, dim);
  mata << 0.24, 0.39, 0.42, -0.16,
          0.39, -0.11, 0.79, 0.63,
          0.42, 0.79, -0.25, 0.48,
         -0.16, 0.63, 0.48, -0.03;

  matb << 4.16, -3.12, 0.56, -0.10,
         -3.12, 5.03, -0.83, 1.09,
          0.56, -0.83, 0.76, 0.34,
         -0.10, 1.09, 0.34, 1.18;
    
  //std::cout << "Frank matrix:\n" << mat << std::endl;

  rokko::localized_vector<double> eigval(dim);
  rokko::localized_matrix<double, matrix_major> eigvec(dim, dim);
  rokko::parameters params;
  //params.set("upper_value", 1.2);
  //params.set("lower_value", 0.1);
  //params.set("upper_index", 5);
  //params.set("lower_index", 3);
  params.set("uplow", 'L');
  //params.set("uplow", 'lower');
  params.set("verbose", true);
  try {
    //solver.diagonalize(routine, mata, matb, eigval, eigvec, params);
    solver.diagonalize(routine, mata, matb, eigval, params);
  }
  catch (const char *e) {
    std::cout << "Exception : " << e << std::endl;
    exit(22);
  }
  rokko::frank_matrix::generate(mata);

  bool sorted = true;
  for (unsigned int i = 1; i < dim; ++i) sorted &= (eigval(i-1) <= eigval(i));
  if (!sorted) std::cout << "Warning: eigenvalues are not sorted in ascending order!\n";

  std::cout << "eigenvalues:\n" << eigval.transpose() << std::endl
            << "eigvectors:\n" << eigvec << std::endl;
  std::cout << "orthogonality of eigenvectors:" << std::endl
            << eigvec.transpose() * eigvec << std::endl;
  std::cout << "residual of the smallest eigenvalue/vector (A x - lambda x):" << std::endl
            << (mata * eigvec.col(0) - eigval(0) * eigvec.col(0)).transpose() << std::endl;

  solver.finalize();
}