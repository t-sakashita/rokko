/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
#include <rokko/utility/frank_matrix.hpp>
#include <boost/foreach.hpp>
#define BOOST_TEST_MODULE test_solver
#ifndef BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

template<typename MATRIX_MAJOR>
void test(int dim, std::string const& name) {
  rokko::serial_dense_ev solver(name);
  solver.initialize(boost::unit_test::framework::master_test_suite().argc,
                    boost::unit_test::framework::master_test_suite().argv);
  rokko::localized_matrix<double, MATRIX_MAJOR> mat(dim, dim);
  rokko::frank_matrix::generate(mat);
  rokko::localized_vector<double> eigval(dim);
  rokko::localized_matrix<double, MATRIX_MAJOR> eigvec(dim, dim);

  solver.diagonalize(mat, eigval, eigvec);
  
  double sum = 0;
  for(int i = 0; i < dim; ++i) sum += eigval[i];
  BOOST_CHECK_CLOSE(sum, dim * (dim+1) * 0.5, 10e-5);
  
  rokko::frank_matrix::generate(mat);
  for (int i = 0; i < dim; ++i) {
    double w = eigvec.col(i).transpose() * mat * eigvec.col(i);
    BOOST_CHECK_CLOSE(w, eigval[i], 10e-5);
  }

  solver.finalize();
}

BOOST_AUTO_TEST_CASE(test_solver) {
  const int dim = 100;
  std::cout << "dimension = " << dim << std::endl;

  std::vector<std::string> names;
  int argc = boost::unit_test::framework::master_test_suite().argc;
  if (argc == 1) {
    names = rokko::serial_dense_ev::solvers();
  } else {
    for (int num=1; num < argc; ++num) {
      names.push_back(boost::unit_test::framework::master_test_suite().argv[num]);
    }
  }

  BOOST_FOREACH(std::string name, names) {
    std::cout << "solver = " << name << std::endl;
    std::cout << "  test for row major" << std::endl;
    test<rokko::matrix_row_major>(dim, name);
    std::cout << "  test for column major" << std::endl;
    test<rokko::matrix_col_major>(dim, name);
  }
}
