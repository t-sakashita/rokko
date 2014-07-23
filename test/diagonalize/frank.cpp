/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2013 by Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>,
*                            Tatsuya Sakashita <t-sakashit@issp.u-tokyo.ac.jp>,
*               2014-2014    Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <boost/foreach.hpp>

#include <rokko/solver.hpp>
#include <rokko/utility/frank_matrix.hpp>

#define BOOST_TEST_MODULE test_solver
#ifndef BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

BOOST_AUTO_TEST_CASE(test_solver) {
  const int dim = 100;

  std::vector<std::string> names;
  int argc = boost::unit_test::framework::master_test_suite().argc;
  if (argc == 1) {
    names = rokko::serial_dense_solver::solvers();
  } else {
    for (int num=1; num < argc; ++num) {
      names.push_back(boost::unit_test::framework::master_test_suite().argv[num]);
    }
  }

  BOOST_FOREACH(std::string name, names) {
    std::cout << "solver=" << name << std::endl;
    rokko::serial_dense_solver solver(name);
    solver.initialize(boost::unit_test::framework::master_test_suite().argc,
                      boost::unit_test::framework::master_test_suite().argv);
    rokko::localized_matrix<rokko::matrix_col_major> mat(dim, dim);
    rokko::frank_matrix::generate(mat);
    rokko::localized_vector w(dim);
    rokko::localized_matrix<rokko::matrix_col_major> Z(dim, dim);

    solver.diagonalize(mat, w, Z);
    
    double sum = 0;
    for(int i=0; i<dim; ++i) {
      sum += w[i];
    }

    BOOST_CHECK_CLOSE(sum, dim * (dim+1) * 0.5, 10e-5);
    
    solver.finalize();
  }
}
