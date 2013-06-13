/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2013 by Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <boost/shared_ptr.hpp>
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
  MPI_Init(&boost::unit_test::framework::master_test_suite().argc,
             &boost::unit_test::framework::master_test_suite().argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  const int dim = 100;

  boost::shared_ptr<rokko::solver_factory> instance(rokko::solver_factory::instance());

  std::vector<std::string> names;
  int argc = boost::unit_test::framework::master_test_suite().argc;
  if (argc == 1) {
    names = instance->solver_names();
  } else {
    for (int num=1; num < argc; ++num) {
      names.push_back(boost::unit_test::framework::master_test_suite().argv[num]);
    }
  }

  BOOST_FOREACH(std::string name, names) {
    std::cout << "solver=" << name << std::endl;
    rokko::solver solver(name);
    std::cout << "line=" << __LINE__ << std::endl;
    solver.initialize(boost::unit_test::framework::master_test_suite().argc,
                      boost::unit_test::framework::master_test_suite().argv);
    std::cout << "line=" << __LINE__ <<std::endl;
    rokko::grid g(comm, rokko::grid_col_major);
    std::cout << "line=" << __LINE__ <<std::endl;
    rokko::distributed_matrix<rokko::matrix_col_major> mat(dim, dim, g, solver);
    std::cout << "line=" << __LINE__ <<std::endl;
    rokko::frank_matrix::generate(mat);
    std::cout << "line=" << __LINE__ <<std::endl;
    rokko::localized_vector w(dim);
    std::cout << "line=" << __LINE__ <<std::endl;

    rokko::distributed_matrix<rokko::matrix_col_major> Z(dim, dim, g, solver);
    std::cout << "line=" << __LINE__ <<std::endl;

    solver.diagonalize(mat, w, Z);

    std::cout << "line=" << __LINE__ <<std::endl;
    double sum = 0;
    std::cout << "line=" << __LINE__ <<std::endl;

    for(int i=0; i<dim; ++i) {
      sum += w[i];
    }
    std::cout << "line=" << __LINE__ <<std::endl;

    BOOST_CHECK_CLOSE(sum, dim * (dim+1) * 0.5, 10e-5);
    std::cout << "line=" << __LINE__ <<std::endl;

    solver.finalize();
    std::cout << "line=" << __LINE__ <<std::endl;
  }
  MPI_Finalize();
}
