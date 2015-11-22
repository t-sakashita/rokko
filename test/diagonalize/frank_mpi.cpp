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

BOOST_AUTO_TEST_CASE(test_solver) {
  int provided;
  MPI_Init_thread(&boost::unit_test::framework::master_test_suite().argc,
                  &boost::unit_test::framework::master_test_suite().argv,
                  MPI_THREAD_MULTIPLE, &provided);

  MPI_Comm comm = MPI_COMM_WORLD;
  const int dim = 10;

  std::vector<std::string> names;
  int argc = boost::unit_test::framework::master_test_suite().argc;
  if (argc == 1) {
    names = rokko::parallel_dense_ev::solvers();
  } else {
    for (int num=1; num < argc; ++num) {
      names.push_back(boost::unit_test::framework::master_test_suite().argv[num]);
    }
  }

  BOOST_FOREACH(std::string name, names) {
    std::cout << "solver=" << name << std::endl;
    rokko::parallel_dense_ev solver(name);
    if (solver.is_available_grid_major(rokko::grid_col_major)) {
      solver.initialize(boost::unit_test::framework::master_test_suite().argc,
                        boost::unit_test::framework::master_test_suite().argv);
      rokko::grid g(comm, rokko::grid_col_major);
      rokko::distributed_matrix<double, rokko::matrix_col_major> mat(dim, dim, g, solver);
      rokko::frank_matrix::generate(mat);
      rokko::localized_vector<double> w(dim);
      rokko::distributed_matrix<double, rokko::matrix_col_major> Z(dim, dim, g, solver);

      solver.diagonalize(mat, w, Z);
      
      double sum = 0;
      for(int i=0; i<dim; ++i) {
        sum += w[i];
      }
      BOOST_CHECK_CLOSE(sum, dim * (dim+1) * 0.5, 10e-5);
      
      solver.finalize();
    }
  }


  BOOST_FOREACH(std::string name, names) {
    std::cout << "solver=" << name << std::endl;
    rokko::parallel_dense_ev solver(name);
    if (solver.is_available_grid_major(rokko::grid_row_major)) {
      solver.initialize(boost::unit_test::framework::master_test_suite().argc,
                        boost::unit_test::framework::master_test_suite().argv);
      rokko::grid g(comm, rokko::grid_row_major);
      rokko::distributed_matrix<double, rokko::matrix_col_major> mat(dim, dim, g, solver);
      rokko::frank_matrix::generate(mat);
      rokko::localized_vector<double> w(dim);
      rokko::distributed_matrix<double, rokko::matrix_col_major> Z(dim, dim, g, solver);
      
      solver.diagonalize(mat, w, Z);
      
      double sum = 0;
      for(int i=0; i<dim; ++i) {
        sum += w[i];
      }
      BOOST_CHECK_CLOSE(sum, dim * (dim+1) * 0.5, 10e-5);

      solver.finalize();
    }
  }
  MPI_Finalize();
}
