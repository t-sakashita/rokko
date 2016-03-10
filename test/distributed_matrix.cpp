/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2013 by Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
*               2014-2014    Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <boost/foreach.hpp>

#include <rokko/grid.hpp>
#include <rokko/solver.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/utility/frank_matrix.hpp>

#define BOOST_TEST_MODULE test_distributed_matrix
#ifndef BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

BOOST_AUTO_TEST_CASE(test_distributed_matrix) {
  int argc = boost::unit_test::framework::master_test_suite().argc;
  char** argv = boost::unit_test::framework::master_test_suite().argv;
  unsigned int dim = 1000;
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid g(comm);
  BOOST_FOREACH(std::string name, rokko::parallel_dense_ev::solvers()) {
    rokko::parallel_dense_ev solver(name);
    solver.initialize(argc, argv);
    rokko::mapping_bc<matrix_major> map = solver.default_mapping(dim, g);
    rokko::distributed_matrix<rokko::matrix_col_major> mat(map);
    rokko::frank_matrix::generate(mat);
    //    BOOST_CHECK_CLOSE(mat.get_global(0, 0), dim, 1e-14);
    solver.finalize();
  }
  MPI_Finalize();
}
