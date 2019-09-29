/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
#include <rokko/collective.hpp>
#include <boost/lexical_cast.hpp>
#include <random>
#define BOOST_TEST_MODULE test_product_v
#ifndef BOOST_TEST_DYN_LINK
# include <boost/test/included/unit_test.hpp>
#else
# include <boost/test/unit_test.hpp>
#endif

BOOST_AUTO_TEST_CASE(test_product_v) {
  MPI_Init(&boost::unit_test::framework::master_test_suite().argc,
           &boost::unit_test::framework::master_test_suite().argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  int dim = 100;
  if (boost::unit_test::framework::master_test_suite().argc > 1) {
    dim = boost::lexical_cast<int>(boost::unit_test::framework::master_test_suite().argv[1]);
  }

  std::mt19937 engine(123lu);
  std::uniform_real_distribution<> dist(-1.0, 1.0);

  if (rank == 0) std::cout << "dimension = " << dim << std::endl;
  rokko::parallel_dense_ev solver(rokko::parallel_dense_ev::default_solver());
  rokko::grid g(comm);
  rokko::mapping_bc<rokko::matrix_col_major> map = solver.default_mapping(dim, g);
  rokko::distributed_matrix<double, rokko::matrix_col_major> matA(map);
  rokko::mapping_bc<rokko::matrix_col_major> mapvec(dim, 1, g, 1, 1);
  rokko::distributed_matrix<double, rokko::matrix_col_major> vecX(mapvec);
  rokko::distributed_matrix<double, rokko::matrix_col_major> vecY(mapvec);
  rokko::localized_matrix<double, rokko::matrix_col_major> locA(dim, dim);
  rokko::localized_matrix<double, rokko::matrix_col_major> locX(dim, 1);
  rokko::localized_matrix<double, rokko::matrix_col_major> locY(dim, 1);
  for (int j = 0; j < dim; ++j) for (int i = 0; i < dim; ++i) locA(i, j) = dist(engine);
  for (int i = 0; i < dim; ++i) locX(i, 0) = dist(engine);
  for (int i = 0; i < dim; ++i) locY(i, 0) = dist(engine);
  rokko::scatter(locA, matA, 0);
  rokko::scatter(locX, vecX, 0);
  rokko::scatter(locY, vecY, 0);

  // local calculation
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      locY(i, 0) += locA(i, j) * locX(j, 0);

  // global calculation
  rokko::product_v(1.0, matA, false, vecX, false, 0, 1, vecY, false, 0);

  int success_local = 1;
  for (int i = 0; i < dim; ++i) {
    if (vecY.is_gindex(i, 0) && std::abs(vecY.get_global(i, 0) - locY(i, 0)) >  10e-12)
      success_local = 0;
  }
  int success;
  MPI_Allreduce(&success_local, &success, 1, MPI_INT, MPI_PROD, comm);
  BOOST_CHECK_EQUAL(success, true);
  MPI_Finalize();
}
