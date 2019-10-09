/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <gtest/gtest.h>
#include <random>
#include <rokko/eigen3.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/blas.hpp>
#include <rokko/lapack.hpp>
#include <rokko/pblas.hpp>
#include <rokko/scalapack.hpp>
#include <rokko/collective.hpp>

TEST(pgemm, pdgemm) {
  rokko::grid grid(MPI_COMM_WORLD);
  std::mt19937 engine(123lu);
  std::uniform_real_distribution<> dist(-1.0, 1.0);

  int n = 64;
  double alpha = 3.3;
  double beta = 2.1;

  Eigen::MatrixXd a(n, n);
  Eigen::MatrixXd b(n, n);
  Eigen::MatrixXd c(n, n);
  if (grid.get_myrank() == 0) {
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        a(i, j) = dist(engine);
        b(i, j) = dist(engine);
        c(i, j) = dist(engine);
      }
    }
  }
  
  rokko::mapping_bc<rokko::matrix_col_major> map(n, 1);
  rokko::distributed_matrix<double> da(map);
  rokko::distributed_matrix<double> db(map);
  rokko::distributed_matrix<double> dc(map);
  rokko::scatter(a, da, 0);
  rokko::scatter(b, db, 0);
  rokko::scatter(c, dc, 0);
  
  rokko::blas::gemm(CblasNoTrans, CblasNoTrans, alpha, a, b, beta, c);
  rokko::pblas::pgemm('N', 'N', alpha, da, db, beta, dc);
  
  double r = rokko::lapack::lange('F', c);
  double dr = rokko::scalapack::plange('F', dc);
  
  if (grid.get_myrank() == 0) {
    EXPECT_NEAR(r, dr, 1.0e-12);
  }
}

int main(int argc, char** argv) {
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
