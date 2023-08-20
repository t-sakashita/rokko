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

constexpr double eps = 1e-12;

TEST(pgemm, pdgemm) {
  constexpr int root_proc = 0;

  const rokko::grid grid(MPI_COMM_WORLD);
  std::mt19937 engine(123lu);
  std::uniform_real_distribution<> dist(-1.0, 1.0);

  constexpr int n = 64;
  constexpr double alpha = 3.3;
  constexpr double beta = 2.1;

  Eigen::MatrixXd a(n, n);
  Eigen::MatrixXd b(n, n);
  Eigen::MatrixXd c(n, n);
  if (grid.get_myrank() == root_proc) {
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        a(i, j) = dist(engine);
        b(i, j) = dist(engine);
        c(i, j) = dist(engine);
      }
    }
  }
  
  const rokko::mapping_bc<rokko::matrix_col_major> map(n, 1);
  rokko::distributed_matrix<double> da(map);
  rokko::distributed_matrix<double> db(map);
  rokko::distributed_matrix<double> dc(map);
  rokko::scatter(a, da, root_proc);
  rokko::scatter(b, db, root_proc);
  rokko::scatter(c, dc, root_proc);
  
  rokko::blas::gemm(CblasNoTrans, CblasNoTrans, alpha, a, b, beta, c);
  rokko::pblas::pgemm('N', 'N', alpha, da, db, beta, dc);
  
  const auto r = rokko::lapack::lange('F', c);
  const auto dr = rokko::scalapack::plange('F', dc);
  
  if (grid.get_myrank() == root_proc) {
    EXPECT_NEAR(r, dr, eps);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  const auto result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
