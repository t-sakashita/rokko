/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2017 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <gtest/gtest.h>
#include <cmath>
#include <rokko/eigen3.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/scalapack.hpp>
#include <rokko/lapack.hpp>
#include <rokko/collective.hpp>

constexpr double eps = 1e-10;

TEST(lange, pdlange) {
  rokko::grid grid(MPI_COMM_WORLD);

  int m = 20;
  int n = 30;

  // generate matrix
  Eigen::MatrixXd a(m, n);
  if (grid.get_myrank() == 0) a = Eigen::MatrixXd::Random(m, n);
  
  // distributed matrix
  rokko::mapping_bc<rokko::matrix_col_major> map({m, n}, grid, {1, 1});
  rokko::distributed_matrix<double, rokko::matrix_col_major> mat(map);
  rokko::scatter(a, mat, 0);

  double norm, expect;
  
  norm = rokko::scalapack::plange('M', mat);
  if (grid.get_myrank() == 0) {
    expect = rokko::lapack::lange('M', a);
    std::cout << "max norm of A = " << norm << " (expect: " << expect << ")" << std::endl;
    EXPECT_NEAR(expect, norm, eps);
  }
  
  norm = rokko::scalapack::plange('1', mat);
  if (grid.get_myrank() == 0) {
    expect = rokko::lapack::lange('1', a);
    std::cout << "1-norm of A = " << norm << " (expect: " << expect << ")" << std::endl;
    EXPECT_NEAR(expect, norm, eps);
  }
  
  norm = rokko::scalapack::plange('I', mat);
  if (grid.get_myrank() == 0) {
    expect = rokko::lapack::lange('I', a);
    std::cout << "infinity norm of A = " << norm << " (expect: " << expect << ")" << std::endl;
    EXPECT_NEAR(expect, norm, eps);
  }
  
  norm = rokko::scalapack::plange('F', mat);
  if (grid.get_myrank() == 0) {
    expect = rokko::lapack::lange('F', a);
    std::cout << "Frobenius norm of A = " << norm << " (expect: " << expect << ")" << std::endl;
    EXPECT_NEAR(expect, norm, eps);
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
