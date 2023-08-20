/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <iostream>

#include <rokko/eigen3.hpp>
#include <rokko/distributed_mfree.hpp>

#include <gtest/gtest.h>

void test_fill_diagonal(rokko::distributed_mfree_default const& mat) {
  const auto myrank = mat.get_mpi_comm().get_myrank();
  constexpr int root = 0;

  if (myrank == root)  std::cout << "dim = " << mat.get_dim() << std::endl;
  const auto N = mat.get_num_local_rows();
  Eigen::VectorXd v(N), w(N), diag(N);
  mat.fill_diagonal(diag.data());

  for (int i=0; i<mat.get_dim(); ++i) {
    const auto local_i = i - mat.start_row();
    v.setZero();
    if ((i >= mat.start_row()) && (i < mat.end_row()))  v(local_i) = 1.;
    w.setZero();
    mat.multiply(v.data(), w.data());
    if ((i >= mat.start_row()) && (i < mat.end_row()))  ASSERT_EQ(diag(local_i), w(local_i));
  }
}
