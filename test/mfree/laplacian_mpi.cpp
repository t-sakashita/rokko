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

#include <rokko/utility/laplacian_mfree.hpp>
#include "test_fill_diagonal.hpp"

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

TEST(laplacian_mfree, fill_diagonal) {
  constexpr std::size_t dim = 20;

  const rokko::laplacian_mfree mat(dim, MPI_COMM_WORLD);
  test_fill_diagonal(mat);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  global_argc = argc;
  global_argv = argv;
  const auto result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
