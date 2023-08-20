/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/mapping_bc.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

template <typename MATRIX_MAJOR, typename GRID_MAJOR>
void run_test(std::array<int,2> const& global_size, std::array<int,2> const& block_size, GRID_MAJOR) {
  rokko::grid g(MPI_COMM_WORLD, GRID_MAJOR{});
  rokko::mapping_bc<MATRIX_MAJOR> map(global_size, block_size, g);
  constexpr int root_proc = 0;

  if (g.get_myrank() == root_proc) {
    std::array<int,2> local_sum{0,0};

    for (int proc=0; proc<g.get_nprow(); ++proc) {
      local_sum[0] += map.template calculate_default_local_size<0>(proc);
    }
    ASSERT_EQ(local_sum[0], global_size[0]);

    for (int proc=0; proc<g.get_npcol(); ++proc) {
      local_sum[1] += map.template calculate_default_local_size<1>(proc);
    }
    ASSERT_EQ(local_sum[1], global_size[1]);
  }
}

template <typename MATRIX_MAJOR>
void test_both_grid_major(std::array<int,2> const& global_size, std::array<int,2> const& block_size) {
  run_test<MATRIX_MAJOR>(global_size, block_size, rokko::grid_row_major);
  run_test<MATRIX_MAJOR>(global_size, block_size, rokko::grid_col_major);
}

TEST(mapping_bc, sum_of_local_size) {
  test_both_grid_major<rokko::matrix_col_major>({20,8}, {4,3});
  test_both_grid_major<rokko::matrix_col_major>({8,20}, {4,3});
  test_both_grid_major<rokko::matrix_col_major>({5,7}, {9,12});
  test_both_grid_major<rokko::matrix_col_major>({7,9}, {9,5});
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
