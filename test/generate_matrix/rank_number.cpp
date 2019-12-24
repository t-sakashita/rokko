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

#include <rokko/distributed_matrix.hpp>
#include <rokko/collective.hpp>

#include <algorithm>
#include <gtest/gtest.h>

int global_argc;
char** global_argv;

template <typename MATRIX_MAJOR, typename GRID_MAJOR>
void run_test(std::array<int,2> const& global_size, std::array<int,2> const& block_size, GRID_MAJOR) {
  rokko::grid g(MPI_COMM_WORLD, GRID_MAJOR{});
  const auto grid_size = g.get_size();
  const int myrank = g.get_myrank();
  constexpr int root_proc = 0;
  rokko::mapping_bc<MATRIX_MAJOR> map(global_size, block_size, g);
  rokko::distributed_matrix<double, MATRIX_MAJOR> mat(map);

  for (int local_i=0; local_i<map.get_m_local(); ++local_i) {
    for (int local_j=0; local_j<map.get_n_local(); ++local_j) {
      mat.set_local(local_i, local_j, myrank);
    }
  }

  auto global_size_proc = (myrank == root_proc) ? global_size : std::array<int,2>({0,0});
  using eigen_matrix = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<MATRIX_MAJOR>>;
  eigen_matrix mat_loc(global_size_proc[0], global_size_proc[1]);
  rokko::gather(mat, mat_loc, root_proc);

  if (myrank == root_proc) {
    for (int i=0; i<global_size[0]; i += block_size[0]) {
      const int i_block = i / block_size[0];
      const int ip = i_block % grid_size[0];
      const int block0 = std::min(global_size[0] - i, block_size[0]);
      for (int j=0; j<global_size[1]; j += block_size[1]) {
        const int j_block = j / block_size[1];
        const int jp = j_block % grid_size[1];
        const int block1 = std::min(global_size[1] - j, block_size[1]);
        const int rank = g.calculate_rank_form_coords(ip, jp);
        ASSERT_TRUE(mat_loc.block(i, j, block0, block1) == eigen_matrix::Constant(block0, block1, rank));
      }
    }
  }
}

template <typename MATRIX_MAJOR>
void test_both_grid_major(std::array<int,2> const& global_size, std::array<int,2> const& block_size) {
  run_test<MATRIX_MAJOR>(global_size, block_size, rokko::grid_row_major);
  run_test<MATRIX_MAJOR>(global_size, block_size, rokko::grid_col_major);
}

TEST(generate_matrix, rank_number) {
  test_both_grid_major<rokko::matrix_col_major>({20,8}, {4,3});
  test_both_grid_major<rokko::matrix_col_major>({8,20}, {4,3});
  test_both_grid_major<rokko::matrix_col_major>({5,7}, {9,12});
  test_both_grid_major<rokko::matrix_col_major>({7,9}, {9,5});
}

int main(int argc, char** argv) {
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  global_argc = argc;
  global_argv = argv;
  result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
