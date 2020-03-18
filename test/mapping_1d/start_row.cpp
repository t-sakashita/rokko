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

#include <rokko/skel/mapping_1d.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

int start_row(rokko::skel::mapping_1d const& map) {
  int index = 0;
  for (int proc=0; proc<map.get_mpi_comm().get_myrank(); ++proc)
    index += map.calculate_num_local_rows(proc);
  return index;
}

int dim_by_sum(rokko::skel::mapping_1d const& map) {
  int sum = 0;
  for (int proc=0; proc<map.get_mpi_comm().get_nprocs(); ++proc)
    sum += map.calculate_num_local_rows(proc);
  return sum;
}

TEST(mapping_1d, start_row) {
  constexpr int dim = 100;
  rokko::skel::mapping_1d map(dim);

  ASSERT_EQ(map.start_row(), start_row(map));
  ASSERT_EQ(dim_by_sum(map), dim);

  const int end_proc = map.get_mpi_comm().get_nprocs() - 1;
  if (map.get_mpi_comm().get_myrank() == end_proc)
    ASSERT_EQ(map.end_row(), dim);
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
