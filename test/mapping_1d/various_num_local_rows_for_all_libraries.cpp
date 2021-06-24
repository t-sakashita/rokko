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

#include <rokko/parallel_sparse_ev.hpp>
#include <rokko/mapping_1d.hpp>
#include <rokko/skel/mapping_1d.hpp>
#include <rokko/distributed_crs_matrix.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

void check_map_mat(rokko::skel::mapping_1d const& map, rokko::distributed_crs_matrix const& mat) {
  ASSERT_EQ(map.get_dim(), mat.get_dim());
  ASSERT_EQ(map.get_num_local_rows(), mat.get_num_local_rows());
  ASSERT_EQ(map.start_row(), mat.start_row());
  ASSERT_EQ(map.end_row(), mat.end_row());
}

TEST(mapping_1d, various_num_local_rows_for_all_libraries) {
  constexpr int dim = 100;

  rokko::mpi_comm comm{MPI_COMM_WORLD};
  const int nprocs = comm.get_nprocs();
  const int myrank = comm.get_myrank();

  const int num_local_rows = (myrank == (nprocs-1)) ? dim - (nprocs-1) : 1;
  rokko::skel::mapping_1d skel_map(dim, num_local_rows, comm);

  for (auto const& name : rokko::parallel_sparse_ev::solvers()) {
    std::cout << "name=" << name << std::endl;
    rokko::parallel_sparse_ev solver(name);

    auto const map = solver.default_mapping(dim, num_local_rows, comm);
    constexpr int num_entries_per_row = 3;
    rokko::distributed_crs_matrix mat(map, num_entries_per_row);
    check_map_mat(skel_map, mat);
  }
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
