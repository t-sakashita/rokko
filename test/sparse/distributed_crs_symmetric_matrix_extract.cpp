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

#include <rokko/rokko.hpp>
#include <rokko/utility/various_mpi_comm.hpp>

#include <gtest/gtest.h>

constexpr double eps = 1e-8;

int global_argc;
char** global_argv;

void run_test(std::string const& library, MPI_Comm comm) {
  constexpr int dim = 8;
  const std::vector<std::vector<int>> cols = {{0, 4}, {3}, {5}, {1, 7}, {0, 5, 6}, {2, 4}, {4, 7}, {3, 6}};
  const std::vector<std::vector<double>> values = {{7.1, 2.8}, {6.4}, {0.5}, {6.4, 3.5}, {2.8, 0.2, 1.4}, {0.5, 0.2}, {1.4, 4.3}, {3.5, 4.3}};

  const auto num_entries_per_row = std::max_element(cols.cbegin(), cols.cend(),
                                              [] (auto const& a, auto const& b) {
                                                return a.size() < b.size();
                                              })->size();

  std::cout << "library = " << library << std::endl;
  rokko::parallel_sparse_ev solver(library);
  if (comm != MPI_COMM_NULL) {
    const auto map = solver.default_mapping(dim, rokko::mpi_comm{comm});
    rokko::distributed_crs_matrix mat(map, num_entries_per_row);
    // storing
    for (int row = map.start_row(); row < map.end_row(); ++row) {
      mat.insert(row, cols[row], values[row]);
    }
    mat.complete();
    // checking
    std::vector<int> cols_check;
    std::vector<double> values_check;
    for (int row = map.start_row(); row < map.end_row(); ++row) {
      mat.extract(row, cols_check, values_check);
      for (size_t i=0; i<cols_check.size(); ++i) {
        ASSERT_EQ(cols_check[i], cols[row][i]);
        ASSERT_EQ(values_check[i], values[row][i]);
      }
    }
  }
  solver.finalize();
}

TEST(laplacian_crs, eigenvalue) {
  for(auto library : rokko::parallel_sparse_ev::solvers()) {
    run_test(library, MPI_COMM_WORLD);
    run_test(library, create_even_odd_comm_by_split());
    run_test(library, create_even_odd_comm());
    run_test(library, create_even_comm());  // MPI_COMM_NULL for odd rank number
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
