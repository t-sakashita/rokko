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
#include <rokko/utility/heisenberg_hamiltonian_mfree.hpp>
#include <rokko/utility/output_matrix_market.hpp>

#include <sstream>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

TEST(heisenberg_mfree, output_market_through_crs) {
  const std::string library = (global_argc >= 2) ? global_argv[1] : rokko::parallel_sparse_ev::default_solver();

  const int L = (global_argc >= 3) ? std::stoi(global_argv[2]) : 4;
  const auto dim = 1 << L;
  std::vector<std::pair<int, int>> lattice;
  for (int i = 0; i < L; ++i) lattice.emplace_back(std::make_pair(i, (i+1) % L));

  const rokko::heisenberg_mfree op(L, lattice);
  std::ostringstream os_direct, os_crs;

  rokko::output_matrix_market(op, os_direct);

  rokko::parallel_sparse_ev solver(library);
  const auto map = solver.default_mapping(dim, rokko::mpi_comm{op.get_comm()});
  const int num_entries_per_row = lattice.size() + 1;
  rokko::distributed_crs_matrix mat(map, num_entries_per_row);
  rokko::output_matrix_market(op, os_crs);

  ASSERT_TRUE(os_direct.str() == os_crs.str());

  solver.finalize();
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
