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

#include <rokko/utility/heisenberg_hamiltonian_mfree.hpp>
#include "test_fill_diagonal.hpp"

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

TEST(heisenberg_mfree, fill_diagonal) {
  constexpr std::size_t L = 8;
  std::vector<std::pair<int, int>> lattice;
  for (std::size_t i=0; i<L; ++i) lattice.emplace_back(std::make_pair(i, (i+1) % L));

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  constexpr int root = 0;

  if (myrank == root) {
    std::cout << "L=" << L << " lattice.size=" << lattice.size() << std::endl;
    for (std::size_t i=0; i < lattice.size(); ++i) {
      std::cout << lattice[i].first << " " << lattice[i].second << std::endl;
    }
  }

  rokko::heisenberg_mfree mat(L, lattice);
  test_fill_diagonal(mat);
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
