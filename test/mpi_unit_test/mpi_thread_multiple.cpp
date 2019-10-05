/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <mpi.h>
#include <iostream>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

TEST(mpi, mpi_thread_multiple) {
  const int required = MPI_THREAD_MULTIPLE;
  int provided;
  MPI_Init_thread(&global_argc, &global_argv, required, &provided);
  if (provided == required) {
     std::cerr << "MPI_THREAD_MULTIPLE is supported\n";
   } else {
     std::cerr << "MPI_THREAD_MULTIPLE is NOT supported\n";
   }
  MPI_Finalize();
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  global_argc = argc;
  global_argv = argv;
  return RUN_ALL_TESTS();
}
