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

#include <rokko/solver.hpp>
#ifdef ROKKO_HAVE_MPI
# include <mpi.h>
#endif

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

TEST(solver, all) {
  for(auto name : rokko::serial_dense_ev::solvers()) {
    std::cerr << name << std::endl;
    rokko::serial_dense_ev solver(name);
    solver.initialize(global_argc, global_argv);
    solver.finalize();
  }

#ifdef ROKKO_HAVE_MPI
  int provided;
  MPI_Init_thread(&global_argc, &global_argv, MPI_THREAD_MULTIPLE, &provided);

#ifdef ROKKO_HAVE_PARALLEL_DENSE_SOLVER
  for(auto name : rokko::parallel_dense_ev::solvers()) {
    std::cerr << name << std::endl;
    rokko::parallel_dense_ev solver(name);
    solver.initialize(global_argc, global_argv);
    solver.finalize();
  }
#endif // ROKKO_HAVE_PARALLEL_DENSE_SOLVER

#ifdef ROKKO_HAVE_PARALLEL_SPARSE_SOLVER
  for(auto name : rokko::parallel_sparse_ev::solvers()) {
    std::cerr << name << std::endl;
    rokko::parallel_sparse_ev solver(name);
    solver.initialize(global_argc, global_argv);
    solver.finalize();
  }
#endif // ROKKO_HAVE_PARALLEL_SPARSE_SOLVER

  MPI_Finalize();
#endif // ROKKO_HAVE_MPI
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  global_argc = argc;
  global_argv = argv;
  return RUN_ALL_TESTS();
}
