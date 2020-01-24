/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
#include <rokko/utility/heisenberg_hamiltonian.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

TEST(diagonalize, heisenberg) {
  std::vector<std::string> names;
  if (global_argc == 1) {
    names = rokko::serial_dense_ev::solvers();
  } else {
    for (int num=1; num < global_argc; ++num) {
      names.emplace_back(global_argv[num]);
    }
  }

  for(auto name : names) {
    int L = 5;
    int dim = 1 << L;
    std::vector<std::pair<int, int>> lattice;
    for (int i=0; i<L-1; ++i) {
      lattice.emplace_back(std::make_pair(i, i+1));
    }

    std::cout << "solver=" << name << std::endl;
    rokko::serial_dense_ev solver(name);
    solver.initialize(global_argc, global_argv);
    Eigen::MatrixXd mat(dim, dim);
    rokko::heisenberg_hamiltonian::generate(L, lattice, mat);
    Eigen::VectorXd w(dim);
    Eigen::MatrixXd Z(dim, dim);
    std::cout << "mat=" << mat << std::endl;
    solver.diagonalize(mat, w, Z);
    
    double sum = 0;
    for(int i=0; i<dim; ++i) {
      sum += w[i];
    }

    std::cout << "w=" << w.transpose() << std::endl;
    
    solver.finalize();
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  global_argc = argc;
  global_argv = argv;
  return RUN_ALL_TESTS();
}
