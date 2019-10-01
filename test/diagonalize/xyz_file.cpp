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

#include <fstream>
#include <rokko/rokko.hpp>
#include <rokko/utility/xyz_hamiltonian.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

TEST(diagonalize, xyz_file) {
  std::vector<std::string> names;
  if (global_argc == 1) {
    names = rokko::serial_dense_ev::solvers();
  } else {
    for (int num=1; num < global_argc; ++num) {
      names.push_back(global_argv[num]);
    }
  }

  for(auto name : names) {
    std::ifstream ifs("./heisenberg.ip"); //str);
    if (!ifs) {
      std::cout << "can't open file" << std::endl;
      exit(1);
    }

    int L, num_bonds;
    std::vector<std::pair<int, int> > lattice;
    std::vector<std::tuple<double, double, double> > coupling;
    ifs >> L >> num_bonds;
    for (int i=0; i<num_bonds; ++i) {
      int j, k;
      ifs >> j >> k;
      lattice.push_back(std::make_pair(j, k));
    }

    for (int i=0; i<num_bonds; ++i) {
      double jx, jy, jz;
      ifs >> jx >> jy >> jz;
      coupling.push_back(std::make_tuple(jx, jy, jz));
    }

    std::cout << "L=" << L << " num_bonds=" << num_bonds << std::endl;
    for (int i=0; i<num_bonds; ++i) {
      std::cout << lattice[i].first << " " << lattice[i].second << " " << std::get<0>(coupling[i]) << " " << std::get<1>(coupling[i]) << " " << std::get<2>(coupling[i]) << std::endl;
    }
    int dim = 1 << L;

    std::cout << "solver=" << name << std::endl;
    rokko::serial_dense_ev solver(name);
    solver.initialize(global_argc, global_argv);
    rokko::localized_matrix<double, rokko::matrix_col_major> mat(dim, dim);
    rokko::xyz_hamiltonian::generate(L, lattice, coupling, mat);
    rokko::localized_vector<double> w(dim);
    rokko::localized_matrix<double, rokko::matrix_col_major> Z(dim, dim);
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
