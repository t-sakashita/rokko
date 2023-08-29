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
#include <rokko/utility/xyz_lattice.hpp>
#include <rokko/utility/xyz_hamiltonian.hpp>
#include <rokko/utility/command_line_parameters.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

TEST(diagonalize, xyz_file) {
  const auto names = global_argc == 1 ? rokko::serial_dense_ev::solvers()
    : rokko::get_command_line_args(global_argc, global_argv);

  const auto [num_sites, lattice, coupling] = rokko::read_lattice_file("./heisenberg.ip");
  const auto dim = 1 << num_sites;

  std::cout << "num_sites=" << num_sites << " num_bonds=" << lattice.size() << std::endl;
  for (int i=0; i<coupling.size(); ++i) {
    std::cout << lattice[i].first << " " << lattice[i].second << " " << std::get<0>(coupling[i]) << " " << std::get<1>(coupling[i]) << " " << std::get<2>(coupling[i]) << std::endl;
  }

  for(auto const& name : names) {
    std::cout << "library=" << name << std::endl;
    rokko::serial_dense_ev solver(name);
    solver.initialize(global_argc, global_argv);
    Eigen::MatrixXd mat(dim, dim);
    rokko::xyz_hamiltonian::generate(num_sites, lattice, coupling, mat);
    Eigen::VectorXd w(dim);
    Eigen::MatrixXd Z(dim, dim);
    std::cout << "mat=" << mat << std::endl;
    solver.diagonalize(mat, w, Z);

    const auto sum = w.trace();

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
