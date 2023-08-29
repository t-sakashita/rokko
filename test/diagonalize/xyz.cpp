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
#include <rokko/utility/xyz_lattice.hpp>
#include <rokko/utility/xyz_hamiltonian.hpp>
#include <rokko/utility/command_line_parameters.hpp>

#include <gtest/gtest.h>

int global_argc;
char** global_argv;

auto create_periodic_1dim_lattice_coupling(int L) {
  std::vector<std::pair<int, int>> lattice;
  std::vector<std::tuple<double, double, double>> coupling;
  const auto dim = 1 << L;

  const auto num_bonds = L-1;
  for (auto i=0; i<num_bonds; ++i) {
    lattice.emplace_back(std::make_pair(i, i+1));
    coupling.emplace_back(std::make_tuple(1, 1, 1));
  }

  assert(lattice.size() == num_bonds);
  assert(coupling.size() == num_bonds);

  rokko::print_lattice_coupling(L, lattice, coupling);

  return std::tuple{lattice, coupling, dim};
}

TEST(diagonalize, xyz) {
  const auto names = global_argc == 1 ? rokko::serial_dense_ev::solvers()
    : rokko::get_command_line_args(global_argc, global_argv);

  constexpr int L = 5;
  const auto [lattice, coupling, dim] = create_periodic_1dim_lattice_coupling(L);

  for(auto const& name : names) {
    std::cout << "library=" << name << std::endl;
    rokko::serial_dense_ev solver(name);
    solver.initialize(global_argc, global_argv);
    Eigen::MatrixXd mat(dim, dim);
    rokko::xyz_hamiltonian::generate(L, lattice, coupling, mat);
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
