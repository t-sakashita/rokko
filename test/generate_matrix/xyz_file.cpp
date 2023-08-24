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

#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>

#include <rokko/utility/xyz_hamiltonian.hpp>
#include <rokko/eigen3.hpp>

int main(int argc, char *argv[]) {
  if (argc <= 1) {
    throw std::invalid_argument("Specify input file name by command line argument");
  }

  std::ifstream ifs(argv[1]);
  if (!ifs) {
    throw std::runtime_error("can't open file \"" + std::string(argv[1]) + "\"");
  }

  std::size_t L, num_bonds;
  std::vector<std::pair<int, int>> lattice;
  std::vector<std::tuple<double, double, double>> coupling;
  ifs >> L >> num_bonds;
  for (int i=0; i<num_bonds; ++i) {
    int j, k;
    ifs >> j >> k;
    lattice.emplace_back(std::make_pair(j, k));
  }

  for (int i=0; i<num_bonds; ++i) {
    double jx, jy, jz;
    ifs >> jx >> jy >> jz;
    coupling.emplace_back(std::make_tuple(jx, jy, jz));
  }

  std::cout << "L=" << L << " num_bonds=" << num_bonds << std::endl;
  for (int i=0; i<num_bonds; ++i) {
    std::cout << lattice[i].first << " " << lattice[i].second << " " << std::get<0>(coupling[i]) << " " << std::get<1>(coupling[i]) << " " << std::get<2>(coupling[i]) << std::endl;
  }
  const auto dim = 1 << L;
  const auto N = dim;
  std::cout << "dim=" << dim << std::endl;

  Eigen::MatrixXd mat1(N, N);
  std::cout << "multiply:" << std::endl;
  for (int i=0; i<N; ++i) {
    Eigen::VectorXd v(N), w(N);
    v.setZero();
    v(i) = 1;
    w.setZero();
    rokko::xyz_hamiltonian::multiply(L, lattice, coupling, v, w);
    mat1.col(i) = w;
    std::cout << w.transpose() << std::endl;
  }

  std::cout << "fill_diagonal:" << std::endl;
  Eigen::VectorXd diagonal(N);
  rokko::xyz_hamiltonian::fill_diagonal(L, lattice, coupling, diagonal);
  std::cout << "diagonal = " << diagonal.transpose() << std::endl;

  std::cout << "fill_matrix:" << std::endl;
  Eigen::MatrixXd mat2(N, N);
  rokko::xyz_hamiltonian::generate(L, lattice, coupling, mat2);
  std::cout << mat2 << std::endl;

  if (mat1 == mat2) {
    std::cout << "OK: matrix by 'multiply' equals to a matrix by 'generate'." << std::endl;
  } else {
    std::cout << "ERROR: matrix by 'multiply' is differnet from a matrix by 'generate'."<< std::endl;
    exit(1);
  }

  if (diagonal == mat2.diagonal()) {
    std::cout << "OK: diagonal by 'fill_diagonal' equals to diagonal elementas of a matrix by 'genertate'."<< std::endl;
  } else {
    std::cout << "ERROR: diagonal by 'fill_diagonal' is differnet from diagonal elementas of a matrix by 'genertate'."<< std::endl;
    exit(1);
  }

}


