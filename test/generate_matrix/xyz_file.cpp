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
#include <rokko/utility/xyz_lattice.hpp>
#include <rokko/utility/xyz_hamiltonian.hpp>
#include <rokko/eigen3.hpp>

int main(int argc, char *argv[]) {
  if (argc <= 1) {
    throw std::invalid_argument("Specify input file name by command line argument");
  }
  const auto [num_sites, lattice, coupling] = rokko::read_lattice_file(argv[1]);
  rokko::print_lattice_coupling(num_sites, lattice, coupling);
  const auto dim = 1 << num_sites;
  const auto N = dim;
  std::cout << "dim=" << dim << std::endl;

  Eigen::MatrixXd mat1(N, N);
  std::cout << "multiply:" << std::endl;
  for (int i=0; i<N; ++i) {
    Eigen::VectorXd v(N), w(N);
    v.setZero();
    v(i) = 1;
    w.setZero();
    rokko::xyz_hamiltonian::multiply(num_sites, lattice, coupling, v, w);
    mat1.col(i) = w;
    std::cout << w.transpose() << std::endl;
  }

  std::cout << "fill_diagonal:" << std::endl;
  Eigen::VectorXd diagonal(N);
  rokko::xyz_hamiltonian::fill_diagonal(num_sites, lattice, coupling, diagonal);
  std::cout << "diagonal = " << diagonal.transpose() << std::endl;

  std::cout << "fill_matrix:" << std::endl;
  Eigen::MatrixXd mat2(N, N);
  rokko::xyz_hamiltonian::generate(num_sites, lattice, coupling, mat2);
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
