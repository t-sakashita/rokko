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

#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>

#include <rokko/utility/xyz_hamiltonian.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>

int main() {
  int L, N, num_bonds;
  L = 4;
  num_bonds = L - 1;
  N = 1 << L;
  std::vector<std::pair<int, int> > lattice;
  std::vector<std::tuple<double, double, double> > coupling;
  for (int i=0; i<L-1; ++i) {
    lattice.push_back(std::make_pair(i, i+1));
    coupling.push_back(std::make_tuple(1, 1, 1));
  }

  std::cout << "dim=" << N << std::endl;  
  std::cout << "L=" << L << std::endl;
  for (int i=0; i<num_bonds; ++i) {
    std::cout << lattice[i].first << " " << lattice[i].second << " " << std::get<0>(coupling[i]) << " " << std::get<1>(coupling[i]) << " " << std::get<2>(coupling[i]) << std::endl;
  }


  rokko::localized_matrix<double, rokko::matrix_col_major> mat1(N, N);
  std::cout << "multiply:" << std::endl;
  for (int i=0; i<N; ++i) {
    std::vector<double> v, w;
    v.assign(N, 0);
    v[i] = 1;
    w.assign(N, 0);
    rokko::xyz_hamiltonian::multiply(L, lattice, coupling, v, w);
    for (int j=0; j<N; ++j) {
      mat1(j,i) = w[j];
      std::cout << w[j] << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "fill_diagonal:" << std::endl;
  rokko::localized_vector<double> diagonal(N);
  std::vector<double> v(N);
  rokko::xyz_hamiltonian::fill_diagonal(L, lattice, coupling, v);
  for (int j=0; j<N; ++j) {
    diagonal(j) = v[j];
    std::cout << v[j] << " ";
  }
  std::cout << std::endl;

  std::cout << "fill_matrix:" << std::endl;
  rokko::localized_matrix<double, rokko::matrix_col_major> mat2(N, N);
  rokko::xyz_hamiltonian::generate(L, lattice, coupling, mat2);
  for (int i=0; i<N; ++i) {
    for (int j=0; j<N; ++j) {
      std::cout << mat2(i,j) << " ";
    }
    std::cout << std::endl;
  }

  if (mat1 == mat2) {
    std::cout << "OK: matrix by 'multiply' equals to a matrix by 'generate'." << std::endl;
  } else {
    std::cout << "ERROR: matrix by 'multiply' is differnet from a matrix by 'generate'."<< std::endl;
    exit(1);
  }

  if (diagonal == mat2.diagonal()) {
    std::cout << "OK: diagonal by 'fill_diagonal' equals to diagonal elementas of a matrix by 'genertate'."<< std::endl;
  } else {
    std::cout << "ERROR: diagonal by 'fill_diagonal' is differnet from diagonal elements of a matrix by 'genertate'."<< std::endl;
    exit(1);
  }

}


