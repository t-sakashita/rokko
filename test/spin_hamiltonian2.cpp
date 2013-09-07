/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <iostream>
#include <vector>
#include <boost/tuple/tuple.hpp>

#include <rokko/utility/spin_hamiltonian2.hpp>
#include <rokko/localized_matrix.hpp>

int main()
{
  int L = 3;
  int N = 1 << L;
  std::vector<std::pair<int, int> > lattice;
  std::vector<boost::tuple<double, double, double> > vec_j;
  for (int i=0; i<L-1; ++i) {
    lattice.push_back(std::make_pair(i, i+1));
    vec_j.push_back(boost::make_tuple(1., 1., 1.));
  }

  std::cout << "multiply:" << std::endl;
  for (int i=0; i<N; ++i) {
    std::vector<double> v, w;
    v.assign(N, 0);
    v[i] = 1;
    w.assign(N, 0);
    rokko::spin_hamiltonian2::multiply(L, lattice, vec_j, v, w);
    for (int j=0; j<N; ++j) {
      std::cout << w[j] << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "fill_diagonal:" << std::endl;
  std::vector<double> v(N);
  rokko::spin_hamiltonian2::fill_diagonal(L, lattice, vec_j, v);
  for (int j=0; j<N; ++j) {
    std::cout << v[j] << " ";
  }
  std::cout << std::endl;

  std::cout << "fill_matrix:" << std::endl;
  rokko::localized_matrix<rokko::matrix_col_major> mat(N, N);
  rokko::spin_hamiltonian2::generate(L, lattice, vec_j, mat);
  for (int i=0; i<N; ++i) {
    for (int j=0; j<N; ++j) {
      std::cout << mat(i,j) << " ";
    }
    std::cout << std::endl;
  }

}


