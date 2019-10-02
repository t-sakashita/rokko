/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_UTILITY_HEISENBERG_HAMILTONIAN_HPP
#define ROKKO_UTILITY_HEISENBERG_HAMILTONIAN_HPP

#include <vector>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>

namespace rokko {

namespace heisenberg_hamiltonian {

void multiply(int L, const std::vector<std::pair<int, int>>& lattice, const double* v, double* w) {
  int N = 1 << L;
  for (std::size_t l = 0; l < lattice.size(); ++l) {
    int i = lattice[l].first;
    int j = lattice[l].second;
    int m1 = 1 << i;
    int m2 = 1 << j;
    int m3 = m1 + m2;
    for (int k = 0; k < N; ++k) {
      if (((k & m3) == m1) || ((k & m3) == m2)) {  // when (bit i == 1, bit j == 0) || (bit i == 0, bit j == 1) 
        w[k] += 0.5 * v[k^m3] - 0.25 * v[k];
      } else {
        w[k] += 0.25 * v[k];
      }
    }
  }
}

void multiply(int L, const std::vector<std::pair<int, int>>& lattice, const std::vector<double>& v, std::vector<double>& w) {
  multiply(L, lattice, &v[0], &w[0]);
}

void fill_diagonal(int L, const std::vector<std::pair<int, int>>& lattice, double* w) {
  int N = 1 << L;
  for (int k = 0; k < N; ++k) {
    w[k] = 0;
  }

  for (std::size_t l = 0; l < lattice.size(); ++l) {
    int i = lattice[l].first;
    int j = lattice[l].second;
    int m1 = 1 << i;
    int m2 = 1 << j;
    int m3 = m1 + m2;
    for (int k = 0; k < N; ++k) {
      if (((k & m3) == m1) || ((k & m3) == m2)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
        w[k] -= 0.25;
      }
      else {
        w[k] += 0.25;
      }
    }
  }
}

void fill_diagonal(int L, const std::vector<std::pair<int, int>>& lattice, std::vector<double>& w) {
  fill_diagonal(L, lattice, &w[0]);
}

template<typename T>
void fill_diagonal(int L, const std::vector<std::pair<int, int>>& lattice, Eigen::Vector<T>& w) {
  fill_diagonal(L, lattice, &w[0]);
}

template<typename T, typename MATRIX_MAJOR>
void generate(int L, const std::vector<std::pair<int, int>>& lattice, rokko::localized_matrix<T, MATRIX_MAJOR>& mat) {
  mat.setZero();
  int N = 1 << L;
  for (std::size_t l = 0; l < lattice.size(); ++l) {
    int i = lattice[l].first;
    int j = lattice[l].second;
    int m1 = 1 << i;
    int m2 = 1 << j;
    int m3 = m1 + m2;
    for (int k = 0; k < N; ++k) {
      if (((k & m3) == m1) || ((k & m3) == m2)) {
        // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
        mat(k^m3, k) += 0.5;
        mat(k, k) -= 0.25;
      } else {
        mat(k, k) += 0.25;
      }
    }
  }
}

} // namespace heisenberg_hamiltonian

} // namespace rokko

#endif // ROKKO_UTILITY_HEISENBERG_HAMILTONIAN_HPP
