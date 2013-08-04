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

#ifndef ROKKO_UTILITY_SPIN_HAMILTONIAN_HPP
#define ROKKO_UTILITY_SPIN_HAMILTONIAN_HPP

#include <vector>
#include <rokko/localized_matrix.hpp>

namespace rokko {

namespace spin_hamiltonian {

void multiply(int L, std::vector<std::pair<int, int> >& lattice, const double* v, double* w) {
  int N = 1 << L;
  for (int k=0; k<N; ++k) {
    for (int l=0; l<lattice.size(); ++l) {
      int i = lattice[l].first;
      int j = lattice[l].second;
      //cout << "k=" << k << " i=" << i << " j=" << j << endl;
      int m1 = 1 << i;
      int m2 = 1 << j;
      int m3 = m1 + m2;
      if (((k & m3) == m1)) {  // when (bit i == 1, bit j == 0)
        w[k] += 0.5 * v[k^m3] - 0.25 * v[k];
        //cout << "if" << endl;
      }
      else if ((k & m3) == m2) { // when (bit i == 0, bit j == 1)
        w[k] += 0.5 * v[k^m3] - 0.25 * v[k];
      }
      else {
        w[k] += 0.25 * v[k];
        //      cout << "else" << endl;
      }
    }
  }
}

void multiply(int L, std::vector<std::pair<int, int> >& lattice, const std::vector<double>& v, std::vector<double>& w) {
  multiply(L, lattice, &v[0], &w[0]);
}

void fill_diagonal(int L, std::vector<std::pair<int, int> >& lattice, double* w) {
  int N = 1 << L;
  for (int k=0; k<N; ++k) {
    w[k] = 0;
    for (int l=0; l<lattice.size(); ++l) {
      int i = lattice[l].first;
      int j = lattice[l].second;
      //cout << "k=" << k << " i=" << i << " j=" << j << endl;
      int m1 = 1 << i;
      int m2 = 1 << j;
      int m3 = m1 + m2;
      if (((k & m3) == m1)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
        w[k] -= 0.25;
        //cout << "if" << endl;
      }
      else if ((k & m3) == m2) {
        w[k] -= 0.25;
      }
      else {
        w[k] += 0.25;
        //      cout << "else" << endl;
      }
    }
  }
}

void fill_diagonal(int L, std::vector<std::pair<int, int> >& lattice, std::vector<double>& w) {
  fill_diagonal(L, lattice, &w[0]);
}

template <typename MATRIX_MAJOR>
void generate(int L, std::vector<std::pair<int, int> >& lattice, rokko::localized_matrix<MATRIX_MAJOR>& mat) {
  mat.setZero();
  int N = 1 << L;
  for (int k1=0; k1<N; ++k1) {
    for (int k2=0; k2<N; ++k2) {
      for (int l=0; l<lattice.size(); ++l) {
        int i = lattice[l].first;
        int j = lattice[l].second;
        //cout << "k=" << k << " i=" << i << " j=" << j << endl;
        int m1 = 1 << i;
        int m2 = 1 << j;
        int m3 = m1 + m2;
        if (((k2 & m3) == m1) || ((k2 & m3) == m2)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
          if (k1 == (k2^m3)) {
            mat(k1, k2) += 0.5;
          }
          if (k1 == k2) {
            mat(k1, k2) -= 0.25; // 0.5 * v[k2^m3 ==] - 0.25 * v[k];
          }
        } else if (k1 == k2) {
          mat(k1,k2) += 0.25;
        }
      }
    }
  }
}

} // namespace spin_hamiltonian

} // namespace rokko

#endif ROKKO_UTILITY_SPIN_HAMILTONIAN_HPP
