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

#ifndef ROKKO_UTILITY_XYZ_HAMILTONIAN_HPP
#define ROKKO_UTILITY_XYZ_HAMILTONIAN_HPP

#include <vector>
#include <boost/tuple/tuple.hpp>

#include <rokko/localized_matrix.hpp>

namespace rokko {

namespace xyz_hamiltonian {

void multiply(int L, const std::vector<std::pair<int, int> >& lattice, const std::vector<boost::tuple<double, double, double> >& coupling, const double* v, double* w) {
  int N = 1 << L;
  for (int l=0; l<lattice.size(); ++l) {
    int i = lattice[l].first;
    int j = lattice[l].second;
    double jx = coupling[l].get<0>();
    double jy = coupling[l].get<1>();
    double jz = coupling[l].get<2>();

    double diag_plus = jz / 4.0;
    double diag_minus = - jz / 4.0;
    double offdiag_plus = (jx + jy) / 4.0;
    double offdiag_minus = (jx - jy) / 4.0;

    int m1 = 1 << i;
    int m2 = 1 << j;
    int m3 = m1 + m2;
    for (int k=0; k<N; ++k) {
      if (((k & m3) == m1) || ((k & m3) == m2)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
        w[k] += diag_minus * v[k] + offdiag_plus * v[k^m3];
      } else {
        w[k] += diag_plus * v[k] + offdiag_minus * v[k^m3];
      }
    }
  }
}

void multiply(int L, const std::vector<std::pair<int, int> >& lattice, const std::vector<boost::tuple<double, double, double> >& coupling, const std::vector<double>& v, std::vector<double>& w) {
  multiply(L, lattice, coupling, &v[0], &w[0]);
}

void fill_diagonal(int L, const std::vector<std::pair<int, int> >& lattice, const std::vector<boost::tuple<double, double, double> >& coupling, double* w) {
  int N = 1 << L;
  for (int k=0; k<N; ++k) {
    w[k] = 0;
  }
  for (int l=0; l<lattice.size(); ++l) {
    int i = lattice[l].first;
    int j = lattice[l].second;
    double jx = coupling[l].get<0>();
    double jy = coupling[l].get<1>();
    double jz = coupling[l].get<2>();

    double diag_plus = jz / 4.0;
    double diag_minus = - jz / 4.0;
    double offdiag_plus = (jx + jy) / 4.0;
    double offdiag_minus = (jx - jy) / 4.0;

    int m1 = 1 << i;
    int m2 = 1 << j;
    int m3 = m1 + m2;
    for (int k=0; k<N; ++k) {
      if (((k & m3) == m1) || ((k & m3) == m2)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
        w[k] += diag_minus;
      } else {
        w[k] += diag_plus;
      }
    }
  }
}

void fill_diagonal(int L, const std::vector<std::pair<int, int> >& lattice, const std::vector<boost::tuple<double, double, double> >& coupling, std::vector<double>& w) {
  fill_diagonal(L, lattice, coupling, &w[0]);
}

template <typename MATRIX_MAJOR>
void generate(int L, const std::vector<std::pair<int, int> >& lattice, const std::vector<boost::tuple<double, double, double> >& coupling, rokko::localized_matrix<MATRIX_MAJOR>& mat) {
  mat.setZero();
  int N = 1 << L;
  for (int l=0; l<lattice.size(); ++l) {
    int i = lattice[l].first;
    int j = lattice[l].second;
    double jx = coupling[l].get<0>();
    double jy = coupling[l].get<1>();
    double jz = coupling[l].get<2>();
    double diag_plus = jz / 4.0;
    double diag_minus = - jz/ 4.0;
    double offdiag_plus = (jx + jy) / 4.0;
    double offdiag_minus = (jx - jy) / 4.0;

    int m1 = 1 << i;
    int m2 = 1 << j;
    int m3 = m1 + m2;
    for (int k1=0; k1<N; ++k1) {
      for (int k2=0; k2<N; ++k2) {
        if (((k2 & m3) == m1) || ((k2 & m3) == m2)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
          if (k1 == (k2^m3)) {
            mat(k1, k2) += offdiag_plus;
          }
          if (k1 == k2) {
            mat(k1, k2) += diag_minus;
          }
        } else {
          if (k1 == (k2^m3)) {
            mat(k1, k2) += offdiag_minus;
          }
          if (k1 == k2) {
            mat(k1, k2) += diag_plus;
          }
        }
      }
    }
  }
}

} // namespace spin_hamiltonian

} // namespace rokko

#endif // ROKKO_UTILITY_XYZ_HAMILTONIAN_HPP
