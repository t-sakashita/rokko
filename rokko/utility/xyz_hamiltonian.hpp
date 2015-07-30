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

#ifndef ROKKO_UTILITY_XYZ_HAMILTONIAN_HPP
#define ROKKO_UTILITY_XYZ_HAMILTONIAN_HPP

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <rokko/localized_matrix.hpp>

namespace rokko {

namespace xyz_hamiltonian {

template<typename T>
void multiply(int L, const std::vector<std::pair<int, int> >& lattice,
  const std::vector<boost::tuple<double, double, double> >& coupling, const T* v, T* w) {
  int N = 1 << L;
  for (int l = 0; l < lattice.size(); ++l) {
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
      if (((k & m3) == m1) || ((k & m3) == m2)) {
        // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
        w[k] += diag_minus * v[k] + offdiag_plus * v[k^m3];
      } else {
        w[k] += diag_plus * v[k] + offdiag_minus * v[k^m3];
      }
    }
  }
}

template<typename T>
void multiply(int L, const std::vector<std::pair<int, int> >& lattice,
  const std::vector<boost::tuple<double, double, double> >& coupling, const std::vector<T>& v,
  std::vector<T>& w) {
  multiply(L, lattice, coupling, &v[0], &w[0]);
}

template<typename T>
void fill_diagonal(int L, const std::vector<std::pair<int, int> >& lattice,
  const std::vector<boost::tuple<double, double, double> >& coupling, T* w) {
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
      if (((k & m3) == m1) || ((k & m3) == m2)) {
        // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
        w[k] += diag_minus;
      } else {
        w[k] += diag_plus;
      }
    }
  }
}

template<typename T>
void fill_diagonal(int L, const std::vector<std::pair<int, int> >& lattice,
  const std::vector<boost::tuple<double, double, double> >& coupling, std::vector<T>& w) {
  fill_diagonal(L, lattice, coupling, &w[0]);
}

template<typename T, typename MATRIX_MAJOR>
void generate(int L, const std::vector<std::pair<int, int> >& lattice,
  const std::vector<boost::tuple<double, double, double> >& coupling,
  rokko::localized_matrix<T, MATRIX_MAJOR>& mat) {
  mat.set_zeros();
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
    for (int k=0; k<N; ++k) {
      if (((k & m3) == m1) || ((k & m3) == m2)) {
        // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
        mat(k^m3, k) += offdiag_plus;
        mat(k, k) += diag_minus;
      } else {
        mat(k^m3, k) += offdiag_minus;
        mat(k, k) += diag_plus;
      }
    }
  }
}

} // namespace xyz_hamiltonian

} // namespace rokko

#endif // ROKKO_UTILITY_XYZ_HAMILTONIAN_HPP
