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

#pragma once

#include <vector>
#include <tuple>
#include <rokko/eigen3.hpp>

namespace rokko {

namespace xyz_hamiltonian {

template<typename T>
void multiply(int L, const std::vector<std::pair<int, int>>& lattice,
  const std::vector<std::tuple<double, double, double>>& coupling, const T* v, T* w) {
  const auto N = 1 << L;
  for (std::size_t l = 0; l < lattice.size(); ++l) {
    const auto i = lattice[l].first;
    const auto j = lattice[l].second;
    const auto jx = std::get<0>(coupling[l]);
    const auto jy = std::get<1>(coupling[l]);
    const auto jz = std::get<2>(coupling[l]);

    const auto diag_plus = jz / 4.0;
    const auto diag_minus = - jz / 4.0;
    const auto offdiag_plus = (jx + jy) / 4.0;
    const auto offdiag_minus = (jx - jy) / 4.0;

    const auto m1 = 1 << i;
    const auto m2 = 1 << j;
    const auto m3 = m1 + m2;
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
void multiply(int L, const std::vector<std::pair<int, int>>& lattice,
  const std::vector<std::tuple<double, double, double>>& coupling, const std::vector<T>& v,
  std::vector<T>& w) {
  multiply(L, lattice, coupling, v.data(), w.data());
}

template<typename T>
void multiply(int L, const std::vector<std::pair<int, int>>& lattice,
  const std::vector<std::tuple<double, double, double>>& coupling, const Eigen::Vector<T>& v,
  Eigen::Vector<T>& w) {
  multiply(L, lattice, coupling, v.data(), w.data());
}

template<typename T>
void fill_diagonal(int L, const std::vector<std::pair<int, int>>& lattice,
  const std::vector<std::tuple<double, double, double>>& coupling, T* w) {
  const auto N = 1 << L;
  for (std::size_t k=0; k<N; ++k) {
    w[k] = 0;
  }
  for (std::size_t l = 0; l < lattice.size(); ++l) {
    const auto i = lattice[l].first;
    const auto j = lattice[l].second;
    const auto jz = std::get<2>(coupling[l]); // jx and jy are unused

    const auto diag_plus = jz / 4.0;
    const auto diag_minus = - jz / 4.0;

    const auto m1 = 1 << i;
    const auto m2 = 1 << j;
    const auto m3 = m1 + m2;
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
void fill_diagonal(int L, const std::vector<std::pair<int, int>>& lattice,
  const std::vector<std::tuple<double, double, double>>& coupling, std::vector<T>& w) {
  fill_diagonal(L, lattice, coupling, w.data());
}

template<typename T>
void fill_diagonal(int L, const std::vector<std::pair<int, int>>& lattice,
  const std::vector<std::tuple<double, double, double>>& coupling, Eigen::Vector<T>& w) {
  fill_diagonal(L, lattice, coupling, w.data());
}

template<typename T, int MATRIX_MAJOR>
void generate(int L, const std::vector<std::pair<int, int>>& lattice,
  const std::vector<std::tuple<T, T, T>>& coupling,
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat) {
  mat.setZero();
  const auto N = 1 << L;
  for (std::size_t l = 0; l < lattice.size(); ++l) {
    const auto i = lattice[l].first;
    const auto j = lattice[l].second;
    const auto jx = std::get<0>(coupling[l]);
    const auto jy = std::get<1>(coupling[l]);
    const auto jz = std::get<2>(coupling[l]);
    const auto diag_plus = jz / 4.0;
    const auto diag_minus = - jz/ 4.0;
    const auto offdiag_plus = (jx + jy) / 4.0;
    const auto offdiag_minus = (jx - jy) / 4.0;

    const auto m1 = 1 << i;
    const auto m2 = 1 << j;
    const auto m3 = m1 + m2;
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
