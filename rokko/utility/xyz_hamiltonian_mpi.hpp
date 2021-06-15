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

#include "mpi.h"
#include <vector>
#include <tuple>

#include <iostream>

#include <rokko/distributed_matrix.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/math.hpp>

namespace rokko {

namespace xyz_hamiltonian {

void multiply(const MPI_Comm& comm, int L, const std::vector<std::pair<int, int>>& lattice, const std::vector<std::tuple<double, double, double>>& coupling, const double* v, double* w, double* buffer) {
  int myrank, nproc;
  MPI_Status status;

  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &myrank);

  const int p = rokko::find_power_of_two(nproc);
  if (nproc != (1 << p)) {
    throw std::invalid_argument("This program can be run only with 2^n MPI processes");
  }
  int N = 1 << (L-p);

  for (std::size_t l = 0; l < lattice.size(); ++l) {
    int i = lattice[l].first;
    int j = lattice[l].second;
    double jx = std::get<0>(coupling[l]);
    double jy = std::get<1>(coupling[l]);
    double jz = std::get<2>(coupling[l]);

    double diag_plus = jz / 4.0;
    double diag_minus = - jz / 4.0;
    double offdiag_plus = (jx + jy) / 4.0;
    double offdiag_minus = (jx - jy) / 4.0;

    if (i < (L-p)) {
      if (j < (L-p)) {
        int m1 = 1 << i;
        int m2 = 1 << j;
        int m3 = m1 + m2;
        for (int k = 0; k < N; ++k) {
          if (((k & m3) == m1) || ((k & m3) == m2)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
            w[k] += diag_minus * v[k] + offdiag_plus * v[k^m3];
          } else {
            w[k] += diag_plus * v[k] + offdiag_minus * v[k^m3];
          }
        }
      } else {
        int m = 1 << (j-(L-p));
        MPI_Sendrecv(v, N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     buffer, N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     comm, &status);
        int m1 = 1 << i;
        if ((myrank & m) == m) { 
          for (int k = 0; k < N; ++k) {
            if ((k & m1) == m1) {
              w[k] += diag_plus * v[k] + offdiag_minus * buffer[k^m1];
            } else {
              w[k] += diag_minus * v[k] + offdiag_plus * buffer[k^m1];
            }
          }
        } else {
          for (int k = 0; k < N; ++k) {
            if ((k & m1) == m1) {
              w[k] += diag_minus * v[k] + offdiag_plus * buffer[k^m1];
            } else {
              w[k] += diag_plus * v[k] + offdiag_minus * buffer[k^m1];
            }
          }
        }
      }
    } else {
      if (j < (L-p)) {
        int m = 1 << (i-(L-p));
        MPI_Sendrecv(v, N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     buffer, N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     comm, &status);
        int m1 = 1 << j;
        if ((myrank & m) == m) {
          for (int k = 0; k < N; ++k) {
            if ((k & m1) == m1) {
              w[k] += diag_plus * v[k] + offdiag_minus * buffer[k^m1];
            } else {
              w[k] += diag_minus * v[k] + offdiag_plus * buffer[k^m1];
            }
          }
        } else {
          for (int k = 0; k < N; ++k) {
            if ((k & m1) == m1) {
              w[k] += diag_minus * v[k] + offdiag_plus * buffer[k^m1];
            } else {
              w[k] += diag_plus * v[k] + offdiag_minus * buffer[k^m1];
            }
          }
        }
      } else {
        int m = (1 << (i-(L-p))) + (1 << (j-(L-p)));
        MPI_Sendrecv(v, N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     buffer, N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     comm, &status);
        if (((myrank & m) != m) && ((myrank & m) != 0)) {
          for (int k = 0; k < N; ++k) {
            w[k] += diag_minus * v[k] + offdiag_plus * buffer[k];
          }
        } else {
          for (int k = 0; k < N; ++k) {
            w[k] += diag_plus * v[k] + offdiag_minus * buffer[k];
          }
        }
      }
    }
  }
}

void multiply(const MPI_Comm& comm, int L, const std::vector<std::pair<int, int>>& lattice, const std::vector<std::tuple<double, double, double>>& coupling, const std::vector<double>& v, std::vector<double>& w, std::vector<double>& buffer) {
  multiply(comm, L, lattice, coupling, v.data(), w.data(), buffer.data());
}

void fill_diagonal(const MPI_Comm& comm, int L, const std::vector<std::pair<int, int>>& lattice, const std::vector<std::tuple<double, double, double>>& coupling, double* w) {
  int myrank, nproc;

  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &myrank);

  const int p = rokko::find_power_of_two(nproc);
  if (nproc != (1 << p)) {
    throw std::invalid_argument("This program can be run only with 2^n MPI processes");
  }

  int N_seq = 1 << L;
  int N = 1 << (L-p);
  int myrank_shift = myrank * N;
  int nproc_shift = (nproc-1) * N;
  int mask = N - 1;

  for (int k = 0; k < N; ++k) {
    w[k] = 0;
  }

  for (std::size_t l = 0; l < lattice.size(); ++l) {
    int i = lattice[l].first;
    int j = lattice[l].second;
    double jz = std::get<2>(coupling[l]); // jx and jy are unused
    double diag_plus = jz / 4.0;
    double diag_minus = - jz / 4.0;

    int m1 = 1 << i;
    int m2 = 1 << j;
    int m3 = m1 + m2;

    for (int k = 0; k < N_seq; ++k) {
      if (myrank_shift == (k & nproc_shift)) {
        if (((k & m3) == m1) || ((k & m3) == m2)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
          w[k & mask] += diag_minus;
        } else {
          w[k & mask] += diag_plus;
        }        
      }
    }  // end for k
  } // end for lattice
}

void fill_diagonal(const MPI_Comm& comm, int L, const std::vector<std::pair<int, int>>& lattice, const std::vector<std::tuple<double, double, double>>& coupling, std::vector<double>& w) {
  fill_diagonal(comm, L, lattice, coupling, w.data());
}

template<typename T, typename MATRIX_MAJOR>
void generate(int L, const std::vector<std::pair<int, int>>& lattice,
  const std::vector<std::tuple<double, double, double>>& coupling,
  rokko::distributed_matrix<T, MATRIX_MAJOR>& mat) {
  const auto& map = mat.get_mapping();

  mat.set_zeros();
  int N = 1 << L;
  for (std::size_t l = 0; l < lattice.size(); ++l) {
    int i = lattice[l].first;
    int j = lattice[l].second;
    double jx = std::get<0>(coupling[l]);
    double jy = std::get<1>(coupling[l]);
    double jz = std::get<2>(coupling[l]);
    double diag_plus = jz / 4.0;
    double diag_minus = - jz/ 4.0;
    double offdiag_plus = (jx + jy) / 4.0;
    double offdiag_minus = (jx - jy) / 4.0;

    int m1 = 1 << i;
    int m2 = 1 << j;
    int m3 = m1 + m2;

    for (int k = 0; k < N; ++k) {
      if (map.has_global_col_index(k)) {
        int local_k = map.translate_g2l_col(k);
        if (((k & m3) == m1) || ((k & m3) == m2)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
          if (map.has_global_row_index(k^m3)) {
            mat.update_local(map.translate_g2l_row(k^m3), local_k, offdiag_plus);
          }
          if (map.has_global_row_index(k)) {
            mat.update_local(map.translate_g2l_row(k), local_k, diag_minus);
          }
        } else {
          if (map.has_global_row_index(k^m3)) {
            mat.update_local(map.translate_g2l_row(k^m3), local_k, offdiag_minus);
          }
          if (map.has_global_row_index(k)) {
            mat.update_local(map.translate_g2l_row(k), local_k, diag_plus);
          }
        }
      }
    }
  }
}

// The following routine uses local indices.  It works correctly.
/*
template <typename MATRIX_MAJOR>
void generate(int L, const std::vector<std::pair<int, int>>& lattice, const std::vector<std::tuple<double, double, double>>& coupling, rokko::distributed_matrix<MATRIX_MAJOR>& mat) {
  const auto& map = mat.get_mapping();

  mat.set_zeros();
  int N = 1 << L;
  for (int l=0; l<lattice.size(); ++l) {
    int i = lattice[l].first;
    int j = lattice[l].second;
    double jx = std::get<0>(coupling[l]);
    double jy = std::get<1>(coupling[l]);
    double jz = std::get<2>(coupling[l]);
    double diag_plus = jz / 4.0;
    double diag_minus = - jz/ 4.0;
    double offdiag_plus = (jx + jy) / 4.0;
    double offdiag_minus = (jx - jy) / 4.0;

    int m1 = 1 << i;
    int m2 = 1 << j;
    int m3 = m1 + m2;
    for(int local_i = 0; local_i < map.get_m_local(); ++local_i) {
      int k1 = map.translate_l2g_row(local_i);
      for(int local_j = 0; local_j < map.get_n_local(); ++local_j) {
        int k2 = map.translate_l2g_col(local_j);
        if (((k2 & m3) == m1) || ((k2 & m3) == m2)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
          if (k1 == (k2^m3)) {
            mat.update_local(local_i, local_j, offdiag_plus);
          }
          if (k1 == k2) {
            mat.update_local(local_i, local_j, diag_minus);
          }
        } else {
          if (k1 == (k2^m3)) {
            mat.update_local(local_i, local_j, offdiag_minus);
          }
          if (k1 == k2) {
            //std::cout << "k1=" << k1 << " k2=" << k2 << std::endl;
            mat.update_local(local_i, local_j, diag_plus);
          }
        }
      }
    }
  }
}
*/

} // namespace xyz_hamiltonian

} // namespace rokko
