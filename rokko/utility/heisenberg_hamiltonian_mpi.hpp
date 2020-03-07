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

#ifndef ROKKO_UTILITY_HEISENBERG_HAMILTONIAN_MPI_HPP
#define ROKKO_UTILITY_HEISENBERG_HAMILTONIAN_MPI_HPP

#include "mpi.h"
#include <vector>
#include <iostream>
#include <rokko/distributed_matrix.hpp>
#include <rokko/utility/math.hpp>

namespace rokko {

namespace heisenberg_hamiltonian {

void multiply(const MPI_Comm& comm, int L, const std::vector<std::pair<int, int>>& lattice, const double* v, double* w, double* buffer) {
  int myrank, nproc;
  MPI_Status status;

  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &myrank);

  const int p = rokko::find_power_of_two(nproc);
  if (nproc != (1 << p)) {
    throw std::invalid_argument("This program can be run only with 2^n MPI processes");
  }

  int N = 1 << (L-p);

  for(int k=0; k<N; ++k) {
    w[k] = 0.;
  }

  for (std::size_t l = 0; l < lattice.size(); ++l) {
    int i = lattice[l].first;
    int j = lattice[l].second;
    if (i < (L-p)) {
      if (j < (L-p)) {
        int m1 = 1 << i;
        int m2 = 1 << j;
        int m3 = m1 + m2;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int k=0; k<N; ++k) {
          if ((k & m3) == m1) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
            w[k] += 0.5 * v[k^m3] - 0.25 * v[k];
          } else if ((k & m3) == m2) {
            w[k] += 0.5 * v[k^m3] - 0.25 * v[k];
          } else {
            w[k] += 0.25 * v[k];
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
#ifdef _OPENMP
#pragma omp parallel for
#endif
          for (int k=0; k<N; ++k) {
            if ((k & m1) == m1) {
              w[k] += 0.25 * v[k];
            } else {
              w[k] += 0.5 * buffer[k^m1] - 0.25 * v[k];
            }
          }
        } else {
#ifdef _OPENMP
#pragma omp parallel for
#endif
          for (int k=0; k<N; ++k) {
            if ((k & m1) == m1) {
              w[k] += 0.5 * buffer[k^m1] - 0.25 * v[k];
            } else {
              w[k] += 0.25 * v[k];
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
#ifdef _OPENMP
#pragma omp parallel for
#endif
          for (int k=0; k<N; ++k) {
            if ((k & m1) == m1) {
              w[k] += 0.25 * v[k];
            } else {
              w[k] += 0.5 * buffer[k^m1] - 0.25 * v[k];
            }
          }
        } else {
#ifdef _OPENMP
#pragma omp parallel for
#endif
          for (int k=0; k<N; ++k) {
            if ((k & m1) == m1) {
              w[k] += 0.5 * buffer[k^m1] - 0.25 * v[k];
            } else {
              w[k] += 0.25 * v[k];
            }
          }
        }
      } else {
        int m = (1 << (i-(L-p))) + (1 << (j-(L-p)));
        if (((myrank & m) != m) && ((myrank & m) != 0)) {
          MPI_Sendrecv(v, N, MPI_DOUBLE,
                       myrank ^ m, 0,
                       buffer, N, MPI_DOUBLE,
                       myrank ^ m, 0,
                       comm, &status);
#ifdef _OPENMP
#pragma omp parallel for
#endif
          for (int k=0; k<N; ++k) {
            w[k] += 0.5 * buffer[k] - 0.25 * v[k];
          }
        } else {
#ifdef _OPENMP
#pragma omp parallel for
#endif
          for (int k=0; k<N; ++k) {
            w[k] += 0.25 * v[k];
          }
        }
      }
    }
  }
}

void multiply(const MPI_Comm& comm, int L, const std::vector<std::pair<int, int>>& lattice, const std::vector<double>& v, std::vector<double>& w, std::vector<double>& buffer) {
  multiply(comm, L, lattice, v.data(), w.data(), buffer.data());
}


void fill_diagonal(const MPI_Comm& comm, int L, std::vector<std::pair<int, int>> const& lattice, double* w) {
  int myrank, nproc;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &myrank);

  const int p = rokko::find_power_of_two(nproc);
  if (nproc != (1 << p)) {
    throw std::invalid_argument("This program can be run only with 2^n MPI processes");
  }
  const int N = 1 << (L-p);
  const int myrank_shift = myrank * N;

  for (int k=0; k<N; ++k) {
    w[k] = 0;
  }

  for (int local_k=0; local_k<N; ++local_k) {
    const int k = local_k + myrank_shift;
    for (std::size_t l = 0; l < lattice.size(); ++l) {
      int i = lattice[l].first;
      int j = lattice[l].second;

      int m1 = 1 << i;
      int m2 = 1 << j;
      int m3 = m1 + m2;
      if (((k & m3) == m1) || ((k & m3) == m2)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
        w[local_k] -= 0.25;
      } else {
        w[local_k] += 0.25;
      }
    } // end for lattice
  } // end for k
}

void fill_diagonal(const MPI_Comm& comm, int L, std::vector<std::pair<int, int>>& lattice, std::vector<double>& w) {
  fill_diagonal(comm, L, lattice, w.data());
}

template<typename T, typename MATRIX_MAJOR>
void generate(int /* L */, std::vector<std::pair<int, int>>& lattice,
  rokko::distributed_matrix<T, MATRIX_MAJOR>& mat) {
  mat.set_zeros();
  for (std::size_t l = 0; l < lattice.size(); ++l) {
    int i = lattice[l].first;
    int j = lattice[l].second;
    int m1 = 1 << i;
    int m2 = 1 << j;
    int m3 = m1 + m2;
    for (int k = 0; k < mat.get_n_global(); ++k) {
      if (mat.is_gindex_mycol(k)) {
        int local_k = mat.translate_g2l_col(k);
        if (((k & m3) == m1) || ((k & m3) == m2)) {
          // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
          if (mat.is_gindex_myrow(k^m3)) {
            mat.update_local(mat.translate_g2l_row(k^m3), local_k, 0.5);
          }
          if (mat.is_gindex_myrow(k)) {
            mat.update_local(mat.translate_g2l_row(k), local_k, -0.25);
          }
        } else if (mat.is_gindex_myrow(k)) {
          mat.update_local(mat.translate_g2l_row(k), local_k, 0.25);
        }
      }
    }
  }
}

} // namespace heisenberg_hamiltonian

} // namespace rokko

#endif // ROKKO_UTILITY_HEISENBERG_HAMILTONIAN_MPI_HPP
