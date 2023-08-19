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

#pragma once

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

  const auto p = rokko::find_power_of_two(nproc);
  if (nproc != (1 << p)) {
    throw std::invalid_argument("This program can be run only with 2^n MPI processes");
  }

  const auto N = 1 << (L-p);

  for(std::size_t k=0; k<N; ++k) {
    w[k] = 0.;
  }

  for (std::size_t l = 0; l < lattice.size(); ++l) {
    const auto i = lattice[l].first;
    const auto j = lattice[l].second;
    if (i < (L-p)) {
      if (j < (L-p)) {
        const auto m1 = 1 << i;
        const auto m2 = 1 << j;
        const auto m3 = m1 + m2;
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
        const auto m = 1 << (j-(L-p));
        MPI_Sendrecv(v, N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     buffer, N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     comm, &status);
        const auto m1 = 1 << i;
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
        const auto m = 1 << (i-(L-p));
        MPI_Sendrecv(v, N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     buffer, N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     comm, &status);
        const auto m1 = 1 << j;
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
        const auto m = (1 << (i-(L-p))) + (1 << (j-(L-p)));
        if (((myrank & m) != m) && ((myrank & m) != 0)) {
          MPI_Sendrecv(v, N, MPI_DOUBLE,
                       myrank ^ m, 0,
                       buffer, N, MPI_DOUBLE,
                       myrank ^ m, 0,
                       comm, &status);
#ifdef _OPENMP
#pragma omp parallel for
#endif
          for (std::size_t k=0; k<N; ++k) {
            w[k] += 0.5 * buffer[k] - 0.25 * v[k];
          }
        } else {
#ifdef _OPENMP
#pragma omp parallel for
#endif
          for (std::size_t k=0; k<N; ++k) {
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

  const auto p = rokko::find_power_of_two(nproc);
  if (nproc != (1 << p)) {
    throw std::invalid_argument("This program can be run only with 2^n MPI processes");
  }
  const auto N = 1 << (L-p);
  const auto myrank_shift = myrank * N;

  for (std::size_t k=0; k<N; ++k) {
    w[k] = 0;
  }

  for (int local_k=0; local_k<N; ++local_k) {
    const auto k = local_k + myrank_shift;
    for (std::size_t l = 0; l < lattice.size(); ++l) {
      const auto i = lattice[l].first;
      const auto j = lattice[l].second;

      const auto m1 = 1 << i;
      const auto m2 = 1 << j;
      const auto m3 = m1 + m2;
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
void generate(int /* L */, std::vector<std::pair<int, int>> const& lattice,
  rokko::distributed_matrix<T, MATRIX_MAJOR>& mat) {
  mat.set_zeros();
  for (std::size_t l = 0; l < lattice.size(); ++l) {
    const auto i = lattice[l].first;
    const auto j = lattice[l].second;
    const auto m1 = 1 << i;
    const auto m2 = 1 << j;
    const auto m3 = m1 + m2;
    for (int k = 0; k < mat.get_n_global(); ++k) {
      if (mat.has_global_col_index(k)) {
        const auto local_k = mat.translate_g2l_col(k);
        if (((k & m3) == m1) || ((k & m3) == m2)) {
          // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
          if (mat.has_global_row_index(k^m3)) {
            mat.update_local(mat.translate_g2l_row(k^m3), local_k, 0.5);
          }
          if (mat.has_global_row_index(k)) {
            mat.update_local(mat.translate_g2l_row(k), local_k, -0.25);
          }
        } else if (mat.has_global_row_index(k)) {
          mat.update_local(mat.translate_g2l_row(k), local_k, 0.25);
        }
      }
    }
  }
}

} // namespace heisenberg_hamiltonian

} // namespace rokko
