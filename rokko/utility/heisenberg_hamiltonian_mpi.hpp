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

#ifndef ROKKO_UTILITY_SPIN_HAMILTONIANPARALLEL_HPP
#define ROKKO_UTILITY_SPIN_HAMILTONIANPARALLEL_HPP

#include "mpi.h"
#include <vector>
#include <iostream>

//#include <rokko/localized_matrix.hpp>
#include <rokko/distributed_matrix.hpp>

namespace rokko {

namespace heisenberg_hamiltonian {

void multiply(const MPI_Comm& comm, int L, const std::vector<std::pair<int, int> >& lattice, const double* v, double* w, double* buffer) {
  int myrank, nproc;
  MPI_Status status;
  int ierr;

  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &myrank);

  int n = nproc;
  int p = -1;
  do {
    n /= 2;
    ++p;
  } while (n > 0);

  if (nproc != (1 << p)) {
    if ( myrank == 0 ) {
      std::cout << "This program can be run only for powers of 2" << std::endl;
    }
    MPI_Abort(comm, 1);
  }

  int N = 1 << (L-p);

  for (int l=0; l<lattice.size(); ++l) {
    int i = lattice[l].first;
    int j = lattice[l].second;
    if (i < (L-p)) {
      if (j < (L-p)) {
        int m1 = 1 << i;
        int m2 = 1 << j;
        int m3 = m1 + m2;
        for (int k=0; k<N; ++k) {
          if (((k & m3) == m1)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
            w[k] += 0.5 * v[k^m3] - 0.25 * v[k];
          } else if ((k & m3) == m2) {
            w[k] += 0.5 * v[k^m3] - 0.25 * v[k];
          } else {
            w[k] += 0.25 * v[k];
          }
        }
      } else {
        int m = 1 << (j-(L-p));
        MPI_Sendrecv((void*)&v[0], N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     &buffer[0], N, MPI_DOUBLE, 
                     myrank ^ m, 0,
                     comm, &status);
        int m1 = 1 << i;
        if ((myrank & m) == m) { 
          for (int k=0; k<N; ++k) {
            if ((k & m1) == m1) {
              w[k] += 0.25 * v[k];
            } else {
              w[k] += 0.5 * buffer[k^m1] - 0.25 * v[k];
            }
          }
        } else {
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
        MPI_Sendrecv((void*)&v[0], N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     &buffer[0], N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     comm, &status);
        int m1 = 1 << j;
        if ((myrank & m) == m) {
          for (int k=0; k<N; ++k) {
            if ((k & m1) == m1) {
              w[k] += 0.25 * v[k];
            } else {
              w[k] += 0.5 * buffer[k^m1] - 0.25 * v[k];
            }
          }
        } else {
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
          MPI_Sendrecv((void*)&v[0], N, MPI_DOUBLE,
                       myrank ^ m, 0,
                       &buffer[0], N, MPI_DOUBLE,
                       myrank ^ m, 0,
                       comm, &status);
          for (int k=0; k<N; ++k) {
            w[k] += 0.5 * buffer[k] - 0.25 * v[k];
          }
        } else {
          for (int k=0; k<N; ++k) {
            w[k] += 0.25 * v[k];
          }
        }
      }
    }
  }
}

void multiply(const MPI_Comm& comm, int L, std::vector<std::pair<int, int> >& lattice, const std::vector<double>& v, std::vector<double>& w, std::vector<double>& buffer) {
  multiply(comm, L, lattice, &v[0], &w[0], &buffer[0]);
}


void fill_diagonal(const MPI_Comm& comm, int L, std::vector<std::pair<int, int> >& lattice, double* w) {
  int myrank, nproc;
  MPI_Status status;
  int ierr;

  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &myrank);

  int n = nproc;
  int p = -1;
  do {
    n /= 2;
    ++p;
  } while (n > 0);

  if (nproc != (1 << p)) {
    if ( myrank == 0 ) {
      std::cout << "This program can be run only for powers of 2" << std::endl;
    }
    MPI_Abort(comm, 1);
  }

  int N = 1 << (L-p);

  for (int l=0; l<lattice.size(); ++l) {
    int i = lattice[l].first;
    int j = lattice[l].second;
    if (i < (L-p)) {
      if (j < (L-p)) {
        int m1 = 1 << i;
        int m2 = 1 << j;
        int m3 = m1 + m2;
        for (int k=0; k<N; ++k) {
          if (((k & m3) == m1) || ((k & m3) == m2)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
            w[k] -= 0.25;
          } else {
            w[k] += 0.25;
          }
        }
      }
    }
  }
}

void fill_diagonal(const MPI_Comm& comm, int L, std::vector<std::pair<int, int> >& lattice, std::vector<double>& w) {
  fill_diagonal(comm, L, lattice, &w[0]);
}

template <typename MATRIX_MAJOR>
void generate(int L, std::vector<std::pair<int, int> >& lattice, rokko::distributed_matrix<MATRIX_MAJOR>& mat) {
  mat.set_zeros();
  for (int l=0; l<lattice.size(); ++l) {
    int i = lattice[l].first;
    int j = lattice[l].second;
    int m1 = 1 << i;
    int m2 = 1 << j;
    int m3 = m1 + m2;
    for (int k=0; k<mat.get_n_global(); ++k) {
      if (mat.is_gindex_mycol(k)) {
        int local_k = mat.translate_g2l_col(k);
        if (((k & m3) == m1) || ((k & m3) == m2)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
          if (mat.is_gindex_myrow(k^m3)) {
            mat.update_local(mat.translate_g2l_row(k^m3), local_k, 0.5);
          }
          if (mat.is_gindex_myrow(k)) {
            mat.update_local(local_k, local_k, -0.25);
          }
        } else if (mat.is_gindex_myrow(k)) {
          mat.update_local(local_k, local_k, 0.25);
        }
      }
    }
  }
}

} // namespace heisenberg_hamiltonian

} // namespace rokko

#endif // ROKKO_UTILITY_HEISENBERG_HAMILTONIANPARALLEL_HPP
