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

void multiply(const MPI_Comm comm, int L, int lattice_size, int lattice_first[], int lattice_second[], const double* v, double* w, double* buffer) {
  int myrank, nproc;
  MPI_Status status;

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
      printf("This program can be run only for powers of 2\n");
    }
    MPI_Abort(comm, 1);
  }

  int N = 1 << (L-p);

  for(size_t k=0; k<N; ++k) {
    w[k] = 0.;
  }

  for (size_t l = 0; l < lattice_size; ++l) {
    int i = lattice_first[l];
    int j = lattice_second[l];
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
