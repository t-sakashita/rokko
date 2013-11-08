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

#ifndef ROKKO_UTILITY_XYZ_HAMILTONIAN_MPI_HPP
#define ROKKO_UTILITY_XYZ_HAMILTONIAN_MPI_HPP

#include "mpi.h"
#include <vector>
#include <boost/tuple/tuple.hpp>

#include <iostream>

#include <rokko/localized_matrix.hpp>

namespace rokko {

namespace xyz_hamiltonian {

void multiply(const MPI_Comm& comm, int L, const std::vector<std::pair<int, int> >& lattice, const std::vector<boost::tuple<double, double, double> >& coupling, const double* v, double* w, double* buffer) {
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
    int jx = coupling[l].get<0>();
    int jy = coupling[l].get<1>();
    int jz = coupling[l].get<2>();

    double diag_plus = jz / 4.0;
    double diag_minus = - jz / 4.0;
    double offdiag_plus = (jx + jy) / 4.0;
    double offdiag_minus = (jx - jy) / 4.0;

    if (i < (L-p)) {
      if (j < (L-p)) {
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
      } else {
        int m = 1 << (j-(L-p));
        MPI_Sendrecv(const_cast<double*>(&v[0]), N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     &buffer[0], N, MPI_DOUBLE, 
                     myrank ^ m, 0,
                     comm, &status);
        int m1 = 1 << i;
        if ((myrank & m) == m) { 
          for (int k=0; k<N; ++k) {
            if ((k & m1) == m1) {
              w[k] += diag_plus * v[k] + offdiag_minus * buffer[k^m1];
            } else {
              w[k] += diag_minus * v[k] + offdiag_plus * buffer[k^m1];
            }
          }
        } else {
          for (int k=0; k<N; ++k) {
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
        MPI_Sendrecv(const_cast<double*>(&v[0]), N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     &buffer[0], N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     comm, &status);
        int m1 = 1 << j;
        if ((myrank & m) == m) {
          for (int k=0; k<N; ++k) {
            if ((k & m1) == m1) {
              w[k] += diag_plus * v[k] + offdiag_minus * buffer[k^m1];
            } else {
              w[k] += diag_minus * v[k] + offdiag_plus * buffer[k^m1];
            }
          }
        } else {
          for (int k=0; k<N; ++k) {
            if ((k & m1) == m1) {
              w[k] += diag_minus * v[k] + offdiag_plus * buffer[k^m1];
            } else {
              w[k] += diag_plus * v[k] + offdiag_minus * buffer[k^m1];
            }
          }
        }
      } else {
        int m = (1 << (i-(L-p))) + (1 << (j-(L-p)));
        MPI_Sendrecv(const_cast<double*>(&v[0]), N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     &buffer[0], N, MPI_DOUBLE,
                     myrank ^ m, 0,
                     comm, &status);
        if (((myrank & m) != m) && ((myrank & m) != 0)) {
          for (int k=0; k<N; ++k) {
            w[k] += diag_minus * v[k] + offdiag_plus * buffer[k];
          }
        } else {
          for (int k=0; k<N; ++k) {
            w[k] += diag_plus * v[k] + offdiag_minus * buffer[k];
          }
        }
      }
    }
  }
}

void multiply(const MPI_Comm& comm, int L, std::vector<std::pair<int, int> >& lattice, const std::vector<boost::tuple<double, double, double> >& coupling, const std::vector<double>& v, std::vector<double>& w, std::vector<double>& buffer) {
  multiply(comm, L, lattice, coupling, &v[0], &w[0], &buffer[0]);
}

/*
void fill_diagonal(const MPI_Comm& comm, int L, std::vector<std::pair<int, int> >& lattice, double* w, double* buffer) {
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
          }
          else if ((k & m3) == m2) {
            w[k] += 0.5 * v[k^m3] - 0.25 * v[k];
          }
          else {
            w[k] += 0.25 * v[k];
          }
        }
      } else {
        int m = 1 << (j-(L-p));
        MPI_Sendrecv(&v[0], N, MPI_DOUBLE,
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
        MPI_Sendrecv(&v[0], N, MPI_DOUBLE,
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
          MPI_Sendrecv(&v[0], N, MPI_DOUBLE,
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

      else {
        w[k] += 0.25;
        //      cout << "else" << endl;
      }
    }
  }
}

void fill_diagonal(const MPI_Comm& comm, int L, std::vector<std::pair<int, int> >& lattice, std::vector<double>& w, std::vector<double>& buffer) {
  fill_diagonal(comm, L, lattice, &w[0], &buffer[0]);
}
*/

} // namespace xyz_hamiltonian

} // namespace rokko

#endif // ROKKO_UTILITY_XYZ_HAMILTONIAN_MPI_HPP
