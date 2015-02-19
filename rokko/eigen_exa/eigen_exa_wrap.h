/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_EIGEN_EXA_WRAP_H
#define ROKKO_EIGEN_EXA_WRAP_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

void ROKKO_eigen_exa_init(MPI_Comm comm, char grid_major);

void ROKKO_eigen_exa_free(int flag);

int ROKKO_eigen_exa_get_optdim(int n_min, int n_unroll, int delta_L1, int delta_L2);

void ROKKO_eigen_exa_get_matdims(int nprow, int npcol, int n, int* nx, int* ny);
void ROKKO_eigen_exa_s(int n, int nvec, double* a, int lda, double* w, double* z, int ldz,
                       int m_forward, int m_backword, char mode);
void ROKKO_eigen_exa_sx(int, int, double*, int, double*, double*, int, int, int, char);

#ifdef __cplusplus
}
#endif

#endif // ROKKO_EIGEN_EXA_WRAP_H
