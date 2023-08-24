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

#pragma once

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

void ceigenexa_init(void);
void ceigenexa_init1(MPI_Comm comm);
void ceigenexa_init2(MPI_Comm comm, char grid_major);

void ceigenexa_free(void);
void ceigenexa_free1(int flag);

// int ceigenexa_get_optdim(int n_min, int n_unroll, int delta_L1, int delta_L2);

void ceigenexa_get_matdims(int nprow, int npcol, int n, int* nx, int* ny);

void ceigenexa_get_procs(int* procs, int* x_procs, int* y_procs);

void ceigenexa_get_id(int* id, int* x_id, int* y_id);

int ceigenexa_loop_start(int istart, int nnod, int inod);

int ceigenexa_loop_end(int iend, int nnod, int inod);

int ceigenexa_translate_l2g(int ictr, int nnod, int inod);

int ceigenexa_translate_g2l(int ictr, int nnod, int inod);

void ceigenexa_eigen_s(int n, int nvec, double* a, int lda, double* w, double* z, int ldz,
                       int m_forward, int m_backword, char mode);

void ceigenexa_eigen_sx(int n, int nvec, double* a, int lda, double* w, double* z, int ldz,
                        int m_forward, int m_backword, char mode);

#ifdef __cplusplus
}
#endif
