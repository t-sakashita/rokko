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

#ifndef ROKKO_ELPA_H
#define ROKKO_ELPA_H

#ifdef __cplusplus
extern "C" {
#endif

int elpa_solve_evp_real_1stage(int na, int nev, double *a, int lda, double *ev, double *q, int ldq, int nblk, int matrixCols,
			       int mpi_comm_rows, int mpi_comm_cols);

int elpa_solve_evp_real_2stage(int na, int nev, double *a, int lda, double *ev, double *q, int ldq, int nblk, int matrixCols,
			       int mpi_comm_rows, int mpi_comm_cols, int mpi_comm_all, int THIS_REAL_ELPA_KERNEL_API, int useQR);

int elpa_get_communicators(int mpi_comm_world, int my_prow, int my_pcol, int *mpi_comm_rows, int *mpi_comm_cols);

#ifdef __cplusplus
}
#endif

/* Fix me: we should include ELPA's header file of its install dir to use the following costants. */
/* However, to do so, we should set include path of ELPA's install dir, which seems diificult. */

#define ELPA2_REAL_KERNEL_GENERIC 1
#define ELPA2_REAL_KERNEL_GENERIC_SIMPLE 2
#define ELPA2_REAL_KERNEL_BGP 3
#define ELPA2_REAL_KERNEL_BGQ 4
#define ELPA2_REAL_KERNEL_SSE 5
#define ELPA2_REAL_KERNEL_AVX_BLOCK2 6
#define ELPA2_REAL_KERNEL_AVX_BLOCK4 7
#define ELPA2_REAL_KERNEL_AVX_BLOCK6 8

#define ELPA2_NUMBER_OF_REAL_KERNELS 8

#endif // ROKKO_ELPA_H
