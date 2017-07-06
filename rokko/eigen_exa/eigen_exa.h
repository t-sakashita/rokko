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

#ifndef ROKKO_EIGEN_EXA_H
#define ROKKO_EIGEN_EXA_H

#include <mpi.h>
#include <rokko/mangling.h>

#ifdef __cplusplus
extern "C" {
#endif

#define EIGEN_EXA_init_wrap ROKKO_GLOBAL(eigen_init_wrap,EIGEN_INIT_WRAP)
void EIGEN_EXA_init_wrap(const MPI_Fint* comm, const char* grid_major);

#define EIGEN_EXA_free_wrap ROKKO_GLOBAL(eigen_free_wrap,EIGEN_FREE_WRAP)
void EIGEN_EXA_free_wrap(const int* flag);

#define EIGEN_EXA_cstab_get_optdim ROKKO_GLOBAL(cstab_get_optdim,CSTAB_GET_OPTDIM)
void EIGEN_EXA_cstab_get_optdim(const int* n_min, const int* n_unroll, const int* delta_L1,
                                const int* delta_L2, int* n_opt);

#define EIGEN_EXA_get_matdims_wrap ROKKO_GLOBAL(eigen_get_matdims_wrap,EIGEN_GET_MATDIMS_WRAP)
void EIGEN_EXA_get_matdims_wrap(const int* nprow, const int* npcol, const int* n, int* nx, int* ny);

#define EIGEN_EXA_eigen_s ROKKO_GLOBAL(eigen_s,EIGEN_S)
void EIGEN_EXA_eigen_s(const int* n, const int* nvec, double* a, const int* lda,
                       double* w, double* z, const int* ldz,
                       const int* m_forward, const int* m_backword, const char* mode);

#define EIGEN_EXA_eigen_sx ROKKO_GLOBAL(eigen_sx,EIGEN_SX)
void EIGEN_EXA_eigen_sx(const int*, const int*, double*, const int*, double*, double*, const int*,
	       const int*, const int*, const char*);

#ifdef __cplusplus
}
#endif

#endif // ROKKO_EIGEN_EXA_H
