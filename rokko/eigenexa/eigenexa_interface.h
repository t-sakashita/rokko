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

#ifndef ROKKO_EIGENEXA_EIGENEXA_INTERFACE_H
#define ROKKO_EIGENEXA_EIGENEXA_INTERFACE_H

#include <mpi.h>
#include <rokko/mangling.h>

#ifdef __cplusplus
extern "C" {
#endif

#define EIGENEXA_init_wrap0 ROKKO_GLOBAL(eigen_init_wrap0,EIGEN_INIT_WRAP0)
void EIGENEXA_init_wrap0();

#define EIGENEXA_init_wrap1 ROKKO_GLOBAL(eigen_init_wrap1,EIGEN_INIT_WRAP1)
void EIGENEXA_init_wrap1(const MPI_Fint* comm);

#define EIGENEXA_init_wrap2 ROKKO_GLOBAL(eigen_init_wrap2,EIGEN_INIT_WRAP2)
void EIGENEXA_init_wrap2(const MPI_Fint* comm, const char* grid_major);

#define EIGENEXA_free_wrap0 ROKKO_GLOBAL(eigen_free_wrap0,EIGEN_FREE_WRAP0)
void EIGENEXA_free_wrap0();

#define EIGENEXA_free_wrap1 ROKKO_GLOBAL(eigen_free_wrap1,EIGEN_FREE_WRAP1)
void EIGENEXA_free_wrap1(const int* flag);

#define EIGENEXA_get_matdims_wrap ROKKO_GLOBAL(eigen_get_matdims_wrap,EIGEN_GET_MATDIMS_WRAP)
void EIGENEXA_get_matdims_wrap(const int* nprow, const int* npcol, const int* n, int* nx, int* ny);

#define EIGENEXA_get_procs_wrap ROKKO_GLOBAL(eigen_get_procs_wrap,EIGEN_GET_PROCS_WRAP)
void EIGENEXA_get_procs_wrap(int* procs, int* x_procs, int* y_procs);

#define EIGENEXA_get_id_wrap ROKKO_GLOBAL(eigen_get_id_wrap,EIGEN_GET_ID_WRAP)
void EIGENEXA_get_id_wrap(int* id, int* x_id, int* y_id);

#define EIGENEXA_loop_start_wrap ROKKO_GLOBAL(eigen_loop_start_wrap,EIGEN_LOOP_START_WRAP)
int EIGENEXA_loop_start_wrap(const int* istart, const int* nnod, const int* inod);

#define EIGENEXA_loop_end_wrap ROKKO_GLOBAL(eigen_loop_end_wrap,EIGEN_LOOP_END_WRAP)
int EIGENEXA_loop_end_wrap(const int* iend, const int* nnod, const int* inod);

#define EIGENEXA_translate_l2g_wrap ROKKO_GLOBAL(eigen_translate_l2g_wrap,EIGEN_TRANSLATE_L2G_WRAP)
int EIGENEXA_translate_l2g_wrap(const int* ictr, const int* nnod, const int* inod);

#define EIGENEXA_translate_g2l_wrap ROKKO_GLOBAL(eigen_translate_g2l_wrap,EIGEN_TRANSLATE_G2L_WRAP)
int EIGENEXA_translate_g2l_wrap(const int* ictr, const int* nnod, const int* inod);

#define EIGENEXA_eigen_s ROKKO_GLOBAL(eigen_s0,EIGEN_S0)
void EIGENEXA_eigen_s(const int* n, const int* nvec, double* a, const int* lda,
                      double* w, double* z, const int* ldz,
                      const int* m_forward, const int* m_backword, const char* mode);
  
#define EIGENEXA_eigen_sx ROKKO_GLOBAL(eigen_sx,EIGEN_SX)
void EIGENEXA_eigen_sx(const int* n, const int* nvec, double* a, const int* lda,
                       double* w, double* z, const int* ldz,
                       const int* m_forward, const int* m_backword, const char* mode);

#ifdef __cplusplus
}
#endif

#endif // ROKKO_EIGENEXA_EIGENEXA_INTERFACE_H
