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

#ifndef ROKKO_EIGEN_EXA_H
#define ROKKO_EIGEN_EXA_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

void eigen_init_wrap_(const MPI_Fint*, const char*);
void eigen_free_wrap_(const int*);
void cstab_get_optdim_(const int*, const int*, const int*, const int*, int*);
void eigen_get_matdims_wrap_(const int*, const int*, const int*, int*, int*);
void eigen_s_(const int*, const int*, double*, const int*, double*, double*, const int*,
	       const int*, const int*, const char*);
void eigen_sx_(const int*, const int*, double*, const int*, double*, double*, const int*,
	       const int*, const int*, const char*);

void ROKKO_eigen_init(MPI_Comm, char);
void ROKKO_eigen_free(int);
void ROKKO_cstab_get_optdim(int, int, int, int, int*);
void ROKKO_eigen_get_matdims(int, int, int, int*, int*);
void ROKKO_eigen_s(int, int, double*, int, double*, double*, int, int, int, char);
void ROKKO_eigen_sx(int, int, double*, int, double*, double*, int, int, int, char);

#ifdef __cplusplus
}
#endif

#endif // ROKKO_EIGEN_EXA_H
