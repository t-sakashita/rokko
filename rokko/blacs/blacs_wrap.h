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

#ifndef ROKKO_BLACS_WRAP_H
#define ROKKO_BLACS_WRAP_H

#ifdef __cplusplus
extern "C" {
#endif
  
void ROKKO_blacs_barrier(int ictxt, char score);

void ROKKO_blacs_exit(int conti);

int ROKKO_blacs_get(int context, int request);

void ROKKO_blacs_gridexit(int* ictxt);

void ROKKO_blacs_gridinfo(int ictxt, int* nprow, int* npcol, int* myrow, int* mycol);

void ROKKO_blacs_gridinit(int* ictxt, char order, int nprow, int npcol);

void ROKKO_blacs_pinfo(int* mypnum, int* nprocs);

int ROKKO_descinit(int* desc, int m, int n, int mb, int nb, int irsrc, int icsrc,
                   int ictxt, int lld);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // ROKKO_BLACS_WRAP_H
