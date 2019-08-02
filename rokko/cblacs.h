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

#ifndef ROKKO_CBLACS_H
#define ROKKO_CBLACS_H

#ifdef __cplusplus
extern "C" {
#endif
  
void cblacs_barrier(int ictxt, char score);

int cblacs_descinit(int* desc, int m, int n, int mb, int nb, int irsrc, int icsrc,
                    int ictxt, int lld);

void cblacs_exit(int conti);
  
void cblacs_free_blacs_system_handle(int *ISysCxt);

int cblacs_get(int context, int request);

void cblacs_gridexit(int* ictxt);

void cblacs_gridinfo(int ictxt, int* nprow, int* npcol, int* myrow, int* mycol);

void cblacs_gridinit(int* ictxt, char order, int nprow, int npcol);

void cblacs_pinfo(int* mypnum, int* nprocs);

int cblacs_sys2blacs_handle(int *SysCtxt);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // ROKKO_CBLACS_H
