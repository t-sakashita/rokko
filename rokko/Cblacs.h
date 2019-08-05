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

#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#endif
  
MPI_Comm Cblacs2sys_handle(int BlacsCtxt);

void Cblacs_barrier(int ConTxt, char *scope);

void Cblacs_exit(int NotDone);

void Cblacs_get(int ConTxt, int what, int *val);

void Cblacs_gridexit(int ConTxt);

void Cblacs_gridinfo(int ConTxt, int *nprow, int *npcol, int *myrow, int *mycol);

void Cblacs_gridinit(int *ConTxt, char *order, int nprow, int npcol);

void Cblacs_gridmap(int *ConTxt, int *usermap, int ldup, int nprow0, int npcol0);

void Cblacs_pinfo(int *mypnum, int *nprocs);

void Cfree_blacs_system_handle(int ISysCtxt);

int Csys2blacs_handle(MPI_Comm SysCtxt);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // ROKKO_CBLACS_H
