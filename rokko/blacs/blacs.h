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

#ifndef ROKKO_BLACS_H
#define ROKKO_BLACS_H

#include <rokko/mangling.h>

#ifdef __cplusplus
extern "C" {
#endif
  
#define BLACS_barrier LAPACK_GLOBAL(blacs_barrier,BLACS_BARRIER)
void BLACS_barrier(const int* ictxt, const char* score);

#define BLACS_exit LAPACK_GLOBAL(blacs_exit,BLACS_EXIT)
void BLACS_exit(const int* conti);
  
#define BLACS_get LAPACK_GLOBAL(blacs_get,BLACS_GET)
void BLACS_get(const int* context, const int* request, int* value);

#define BLACS_gridexit LAPACK_GLOBAL(blacs_gridexit,BLACS_GRIDEXIT)
void BLACS_gridexit(int* ictxt);

#define BLACS_gridinfo LAPACK_GLOBAL(blacs_gridinfo,BLACS_GRIDINFO)
void BLACS_gridinfo(const int* ictxt, int* nprow, int* npcol, int* myrow, int* mycol);

#define BLACS_gridinit LAPACK_GLOBAL(blacs_gridinit,BLACS_GRIDINIT)
void BLACS_gridinit(int* ictxt, const char* order, const int* nprow, const int* npcol);

#define BLACS_pinfo LAPACK_GLOBAL(blacs_pinfo,BLACS_PINFO)
void BLACS_pinfo(int* mypnum, int* nprocs);

#define BLACS_descinit LAPACK_GLOBAL(descinit,DESCINIT)
void BLACS_descinit(int* desc, const int* m, const int* n, const int* mb, const int* nb,
                    const int* irsrc, const int* icsrc, const int* ixtxt, const int* lld,
                    int* info);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // ROKKO_BLACS_H
