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

#ifndef ROKKO_BLACS_H
#define ROKKO_BLACS_H

#ifdef __cplusplus
extern "C" {
#endif
  
void blacs_barrier_(const int* ictxt, const char* score);
void blacs_exit_(const int* conti);
void blacs_get_(const int* context, const int* request, int* value);
void blacs_gridexit_(int* ictxt);
void blacs_gridinfo_(const int* ictxt, int* nprow, int* npcol, int* myrow, int* mycol);
void blacs_gridinit_(const int* ictxt, const char* order, const int* nprow, const int* npcol);
void blacs_pinfo_(int* mypnum, int* nprocs);
void descinit_(int* desc, const int* m, const int* n, const int* mb, const int* nb,
               const int* irsrc, const int* icsrc, const int* ixtxt, const int* lld, int* info);

void ROKKO_blacs_barrier(int ictxt, char score);
void ROKKO_blacs_exit(int conti);
void ROKKO_blacs_get(int context, int request, int* value);
void ROKKO_blacs_gridexit(int* ictxt);
void ROKKO_blacs_gridinfo(int ictxt, int* nprow, int* npcol, int* myrow, int* mycol);
void ROKKO_blacs_gridinit(int ictxt, char order, int nprow, int npcol);
void ROKKO_blacs_pinfo(int* mypnum, int* nprocs);
void ROKKO_descinit(int* desc, int m, int n, int mb, int nb, int irsrc, int icsrc,
                    int ixtxt, int lld, int* info);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // ROKKO_BLACS_H
