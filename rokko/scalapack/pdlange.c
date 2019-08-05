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

#include <rokko/cscalapack.h>
#include <stdlib.h>
#include <rokko/Cblacs.h>
#include <rokko/scalapack/scalapack_interface.h>

double cscalapack_pdlange(char norm, int m, int n, const double* A, const int* descA) {
  int ia = 0;
  int ja = 0;
  int lwork = 0;
  if (norm == '1' || norm == 'O' || norm == 'o' || norm == 'I' || norm == 'i') {
    int context = descA[1];
    int nprow, npcol, myrow, mycol;
    Cblacs_gridinfo(context, &nprow, &npcol, &myrow, &mycol);
    if (norm == '1' || norm == 'O' || norm == 'o') {
      int nb_a = descA[5];
      int csrc_a = descA[7];
      int icoffa = ja % nb_a;
      int iacol = cscalapack_indxg2p(ja, nb_a, mycol, csrc_a, npcol);
      lwork = cscalapack_numroc(n + icoffa, nb_a, mycol, iacol, npcol);
    } else if (norm == 'I' || norm == 'i') {
      int mb_a = descA[4];
      int rsrc_a = descA[6];
      int ifoffa = ia % mb_a;
      int iarow = cscalapack_indxg2p(ia, mb_a, myrow, rsrc_a, nprow);
      lwork = cscalapack_numroc(m + ifoffa, mb_a, myrow, iarow, nprow);
    }
  }
  double* work = NULL;
  if (lwork) work = (double*)malloc(sizeof(double) * lwork);
  double result = cscalapack_pdlange_work(norm, m, n, A, descA, work);
  if (lwork) free(work);
  return result;
}
