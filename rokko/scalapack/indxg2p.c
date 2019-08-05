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
#include <rokko/scalapack/scalapack_interface.h>

int cscalapack_indxg2p(int indxglob, int nb, int iproc, int isrcproc, int nprocs) {
  int indxglob_f = indxglob + 1;
  int iproc_f = iproc + 1;
  int isrcproc_f = isrcproc + 1;
  return SCALAPACK_indxg2p(&indxglob_f, &nb, &iproc_f, &isrcproc_f, &nprocs);
}
