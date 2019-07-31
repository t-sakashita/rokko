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

#include <rokko/cblacs.h>
#include <rokko/blacs/blacs.h>

int CBLACS_descinit(int* desc, int m, int n, int mb, int nb, int irsrc, int icsrc,
                    int ictxt, int lld) {
  int info;
  BLACS_descinit(desc, &m, &n, &mb, &nb, &irsrc, &icsrc, &ictxt, &lld, &info);
  return info;
}
