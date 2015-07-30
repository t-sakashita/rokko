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

#include <rokko/blacs/blacs.h>
#include <rokko/blacs/blacs_wrap.h>

void ROKKO_blacs_gridinit(int* ictxt, char order, int nprow, int npcol) {
  BLACS_gridinit(ictxt, &order, &nprow, &npcol);
}
