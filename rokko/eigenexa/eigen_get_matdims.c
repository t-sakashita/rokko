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

#include <rokko/ceigenexa.h>
#include <rokko/eigenexa/eigenexa_interface.h>

void ceigenexa_get_matdims(int nprow, int npcol, int n, int *nx, int *ny) {
  EIGENEXA_get_matdims_wrap(&nprow, &npcol, &n, nx, ny);
}
