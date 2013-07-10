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

#include "rokko/blacs.h"

void ROKKO_blacs_gridinfo(int ictxt, int* nprow, int* npcol, int* myrow, int* mycol) {
  blacs_gridinfo_(&ictxt, nprow, npcol, myrow, mycol);
}
