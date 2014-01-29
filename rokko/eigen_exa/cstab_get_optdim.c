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

#include <rokko/eigen_exa.h>

void ROKKO_cstab_get_optdim(int n_min, int n_unroll, int delta_L1, int delta_L2, int *n_opt) {
  cstab_get_optdim_(&n_min, &n_unroll, &delta_L1, &delta_L2, n_opt);
}
