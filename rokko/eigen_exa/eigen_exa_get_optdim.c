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

#include <rokko/eigen_exa/eigen_exa.h>
#include <rokko/eigen_exa/eigen_exa_wrap.h>

int ROKKO_EIGEN_EXA_get_optdim(int n_min, int n_unroll, int delta_L1, int delta_L2) {
  int n_opt;
  EIGEN_EXA_cstab_get_optdim(&n_min, &n_unroll, &delta_L1, &delta_L2, &n_opt);
  return n_opt;
}
