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

void ROKKO_eigen_exa_free(int flag) {
  EIGEN_EXA_free_wrap(&flag);
}
