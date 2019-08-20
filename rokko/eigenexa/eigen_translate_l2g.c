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

int ceigenexa_translate_l2g(int ictr, int nnod, int inod) {
  ictr += 1;
  inod += 1;
  return EIGENEXA_translate_l2g_wrap(&ictr, &nnod, &inod) - 1;
}
