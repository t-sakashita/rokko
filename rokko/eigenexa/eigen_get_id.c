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

void ceigenexa_get_id(int* id, int* x_id, int* y_id) {
  EIGENEXA_get_id_wrap(id, x_id, y_id);
  *id -= 1;
  *x_id -= 1;
  *y_id -= 1;
}
