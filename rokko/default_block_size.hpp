/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2023 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once


namespace rokko {

int get_default_block_size(int dim, int num_proc_axis) {
  int num_block = dim / num_proc_axis;
  if (num_block == 0)  num_block = 1;
  return num_block;
}

} // end namespace rokko
