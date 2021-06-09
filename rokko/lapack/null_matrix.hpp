/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

namespace rokko {
namespace lapack {

template<typename T>
class null_matrix {};

template<typename T>
int ld(null_matrix<T>) {
  return 1;
}

template<typename T>
T* storage(null_matrix<T>) {
  return NULL;
}

} // end namespace lapack
} // end namespace rokko
