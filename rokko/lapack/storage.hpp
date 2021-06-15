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

#include <rokko/lapack/null_matrix.hpp>
#include <vector>

namespace rokko {
namespace lapack {

template <typename T>
T* storage(T* ptr) { return ptr; }

template <typename T>
const T* storage(const T* ptr) { return ptr; }

template <typename T>
T* storage(std::vector<T>& vec) { return vec.data(); }

} // end namespace lapack
} // end namespace rokko
