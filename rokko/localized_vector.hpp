/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_LOCALIZED_VECTOR_HPP
#define ROKKO_LOCALIZED_VECTOR_HPP

#include <rokko/eigen3.hpp>
#include <iostream>

namespace rokko {

template<typename T, int ROWS = Eigen::Dynamic, int MAJOR = Eigen::ColMajor>
using Vector = Eigen::Matrix<T, ROWS, 1, MAJOR>;

template<typename T>
using RefColVec = Eigen::Ref<Vector<T>>;

} // namespace rokko

#endif // ROKKO_LOCALIZED_VECTOR_HPP
