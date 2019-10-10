/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef PYROKKO_GENERATE_MATRIX_HPP
#define PYROKKO_GENERATE_MATRIX_HPP

#include <rokko/eigen3.hpp>

namespace rokko {

template<int MATRIX_MAJOR>
void pyrokko_generate(Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>> mat, std::function<double(int, int)> const& func) {
  for(int i = 0; i < mat.rows(); ++i) {
    for(int j = 0; j < mat.cols(); ++j) {
      mat(i, j) = func(i, j);
    }
  }
}

} // namespace rokko

#endif // PYROKKO_GENERATE_MATRIX_HPP
