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

#pragma once

#include <rokko/utility/frank_matrix.hpp>

#include <rokko/pyrokko_distributed_matrix.hpp>

namespace rokko {

class wrap_frank_matrix {
public:
  template <typename T, int MAJOR>
  static void generate(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>> mat_in) {
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR> mat;
    new (&mat) Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>>(mat_in);
    frank_matrix::generate(mat);
    new (&mat) Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>();
  }

  template<typename T, typename MATRIX_MAJOR>
  static void generate(wrap_distributed_matrix<T,MATRIX_MAJOR>& mat) {
    frank_matrix::generate(mat);
  }
};

} // namespace rokko
