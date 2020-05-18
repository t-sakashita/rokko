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

#ifndef PYROKKO_MINIJ_MATRIX_HPP
#define PYROKKO_MINIJ_MATRIX_HPP

#include <rokko/utility/minij_matrix.hpp>

#include <rokko/pyrokko_distributed_matrix.hpp>

namespace rokko {

class wrap_minij_matrix {
public:
  static void generate_row(Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>> mat_in) {
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> mat(mat_in.rows(), mat_in.cols());
    new (&mat) Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>(mat_in);
    minij_matrix::generate(mat);
    new (&mat) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>();
  }

  static void generate_col(Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> mat_in) {
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> mat(mat_in.rows(), mat_in.cols());
    new (&mat) Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(mat_in);
    minij_matrix::generate(mat);
    new (&mat) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>();
  }

  template <typename T, typename MATRIX_MAJOR>
  static void generate(wrap_distributed_matrix<T,MATRIX_MAJOR>& mat) {
    minij_matrix::generate(mat);
  }
};

} // namespace rokko

#endif // PYROKKO_MINIJ_MATRIX_HPP
