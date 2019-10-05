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

#ifndef PYROKKO_LAPLACIAN_MATRIX_HPP
#define PYROKKO_LAPLACIAN_MATRIX_HPP

#include <rokko/utility/laplacian_matrix.hpp>

#include <rokko/pyrokko_distributed_matrix.hpp>

namespace rokko {

class wrap_laplacian_matrix {
public:
  static void generate_row(Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>> mat_in) {
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> mat;
    new (&mat) Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>(mat_in);
    laplacian_matrix::generate(mat);
    new (&mat) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>();
  }

  static void generate_col(Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> mat_in) {
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> mat;
    new (&mat) Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(mat_in);
    laplacian_matrix::generate(mat);
    new (&mat) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>();
  }

  static void generate(wrap_distributed_matrix& mat) {
    if (mat.is_major_col())
      laplacian_matrix::generate(mat.col_ver());
    else
      laplacian_matrix::generate(mat.row_ver());
  }
};

} // namespace rokko

#endif // PYROKKO_LAPLACIAN_MATRIX_HPP
