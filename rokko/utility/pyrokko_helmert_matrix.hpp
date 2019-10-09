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

#ifndef PYROKKO_HELMERT_MATRIX_HPP
#define PYROKKO_HELMERT_MATRIX_HPP

#include <rokko/utility/helmert_matrix.hpp>

#include <rokko/pyrokko_distributed_matrix.hpp>

namespace rokko {

class wrap_helmert_matrix {
public:

  template <int MAJOR>
  static void generate_eigen(Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>> mat_in) {
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR> mat(mat_in.rows(), mat_in.cols());
    new (&mat) Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>>(mat_in);
    helmert_matrix::generate(mat);
    new (&mat) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>();
  }

  template <int MAJOR>
  static void generate_for_given_eigenvalues_eigen(Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>> mat_in, Eigen::Ref<Eigen::VectorXd> diag_in) {
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR> mat(mat_in.rows(), mat_in.cols());
    Eigen::VectorXd diag;

    new (&mat) Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>>(mat_in);
    new (&diag) Eigen::Ref<Eigen::VectorXd>(diag_in);

    helmert_matrix::generate_for_given_eigenvalues(mat, diag);

    new (&mat) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>();
    new (&diag) Eigen::VectorXd();
  }

  static void generate(wrap_distributed_matrix& mat) {
    if (mat.is_major_col())
      helmert_matrix::generate(mat.col_ver());
    else
      helmert_matrix::generate(mat.row_ver());
  }

  static void generate_for_given_eigenvalues(wrap_distributed_matrix& mat, Eigen::Ref<Eigen::VectorXd> diag_in) {
    Eigen::VectorXd diag;
    new (&diag) Eigen::Ref<Eigen::VectorXd>(diag_in);

    if (mat.is_major_col())
      helmert_matrix::generate_for_given_eigenvalues(mat.col_ver(), diag);
    else
      helmert_matrix::generate_for_given_eigenvalues(mat.row_ver(), diag);

    new (&diag) Eigen::VectorXd();
  }
};

} // namespace rokko

#endif // PYROKKO_HELMERT_MATRIX_HPP
