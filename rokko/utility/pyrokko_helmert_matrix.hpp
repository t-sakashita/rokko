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

  template <typename T, int MAJOR>
  static void generate_eigen(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>> mat_in) {
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR> mat(mat_in.rows(), mat_in.cols());
    new (&mat) Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>>(mat_in);
    helmert_matrix::generate(mat);
    new (&mat) Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>();
  }

  template <typename T, int MAJOR>
  static void generate_for_given_eigenvalues_eigen(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>> mat_in, Eigen::Ref<Eigen::Vector<T>> diag_in) {
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR> mat(mat_in.rows(), mat_in.cols());
    Eigen::Vector<T> diag;

    new (&mat) Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>>(mat_in);
    new (&diag) Eigen::Ref<Eigen::Vector<T>>(diag_in);

    helmert_matrix::generate_for_given_eigenvalues(mat, diag);

    new (&mat) Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MAJOR>();
    new (&diag) Eigen::Vector<T>();
  }

  template <typename T, typename MATRIX_MAJOR>
  static void generate(wrap_distributed_matrix<T,MATRIX_MAJOR>& mat) {
    helmert_matrix::generate(mat);
  }

  template <typename T, typename MATRIX_MAJOR>
  static void generate_for_given_eigenvalues(wrap_distributed_matrix<T,MATRIX_MAJOR>& mat, Eigen::Ref<Eigen::Vector<T>> diag_in) {
    Eigen::Vector<T> diag;
    new (&diag) Eigen::Ref<Eigen::Vector<T>>(diag_in);

    helmert_matrix::generate_for_given_eigenvalues(mat, diag);

    new (&diag) Eigen::Vector<T>();
  }
};

} // namespace rokko

#endif // PYROKKO_HELMERT_MATRIX_HPP
