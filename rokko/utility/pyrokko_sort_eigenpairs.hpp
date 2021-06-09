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

#include <rokko/eigen3.hpp>

#include <rokko/utility/sort_eigenpairs.hpp>

namespace rokko {

template <int MAJOR>
void pyrokko_sort_eigenpairs(Eigen::Ref<Eigen::VectorXd> eigval_in,
                             Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>> eigvec_in,
                             Eigen::Ref<Eigen::VectorXd> eigval_sorted_in,
                             Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>> eigvec_sorted_in,
                             bool ascending = true) {
  Eigen::VectorXd eigval, eigval_sorted;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR> eigvec(eigvec_in.rows(),eigvec_in.cols()), eigvec_sorted(eigvec_in.rows(),eigvec_in.cols());

  new (&eigval) Eigen::Ref<Eigen::VectorXd>(eigval_in);
  new (&eigvec) Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>>(eigvec_in);
  new (&eigval_sorted) Eigen::Ref<Eigen::VectorXd>(eigval_sorted_in);
  new (&eigvec_sorted) Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>>(eigvec_sorted_in);

  sort_eigenpairs(eigval, eigvec, eigval_sorted, eigvec_sorted, ascending);

  new (&eigval) Eigen::VectorXd();
  new (&eigvec) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>();
  new (&eigval_sorted) Eigen::VectorXd();
  new (&eigvec_sorted) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,MAJOR>();
}

} // namespace rokko
