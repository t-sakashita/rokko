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

#ifndef PYROKKO_COLLECTIVE_HPP
#define PYROKKO_COLLECTIVE_HPP

#include <rokko/pyrokko_distributed_matrix.hpp>

#include <rokko/collective.hpp>

namespace rokko {

void pyrokko_gather(wrap_distributed_matrix<matrix_row_major> const& from, Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>> to_in, int root) {
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> to;
  new (&to) Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>(to_in.data(), to_in.rows(), to_in.cols());

  gather(from, to, root);

  new (&to) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>();
}

void pyrokko_gather(wrap_distributed_matrix<matrix_col_major> const& from, Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> to_in, int root) {
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> to;
  new (&to) Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(to_in.data(), to_in.rows(), to_in.cols());

  gather(from, to, root);

  new (&to) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>();
}

void pyrokko_scatter(Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>> from_in, wrap_distributed_matrix<matrix_row_major>& to, int root) {
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> from;
  new (&from) Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>(from_in.data(), from_in.rows(), from_in.cols());

  scatter(from, to, root);

  new (&from) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>();
}

void pyrokko_scatter(Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> from_in, wrap_distributed_matrix<matrix_col_major>& to, int root) {
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> from;
  new (&from) Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>(from_in.data(), from_in.rows(), from_in.cols());
  
  scatter(from, to, root);

  new (&from) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>();
}

} // end namespace rokko

#endif // PYROKKO_COLLECTIVE_HPP
