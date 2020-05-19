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

#ifndef PYROKKO_COLLECTIVE_HPP
#define PYROKKO_COLLECTIVE_HPP

#include <rokko/pyrokko_distributed_matrix.hpp>

#include <rokko/collective.hpp>

namespace rokko {

template<typename T, typename MATRIX_MAJOR>
void pyrokko_gather(wrap_distributed_matrix<double,MATRIX_MAJOR> const& from, Eigen::Ref<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<MATRIX_MAJOR>>> to_in, int root) {
  constexpr int major = rokko::eigen3_major<MATRIX_MAJOR>;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,major> to;
  new (&to) Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,major>>(to_in.data(), to_in.rows(), to_in.cols());

  gather(from, to, root);

  new (&to) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,major>();
}

template<typename T, typename MATRIX_MAJOR>
void pyrokko_scatter(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,rokko::eigen3_major<MATRIX_MAJOR>>> from_in, wrap_distributed_matrix<double,MATRIX_MAJOR>& to, int root) {
  constexpr int major = rokko::eigen3_major<MATRIX_MAJOR>;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,major> from;
  new (&from) Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,major>>(from_in.data(), from_in.rows(), from_in.cols());

  scatter(from, to, root);

  new (&from) Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,major>();
}

} // end namespace rokko

#endif // PYROKKO_COLLECTIVE_HPP
