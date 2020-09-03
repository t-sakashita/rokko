/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_LAPACK_DIAGONALIZE_TRIDIAGONALIZED_MATRIX_BISECTION_HPP
#define ROKKO_LAPACK_DIAGONALIZE_TRIDIAGONALIZED_MATRIX_BISECTION_HPP

#include <rokko/parameters.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_get_parameters.hpp>
#include <rokko/lapack/stebz.hpp>
#include <rokko/lapack/stein.hpp>

namespace rokko {
namespace lapack {

// eigenvalues / eigenvectors
template<typename VEC, typename T, int MATRIX_MAJOR>
parameters diagonalize_bisection(VEC& alpha, VEC& beta,
			      VEC& eigvals,
			      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			      parameters const& params) {
  // calculating eigenvalues
  int num_conv;
  const int dim = alpha.size();
  eigvals.resize(dim);
  constexpr T abstol = 0.;
  int il = dim, iu = dim;
  int nsplit;
  Eigen::VectorXi iblock(dim), isplit(dim);
  int info = rokko::lapack::stebz('E', il, iu, abstol, alpha, beta, num_conv, nsplit, eigvals, iblock, isplit);

  // calculating eigenvectors
  eigvecs.resize(dim, num_conv);
  Eigen::VectorXi ifailv(num_conv);
  info = rokko::lapack::stein(alpha, beta, num_conv, eigvals, eigvecs, iblock, isplit, ifailv);

  rokko::parameters params_out;
  params_out.set("num_conv", num_conv);
  params_out.set("info", info);
  if (info) {
    std::cerr << "error at stein function. info=" << info  << std::endl;
  }

  return params_out;
}

} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_TRIDIAGONALIZED_MATRIX_BISECTION_HPP
