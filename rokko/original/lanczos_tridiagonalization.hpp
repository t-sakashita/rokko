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

#ifndef ROKKO_ORIGINAL_LANCZOS_TRIDIAGONALIZATION_HPP
#define ROKKO_ORIGINAL_LANCZOS_TRIDIAGONALIZATION_HPP

#include <rokko/mpi/reduction.hpp>

namespace rokko {

namespace original {

class lanczos {
public:
  template<typename VEC, typename MAT>
  static void tridiagonalization(const rokko::distributed_mfree& mat, VEC& alpha, VEC& beta, MAT& u) {
    constexpr double eps_lanczos = 1e-8;

    const int dim = mat.get_dim();
    MPI_Comm comm = mat.get_comm();
    Eigen::VectorXd tmp(dim), y(dim);

    u.col(0).setRandom();
    u.col(0) /= sqrt(rokko::all_reduce(u.col(0).squaredNorm(), comm));

    y.setZero();
    mat.multiply(u.col(0).data(), y.data());
    alpha(0) = rokko::all_reduce(y.dot(u.col(0)), comm);
    tmp = y - alpha(0) * u.col(0);
    for (int i=1; i<dim; ++i) {
      beta(i-1) = sqrt(rokko::all_reduce(tmp.squaredNorm(), comm));
      if (beta(i-1) < eps_lanczos)  throw std::runtime_error("beta is too small");
      u.col(i) = tmp / beta(i-1);

      y.setZero();
      mat.multiply(u.col(i).data(), y.data());
      alpha(i) = rokko::all_reduce(y.dot(u.col(i)), comm);
      tmp = y - alpha(i) * u.col(i) - beta(i-1) * u.col(i-1);
    }
  }
};

} // namespace original

} // namespace rokko

#endif // ROKKO_ORIGINAL_LANCZOS_TRIDIAGONALIZATION_HPP
