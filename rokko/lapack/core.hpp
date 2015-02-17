/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2014 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_LAPACK_CORE_HPP
#define ROKKO_LAPACK_CORE_HPP

#include <rokko/lapack/diagonalize.hpp>
// #include <rokko/lapack/diagonalize_dsyevd.hpp>
// #include <rokko/lapack/diagonalize_dsyevx.hpp>

namespace rokko {
namespace lapack {

struct dsyev {};
// struct dsyevd {};
// struct dsyevx {};

template<typename ROUTINE>
class solver {
public:
  void initialize(int& argc, char**& argv) {}
  void finalize() {}
  template<typename MATRIX_MAJOR, typename VEC>
  void diagonalize(localized_matrix<MATRIX_MAJOR>& mat, VEC& eigvals,
                   localized_matrix<MATRIX_MAJOR>& eigvecs, timer& timer);
};

template<>
template<typename MATRIX_MAJOR, typename VEC>
void solver<rokko::lapack::dsyev>::diagonalize(localized_matrix<MATRIX_MAJOR>& mat,
  VEC& eigvals, localized_matrix<MATRIX_MAJOR>& eigvecs, timer& timer) {
  rokko::lapack::diagonalize(mat, eigvals, eigvecs, timer);
}

// template<>
// template<typename MATRIX_MAJOR, typename VEC>
// void solver<rokko::lapack::dsyevd>::diagonalize(localized_matrix<MATRIX_MAJOR>& mat,
//   VEC& eigvals, localized_matrix<MATRIX_MAJOR>& eigvecs, timer& timer) {
//   rokko::lapack::diagonalize_d(mat, eigvals, eigvecs, timer);
// }

// template<>
// template<typename MATRIX_MAJOR, typename VEC>
// void solver<rokko::lapack::dsyevx>::diagonalize(localized_matrix<MATRIX_MAJOR>& mat,
//   VEC& eigvals, localized_matrix<MATRIX_MAJOR>& eigvecs, timer& timer) {
//   rokko::lapack::diagonalize_x(mat, eigvals, eigvecs, timer);
// }

} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_CORE_HPP
