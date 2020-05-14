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

#ifndef ROKKO_LAPACK_DIAGONALIZE_SYEVX_HPP
#define ROKKO_LAPACK_DIAGONALIZE_SYEVX_HPP

#include <rokko/parameters.hpp>
#include <rokko/eigen3.hpp>
#include <rokko/utility/timer.hpp>
#include <rokko/lapack.hpp>
#include <rokko/lapack/diagonalize_get_parameters.hpp>
#include <rokko/lapack/syevx.hpp>

namespace rokko {
namespace lapack {

// only eigenvalues
template<typename T, int MATRIX_MAJOR, typename VEC>
parameters diagonalize_syevx(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, VEC& eigvals,
			      parameters const& params) {
  parameters params_out;

  auto abstol = params.defined("abstol") ? params.get<norm_t<T>>("abstol") : 2*LAPACKE_dlamch('S');
  params_out.set("abstol", abstol);

  lapack_int il, iu;
  norm_t<T> vl, vu;
  const char range = get_eigenvalues_range(params, vl, vu, il, iu);
  const char uplow = get_matrix_part(params);

  lapack_int m;  // output: found eigenvalues
  const int dim = mat.outerSize();
  std::vector<lapack_int> ifail(dim);

  int info = syevx(range, uplow, mat,
                   vl, vu, il, iu, abstol,
                   m, eigvals, ifail);

  if (info) {
    std::cerr << "error at syevx function. info=" << info << std::endl;
    if (info < 0) {
      std::cerr << "This means that ";
      std::cerr << "the " << abs(info) << "-th argument had an illegal value." << std::endl;
    }
  }
  params_out.set("info", info);
  params_out.set("m", m);
  params_out.set("ifail", ifail);
  
  if (params.get_bool("verbose")) {
    print_verbose("syevx", 'N', range, uplow, vl, vu, il, iu, params_out);
  }

  return params_out;
}


// eigenvalues / eigenvectors
template<typename T, int MATRIX_MAJOR, typename VEC>
parameters diagonalize_syevx(Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& mat, VEC& eigvals,
			      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,MATRIX_MAJOR>& eigvecs,
			      parameters const& params) {
  rokko::parameters params_out;

  auto abstol = params.defined("abstol") ? params.get<norm_t<T>>("abstol") : 2*LAPACKE_dlamch('S');
  params_out.set("abstol", abstol);

  lapack_int il, iu;
  norm_t<T> vl, vu;
  const char range = get_eigenvalues_range(params, vl, vu, il, iu);
  const char uplow = get_matrix_part(params);

  lapack_int m;  // output: found eigenvalues
  const int dim = mat.outerSize();
  std::vector<lapack_int> ifail(dim);

  int info = syevx(range, uplow, mat,
                   vl, vu, il, iu, abstol,
                   m, eigvals, eigvecs, ifail);

  if (info) {
    std::cerr << "Error at syevx function. info=" << info << std::endl;
    if (params.get_bool("verbose")) {
      std::cerr << "This means that ";
      if (info < 0) {
        std::cerr << "the " << abs(info) << "-th argument had an illegal value." << std::endl;
      } else {
        std::cerr << "This means that "	<< info << " eigenvectors failed to converge." << std::endl;
        std::cerr << "The indices of the eigenvectors that failed to converge:" << std::endl;
        for (std::size_t i = 0; i < ifail.size(); ++i) {
          if (ifail[i] == 0) break;
          std::cerr << ifail[i] << " ";
        }
        std::cerr << std::endl;
      }
    }
  }
  params_out.set("info", info);
  params_out.set("m", m);
  params_out.set("ifail", ifail);
  
  if (params.get_bool("verbose")) {
    print_verbose("syevx", 'V', range, uplow, vl, vu, il, iu, params_out);
  }

  return params_out;
}


} // namespace lapack
} // namespace rokko

#endif // ROKKO_LAPACK_DIAGONALIZE_SYEVX_HPP
