/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <rokko/distributed_matrix.hpp>
#include <El.hpp>
#include <rokko/utility/timer.hpp>

namespace rokko {
namespace elemental {

template<typename T>
struct complex_type_traits {
  using type = T;
};

template<typename T>
struct complex_type_traits<std::complex<T>> {
  using type = El::Complex<T>;
};


template<typename T>
using complex_t = typename complex_type_traits<T>::type;


El::UpperOrLower get_matrix_part(parameters const& params) {
  std::string matrix_part;
  if (params.defined("uplow"))
    matrix_part = params.get_string("uplow");
  if (params.defined("matrix_part"))
    matrix_part = params.get_string("matrix_part");
  if (!matrix_part.empty()) {
    const char matrix_part_letter = matrix_part[0];
    if ((matrix_part_letter == 'u') || (matrix_part_letter == 'U'))
      return El::UPPER;
    if ((matrix_part_letter == 'l') || (matrix_part_letter == 'L'))
      return El::LOWER;
  }
  return El::LOWER;  // default
}

template <typename T>
El::HermitianEigSubset<T> get_subset(parameters const& params) {
  El::HermitianEigSubset<T> subset;
  const bool is_lower_value = get_key(params, "lower_value", subset.lowerBound);
  const bool is_lower_index = get_key(params, "lower_index", subset.lowerIndex);
  const bool is_upper_value = get_key(params, "upper_value", subset.upperBound);
  const bool is_upper_index = get_key(params, "upper_index", subset.upperIndex);
  if (is_lower_value && is_upper_value) {
    subset.rangeSubset = true;
  } else if (is_lower_index && is_upper_index) {
    subset.indexSubset = true;
  } else if (!(is_lower_index && is_lower_value && is_upper_index && is_upper_value)) {  
  } else {
    throw std::invalid_argument("elemental::get_subset() : sepcify either of a pair of upper_value and lower_value or a pair of upper_index and lower_index");
  }
  return subset;
}

template <typename T>
El::HermitianEigCtrl<T> get_ctrl(parameters const& params) {
  El::HermitianEigCtrl<T> ctrl;
  if (params.defined("useSDC")) {
    ctrl.useSDC = params.get<bool>("useSDC");
  }
  if (params.defined("cutoff")) {
    ctrl.sdcCtrl.cutoff = params.get<int>("cutoff");
  }
  if (params.defined("maxInnerIts")) {
    ctrl.sdcCtrl.maxInnerIts = params.get<int>("maxInnerIts");
  }
  if (params.defined("maxOuterIts")) {
    ctrl.sdcCtrl.maxOuterIts = params.get<int>("maxOuterIts");
  }
  if (params.defined("tol")) {
    ctrl.sdcCtrl.tol = params.get<real_t<T>>("tol");
  }
  if (params.defined("spreadFactor")) {
    ctrl.sdcCtrl.spreadFactor = params.get<real_t<T>>("spreadFactor");
  }
  if (params.defined("progress")) {
    ctrl.sdcCtrl.progress = params.get<bool>("progress");
  }
  //if (params.defined("verbose")) {
  //  ctrl.sdcCtrl.progress = true;
  //}
  return ctrl;
}

El::SortType get_sort(parameters const& params) {
  El::SortType elem_sort = El::ASCENDING;
  if (params.defined("sort")) {
    const auto str_sort = params.get<std::string>("sort");
    if (str_sort == "ascending") elem_sort =  El::ASCENDING;
    else if (str_sort == "descending") elem_sort =  El::DESCENDING;
    else if (str_sort == "unsorted") elem_sort =  El::UNSORTED;
  }
  return elem_sort;
}

// eigenvalues / eigenvectors
template<typename T, typename VEC>
parameters diagonalize(distributed_matrix<T, rokko::matrix_col_major>& mat,
		       VEC& eigvals, distributed_matrix<T, rokko::matrix_col_major>& eigvecs,
		       parameters const& params) {
  if((mat.get_mb() != 1) || (mat.get_nb() != 1))
    throw std::invalid_argument("elemental::diagonalize() : elemental supports only 1x1 block size.");
  parameters params_out;
  const MPI_Comm comm = mat.get_grid().get_comm();
  const enum El::GridOrder elemental_grid_order = (mat.get_grid().is_row_major()) ?
    El::ROW_MAJOR : El::COLUMN_MAJOR;
  const El::Grid elem_grid(comm, mat.get_grid().get_nprow(), elemental_grid_order);
  El::DistMatrix<complex_t<T>> elem_mat;
  elem_mat.Attach(mat.get_m_global(), mat.get_n_global(), elem_grid, 0, 0,
                  static_cast<complex_t<T>*>(mat.get_array_pointer()), mat.get_lld());
  El::DistMatrix<complex_t<T>> elem_eigvecs(0, 0, elem_grid);
  El::DistMatrix<real_t<T>, El::VR, El::STAR> elem_w(elem_grid);

  const El::UpperOrLower elem_uplow = get_matrix_part(params);
  El::HermitianEigCtrl<complex_t<T>> ctrl = get_ctrl<complex_t<T>>(params);
  ctrl.tridiagEigCtrl.subset = get_subset<real_t<T>>(params);
  ctrl.tridiagEigCtrl.sort = get_sort(params);

  El::HermitianEig(elem_uplow, elem_mat, elem_w, elem_eigvecs, ctrl);

  for (int i = 0; i < elem_w.Height(); ++i) eigvals(i) = elem_w.Get(i, 0);
  T* result_mat = elem_eigvecs.Buffer();
  for(int local_j=0; local_j<elem_w.LocalHeight(); ++local_j) {
    for(int local_i=0; local_i<mat.get_m_local(); ++local_i) {
      eigvecs.set_local(local_i, local_j, result_mat[local_j * mat.get_lld() + local_i]);
    }
  }
  return params_out;
}

template<typename T, typename VEC>
parameters diagonalize(distributed_matrix<T, rokko::matrix_row_major>& mat,
		       VEC& eigvals, distributed_matrix<T, rokko::matrix_row_major>& eigvecs,
		       parameters const& params) {
  throw std::invalid_argument("elemental::diagonalize() : elemental doesn't support matrix_row_major.  Use it with matrix_col_major.");
}

// only eigenvalues
template<typename T, typename VEC>
parameters diagonalize(distributed_matrix<T, rokko::matrix_col_major>& mat,
		       VEC& eigvals,
		       parameters const& params) {
  if((mat.get_mb() != 1) || (mat.get_nb() != 1))
    throw std::invalid_argument("elemental::diagonalize() : elemental supports only 1x1 block size.");
  parameters params_out;
  const MPI_Comm comm = mat.get_grid().get_comm();
  const enum El::GridOrder elemental_grid_order = (mat.get_grid().is_row_major()) ?
    El::ROW_MAJOR : El::COLUMN_MAJOR;
  const El::Grid elem_grid(comm, mat.get_grid().get_nprow(), elemental_grid_order);
  El::DistMatrix<complex_t<T>> elem_mat;
  elem_mat.Attach(mat.get_m_global(), mat.get_n_global(), elem_grid, 0, 0,
                  static_cast<complex_t<T>*>(mat.get_array_pointer()), mat.get_lld());
  El::DistMatrix<complex_t<T>> elem_eigvecs(0, 0, elem_grid);
  El::DistMatrix<real_t<T>, El::VR, El::STAR> elem_w(elem_grid);

  const El::UpperOrLower elem_uplow = get_matrix_part(params);
  El::HermitianEigCtrl<complex_t<T>> ctrl = get_ctrl<complex_t<T>>(params);
  ctrl.tridiagEigCtrl.subset = get_subset<real_t<T>>(params);
  ctrl.tridiagEigCtrl.sort = get_sort(params);

  El::HermitianEig(elem_uplow, elem_mat, elem_w, ctrl);
  for (int i = 0; i < elem_w.Height(); ++i) eigvals(i) = elem_w.Get(i, 0);

  return params_out;
}

template<typename T, typename VEC>
parameters diagonalize(distributed_matrix<T, rokko::matrix_row_major>& mat,
		       VEC& eigvals,
		       parameters const& params) {
  throw std::invalid_argument("elemental::diagonalize() : elemental doesn't support matrix_row_major.  Use it with matrix_col_major.");
}

} // namespace elemental
} // namespace rokko
