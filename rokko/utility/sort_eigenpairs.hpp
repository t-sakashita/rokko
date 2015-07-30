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

#ifndef ROKKO_UTILITY_SORT_EIGENPAIRS_HPP
#define ROKKO_UTILITY_SORT_EIGENPAIRS_HPP

#include <vector>
#include <utility>
#include <algorithm>
#include <rokko/localized_vector.hpp>
#include <rokko/localized_matrix.hpp>

namespace {

template<typename T>
struct less {
  bool operator()(std::pair<T, int> const& lhs, std::pair<T, int> const& rhs) const {
    return lhs < rhs;
  }
};

template<typename T>
struct more {
  bool operator()(std::pair<T, int> const& lhs, std::pair<T, int> const& rhs) const {
  return lhs > rhs;
  }
};

} // end of unnamed namespace

namespace rokko {

// sort eivanvalue and eigenvectors 


template<typename T, typename MATRIX_MAJOR>  
void sort_eigenpairs(const localized_vector<T>& eigval,
                     const localized_matrix<T, MATRIX_MAJOR>& eigvec,
                     localized_vector<T>& eigval_sorted,
                     localized_matrix<T, MATRIX_MAJOR>& eigvec_sorted,
                     bool ascending = true) {
  int dim = eigval.size();
  std::vector<std::pair<T, std::size_t> > entries;
  entries.reserve(dim);
  for (std::size_t i = 0; i < dim; ++i) {
    entries.push_back(std::make_pair(eigval(i), i));
  }
  if (ascending)
    std::sort(entries.begin(), entries.end(), ::less<T>());
  else
    std::sort(entries.begin(), entries.end(), ::more<T>());

  if (eigval_sorted.size() != dim) eigval_sorted.resize(dim);
  if (eigvec_sorted.is_col_major()) {
    if (eigvec_sorted.cols() != dim) eigvec_sorted.resize(Eigen::NoChange, dim);
    for (int i = 0; i < dim; ++i) {
      eigval_sorted(i) = entries[i].first;
      eigvec_sorted.col(i) = eigvec.col(entries[i].second);
    }
  } else {
    if (eigvec_sorted.rows() != dim) eigvec_sorted.resize(dim, Eigen::NoChange);
    for (int i = 0; i < dim; ++i) {
      eigval_sorted(i) = entries[i].first;
      eigvec_sorted.row(i) = eigvec.row(entries[i].second);
    }
  }
}

} // namespace rokko

#endif // ROKKO_UTILITY_SORT_EIGENPAIRS_HPP
