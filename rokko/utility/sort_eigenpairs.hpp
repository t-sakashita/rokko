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
#include <rokko/eigen3.hpp>

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

// sort eigenvalues and eigenvectors
template<typename T, int SIZE, int ROWS, int COLS, int MATRIX_MAJOR>
void sort_eigenpairs(const Eigen::Vector<T, SIZE>& eigval,
                     const Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>& eigvec,
                     Eigen::Vector<T, SIZE>& eigval_sorted,
                     Eigen::Matrix<T,ROWS,COLS,MATRIX_MAJOR>& eigvec_sorted,
                     bool ascending = true) {
  int dim = eigval.size();
  std::vector<std::pair<T, std::size_t>> entries;
  entries.reserve(dim);
  for (int i = 0; i < dim; ++i) {
    entries.emplace_back(std::make_pair(eigval(i), i));
  }
  if (ascending)
    std::sort(entries.begin(), entries.end(), ::less<T>());
  else
    std::sort(entries.begin(), entries.end(), ::more<T>());

  if (eigval_sorted.size() != dim) eigval_sorted.resize(dim);
  if (MATRIX_MAJOR == Eigen::ColMajor) {
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
