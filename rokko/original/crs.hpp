/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2021 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#pragma once

#include <rokko/skel/mapping_1d.hpp>

#include <algorithm>
#include <iostream>

namespace rokko {
namespace original {

class crs : public skel::mapping_1d {
public:
  crs() = default;
  ~crs() = default;

  explicit crs(std::shared_ptr<const rokko::skel::mapping_1d> map, int num_entries_per_row) : skel::mapping_1d(*map) {
    rows.reserve(num_entries_per_row);
    cols_vector.reserve(num_entries_per_row);
    values_vector.reserve(num_entries_per_row);
  }

  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) {
    rows.emplace_back(row);
    cols_vector.push_back(cols);
    values_vector.push_back(values);
  }

  void insert(int row, int col_size, int const*const cols, double const*const values) {
    rows.emplace_back(row);
    std::vector<int> tmp(cols, cols+col_size);
    cols_vector.emplace_back(tmp);
    std::vector<double> tmp2(values, values+col_size);
    values_vector.emplace_back(tmp2);
  }

  void complete() {
    // to be implemented
  }

  int find_row_index(int row) const {
    auto itr = std::find(rows.begin(), rows.end(), row);
    if (itr == rows.end())
      throw std::out_of_range("find_row_index failed");
    return std::distance(rows.begin(), itr);
  }

  void extract(int row, std::vector<int>& cols, std::vector<double>& values) const {
    int row_index = find_row_index(row);
    cols = cols_vector[row_index];
    values = values_vector[row_index];
  }

  int get_nnz() const {
    return 0; // to be implemented
  }

  void print() const {
    for(size_t i=0; i<rows.size(); ++i) {
      int row = rows[i];
      std::cout << "row=" << row << std::endl;
      for (size_t j=0; j<cols_vector[i].size(); ++j) {
        std::cout << "  col=" << cols_vector[i][j] << "  value=" << values_vector[i][j] << std::endl;
      }
    }
  }

  void multiply(const double *const x, double *const y) const {
    for(size_t i=0; i<rows.size(); ++i) {
      int row = rows[i];
      for (size_t j=0; j<cols_vector[i].size(); ++j) {
        int col = cols_vector[i][j];
        y[row] += values_vector[i][j] * x[col];
      }
    }
  }

private:
  std::vector<int> rows;
  std::vector<std::vector<int>> cols_vector;
  std::vector<std::vector<double>> values_vector;
};

} // namespace original
} // namespace rokko
