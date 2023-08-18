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
#include <rokko/original/crs.hpp>
#include <rokko/distributed_crs_matrix.hpp>

#include <iostream>

namespace rokko {
namespace original {

class distributed_crs_matrix : public rokko::detail::ps_crs_base {
public:
  distributed_crs_matrix() = default;
  ~distributed_crs_matrix() = default;

  explicit distributed_crs_matrix(rokko::mapping_1d const& map, int num_entries_per_row)
    : distributed_crs_matrix(cast_map(map), num_entries_per_row) {}

  explicit distributed_crs_matrix(rokko::skel::mapping_1d const& map, int num_entries_per_row)
    : distributed_crs_matrix(std::make_shared<const rokko::skel::mapping_1d>(map), num_entries_per_row) {}

  explicit distributed_crs_matrix(std::shared_ptr<const rokko::skel::mapping_1d> map, int num_entries_per_row)
    : map_(map), matrix_(std::make_shared<crs>(map, num_entries_per_row)) {}

  static std::shared_ptr<const rokko::skel::mapping_1d> cast_map(rokko::mapping_1d const& map) {
    if (map.get_solver_name() != "original") {
      throw std::invalid_argument("original's distributed_crs_matrix() : " + map.get_solver_name() + "'s mapping_1d is given.");
    }
    return std::static_pointer_cast<const rokko::skel::mapping_1d>(map.get_ptr());
  }

  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) override {
    matrix_->insert(row, cols, values);
  }
  void insert(int row, int col_size, int const*const cols, double const*const values) override {
    matrix_->insert(row, col_size, cols, values);
  }
  void complete() override {
    matrix_->complete();
  }
  void extract(int row, std::vector<int>& cols, std::vector<double>& values) const override {
    matrix_->extract(row, cols, values);
  }
  auto get_matrix() const {
    return matrix_;
  }
  int get_dim() const override {
    return map_->get_dim();
  }
  int get_num_local_rows() const override {
    return map_->get_num_local_rows();
  }
  int start_row() const override {
    return map_->start_row();
  }
  int end_row() const override {
    return map_->end_row();
  }
  int get_nnz() const override {
    return matrix_->get_nnz();
  }
  void print() const override {
    matrix_->print();
  }

  void output_matrix_market(std::ostream& os = std::cout) const override {
    const auto& comm = get_map().get_mpi_comm();
    constexpr int root_proc = 0;
    std::vector<int> cols;
    std::vector<double> values;

    const auto nnz = get_nnz();
    if (comm.get_myrank() == root_proc) {
      os << "%%MatrixMarket matrix coordinate real general" << std::endl
         << get_dim() << " " << get_dim() << " " << nnz << std::endl;
    }
    comm.barrier();
    for (int global_row=0; global_row<get_dim(); ++global_row) {
      if ((global_row >= start_row()) && (global_row < end_row())) {
        extract(global_row, cols, values);
        for (size_t i=0; i<cols.size(); ++i) {
          os << global_row + 1 << " " << cols[i] + 1 << " " << values[i] << std::endl;
        }
      }
      comm.barrier();
    }
  }

  std::shared_ptr<const rokko::skel::mapping_1d> get_map_ptr() const { return map_; }

  const rokko::skel::mapping_1d& get_map() const override { return *map_; }

private:
  std::shared_ptr<const rokko::skel::mapping_1d> map_;
  std::shared_ptr<crs> matrix_;
};

} // namespace original
} // namespace rokko
