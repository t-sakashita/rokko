/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_ANASAZI_DISTRIBUTED_CRS_MATRIX_HPP
#define ROKKO_ANASAZI_DISTRIBUTED_CRS_MATRIX_HPP

#include <rokko/anasazi/mapping_1d.hpp>
#include <rokko/distributed_crs_matrix.hpp>

#include <Epetra_CrsMatrix.h>
#include <AnasaziEpetraAdapter.hpp>
#include <Teuchos_RCPDecl.hpp>

#include <numeric>
#include <iostream>

namespace rokko {
namespace anasazi {

class distributed_crs_matrix : public rokko::detail::ps_crs_base {
public:
  distributed_crs_matrix() = default;
  ~distributed_crs_matrix() = default;

  explicit distributed_crs_matrix(rokko::mapping_1d const& map, int num_entries_per_row)
    : distributed_crs_matrix(cast_map(map), num_entries_per_row) {}

  explicit distributed_crs_matrix(rokko::anasazi::mapping_1d const& map, int num_entries_per_row)
    : distributed_crs_matrix(std::make_shared<const rokko::anasazi::mapping_1d>(map), num_entries_per_row) {}

  explicit distributed_crs_matrix(std::shared_ptr<const rokko::anasazi::mapping_1d> map, int num_entries_per_row)
    : map_(map), matrix_(Teuchos::rcp(new Epetra_CrsMatrix(Copy, map_->get_epetra_map(), num_entries_per_row))) {}

  static std::shared_ptr<const rokko::anasazi::mapping_1d> cast_map(rokko::mapping_1d const& map) {
    if (map.get_solver_name() != "anasazi") {
      throw std::invalid_argument("Anasazi's distributed_crs_matrix() : " + map.get_solver_name() + "'s mapping_1d is given.");
    }
    return std::static_pointer_cast<const rokko::anasazi::mapping_1d>(map.get_ptr());
  }

  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) override {
    matrix_->InsertGlobalValues(row, cols.size(), values.data(), cols.data());
  }
  void insert(int row, int col_size, int const*const cols, double const*const values) override {
    matrix_->InsertGlobalValues(row, col_size, values, cols);
  }
  void complete() override {
    matrix_->FillComplete();
    matrix_->SetTracebackMode(1);
  }
  void extract(int row, std::vector<int>& cols, std::vector<double>& values) const override {
    int num_cols;
    int* cols_tmp;
    double* values_tmp;
    matrix_->ExtractMyRowView(matrix_->LRID(row), num_cols, values_tmp, cols_tmp);
    cols.clear();
    cols.reserve(num_cols);
    values.clear();
    values.reserve(num_cols);
    std::vector<int> idx(num_cols);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [this,cols_tmp](auto i, auto j) { return matrix_->GCID(cols_tmp[i]) < matrix_->GCID(cols_tmp[j]); });
    for (int i=0; i<num_cols; ++i) {
      cols.emplace_back(matrix_->GCID(cols_tmp[idx[i]]));
      values.emplace_back(values_tmp[idx[i]]);
    }
  }
  Teuchos::RCP<Epetra_CrsMatrix> get_matrix() const {
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
    return matrix_->NumGlobalNonzeros();
  }
  void print() const override {
    std::cout << *matrix_ << std::endl;
  }
  void output_matrix_market(std::ostream& os = std::cout) const override {
    constexpr int root_proc = 0;
    std::vector<int> cols;
    std::vector<double> values;
    if (map_->get_epetra_comm().MyPID() == root_proc) {
      os << "%%MatrixMarket matrix coordinate real general" << std::endl
         << get_dim() << " " << get_dim() << " " << get_nnz() << std::endl;
    }
    map_->get_epetra_comm().Barrier();
    for (int global_row=0; global_row<get_dim(); ++global_row) {
      if ((global_row >= start_row()) && (global_row < end_row())) {
        extract(global_row, cols, values);
        for (int i=0; i<cols.size(); ++i) {
          os << global_row + 1 << " " << cols[i] + 1 << " " << values[i] << std::endl;
        }
      }
      map_->get_epetra_comm().Barrier();
    }
  }

  std::shared_ptr<const rokko::anasazi::mapping_1d> get_map_ptr() const { return map_; }

  const rokko::anasazi::mapping_1d& get_map() const override { return *map_; }

private:
  std::shared_ptr<const rokko::anasazi::mapping_1d> map_;
  Teuchos::RCP<Epetra_CrsMatrix> matrix_;
};

} // namespace anasazi
} // namespace rokko

#endif // ROKKO_ANASAZI_DISTRIBUTED_CRS_MATRIX_HPP
