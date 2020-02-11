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

  explicit distributed_crs_matrix(rokko::mapping_1d const& map, int num_entries_per_row) {
    initialize(map, num_entries_per_row);
  }

  explicit distributed_crs_matrix(rokko::anasazi::mapping_1d const& map, int num_entries_per_row) {
    initialize(map, num_entries_per_row);
  }

  explicit distributed_crs_matrix(std::array<int,2> const& dims) {
    initialize(dims);
  }
  explicit distributed_crs_matrix(std::array<int,2> const& dims, int num_entries_per_row) {
    initialize(dims, num_entries_per_row);
  }
  void initialize(rokko::mapping_1d const& map, int num_entries_per_row) {
    if (map.get_solver_name() != "anasazi") {
      throw std::invalid_argument("Anasazi's distributed_crs_matrix() : " + map.get_solver_name() + "'s mapping_1d is given.");
    }
    map_ = static_cast<const rokko::anasazi::mapping_1d*>(map.get_ptr()->get_impl());
    matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map_->get_epetra_map(), num_entries_per_row));
  }
  void initialize(rokko::anasazi::mapping_1d const& map, int num_entries_per_row) {
    map_ = &map;
    matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map_->get_epetra_map(), num_entries_per_row));
  }
  void initialize(std::array<int,2> const& dims) {
    map_ = new rokko::anasazi::mapping_1d(dims[0]);
    matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map_->get_epetra_map(), dims[1]));
  }
  void initialize(std::array<int,2> const& dims, int num_entries_per_row) { // ignoring dims[1]
    map_ = new rokko::anasazi::mapping_1d(dims[0]);
    matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map_->get_epetra_map(), num_entries_per_row));
  }

  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) {
    matrix_->InsertGlobalValues(row, cols.size(), values.data(), cols.data());
  }
  void insert(int row, int col_size, int const*const cols, double const*const values) {
    matrix_->InsertGlobalValues(row, col_size, values, cols);
  }
  void complete() {
    matrix_->FillComplete();
    matrix_->SetTracebackMode(1);
  }
  void extract(int row, std::vector<int>& cols, std::vector<double>& values) const {
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
  int get_dim() const {
    return map_->get_dim();
  }
  int num_local_rows() const {
    return map_->get_num_local_rows();
  }
  int start_row() const {
    return map_->start_row();
  }
  int end_row() const {
    return map_->end_row();
  }
  int get_nnz() const {
    return matrix_->NumGlobalNonzeros();
  }
  void print() const {
    std::cout << *matrix_ << std::endl;
  }
  void output_matrix_market(std::ostream& os = std::cout) const {
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

  const rokko::anasazi::mapping_1d* get_map() const { return map_; }

  ps_crs_base* get_impl() { return this; }
  const ps_crs_base* get_impl() const { return this; }

private:
  const rokko::anasazi::mapping_1d* map_;
  Teuchos::RCP<Epetra_CrsMatrix> matrix_;
};

} // namespace anasazi
} // namespace rokko

#endif // ROKKO_ANASAZI_DISTRIBUTED_CRS_MATRIX_HPP
