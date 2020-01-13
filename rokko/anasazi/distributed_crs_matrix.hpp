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

namespace rokko {
namespace anasazi {

struct comp{
  bool operator()(int a, int b) const {
    return matrix_->GCID(v[a]) < matrix_->GCID(v[b]);
  }
  comp(const int *p, Teuchos::RCP<Epetra_CrsMatrix> const& matrix) : v(p), matrix_(matrix) {}

private:
  const int *const v;
  Teuchos::RCP<Epetra_CrsMatrix> matrix_;
};

class distributed_crs_matrix : public rokko::detail::ps_crs_base {
public:
  distributed_crs_matrix() {}
  ~distributed_crs_matrix() {}

  explicit distributed_crs_matrix(rokko::mapping_1d const& map, int num_entries_per_row) {
    initialize(map, num_entries_per_row);
  }

  explicit distributed_crs_matrix(rokko::anasazi::mapping_1d const& map, int num_entries_per_row) {
    initialize(map, num_entries_per_row);
  }

  explicit distributed_crs_matrix(int row_dim, int col_dim) {
    initialize(row_dim, col_dim);
  }
  explicit distributed_crs_matrix(int row_dim, int col_dim, int num_entries_per_row) {
    initialize(row_dim, col_dim, num_entries_per_row);
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
  void initialize(int row_dim, int col_dim) {
    map_ = new rokko::anasazi::mapping_1d(row_dim);
    matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map_->get_epetra_map(), col_dim));
  }
  void initialize(int row_dim, int col_dim, int num_entries_per_row) {
    map_ = new rokko::anasazi::mapping_1d(row_dim);
    matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map_->get_epetra_map(), num_entries_per_row));
  }

  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) {
    matrix_->InsertGlobalValues(row, cols.size(), values.data(), cols.data());
  }
  void insert(int row, int col_size, int* cols, double* const values) {
    matrix_->InsertGlobalValues(row, col_size, values, cols);
  }
  void complete() {
    matrix_->FillComplete();
    matrix_->SetTracebackMode(1);
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
  void output_matrix_market() const {
    const int MaxNumIndices = matrix_->MaxNumEntries();
    std::vector<int> idx(MaxNumIndices);
    int num_cols;
    int* cols;
    double* values;
    const int NumMyElements = map_->get_epetra_map().NumMyElements();
    std::vector<int> MyGlobalElements(NumMyElements);
    map_->get_epetra_map().MyGlobalElements(MyGlobalElements.data());
    int local_row = 0;
    if (map_->get_epetra_comm().MyPID() == 0) {
      std::cout << "%%MatrixMarket matrix coordinate real general" << std::endl;
      std::cout << get_dim() << " " << get_dim() << " " << get_nnz() << std::endl;
    }
    map_->get_epetra_comm().Barrier();
    for (int global_row=0; global_row<get_dim(); ++global_row) {
      if (local_row < NumMyElements) {
        if (global_row == MyGlobalElements[local_row]) {
          matrix_->ExtractMyRowView(local_row, num_cols, values, cols);
          idx.resize(num_cols);
          std::iota(idx.begin(), idx.end(), 0);
          std::sort(idx.begin(), idx.end(), comp(cols, matrix_));
          for (int i=0; i<num_cols; ++i) {
            std::cout << global_row + 1 << " " << matrix_->GCID(cols[idx[i]]) + 1 << " " << values[idx[i]] << std::endl;
          }
          ++local_row;
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
