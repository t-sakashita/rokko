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

class distributed_crs_matrix : public rokko::detail::distributed_crs_matrix_base {
public:
  explicit distributed_crs_matrix(int row_dim, int col_dim) : dim_(row_dim) {
    initialize(row_dim, col_dim);
    set_sizes();
  }
  explicit distributed_crs_matrix(int row_dim, int col_dim, int num_entries_per_row) : dim_(row_dim) {
    initialize(row_dim, col_dim, num_entries_per_row);
    set_sizes();
  }
  void initialize(int row_dim, int col_dim) {
    map_ = new mapping_1d(row_dim);
    matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map_->get_epetra_map(), col_dim));
  }
  void initialize(int row_dim, int col_dim, int num_entries_per_row) {
    map_ = new mapping_1d(row_dim);
    matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map_->get_epetra_map(), num_entries_per_row));
  }
  void set_sizes() {
    num_local_rows_ = map_->get_epetra_map().NumMyElements();
    start_row_ = map_->get_epetra_map().MinMyGID();
    end_row_ = map_->get_epetra_map().MaxMyGID() + 1; // to follow C++ convention
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
  Teuchos::RCP<Epetra_CrsMatrix> get_matrix() {
    return matrix_;
  }
  int get_dim() const {
    return dim_;
  }
  int num_local_rows() const {
    return num_local_rows_;
  }
  int start_row() const {
    return start_row_;
  }
  int end_row() const {
    return end_row_;
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

private:
  mapping_1d* map_;
  Teuchos::RCP<Epetra_CrsMatrix> matrix_;
  int dim_;
  int num_local_rows_, start_row_, end_row_;
};

} // namespace anasazi
} // namespace rokko

#endif // ROKKO_ANASAZI_DISTRIBUTED_CRS_MATRIX_HPP
