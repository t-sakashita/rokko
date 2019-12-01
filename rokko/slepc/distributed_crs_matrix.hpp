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

#ifndef ROKKO_SLEPC_DISTRIBUTED_CRS_MATRIX_HPP
#define ROKKO_SLEPC_DISTRIBUTED_CRS_MATRIX_HPP

#include <rokko/distributed_crs_matrix.hpp>
#include <rokko/distributed_mfree.hpp>
#include <rokko/utility/timer.hpp>

#include <slepceps.h>
#include <petscblaslapack.h>

namespace rokko {
namespace slepc {

struct comp{
  bool operator()(const int& a, const int& b) const {
    return v[a] < v[b];
  }
  comp(const int *p) : v(p) {}

private:
  const int *v;
};

class distributed_crs_matrix : public rokko::detail::distributed_crs_matrix_base {
public:
  distributed_crs_matrix() {}
  ~distributed_crs_matrix() {}

  explicit distributed_crs_matrix(int row_dim, int col_dim) {
    initialize(row_dim, col_dim);
  }
  explicit distributed_crs_matrix(int row_dim, int col_dim, int num_entries_per_row) {
    initialize(row_dim, col_dim, num_entries_per_row);
  }
  #undef __FUNCT__
  #define __FUNCT__ "distributed_crs_matrix/initialize"
  void initialize(int row_dim, int col_dim) {
    ierr = MatCreate(PETSC_COMM_WORLD, &matrix_);  //CHKERRQ(ierr);
    ierr = MatSetSizes(matrix_, PETSC_DECIDE, PETSC_DECIDE, row_dim, col_dim);  //CHKERRQ(ierr);
    ierr = MatSetFromOptions(matrix_);  //CHKERRQ(ierr);
    ierr = MatSetUp(matrix_);  //CHKERRQ(ierr);
    dim_ = row_dim;
    ierr = MatGetOwnershipRange(matrix_, &start_row_, &end_row_); //CHKERRQ(ierr);
    num_local_rows_ = end_row_ - start_row_;// + 1;
  }
  #undef __FUNCT__
  #define __FUNCT__ "distributed_crs_matrix/initialize with num_entries_per_row"
  void initialize(int row_dim, int col_dim, int num_entries_per_row) {
    ierr = MatCreate(PETSC_COMM_WORLD, &matrix_);  //CHKERRQ(ierr);
    ierr = MatSetSizes(matrix_, PETSC_DECIDE, PETSC_DECIDE, row_dim, col_dim);  //CHKERRQ(ierr);
    ierr = MatSetFromOptions(matrix_);  //CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(matrix_, num_entries_per_row, NULL);
    ierr = MatMPIAIJSetPreallocation(matrix_, num_entries_per_row, NULL, num_entries_per_row, NULL);
    dim_ = row_dim;
    ierr = MatGetOwnershipRange(matrix_, &start_row_, &end_row_); //CHKERRQ(ierr);
    num_local_rows_ = end_row_ - start_row_;// + 1;
  }
  #undef __FUNCT__
  #define __FUNCT__ "distributed_crs_matrix/insert"
  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) {
    ierr = MatSetValues(matrix_, 1, &row, cols.size(), &cols[0], &values[0], INSERT_VALUES);  //CHKERRQ(ierr);
  }
  void insert(int row, int col_size, int* cols, double* const values) {
    ierr = MatSetValues(matrix_, 1, &row, col_size, cols, values, INSERT_VALUES);  //CHKERRQ(ierr);
  }
  #undef __FUNCT__
  #define __FUNCT__ "distributed_crs_matrix/complete"
  void complete() {
    ierr = MatAssemblyBegin(matrix_, MAT_FINAL_ASSEMBLY);  //CHKERRQ(ierr);
    ierr = MatAssemblyEnd(matrix_, MAT_FINAL_ASSEMBLY);  //CHKERRQ(ierr);
  }
  Mat* get_matrix() {
    return &matrix_;
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
    MatInfo info;
    MatGetInfo(matrix_, MAT_GLOBAL_SUM, &info);
    return (int)info.nz_used;
  }
  void print() const {
    MatView(matrix_, PETSC_VIEWER_STDOUT_WORLD);
  }
  void output_matrix_market() const {
    std::vector<int> idx;

    PetscInt num_cols;
    const PetscInt * cols;
    const PetscScalar * values;
    int myrank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
    int nnz = get_nnz();
    if (myrank == 0) {
      std::cout << "%%MatrixMarket matrix coordinate real general" << std::endl;
      std::cout << get_dim() << " " << get_dim() << " " << nnz << std::endl;
    }
    MPI_Barrier(PETSC_COMM_WORLD);
    for (int global_row=0; global_row<get_dim(); ++global_row) {
      if ((global_row >= start_row()) && (global_row < end_row())) {
        MatGetRow(matrix_, global_row, &num_cols, &cols, &values);
        idx.resize(num_cols);
        for (int i=0; i<num_cols; ++i) idx[i] = i;
        std::sort(idx.begin(), idx.end(), comp(cols));
        for (int i=0; i<num_cols; ++i) {
          std::cout << global_row + 1 << " " << cols[idx[i]] + 1 << " " << values[idx[i]] << std::endl;
        }
        MatRestoreRow(matrix_, global_row, &num_cols, &cols, &values);
      }
      MPI_Barrier(PETSC_COMM_WORLD);
    }
  }

private:
  int dim_;
  int num_local_rows_;
  int start_row_, end_row_;
  PetscErrorCode ierr;
  Mat matrix_;
};

} // namespace slepc
} // namespace rokko

#endif // ROKKO_SLEPC_DISTRIBUTED_CRS_MATRIX_HPP
