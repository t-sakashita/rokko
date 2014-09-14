/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014 by Synge Todo <wistaria@comp-phys.org>,
*                       Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>
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
    
class distributed_crs_matrix : public rokko::detail::distributed_crs_matrix_base {
public:
  distributed_crs_matrix() {}
  ~distributed_crs_matrix() {}

  distributed_crs_matrix(int row_dim, int col_dim) {
    initialize(row_dim, col_dim);
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
    num_local_rows_ = end_row_ - start_row_ + 1;
  }
  #undef __FUNCT__
  #define __FUNCT__ "distributed_crs_matrix/insert"
  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) {
    ierr = MatSetValues(matrix_, 1, &row, cols.size(), &cols[0], &values[0], ADD_VALUES);  //CHKERRQ(ierr);
  }
  void insert(int row, int col_size, int* cols, double* const values) {
    ierr = MatSetValues(matrix_, 1, &row, col_size, cols, values, ADD_VALUES);  //CHKERRQ(ierr);
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
  int get_dim() {
    return dim_;
  }
  int num_local_rows() {
    return num_local_rows_;
  }
  int start_row() {
    return start_row_;
  }
  int end_row() {
    return end_row_;
  }

private:
  int dim_;
  int num_local_rows_;
  int start_row_, end_row_;
  PetscErrorCode ierr;
public:
  Mat matrix_;
};

} // namespace slepc
} // namespace rokko

#endif // ROKKO_SLEPC_DISTRIBUTED_CRS_MATRIX_HPP
