/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SLEPC_DISTRIBUTED_CRS_MATRIX_HPP
#define ROKKO_SLEPC_DISTRIBUTED_CRS_MATRIX_HPP

#include <rokko/slepc/mapping_1d.hpp>
#include <rokko/distributed_crs_matrix.hpp>
#include <rokko/utility/timer.hpp>

#include <rokko/slepc.hpp>

#include <numeric>
#include <iostream>

namespace rokko {
namespace slepc {

class distributed_crs_matrix : public rokko::detail::ps_crs_base {
public:
  distributed_crs_matrix() = default;
  ~distributed_crs_matrix() = default;

  explicit distributed_crs_matrix(rokko::mapping_1d const& map, int num_entries_per_row)
    : distributed_crs_matrix(cast_map(map), num_entries_per_row) {}

  explicit distributed_crs_matrix(rokko::slepc::mapping_1d const& map, int num_entries_per_row)
    : distributed_crs_matrix(std::make_shared<const rokko::slepc::mapping_1d>(map), num_entries_per_row) {}

  explicit distributed_crs_matrix(std::shared_ptr<const rokko::slepc::mapping_1d> map, int num_entries_per_row): map_(map) {
    int dim = map_->get_dim();
    PetscErrorCode ierr;
    ierr = MatCreate(map_->get_mpi_comm().get_comm(), &matrix_);  //CHKERRQ(ierr);
    ierr = MatSetSizes(matrix_, map_->get_num_local_rows(), map_->get_num_local_rows(), dim, dim);  //CHKERRQ(ierr);
    ierr = MatSetFromOptions(matrix_);  //CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(matrix_, num_entries_per_row, NULL);
    ierr = MatMPIAIJSetPreallocation(matrix_, num_entries_per_row, NULL, num_entries_per_row, NULL);
  }

  static std::shared_ptr<const rokko::slepc::mapping_1d> cast_map(rokko::mapping_1d const& map) {
    if (map.get_solver_name() != "slepc") {
      throw std::invalid_argument("SLEPc's distributed_crs_matrix() : " + map.get_solver_name() + "'s mapping_1d is given.");
    }
    return std::static_pointer_cast<const rokko::slepc::mapping_1d>(map.get_ptr());
  }

  #undef __FUNCT__
  #define __FUNCT__ "distributed_crs_matrix/insert"
  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) override {
    PetscErrorCode ierr = MatSetValues(matrix_, 1, &row, cols.size(), cols.data(), values.data(), INSERT_VALUES);  //CHKERRQ(ierr);
  }
  void insert(int row, int col_size, int const*const cols, double const*const values) override {
    PetscErrorCode ierr = MatSetValues(matrix_, 1, &row, col_size, cols, values, INSERT_VALUES);  //CHKERRQ(ierr);
  }
  #undef __FUNCT__
  #define __FUNCT__ "distributed_crs_matrix/complete"
  void complete() override {
    PetscErrorCode ierr;
    ierr = MatAssemblyBegin(matrix_, MAT_FINAL_ASSEMBLY);  //CHKERRQ(ierr);
    ierr = MatAssemblyEnd(matrix_, MAT_FINAL_ASSEMBLY);  //CHKERRQ(ierr);
  }
  Mat& get_matrix() {
    return matrix_;
  }
  const Mat& get_matrix() const {
    return matrix_;
  }
  int get_dim() const override {
    return map_->get_dim();
  }
  int get_num_local_rows() const override {
    int start_row, end_row;
    PetscErrorCode ierr = MatGetOwnershipRange(matrix_, &start_row, &end_row); //CHKERRQ(ierr);
    return end_row - start_row;
  }
  int start_row() const override {
    int start_row;
    PetscErrorCode ierr = MatGetOwnershipRange(matrix_, &start_row, NULL); //CHKERRQ(ierr);
    return start_row;
  }
  int end_row() const override {
    int end_row;
    PetscErrorCode ierr = MatGetOwnershipRange(matrix_, NULL, &end_row); //CHKERRQ(ierr);
    return end_row;
  }
  int get_nnz() const override {
    MatInfo info;
    MatGetInfo(matrix_, MAT_GLOBAL_SUM, &info);
    return static_cast<int>(info.nz_used);
  }
  void print() const override {
    MatView(matrix_, PETSC_VIEWER_STDOUT_(map_->get_mpi_comm().get_comm()));
  }
  void extract(int row, std::vector<int>& cols, std::vector<double>& values) const override {
    PetscInt num_cols;
    const PetscInt * cols_tmp;
    const PetscScalar * values_tmp;
    MatGetRow(matrix_, row, &num_cols, &cols_tmp, &values_tmp);
    cols.clear();
    cols.reserve(num_cols);
    values.clear();
    values.reserve(num_cols);
    std::vector<int> idx(num_cols);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&cols_tmp](auto i, auto j) { return cols_tmp[i] < cols_tmp[j]; });
    for (int i=0; i<num_cols; ++i) {
      cols.emplace_back(cols_tmp[idx[i]]);
      values.emplace_back(values_tmp[idx[i]]);
    }
    PetscErrorCode ierr = MatRestoreRow(matrix_, row, &num_cols, &cols_tmp, &values_tmp);
  }
  void output_matrix_market(std::ostream& os = std::cout) const override {
    const auto& comm = get_map().get_mpi_comm();
    constexpr int root_proc = 0;
    std::vector<int> cols;
    std::vector<double> values;

    const int nnz = get_nnz();
    if (comm.get_myrank() == root_proc) {
      os << "%%MatrixMarket matrix coordinate real general" << std::endl
         << get_dim() << " " << get_dim() << " " << nnz << std::endl;
    }
    comm.barrier();
    for (int global_row=0; global_row<get_dim(); ++global_row) {
      if ((global_row >= start_row()) && (global_row < end_row())) {
        extract(global_row, cols, values);
        for (int i=0; i<cols.size(); ++i) {
          os << global_row + 1 << " " << cols[i] + 1 << " " << values[i] << std::endl;
        }
      }
      comm.barrier();
    }
  }

  const rokko::slepc::mapping_1d& get_map() const override { return *map_; }

private:
  std::shared_ptr<const rokko::slepc::mapping_1d> map_;
  Mat matrix_;
};

} // namespace slepc
} // namespace rokko

#endif // ROKKO_SLEPC_DISTRIBUTED_CRS_MATRIX_HPP
